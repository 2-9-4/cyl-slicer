from __future__ import annotations

from dataclasses import replace
import math
import re
import unittest

from cyl_slicer_v2.cli import build_parser, _profile_value
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.cooling_timeline import annotate_revisit_timeline
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def layer(index, angle=.1, support=False):
    role="SUPPORT" if support else "SOLID_INFILL"
    radius=30.2+index*.16
    path=Toolpath(role, radius, [(angle,0),(angle,10)])
    return LayerToolpaths(index,radius,{role:[path]},island_paths=[{role:[path]}])


def motion_commands(code):
    return [line for line in code.splitlines()
            if re.match(r"^(G[01]|G92(?:\.9)?|M106|M204|M201|M900)\b",line)]


class CoolingRevisitTimelineTests(unittest.TestCase):
    def cfg(self, **k):
        return replace(SlicerConfig(
            bed_radius=30, retract_mm=0, cooling_min_layer_time_s=0,
            cooling_min_island_print_time_s=0, cooling_revisit_report=True
        ),**k)

    def test_exact_native_marlin_feed_and_coordinate_resets(self):
        code = (
            ";LAYER:0\n;SAMPLE_RADIUS:30\n"
            ";COOLING_REVISIT_GROUP:0 THETA_RAD:0.1 AXIAL_MM:3 SUPPORT:0\n"
            "G0 X0 Y0 F600\nG1 X10 Y0 E1 F600\n"
            ";LAYER:1\n;SAMPLE_RADIUS:30.16\n"
            "G92.9 Y45\nG92 Y0\nG92 E0\n"
            "G0 Y90 F5400\n"
            ";COOLING_REVISIT_GROUP:1 THETA_RAD:6.383185307180 AXIAL_MM:3 SUPPORT:0\n"
            "G1 X20 Y90 E1 F600\n"
        )
        result=annotate_revisit_timeline(code,match_mm=2,target_s=.8)
        # 10 physical X mm at F600 = 1s; rotary 90 *native
        # degrees* at F5400 = 1s, not 90 *radius radians*.
        self.assertIn(";COOLING_REVISIT_LOWER_BOUND_S:1.0000 GROUP:1 PRIOR_LAYER:0",result)
        self.assertIn(";COOLING_REVISIT_TARGET_CERTIFIED_LOWER_BOUND",result)
        self.assertEqual(motion_commands(code),motion_commands(result))

    def test_retractions_count_as_ideal_e_motion_but_no_material(self):
        code = (
            ";LAYER:0\n;SAMPLE_RADIUS:30\n"
            ";COOLING_REVISIT_GROUP:0 THETA_RAD:.1 AXIAL_MM:0 SUPPORT:0\n"
            "G0 X0 Y0 F600\nG1 X10 E1 F600\n"
            "G1 E.5 F60\nG0 X20 F600\nG1 E1 F60\n"
            ";LAYER:1\n;SAMPLE_RADIUS:30.16\nG92 E0\n"
            ";COOLING_REVISIT_GROUP:1 THETA_RAD:6.38318530718 AXIAL_MM:0 SUPPORT:0\n"
            "G1 X30 E1 F600\n"
        )
        result=annotate_revisit_timeline(code,match_mm=2,target_s=3)
        # 0.5 mm retract at 60 mm/min, 10 mm travel at 600 and
        # 0.5 mm prime at 60 = 2 seconds between deposits.
        self.assertIn("LOWER_BOUND_S:2.0000",result)
        self.assertIn("TARGET_NOT_PROVEN",result)
        self.assertNotIn("TARGET_CERTIFIED_LOWER_BOUND",result)

    def test_support_island_separate_from_model_and_different_phase(self):
        code = emit_gcode([layer(0),layer(1,support=True),layer(2,angle=.1+2*math.pi)],
                          self.cfg(cooling_min_revisit_time_s=2))
        annotations = [x for x in code.splitlines()
                       if x.startswith(";COOLING_REVISIT_LOWER_BOUND_S:")]
        self.assertEqual(len(annotations),1)
        self.assertIn("PRIOR_LAYER:0",annotations[0])
        self.assertNotIn("PRIOR_LAYER:1",annotations[0])
        self.assertNotIn("\nG4",code)

    def test_opt_in_never_changes_motion_or_deposited_e(self):
        layers=[layer(0),layer(1)]
        default=emit_gcode(layers,self.cfg(cooling_revisit_report=False))
        report=emit_gcode(layers,self.cfg(cooling_min_revisit_time_s=5))
        self.assertEqual(motion_commands(default),motion_commands(report))
        self.assertIn(";COOLING_REVISIT_LOWER_BOUND_S:",report)
        self.assertIn(";COOLING_REVISIT_TARGET_NOT_PROVEN:",report)

    def test_new_or_full_wrap_ambiguous_islands_not_falsely_matched(self):
        code=emit_gcode([layer(0),layer(1,angle=2.0)],
                        self.cfg(cooling_revisit_match_mm=1))
        self.assertNotIn(";COOLING_REVISIT_LOWER_BOUND_S:",code)

    def test_cli_and_invalid_settings(self):
        for profile in ("quality","legacy-cli","library"):
            parser=build_parser()
            default=parser.parse_args(["--profile",profile])
            self.assertFalse(_profile_value(default,"cooling_revisit_report"))
            self.assertEqual(_profile_value(default,"cooling_min_revisit_time_s"),0)
            args=parser.parse_args(["--profile",profile,
                "--cooling-revisit-report","--cooling-min-revisit-time-s","7"])
            self.assertTrue(_profile_value(args,"cooling_revisit_report"))
            self.assertEqual(_profile_value(args,"cooling_min_revisit_time_s"),7)
        for fields in (
            {"cooling_revisit_report":1},
            {"cooling_min_revisit_time_s":-1},
            {"cooling_min_revisit_time_s":math.nan},
        ):
            with self.subTest(fields=fields),self.assertRaises(ValueError):
                self.cfg(**fields).validate()
        with self.assertRaisesRegex(ValueError,"relative coordinates"):
            annotate_revisit_timeline("G91\nG1 X10 F600\n",match_mm=5)


if __name__=="__main__":
    unittest.main()
