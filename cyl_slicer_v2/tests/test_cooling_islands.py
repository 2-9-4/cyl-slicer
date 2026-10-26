from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (
    _cooling_print_speed, _island_cooling_print_speed,
    _ordered_layer_paths, _ordered_layer_paths_with_islands, emit_gcode,
)
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


def group(role: str, radius: float, theta: float, mm: float):
    return {role: [Toolpath(role, radius, [(theta, 0), (theta, mm)])]}


def layer_with_large_and_small():
    radius = 31.0
    a = group("EXTERNAL_PERIMETER", radius, .1, 400)
    b = group("EXTERNAL_PERIMETER", radius, 1.1, 14)
    return LayerToolpaths(0, radius,
        {"EXTERNAL_PERIMETER": [*a["EXTERNAL_PERIMETER"], *b["EXTERNAL_PERIMETER"]]},
        island_paths=[a, b])


def extrusion_moves(code):
    return [line for line in code.splitlines()
            if line.startswith("G1 X") and " E" in line]


class IslandCoolingTests(unittest.TestCase):
    def test_default_quality_and_legacy_opt_out(self):
        self.assertEqual(SlicerConfig().cooling_min_island_print_time_s, 2)
        self.assertEqual(_profile_value(build_parser().parse_args([]),
                         "cooling_min_island_print_time_s"), 2)
        self.assertEqual(_profile_value(build_parser().parse_args(
            ["--profile", "legacy-cli"]), "cooling_min_island_print_time_s"), 0)
        self.assertEqual(_profile_value(build_parser().parse_args(
            ["--cooling-min-island-print-time-s", "3.5"]),
            "cooling_min_island_print_time_s"), 3.5)

    def test_short_island_gets_own_speed_when_long_island_sets_layer_time(self):
        c = SlicerConfig(bed_radius=30, cooling_min_layer_time_s=8,
                         cooling_min_island_print_time_s=2,
                         cooling_min_print_speed_mm_per_min=300, retract_mm=0)
        layer = layer_with_large_and_small()
        layer_feed, _ = _cooling_print_speed(layer, c)
        self.assertEqual(layer_feed, c.feed_mm_per_min)
        paths = _ordered_layer_paths_with_islands(layer, GcodeState(), c)
        self.assertEqual([p for p, _ in paths], _ordered_layer_paths(layer,GcodeState(),c))
        code = emit_gcode([layer], c)
        moves = extrusion_moves(code)
        self.assertEqual(len(moves), 2)
        by_len = {round(float(re.search(r"X([\d.]+)", line).group(1)), 1): line
                  for line in moves}
        self.assertIn(" F2700.0", by_len[415.0])
        self.assertIn(" F420.0", by_len[29.0])
        self.assertEqual(code.count(";COOLING_ISLAND:"), 2)
        self.assertIn("SPEED_MM_MIN:420.0000", code)

    def test_lower_bound_duration_uses_cylindrical_mm_and_role_speed_caps(self):
        radius = 60.
        linear = [Toolpath("SOLID_INFILL",radius,[(0,0),(0,12)])]
        rotary = [Toolpath("SOLID_INFILL",radius,[(0,0),(.2,0)])]
        c=SlicerConfig(cooling_min_island_print_time_s=2,
            cooling_min_layer_time_s=0,
            cooling_min_print_speed_mm_per_min=300)
        for paths in (linear, rotary):
            speed,duration=_island_cooling_print_speed(paths,c)
            self.assertAlmostEqual(speed,360,places=4)
            self.assertAlmostEqual(duration,2,places=6)
        support=[Toolpath("SUPPORT_INTERFACE",radius,[(0,0),(0,30)])]
        c2=SlicerConfig(cooling_min_island_print_time_s=2,
            support_interface_print_feed_mm_per_min=600)
        speed,duration=_island_cooling_print_speed(support,c2)
        self.assertEqual(speed, c2.feed_mm_per_min)
        self.assertAlmostEqual(duration,3)

    def test_cannot_force_impossible_time_below_configured_floor(self):
        c=SlicerConfig(cooling_min_island_print_time_s=10,
            cooling_min_print_speed_mm_per_min=300)
        tiny=[Toolpath("TOP_SKIN",31,[(0,0),(0,2)])]
        speed,duration=_island_cooling_print_speed(tiny,c)
        self.assertEqual(speed,300)
        self.assertAlmostEqual(duration,.4)
        code=emit_gcode([LayerToolpaths(0,31,{"TOP_SKIN":tiny})],
                        SlicerConfig(bed_radius=30,retract_mm=0,
                            cooling_min_island_print_time_s=10))
        self.assertNotIn("\nG4",code)
        self.assertIn(" F300.0",extrusion_moves(code)[0])

    def test_no_island_mode_keeps_existing_gcode_motion(self):
        layer=layer_with_large_and_small()
        base=SlicerConfig(bed_radius=30,retract_mm=0,
            cooling_min_island_print_time_s=0,cooling_min_layer_time_s=0)
        alternate=SlicerConfig(bed_radius=30,retract_mm=0,
            cooling_min_island_print_time_s=2,cooling_min_layer_time_s=0)
        a=emit_gcode([layer],base)
        b=emit_gcode([layer],alternate)
        self.assertNotIn(";COOLING_ISLAND:",a)
        a_e=[line.split(" F")[0] for line in extrusion_moves(a)]
        b_e=[line.split(" F")[0] for line in extrusion_moves(b)]
        self.assertEqual(a_e,b_e,"cooling must not change path geometry or filament")
        self.assertEqual(a.count("M106 S"),b.count("M106 S"))

    def test_bridge_fan_and_feed_caps_survive_island_slowdown(self):
        p=Toolpath("BRIDGE",31,[(0,0),(0,24)])
        layer=LayerToolpaths(0,31,{"BRIDGE":[p]})
        c=SlicerConfig(bed_radius=30,retract_mm=0,
            bridge_mode="print",bridge_fan_pwm=170,
            bridge_print_feed_mm_per_min=900,
            cooling_min_layer_time_s=0,
            cooling_min_island_print_time_s=2)
        speed,duration=_island_cooling_print_speed([p],c)
        self.assertAlmostEqual(speed,720,places=5)
        self.assertAlmostEqual(duration,2)
        code=emit_gcode([layer],c)
        self.assertIn("M106 S170 ; bridge fan override",code)
        self.assertIn("restore layer fan after bridge",code)
        self.assertIn(" F720.0",extrusion_moves(code)[0])

    def test_invalid_input(self):
        for value in (-1,math.nan,math.inf):
            with self.subTest(value=value),self.assertRaises(ValueError):
                SlicerConfig(cooling_min_island_print_time_s=value).validate()


if __name__=="__main__":
    unittest.main()
