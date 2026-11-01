from __future__ import annotations

from dataclasses import replace
import math
import re
import unittest

from shapely.geometry import box

from cyl_slicer_v2.cli import build_parser, _profile_value
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode, _surface_feed_for_marlin
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def layer(*, supported=True):
    radius = 31.0
    one = Toolpath("EXTERNAL_PERIMETER", radius, [(.1, 0), (.1, 10)])
    two = Toolpath("INTERNAL_PERIMETER", radius, [(.3, 0), (.3, 10)])
    footprint = (box(.09, 0, .11, 5).union(box(.29, 0, .31, 10))
                 if supported else box(1, 0, 1.1, 10))
    return LayerToolpaths(1, radius, {
        one.role: [one], two.role: [two]
    }, previous_occupied=footprint)


def deposition_lines(code):
    return [s for s in code.splitlines() if s.startswith("G1 X") and " E" in s]


def e_and_motion_without_speed(code):
    return [s.split(" F")[0] for s in deposition_lines(code)]


class UnsupportedWallFeedTests(unittest.TestCase):
    def cfg(self, **changes):
        return replace(SlicerConfig(
            bed_radius=30, retract_mm=0,
            cooling_min_layer_time_s=0,
            cooling_min_island_print_time_s=0,
            wall_unsupported_print_feed_mm_per_min=600,
            wall_unsupported_min_span_mm=2,
        ), **changes)

    def test_only_partially_supported_external_wall_gets_cap(self):
        original = emit_gcode([layer()], self.cfg(
            wall_unsupported_print_feed_mm_per_min=None))
        result = emit_gcode([layer()], self.cfg())
        self.assertEqual(e_and_motion_without_speed(original),
                         e_and_motion_without_speed(result))
        # Internal wall is fully supported by the previous layer.
        self.assertEqual(result.count(";WALL_UNSUPPORTED_FEED_CAP:"), 1)
        self.assertIn("FEED_MM_MIN:600", result)
        lines = result.split(";TYPE:")
        per_role = {piece.splitlines()[0]:
                    [m for m in piece.splitlines() if m.startswith("G1 X") and " E" in m]
                    for piece in lines[1:]}
        self.assertEqual(len(per_role["EXTERNAL_PERIMETER"]), 1)
        self.assertIn(" F600.0", per_role["EXTERNAL_PERIMETER"][0])
        self.assertIn(" F2700.0", per_role["INTERNAL_PERIMETER"][0])
        self.assertEqual([s for s in original.splitlines() if s.startswith("M106 ")],
                         [s for s in result.splitlines() if s.startswith("M106 ")])

    def test_min_span_and_existing_slower_feed_take_precedence(self):
        normal = emit_gcode([layer()], self.cfg(wall_unsupported_min_span_mm=6))
        self.assertNotIn(";WALL_UNSUPPORTED_FEED_CAP:",normal)
        slower = emit_gcode([layer()], self.cfg(feed_mm_per_min=350))
        self.assertNotIn(";WALL_UNSUPPORTED_FEED_CAP:",slower)
        self.assertIn(" F350.0",slower)
        missing = layer()
        missing.previous_occupied = None
        self.assertNotIn(";WALL_UNSUPPORTED_FEED_CAP:",emit_gcode([missing],self.cfg()))

    def test_rotary_distance_uses_physical_arc_and_cap_preserves_e(self):
        radius = 60.
        p = Toolpath("EXTERNAL_PERIMETER",radius,[(0,0),(.2,0)])
        lp = LayerToolpaths(1,radius,{"EXTERNAL_PERIMETER":[p]},
                           previous_occupied=box(-.01,-1,.1,1))
        cfg = self.cfg(bed_radius=59,wall_unsupported_print_feed_mm_per_min=480,
                       wall_unsupported_min_span_mm=3)
        g=emit_gcode([lp],cfg)
        self.assertIn(";WALL_UNSUPPORTED_FEED_CAP:",g)
        output=deposition_lines(g)[0]
        # Planner feed is corrected for native rotary degrees, not left
        # at an apparent machine speed of 480 degrees/min.
        actual=float(re.search(r" F([0-9.]+)",output).group(1))
        self.assertAlmostEqual(actual,
            _surface_feed_for_marlin(0,.2,radius,480),delta=.05)
        base=emit_gcode([lp],self.cfg(bed_radius=59,
                       wall_unsupported_print_feed_mm_per_min=None))
        self.assertEqual(e_and_motion_without_speed(base),e_and_motion_without_speed(g))

    def test_invalid_inputs_and_cli_profiles(self):
        for profile in ("quality","legacy-cli","library"):
            args=build_parser().parse_args(["--profile",profile])
            self.assertIsNone(_profile_value(args,"wall_unsupported_print_feed_mm_per_min"))
            args=build_parser().parse_args(["--profile",profile,
                "--wall-unsupported-print-feed","900",
                "--wall-unsupported-min-span-mm","3"])
            self.assertEqual(_profile_value(args,"wall_unsupported_print_feed_mm_per_min"),900)
            self.assertEqual(_profile_value(args,"wall_unsupported_min_span_mm"),3)
        for bad in ({"wall_unsupported_print_feed_mm_per_min":0},
                    {"wall_unsupported_print_feed_mm_per_min":math.nan},
                    {"wall_unsupported_min_span_mm":-1},
                    {"wall_unsupported_min_span_mm":math.inf}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.cfg(**bad).validate()


if __name__=="__main__":
    unittest.main()
