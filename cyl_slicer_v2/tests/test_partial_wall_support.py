from __future__ import annotations

import math
import re
import unittest

from shapely.geometry import box, GeometryCollection

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import Toolpath, LayerToolpaths
from cyl_slicer_v2.wall_support import measure_wall_support


class PartialWallSupportTests(unittest.TestCase):
    def test_axial_partly_supported_and_centerline_mm(self):
        p = Toolpath("EXTERNAL_PERIMETER", 30, [(.1, 0), (.1, 10)])
        previous = box(.09, 0, .11, 5)
        stats = measure_wall_support(p, previous, .4)
        self.assertEqual(stats.category, "partial")
        self.assertAlmostEqual(stats.total_length_mm, 10)
        self.assertAlmostEqual(stats.covered_fraction, .52, delta=.01)
        self.assertAlmostEqual(stats.max_unsupported_run_mm, 4.8, delta=.01)

    def test_rotary_segment_uses_arc_not_radians(self):
        radius = 60.
        p = Toolpath("EXTERNAL_PERIMETER", radius, [(0., 2.), (.2, 2.)])
        prior = box(-.01, 1, .10, 3)
        stats = measure_wall_support(p, prior, .4)
        self.assertAlmostEqual(stats.total_length_mm, 12.)
        self.assertAlmostEqual(stats.covered_fraction, .1 / .2 + .2/12., delta=.003)
        self.assertEqual(stats.category, "partial")

    def test_seam_and_rotary_turn_shift_are_periodic(self):
        prior = box(0., -1., .07, 1.)
        p = Toolpath("EXTERNAL_PERIMETER", 30,
                     [(2*math.pi + .01, 0), (2*math.pi+.05, 0)])
        for turns in (-5, 0, 8):
            shifted = Toolpath(p.role, p.radius,
                [(theta+turns*2*math.pi, axial) for theta,axial in p.points])
            stats = measure_wall_support(shifted, prior, .4)
            self.assertEqual(stats.category, "supported")
            self.assertAlmostEqual(stats.covered_fraction, 1.)

    def test_closed_loop_joins_unsupported_runs_across_start(self):
        p = Toolpath("EXTERNAL_PERIMETER", 30,
                     [(0, 0), (0, 10), (.1, 10), (.1, 0), (0, 0)], closed=True)
        # Only near the top edge is supported, so the remaining run crosses
        # the loop's arbitrary start/end and is measured as a single run.
        prior = box(-.02, 7.5, .12, 11)
        stats = measure_wall_support(p, prior, .4)
        self.assertEqual(stats.category, "partial")
        self.assertGreater(stats.max_unsupported_run_mm, 14.)

    def test_no_prior_is_entirely_unsupported(self):
        p = Toolpath("EXTERNAL_PERIMETER", 30, [(0, 0), (0, 5)])
        stats = measure_wall_support(p, GeometryCollection(), .4)
        self.assertEqual(stats.category, "unsupported")
        self.assertEqual(stats.max_unsupported_run_mm, 5.)

    def test_gcode_only_comments_opt_in_preserves_geometry_and_e(self):
        p = Toolpath("EXTERNAL_PERIMETER", 30.5, [(.1, 0), (.1, 10)])
        layer = LayerToolpaths(1, 30.5, {"EXTERNAL_PERIMETER": [p]},
                               previous_occupied=box(.09, 0, .11, 5))
        base = SlicerConfig(retract_mm=0, cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0)
        original = emit_gcode([layer], base)
        annotated = emit_gcode([layer], SlicerConfig(
            retract_mm=0, cooling_min_layer_time_s=0,
            cooling_min_island_print_time_s=0, wall_support_report=True))
        self.assertIn(";WALL_SUPPORT:partial COVERED_FRACTION:", annotated)
        only_commands = lambda code: [
            line for line in code.splitlines()
            if line and not line.startswith(";")
        ]
        self.assertEqual(only_commands(original), only_commands(annotated))

    def test_default_and_validation(self):
        self.assertFalse(SlicerConfig().wall_support_report)
        args = build_parser().parse_args(["--wall-support-report"])
        self.assertTrue(_profile_value(args, "wall_support_report"))
        args = build_parser().parse_args(["--no-wall-support-report"])
        self.assertFalse(_profile_value(args, "wall_support_report"))
        with self.assertRaises(ValueError):
            SlicerConfig(wall_support_report="yes").validate()


if __name__ == "__main__":
    unittest.main()
