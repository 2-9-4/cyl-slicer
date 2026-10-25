from __future__ import annotations

import math
import unittest

from shapely.geometry import box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _ordered_layer_paths, emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


def rectangle(role, theta_start, x_start, radius=30):
    t, x = theta_start, x_start
    return Toolpath(role, radius, [(t, x), (t + .1, x), (t + .1, x + 2),
                                  (t, x + 2), (t, x)], closed=True)


class WallOrderTests(unittest.TestCase):
    def setUp(self):
        self.e0 = rectangle("EXTERNAL_PERIMETER", 0.1, 0)
        self.i0 = rectangle("INTERNAL_PERIMETER", 0.11, .2)
        self.e1 = rectangle("EXTERNAL_PERIMETER", 1.1, 20)
        self.i1 = rectangle("INTERNAL_PERIMETER", 1.11, 20.2)
        self.s0 = Toolpath("SOLID_INFILL", 30, [(0.15, .5), (0.15, 1.5)])
        self.support = Toolpath("SUPPORT", 30, [(2.0, 5), (2.0, 7)])
        self.layer = LayerToolpaths(
            0, 30,
            {"EXTERNAL_PERIMETER": [self.e0, self.e1],
             "INTERNAL_PERIMETER": [self.i0, self.i1],
             "SOLID_INFILL": [self.s0], "SUPPORT": [self.support]},
            island_paths=[
                {"EXTERNAL_PERIMETER": [self.e0], "INTERNAL_PERIMETER": [self.i0],
                 "SOLID_INFILL": [self.s0]},
                {"EXTERNAL_PERIMETER": [self.e1], "INTERNAL_PERIMETER": [self.i1]},
                {"SUPPORT": [self.support]},
            ],
            travel_region=box(0, -1, 3, 25),
        )

    def test_inner_first_stays_inside_island(self):
        for mode, first in (("outer-first", self.e0), ("inner-first", self.i0)):
            with self.subTest(mode=mode):
                seq = _ordered_layer_paths(self.layer, GcodeState(theta_accum=0, axial=0),
                                           SlicerConfig(wall_order=mode))
                roles = [p.role for p in seq]
                self.assertEqual(roles[:3], (
                    ["EXTERNAL_PERIMETER", "INTERNAL_PERIMETER", "SOLID_INFILL"]
                    if mode == "outer-first"
                    else ["INTERNAL_PERIMETER", "EXTERNAL_PERIMETER", "SOLID_INFILL"]
                ))
                self.assertEqual(roles.count("SUPPORT"), 1)
                self.assertEqual(roles.count("EXTERNAL_PERIMETER"), 2)
                self.assertEqual(roles.count("INTERNAL_PERIMETER"), 2)

    def test_full_gcode_keeps_all_extrusion_and_support(self):
        totals = []
        for mode in ("outer-first", "inner-first"):
            code = emit_gcode([self.layer],
                              SlicerConfig(wall_order=mode, bed_radius=29,
                                           layer_height=.5, retract_mm=0,
                                           axial_min=0))
            self.assertIn(f";WALL_ORDER:{mode}", code)
            sequence = [s for s in code.splitlines() if s.startswith(";TYPE:")]
            self.assertEqual(sequence.count(";TYPE:SUPPORT"), 1)
            self.assertEqual(sequence.count(";TYPE:INTERNAL_PERIMETER"), 2)
            self.assertEqual(sequence.count(";TYPE:EXTERNAL_PERIMETER"), 2)
            self.assertEqual(sequence.count(";TYPE:SOLID_INFILL"), 1)
            self.assertEqual(sequence[0], ";TYPE:EXTERNAL_PERIMETER"
                             if mode == "outer-first" else ";TYPE:INTERNAL_PERIMETER")
            # Motion ordering must not change total deposited physical length.
            segments = parse_extrusion_segments(code)
            lengths = [math.hypot(x1 - x0, 30 * math.radians(y1 - y0))
                       for _, _, x0, y0, x1, y1 in segments]
            totals.append(sum(lengths))
        self.assertAlmostEqual(totals[0], totals[1], places=5)

    def test_only_one_wall_in_legacy_default_and_two_policy_modes(self):
        p = self.e0
        one = LayerToolpaths(0, 30, {"EXTERNAL_PERIMETER": [p]})
        a = _ordered_layer_paths(one, GcodeState(), SlicerConfig(wall_order="outer-first"))
        b = _ordered_layer_paths(one, GcodeState(), SlicerConfig(wall_order="inner-first"))
        self.assertEqual([x.points for x in a], [x.points for x in b])

    def test_validation(self):
        with self.assertRaisesRegex(ValueError, "wall_order"):
            SlicerConfig(wall_order="random").validate()
