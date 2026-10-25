from __future__ import annotations

import math
import re
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import cylindrical_distance
from cyl_slicer_v2.gcode import emit_gcode, emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerSlice, LayerToolpaths, Toolpath
from cyl_slicer_v2.toolpaths import _make_path


class SegmentConservationTests(unittest.TestCase):
    def setUp(self):
        self.config = SlicerConfig(min_path_length=0.6, retract_mm=0.0, axial_min=0.0)

    def test_short_corner_survives_toolpath_construction(self):
        # 0.5 mm rotary leg followed by 1 mm axial leg; a previous
        # distance-only filter removed the right-angle vertex.
        p = _make_path("INFILL_RADIAL", 30.0, [(0, 0), (0.5 / 30, 0), (0.5 / 30, 1)],
                       False, self.config)
        self.assertIsNotNone(p)
        self.assertEqual(len(p.points), 3)
        self.assertAlmostEqual(sum(cylindrical_distance(a, b, p.radius)
                                   for a, b in zip(p.points, p.points[1:])), 1.5)

    def test_very_short_whole_support_feature_is_filtered(self):
        self.assertIsNone(_make_path("SUPPORT", 30, [(0, 0), (0.5 / 30, 0)],
                                     False, self.config))
        self.assertIsNotNone(_make_path("EXTERNAL_PERIMETER", 30,
                                        [(0, 0), (0.5 / 30, 0)], False, self.config))

    def test_direct_gcode_keeps_short_leg_and_correct_e(self):
        p = Toolpath("SUPPORT", 30, [(0, 0), (0.5 / 30, 0), (0.5 / 30, 1.0)])
        state = GcodeState(radial=30)
        lines = emit_toolpath(p, state, self.config, 0.0, retract_after=False)
        motions = [line for line in lines if line.startswith("G1 X")]
        self.assertEqual(len(motions), 2, lines)
        self.assertIn("Y0.955", motions[0])  # 0.5 mm arc at R30.
        self.assertIn("X1.000", motions[1])
        self.assertAlmostEqual(state.extrusion, 1.5 * self.config.extrusion_per_mm, places=9)

    def test_short_corner_is_not_skipped_within_long_diagonal(self):
        p = Toolpath("INFILL_AXIAL", 30, [(0, 0), (0.5 / 30, 0),
                                         (0.5 / 30, 1), (0.5 / 30, 2)])
        state = GcodeState(radial=30)
        code = emit_toolpath(p, state, self.config, 0, retract_after=False)
        moves = [s for s in code if s.startswith("G1 X")]
        self.assertEqual(len(moves), 2)  # Second + third axial legs coalesce.
        self.assertAlmostEqual(state.extrusion, 2.5 * self.config.extrusion_per_mm, places=8)

    def test_collinear_submillimeter_segments_coalesce_without_volume_loss(self):
        points = [(0, 0), (0.005, 0), (0.01, 0), (0.03, 0)]
        p = _make_path("SUPPORT", 30, points, False, self.config)
        self.assertIsNotNone(p)
        self.assertEqual(len(p.points), 2)
        self.assertAlmostEqual(cylindrical_distance(*p.points, p.radius), 0.9)

    def test_short_closed_reverse_path_preserves_total_material(self):
        # A tiny corner followed by a reversal must not simplify to a diagonal.
        points = [(0, 0), (0.01, 0), (0.01, 0.05), (0.01, 0)]
        p = _make_path("INFILL_AXIAL", 30, points, False, self.config)
        self.assertIsNotNone(p)
        self.assertEqual(len(p.points), 4)
        expected = sum(cylindrical_distance(a, b, 30) for a, b in zip(points, points[1:]))
        s = GcodeState(radial=30)
        emit_toolpath(p, s, self.config, 0, retract_after=False)
        self.assertAlmostEqual(s.extrusion, expected * self.config.extrusion_per_mm, places=7)

    def test_direct_short_support_path_is_not_emitted(self):
        p = Toolpath("SUPPORT", 30, [(0, 0), (0.5 / 30, 0)])
        lines = emit_toolpath(p, GcodeState(radial=30), self.config, 0,
                              retract_after=False)
        self.assertEqual(lines, [])
