from __future__ import annotations

import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


def part_layer(index, solid=False):
    source = box(0, 0, 1, 10)
    blank = GeometryCollection()
    return LayerSlice(index, 30 + index, source, source if solid else blank,
                      blank, blank, blank)


class AlternatingInfillTests(unittest.TestCase):
    def test_legacy_both_axes_in_every_sparse_layer(self):
        c = SlicerConfig(wall_count=1, min_path_length=0)
        for index in (0, 1, 2):
            paths = build_layer_toolpaths(part_layer(index), c).paths_by_role
            self.assertTrue(paths.get("INFILL_RADIAL"))
            self.assertTrue(paths.get("INFILL_AXIAL"))

    def test_sparse_alternates_model_layers_not_support_roles(self):
        c = SlicerConfig(wall_count=1, min_path_length=0, sparse_infill_pattern="alternate")
        for index in (0, 1, 2, 3):
            paths = build_layer_toolpaths(part_layer(index), c).paths_by_role
            self.assertEqual(bool(paths.get("INFILL_RADIAL")), index % 2 == 0)
            self.assertEqual(bool(paths.get("INFILL_AXIAL")), index % 2 == 1)
            self.assertTrue(paths.get("EXTERNAL_PERIMETER"))
            self.assertEqual(bool(paths.get("INTERNAL_PERIMETER")), False)

    def test_solid_lines_alternate_and_keep_spacing(self):
        c = SlicerConfig(wall_count=1, min_path_length=0, solid_infill_pattern="alternate")
        for index in (0, 1):
            paths = build_layer_toolpaths(part_layer(index, solid=True), c).paths_by_role
            solid = paths.get("SOLID_INFILL", [])
            self.assertTrue(solid)
            for path in solid:
                if index == 0:  # Lines along axial X; theta fixed.
                    self.assertAlmostEqual(path.points[0][0], path.points[-1][0], places=6)
                    self.assertNotAlmostEqual(path.points[0][1], path.points[-1][1], places=6)
                else:  # Lines along rotary direction; X fixed.
                    self.assertAlmostEqual(path.points[0][1], path.points[-1][1], places=6)
                    self.assertNotAlmostEqual(path.points[0][0], path.points[-1][0], places=6)
            self.assertFalse(paths.get("INFILL_RADIAL"))
            self.assertFalse(paths.get("INFILL_AXIAL"))

    def test_default_solid_skin_remains_axial(self):
        c = SlicerConfig(wall_count=1, min_path_length=0)
        paths = build_layer_toolpaths(part_layer(1, solid=True), c).paths_by_role
        self.assertTrue(paths.get("SOLID_INFILL"))
        for path in paths["SOLID_INFILL"]:
            self.assertAlmostEqual(path.points[0][0], path.points[-1][0], places=6)

    def test_gcode_roles_follow_alternating_patterns(self):
        c = SlicerConfig(wall_count=1, min_path_length=0,
                         sparse_infill_pattern="alternate", bed_radius=29)
        layers = [build_layer_toolpaths(part_layer(i), c) for i in range(2)]
        code = emit_gcode(layers, c)
        first = code.split(";LAYER:0")[1].split(";LAYER:1")[0]
        second = code.split(";LAYER:1")[1]
        self.assertIn(";TYPE:INFILL_RADIAL", first)
        self.assertNotIn(";TYPE:INFILL_AXIAL", first)
        self.assertIn(";TYPE:INFILL_AXIAL", second)
        self.assertNotIn(";TYPE:INFILL_RADIAL", second)

    def test_validation(self):
        for k in ({"sparse_infill_pattern": "zigzag"},
                  {"solid_infill_pattern": "random"}):
            with self.subTest(k=k), self.assertRaises(ValueError):
                SlicerConfig(**k).validate()
