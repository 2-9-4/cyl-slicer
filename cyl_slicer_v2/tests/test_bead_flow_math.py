from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, LineString, box
from shapely.ops import unary_union

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import _to_metric, build_layer_toolpaths, inset_geometry


def rectangular_layer(occupied):
    return LayerSlice(0, 30.1, occupied, occupied, GeometryCollection(),
                      GeometryCollection(), GeometryCollection())


def boundary_lines(paths, radius):
    return unary_union([
        LineString([(theta * radius, axial) for theta, axial in p.points])
        for p in paths
    ])


class RoundedBeadFlowTests(unittest.TestCase):
    def setUp(self):
        self.config = SlicerConfig(
            layer_height=0.2, nozzle_width=0.4,
            bead_model="rounded", wall_count=3, min_path_length=0
        )

    def test_area_and_pitch_are_coupled(self):
        c = self.config
        area = (0.4 - 0.2) * 0.2 + math.pi * 0.2 ** 2 / 4
        pitch = area / 0.2
        self.assertAlmostEqual(area, 0.0714159265359, places=11)
        self.assertAlmostEqual(pitch, 0.3570796326795, places=11)
        self.assertAlmostEqual(c.bead_area, area)
        self.assertAlmostEqual(c.dense_line_spacing, pitch)
        filament_area = math.pi * (c.filament_diameter / 2) ** 2
        self.assertAlmostEqual(c.extrusion_per_mm, area / filament_area)
        self.assertAlmostEqual(c.bead_area / c.dense_line_spacing, c.layer_height)

    def test_legacy_rectangular_area_and_pitch_unchanged(self):
        c = SlicerConfig(layer_height=0.2, nozzle_width=0.4)
        self.assertEqual(c.bead_model, "rectangular")
        self.assertAlmostEqual(c.bead_area, 0.08)
        self.assertAlmostEqual(c.dense_line_spacing, 0.4)
        self.assertAlmostEqual(c.wall_coverage_depth(3), 1.2)

    def test_external_wall_unchanged_and_internal_wall_center_pitch(self):
        c = self.config
        region = box(0.0, 0.0, 1.0, 10.0)
        result = build_layer_toolpaths(rectangular_layer(region), c).paths_by_role
        outer = boundary_lines(result["EXTERNAL_PERIMETER"], 30.1)
        internal = boundary_lines(result["INTERNAL_PERIMETER"], 30.1)
        self.assertLess(outer.hausdorff_distance(
            _to_metric(inset_geometry(region, 30.1, c.nozzle_width / 2), 30.1).boundary
        ), 1e-5)
        want_internal = unary_union([
            _to_metric(inset_geometry(
                region, 30.1, c.nozzle_width / 2 + k * c.dense_line_spacing
            ), 30.1).boundary
            for k in (1, 2)
        ])
        self.assertLess(internal.hausdorff_distance(want_internal), 1e-5)

    def test_skin_clipping_and_pitches_use_same_area(self):
        c = self.config
        region = box(0.0, 0.0, 1.0, 10.0)
        result = build_layer_toolpaths(rectangular_layer(region), c).paths_by_role
        skin = result["SOLID_INFILL"]
        self.assertGreater(len(skin), 3)
        limit = c.wall_coverage_depth(c.wall_count) - (
            c.nozzle_width * c.infill_wall_overlap_fraction
        )
        for path in skin:
            for theta, axial in path.points:
                self.assertGreaterEqual(
                    min(theta * 30.1, (1 - theta) * 30.1, axial, 10 - axial),
                    limit - 1e-6
                )
        # The source generation grid has the same pitch as the extrusion
        # bead's volume-per-unit-height.
        centers = sorted({round(p.points[0][0] * 30.1, 7) for p in skin})
        self.assertGreater(len(centers), 2)
        for a, b in zip(centers, centers[1:]):
            self.assertAlmostEqual(b - a, c.dense_line_spacing, delta=1e-6)

    def test_perimeter_hole_is_offset_into_occupied_material(self):
        c = self.config
        model = box(.1, 0, 1.0, 10).difference(box(.3, 2, .7, 8))
        result = build_layer_toolpaths(rectangular_layer(model), c).paths_by_role
        outer = boundary_lines(result["EXTERNAL_PERIMETER"], 30.1)
        expected = _to_metric(inset_geometry(model, 30.1, .2), 30.1).boundary
        self.assertLess(outer.hausdorff_distance(expected), 1e-5)

    def test_invalid_rounded_dimensions_or_unknown_model_rejected(self):
        for k in (
            {"bead_model": "unknown", "layer_height": .2},
            {"bead_model": "rounded", "layer_height": .5},
            {"bead_model": "rounded", "layer_height": .2, "nozzle_width": .19},
        ):
            with self.subTest(k=k), self.assertRaises(ValueError):
                SlicerConfig(**k).validate()

    def test_gcode_metadata_contains_derived_bead_geometry(self):
        from cyl_slicer_v2.models import LayerToolpaths
        code = emit_gcode([LayerToolpaths(0, 30.1)], self.config)
        self.assertIn(";BEAD_MODEL:rounded", code)
        self.assertIn(";BEAD_AREA_MM2:0.071416", code)
        self.assertIn(";DENSE_SPACING_MM:0.357080", code)

    def test_one_mm_of_dense_line_has_correct_filament_volume(self):
        c = self.config
        # An actual G1 segment should carry material proportional to rounded
        # area, not the old rectangular cross-section.
        from cyl_slicer_v2.models import LayerToolpaths, Toolpath
        layer = LayerToolpaths(0, 30.1, {
            "EXTERNAL_PERIMETER": [Toolpath("EXTERNAL_PERIMETER", 30.1, [(0, 0), (0, 1)])]
        })
        gcode = emit_gcode([layer], c)
        import re
        e = re.findall(r"(?m)^G1 X[^\n]*? E([0-9.]+)", gcode)
        self.assertTrue(e)
        self.assertAlmostEqual(float(e[-1]), c.extrusion_per_mm, delta=0.00001)
