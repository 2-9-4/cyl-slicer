from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, LineString, box
from shapely.ops import unary_union

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.coordinates import TWO_PI
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, inset_geometry, _to_metric


def _layer(geom, radius=30, index=0, periodic=False):
    return LayerSlice(index, radius, geom, GeometryCollection(), GeometryCollection(),
                      GeometryCollection(), GeometryCollection(),
                      theta_window=(0, TWO_PI) if periodic else (0, 0),
                      is_full_wrap=periodic)


def _linework(paths, radius):
    return unary_union([LineString([(theta * radius, x) for theta, x in path.points])
                        for path in paths])


class BeadContourTests(unittest.TestCase):
    def test_outer_and_hole_boundaries_are_both_shifted_half_width_into_material(self):
        model = box(0.1, 0, 1.0, 10).difference(box(0.3, 2, 0.7, 8))
        config = SlicerConfig(nozzle_width=0.4, wall_count=1)
        result = build_layer_toolpaths(_layer(model), config)
        actual = _linework(result.paths_by_role["EXTERNAL_PERIMETER"], 30)
        expected = _to_metric(inset_geometry(model, 30, 0.2), 30).boundary
        self.assertLess(actual.hausdorff_distance(expected), 1e-5)
        # Swept nominal width must not extend beyond the intended occupied
        # patch or enter its hole.
        bead = actual.buffer(0.2, cap_style="flat", join_style="mitre")
        model_metric = _to_metric(model, 30)
        self.assertLess(bead.difference(model_metric).area, 1e-4)

    def test_multiple_walls_have_half_then_full_width_offsets(self):
        model = box(0, 0, 1.0, 10)
        config = SlicerConfig(nozzle_width=0.4, wall_count=3)
        result = build_layer_toolpaths(_layer(model), config)
        paths = result.paths_by_role
        actual = [_linework(paths["EXTERNAL_PERIMETER"], 30)]
        # Internal role carries both inner loops; test per wall through their
        # union rather than assuming their ordering across disconnected rings.
        actual.append(_linework(paths["INTERNAL_PERIMETER"], 30))
        expected_outer = _to_metric(inset_geometry(model, 30, 0.2), 30).boundary
        expected_inner = unary_union([
            _to_metric(inset_geometry(model, 30, offset), 30).boundary
            for offset in (0.6, 1.0)])
        self.assertLess(actual[0].hausdorff_distance(expected_outer), 1e-5)
        self.assertLess(actual[1].hausdorff_distance(expected_inner), 1e-5)

    def test_solid_fill_avoids_wall_centerline_retrace(self):
        model = box(0, 0, 1.0, 10)
        config = SlicerConfig(nozzle_width=0.4, wall_count=3,
                              infill_wall_overlap_fraction=0.15, min_path_length=0)
        layer = _layer(model)
        layer.top_surfaces = model
        paths = build_layer_toolpaths(layer, config).paths_by_role
        fill = paths["SOLID_INFILL"]
        self.assertTrue(fill)
        # The actual fill endpoint limit is 3*0.4 - 0.06 from original outline,
        # not the innermost wall CENTERLINE at 2.5*0.4.
        limit = 3 * 0.4 - 0.06
        for p in fill:
            for theta, axial in p.points:
                self.assertGreaterEqual(min(theta * 30, (1-theta) * 30, axial, 10-axial),
                                        limit - 1e-6)

    def test_overlap_is_only_explicit_fraction_of_width(self):
        model = box(0, 0, 1, 10)
        config = SlicerConfig(nozzle_width=0.4, wall_count=3,
                              infill_wall_overlap_fraction=0.25)
        region = inset_geometry(model, 30,
                                config.wall_count*config.nozzle_width
                                - config.infill_wall_overlap_fraction*config.nozzle_width)
        self.assertAlmostEqual(_to_metric(region,30).bounds[0], 1.1)
        self.assertAlmostEqual(_to_metric(region,30).bounds[1], 1.1)

    def test_periodic_full_wrap_has_no_fake_seam_perimeter(self):
        model = box(0, 0, TWO_PI, 10)
        config = SlicerConfig(wall_count=1, nozzle_width=0.4)
        result = build_layer_toolpaths(_layer(model, periodic=True), config)
        paths = result.paths_by_role["EXTERNAL_PERIMETER"]
        self.assertEqual(len(paths), 2)
        for path in paths:
            self.assertAlmostEqual(path.points[0][0], 0, places=8)
            self.assertAlmostEqual(path.points[-1][0], TWO_PI, places=8)
            self.assertAlmostEqual(path.points[0][1], path.points[-1][1], places=8)
            self.assertTrue(abs(path.points[0][1] - 0.2) < 1e-6
                            or abs(path.points[0][1] - 9.8) < 1e-6)

    def test_bad_overlap_fraction_is_rejected(self):
        for fraction in (-0.01, 0.501, math.nan, math.inf):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                SlicerConfig(infill_wall_overlap_fraction=fraction).validate()
