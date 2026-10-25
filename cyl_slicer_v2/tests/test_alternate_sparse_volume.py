from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


def synthetic_layer(index=0, radius=30.0):
    empty = GeometryCollection()
    return LayerSlice(index, radius, box(0, 0, 2, 60), empty, empty, empty, empty)


def physical_length(paths):
    return sum(
        math.hypot((b[0] - a[0]) * path.radius, b[1] - a[1])
        for path in paths for a, b in zip(path.points, path.points[1:])
    )


class AlternatingDensityTests(unittest.TestCase):
    def test_equal_volume_against_two_direction_crosshatch_per_layer(self):
        # A large, nearly square patch suppresses phase/grid-edge quantization
        # relative to the centerline metric length of the two patterns.
        base = SlicerConfig(wall_count=1, infill_spacing=2,
                            sparse_infill_pattern="both", min_path_length=0)
        alt = SlicerConfig(wall_count=1, infill_spacing=2,
                           sparse_infill_pattern="alternate", min_path_length=0)
        for index in (0, 1, 2, 3):
            both = build_layer_toolpaths(synthetic_layer(index), base).paths_by_role
            only = build_layer_toolpaths(synthetic_layer(index), alt).paths_by_role
            reference = physical_length(both.get("INFILL_RADIAL", [])) + physical_length(
                both.get("INFILL_AXIAL", []))
            observed = physical_length(only.get("INFILL_RADIAL", [])) + physical_length(
                only.get("INFILL_AXIAL", []))
            self.assertGreater(reference, 500)
            self.assertAlmostEqual(observed / reference, 1, delta=.08)
            self.assertEqual(bool(only.get("INFILL_RADIAL")), index % 2 == 0)
            self.assertEqual(bool(only.get("INFILL_AXIAL")), index % 2 == 1)
            self.assertLess(
                physical_length(only.get("INFILL_RADIAL", [])) *
                physical_length(only.get("INFILL_AXIAL", [])), 1e-12
            )

    def test_crosshatch_default_unchanged_and_thin_bead_flow_matched(self):
        c = SlicerConfig(
            wall_count=1, bead_model="rounded", layer_height=.2,
            nozzle_width=.4, infill_spacing=1.8, min_path_length=0,
            sparse_infill_pattern="alternate")
        c.validate()
        ref = SlicerConfig(
            wall_count=1, bead_model="rounded", layer_height=.2,
            nozzle_width=.4, infill_spacing=1.8, min_path_length=0,
            sparse_infill_pattern="both")
        alt = build_layer_toolpaths(synthetic_layer(0), c).paths_by_role
        both = build_layer_toolpaths(synthetic_layer(0), ref).paths_by_role
        a = physical_length(alt.get("INFILL_RADIAL", [])) * c.extrusion_per_mm
        b = sum(physical_length(both.get(role, [])) for role in
                ("INFILL_RADIAL", "INFILL_AXIAL")) * ref.extrusion_per_mm
        self.assertAlmostEqual(a / b, 1, delta=.10)
        self.assertAlmostEqual(c.dense_line_spacing, .357079632679, places=8)

    def test_reject_overpacked_alternation_but_not_legacy_crosshatch(self):
        c = SlicerConfig(infill_spacing=.6, sparse_infill_pattern="alternate")
        with self.assertRaisesRegex(ValueError, "requires infill_spacing"):
            c.validate()
        SlicerConfig(infill_spacing=.6, sparse_infill_pattern="both").validate()

    def test_short_features_and_support_roles_unchanged(self):
        base = SlicerConfig(wall_count=1, infill_spacing=3, min_path_length=0)
        alt = SlicerConfig(wall_count=1, infill_spacing=3, min_path_length=0,
                           sparse_infill_pattern="alternate")
        layer = synthetic_layer()
        a = build_layer_toolpaths(layer, base).paths_by_role
        b = build_layer_toolpaths(layer, alt).paths_by_role
        for role in ("EXTERNAL_PERIMETER", "INTERNAL_PERIMETER",
                     "SUPPORT", "SUPPORT_INTERFACE", "TOP_SKIN"):
            self.assertEqual(len(a.get(role, [])), len(b.get(role, [])))


if __name__ == "__main__":
    unittest.main()
