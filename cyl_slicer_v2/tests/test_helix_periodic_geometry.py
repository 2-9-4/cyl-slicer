from __future__ import annotations

import math
from pathlib import Path
import unittest

import numpy as np
import shapely
from shapely.geometry import GeometryCollection

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import TWO_PI, cylindrical_to_cartesian
from cyl_slicer_v2.gcode import emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.intersection_backend import sample_intersection_layer
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


class HelixPeriodicGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh = load_mesh(Path(__file__).resolve().parents[2] / "helix.stl")

    def test_every_turn_matches_mesh_across_radial_slabs_and_handedness(self):
        config = SlicerConfig(layer_height=0.12, nozzle_width=0.48)
        theta, axial = np.meshgrid(
            np.linspace(0.013, TWO_PI - 0.013, 48),
            np.linspace(-2, 52, 72), indexing="ij",
        )
        for handedness in (1, -1):
            mesh = self.mesh.copy()
            if handedness < 0:
                mesh.apply_transform(np.diag([1., -1., 1., 1.]))
            for radius in (32.8, 35.0, 37.2):
                with self.subTest(handedness=handedness, radius=radius):
                    geom, window, full = sample_intersection_layer(mesh, radius, config)
                    self.assertFalse(full)
                    self.assertEqual(window, (0.0, TWO_PI))
                    self.assertLess(geom.bounds[1], 0)
                    self.assertGreater(geom.bounds[3], 50)
                    expected = mesh.contains(cylindrical_to_cartesian(
                        theta.ravel(), axial.ravel(), radius,
                    ))
                    actual = shapely.intersects_xy(geom, theta.ravel(), axial.ravel())
                    # Check missing material AND incorrectly filled inter-turn
                    # gaps, allowing only boundary tessellation error.
                    self.assertLess(float(np.mean(expected != actual)), 0.003)

    def test_slanted_helix_walls_meet_at_every_seam_crossing(self):
        config = SlicerConfig(layer_height=0.12, nozzle_width=0.48, wall_count=2)
        radius = 35.0
        occupied, window, full = sample_intersection_layer(self.mesh, radius, config)
        empty = GeometryCollection()
        layer = LayerSlice(0, radius, occupied, occupied, empty, empty, empty,
                           theta_window=window, is_full_wrap=full)
        self.assertFalse(full)  # No horizontal band goes all the way around.
        result = build_layer_toolpaths(layer, config)
        self.assertEqual(len(result.island_paths), 1)
        segments = parse_extrusion_segments(emit_gcode([result], config))
        for role in ("EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"):
            with self.subTest(role=role):
                endpoints = [point for path in result.paths_by_role[role]
                             for point in (path.points[0], path.points[-1])]
                lower = sorted(x for t, x in endpoints if abs(t) < 1e-8)
                upper = sorted(x for t, x in endpoints if abs(t - TWO_PI) < 1e-8)
                # Four crossings, two edges each, including the crossings
                # beyond the first turn that normalization used to discard.
                self.assertEqual(len(lower), 8)
                np.testing.assert_allclose(lower, upper, atol=1e-7, rtol=0)
                # No artificial cap runs axially along the coordinate cut.
                for path in result.paths_by_role[role]:
                    for a, b in zip(path.points, path.points[1:]):
                        self.assertFalse(abs(a[0] - b[0]) < 1e-8 and
                                         (abs(a[0]) < 1e-8 or abs(a[0] - TWO_PI) < 1e-8))
                planned_length = sum(
                    math.hypot((b[0] - a[0]) * radius, b[1] - a[1])
                    for path in result.paths_by_role[role]
                    for a, b in zip(path.points, path.points[1:])
                )
                emitted_length = sum(
                    math.hypot(math.radians(y1 - y0) * radius, x1 - x0)
                    for _, kind, x0, y0, x1, y1 in segments if kind == role
                )
                self.assertAlmostEqual(emitted_length, planned_length, delta=0.05)


if __name__ == "__main__":
    unittest.main()
