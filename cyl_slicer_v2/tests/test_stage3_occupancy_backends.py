from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.occupancy_backends import GridOccupancyBackend, MeshContainsBackend, build_occupancy_backend


class OccupancyBackendTests(unittest.TestCase):
    def test_grid_backend_answers_basic_inside_outside_points(self):
        backend = GridOccupancyBackend.build(make_fixture("solid_box"), voxel_mm=1.0)
        points = np.array(
            [
                [20.0, 0.0, 5.0],
                [0.0, 0.0, 0.0],
            ]
        )
        self.assertEqual(backend.contains(points).tolist(), [True, False])

    def test_grid_backend_cache_round_trip(self):
        mesh = make_fixture("solid_box")
        points = np.array([[20.0, 0.0, 5.0], [0.0, 0.0, 0.0]])
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp) / "grid.npz"
            built = build_occupancy_backend(mesh, kind="grid", grid_voxel_mm=1.0, grid_cache=cache)
            loaded = build_occupancy_backend(mesh, kind="grid", grid_voxel_mm=1.0, grid_cache=cache)
            self.assertTrue(cache.exists())
            self.assertEqual(built.contains(points).tolist(), loaded.contains(points).tolist())

    def test_grid_backend_cache_rejects_mismatched_voxel_size(self):
        mesh = make_fixture("solid_box")
        with tempfile.TemporaryDirectory() as tmp:
            cache = Path(tmp) / "grid.npz"
            build_occupancy_backend(mesh, kind="grid", grid_voxel_mm=1.0, grid_cache=cache)
            with self.assertRaises(ValueError):
                build_occupancy_backend(mesh, kind="grid", grid_voxel_mm=0.5, grid_cache=cache)

    def test_grid_backend_slicing_smoke_matches_exact_nonempty_layers(self):
        mesh = make_fixture("tube")
        config = SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0, slice_workers=1)
        exact = slice_mesh(mesh, config, max_layers=3, occupancy_backend=MeshContainsBackend(mesh))
        grid = slice_mesh(mesh, config, max_layers=3, occupancy_backend=GridOccupancyBackend.build(mesh, voxel_mm=0.75))
        self.assertEqual([not layer.occupied.is_empty for layer in exact], [not layer.occupied.is_empty for layer in grid])


if __name__ == "__main__":
    unittest.main()
