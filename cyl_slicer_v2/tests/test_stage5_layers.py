from __future__ import annotations

import math
from pathlib import Path
import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.intersection_backend import _polygon_parts
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.mesh_io import load_mesh


def _corner_count(geom) -> int:
    return sum(max(0, len(poly.exterior.coords) - 1) for poly in _polygon_parts(geom)) + sum(
        max(0, len(ring.coords) - 1) for poly in _polygon_parts(geom) for ring in poly.interiors
    )


class LayerClassificationTests(unittest.TestCase):
    def test_first_occupied_layer_is_bottom_surface(self):
        layers = slice_mesh(make_fixture("solid_box"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0), max_layers=3)
        first = next(layer for layer in layers if not layer.occupied.is_empty)
        self.assertAlmostEqual(first.bottom_surfaces.area, first.occupied.area)

    def test_tube_middle_layers_have_no_false_surface_churn(self):
        layers = slice_mesh(make_fixture("tube"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=0.5))
        occupied = [layer for layer in layers if not layer.occupied.is_empty]
        middle = occupied[len(occupied) // 2]
        self.assertLess(middle.top_surfaces.area, middle.occupied.area * 0.05)
        self.assertLess(middle.bottom_surfaces.area, middle.occupied.area * 0.05)

    def test_overhang_step_creates_support(self):
        layers = slice_mesh(make_fixture("overhang_step"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0))
        self.assertTrue(any(not layer.needs_support.is_empty for layer in layers))
        self.assertTrue(any(not layer.support.is_empty for layer in layers))

    def test_support_regions_have_valid_bounds(self):
        layers = slice_mesh(make_fixture("overhang_step"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0))
        for layer in layers:
            if layer.support.is_empty:
                continue
            t_min, z_min, t_max, z_max = layer.support.bounds
            self.assertGreaterEqual(t_min, 0.0)
            self.assertGreater(t_max, t_min)
            self.assertGreater(z_max, z_min)

    def test_cubeandring2_layers_do_not_seed_full_frame_geometry(self):
        mesh_path = Path(__file__).resolve().parent / "fixtures" / "cubeandring2.stl"
        mesh = load_mesh(mesh_path)
        layers = slice_mesh(
            mesh,
            SlicerConfig(layer_backend="intersection", layer_height=0.25, slice_workers=4),
        )
        frame_area = 2.0 * math.pi * (mesh.bounds[1, 2] - mesh.bounds[0, 2])

        layer25 = layers[25]
        self.assertGreater(layer25.occupied.area, 20.0)
        self.assertTrue(layer25.support.is_empty)

        layer46 = layers[46]
        self.assertFalse(layer46.is_full_wrap)
        self.assertLess(layer46.occupied.area, frame_area * 0.75)

    def test_cubeandring2_layers_79_and_80_do_not_spike_corner_count(self):
        mesh_path = Path(__file__).resolve().parent / "fixtures" / "cubeandring2.stl"
        layers = slice_mesh(
            load_mesh(mesh_path),
            SlicerConfig(layer_backend="intersection", layer_height=0.25, slice_workers=4),
        )
        counts = {index: _corner_count(layers[index].occupied) for index in range(78, 82)}
        neighbor_limit = max(counts[78], counts[81]) * 2

        self.assertLessEqual(counts[79], neighbor_limit)
        self.assertLessEqual(counts[80], neighbor_limit)

    def test_parallel_slicing_matches_serial_slicing(self):
        mesh = make_fixture("tube")
        serial = slice_mesh(mesh, SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0, slice_workers=1), max_layers=3)
        parallel = slice_mesh(mesh, SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0, slice_workers=2), max_layers=3)
        self.assertEqual(len(serial), len(parallel))
        for a, b in zip(serial, parallel):
            self.assertAlmostEqual(a.occupied.area, b.occupied.area)
            self.assertAlmostEqual(a.top_surfaces.area, b.top_surfaces.area)
            self.assertAlmostEqual(a.bottom_surfaces.area, b.bottom_surfaces.area)
            self.assertEqual(a.is_full_wrap, b.is_full_wrap)


if __name__ == "__main__":
    unittest.main()
