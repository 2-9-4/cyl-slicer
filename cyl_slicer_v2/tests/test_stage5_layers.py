from __future__ import annotations

import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.layers import slice_mesh


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
