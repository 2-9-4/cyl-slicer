from __future__ import annotations

import math
import unittest

from shapely.geometry import box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.occupancy import sample_layer
from cyl_slicer_v2.periodic import normalize_periodic


class OccupancyTests(unittest.TestCase):
    def test_tube_produces_full_wrap_layer(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0)
        geom, window, full = sample_layer(make_fixture("tube"), 11.0, config)
        self.assertFalse(geom.is_empty)
        self.assertTrue(full)
        self.assertAlmostEqual(window[1] - window[0], 2 * math.pi)

    def test_offset_boxes_produce_multiple_islands(self):
        config = SlicerConfig(theta_step_deg=5, axial_step=1.0)
        geom, _, _ = sample_layer(make_fixture("offset_box_pair"), 20.0, config)
        parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        self.assertGreaterEqual(len([p for p in parts if not p.is_empty]), 2)

    def test_empty_radius_produces_empty_slice(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0)
        geom, _, full = sample_layer(make_fixture("solid_box"), 2.0, config)
        self.assertTrue(geom.is_empty)
        self.assertFalse(full)

    def test_seam_crossing_cells_remain_connected_after_normalization(self):
        raw = box(-0.12, 0.0, 0.18, 2.0)
        geom, _, full = normalize_periodic(raw, 0.0, 2.0)
        parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        self.assertFalse(full)
        self.assertEqual(len([p for p in parts if not p.is_empty]), 1)
        self.assertGreater(geom.bounds[0], 0.0)


if __name__ == "__main__":
    unittest.main()
