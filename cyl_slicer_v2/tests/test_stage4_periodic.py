from __future__ import annotations

import math
import unittest

from shapely.geometry import box

from cyl_slicer_v2.periodic import fold_to_canonical, is_full_wrap, normalize_periodic


class PeriodicTests(unittest.TestCase):
    def test_island_crossing_zero_is_one_component(self):
        geom, _, _ = normalize_periodic(box(-0.2, 0.0, 0.2, 1.0), 0.0, 1.0)
        parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        self.assertEqual(len([p for p in parts if not p.is_empty]), 1)

    def test_opposite_seam_islands_do_not_merge_when_not_connected(self):
        raw = box(0.1, 0.0, 0.3, 1.0).union(box(2 * math.pi - 0.7, 0.0, 2 * math.pi - 0.5, 1.0))
        geom, _, _ = normalize_periodic(raw, 0.0, 1.0)
        parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        self.assertEqual(len([p for p in parts if not p.is_empty]), 2)

    def test_full_wrap_strip(self):
        raw = box(0.0, 0.0, 2 * math.pi, 1.0)
        geom, window, full = normalize_periodic(raw, 0.0, 1.0)
        self.assertTrue(full)
        self.assertTrue(is_full_wrap(geom))
        self.assertAlmostEqual(window[1] - window[0], 2 * math.pi)

    def test_hole_survives_fold_to_canonical(self):
        outer = box(0.0, 0.0, 2 * math.pi, 4.0)
        hole = box(1.0, 1.0, 2.0, 2.0)
        folded = fold_to_canonical(outer.difference(hole), 0.0, 4.0)
        poly = folded.geoms[0] if hasattr(folded, "geoms") else folded
        self.assertGreater(len(poly.interiors), 0)


if __name__ == "__main__":
    unittest.main()
