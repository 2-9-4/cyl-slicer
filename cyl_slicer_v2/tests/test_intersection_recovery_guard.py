from __future__ import annotations

import math
import unittest
import warnings
from unittest.mock import patch

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.intersection_backend import _recovery_compatible, sample_intersection_layer


class IntersectionRecoveryGuardTests(unittest.TestCase):
    def setUp(self):
        self.config = SlicerConfig(layer_height=0.2, intersection_curve_error_mm=0.025)
        self.center = (box(0, 0, 2 * math.pi, 2), (0, 2 * math.pi), True)
        self.positive = (box(0, 0, 1, 1), (0, 1), False)
        self.negative = (box(0, 0.001, 1, 1.001), (0, 1), False)

    def _test_recovery(self, opposite=None, offsets=(0.001, -0.001), strict=False):
        opposite = self.negative if opposite is None else opposite
        self.config = SlicerConfig(layer_height=0.2, intersection_curve_error_mm=0.025,
                                   intersection_strict_recovery=strict)

        def sample(_mesh, radius, _config, _z_bounds):
            if abs(radius - 30) < 1e-12:
                return self.center
            return self.positive if radius > 30 else opposite

        with patch("cyl_slicer_v2.intersection_backend._sample_intersection_layer_once",
                   side_effect=sample), patch(
            "cyl_slicer_v2.intersection_backend._looks_like_degenerate_full_span",
            side_effect=lambda _m, r, _g, _f, _c: abs(r - 30) < 1e-12
        ), patch(
            "cyl_slicer_v2.intersection_backend._degenerate_radius_offsets",
            return_value=offsets
        ):
            return sample_intersection_layer(object(), 30, self.config)

    def test_consistent_two_sided_recovery_is_reported(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = self._test_recovery()
        self.assertTrue(result[0].equals(self.positive[0]))
        self.assertEqual(len(caught), 1)
        self.assertIn("radial shift=+0.001", str(caught[0].message))
        self.assertIn("boundary deviation=", str(caught[0].message))
        self.assertIn("two-sided verified", str(caught[0].message))

    def test_topology_change_between_sides_is_rejected(self):
        hole = self.positive[0].difference(box(.2, .2, .5, .5))
        with self.assertRaisesRegex(RuntimeError, "degenerate"):
            self._test_recovery(opposite=(hole, (0, 1), False), strict=True)

    def test_large_boundary_shift_rejected_even_with_same_topology(self):
        moved = (box(0, .5, 1, 1.5), (0, 1), False)
        with self.assertRaises(RuntimeError):
            self._test_recovery(opposite=moved, strict=True)

    def test_non_strict_recovery_warns_on_topological_transition(self):
        hole = self.positive[0].difference(box(.2, .2, .5, .5))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = self._test_recovery(opposite=(hole, (0, 1), False))
        self.assertTrue(result[0].equals(self.positive[0]))
        self.assertTrue(any("ONE-SIDED" in str(w.message) for w in caught))

    def test_recovery_cannot_exceed_layer_or_curve_error_budget(self):
        with self.assertRaises(RuntimeError):
            self._test_recovery(offsets=(0.04, -0.04))

    def test_hole_count_and_full_wrap_are_topological_invariants(self):
        self.assertFalse(_recovery_compatible(
            self.positive, (self.negative[0], (0, 1), True), 30, 0.025
        )[0])
        self.assertTrue(_recovery_compatible(
            self.positive, self.negative, 30, 0.025
        )[0])
