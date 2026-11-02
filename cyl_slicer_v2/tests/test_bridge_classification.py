from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, Point, box
from shapely.ops import unary_union

from cyl_slicer_v2.bridges import classify_bridge_regions


class BridgeClassificationTests(unittest.TestCase):
    def bridge(self, current, previous, r=10, max_span=8):
        return classify_bridge_regions(current,previous,r,.42,max_span)

    def test_axial_two_bearings_gap_detected(self):
        current=box(0,0,1,10)
        previous=unary_union([box(0,0,1,4),box(0,6,1,10)])
        regions=self.bridge(current,previous)
        self.assertGreater(regions["axial"].area,0)
        self.assertTrue(regions["radial"].is_empty)
        self.assertTrue(regions["axial"].covers(Point(.5,5)))

    def test_single_ended_cantilever_is_not_a_bridge(self):
        current=box(0,0,1,10)
        previous=box(0,0,1,4)
        regions=self.bridge(current,previous)
        self.assertTrue(regions["axial"].is_empty)
        self.assertTrue(regions["radial"].is_empty)

    def test_radial_gap_across_theta_two_bearings(self):
        current=box(0,0,1,10)
        previous=unary_union([box(0,0,.35,10),box(.65,0,1,10)])
        regions=self.bridge(current,previous)
        self.assertGreater(regions["radial"].area,0)
        self.assertTrue(regions["axial"].is_empty)
        self.assertTrue(regions["radial"].covers(Point(.5,5)))

    def test_unsupported_span_too_long_gets_support_not_bridge(self):
        current=box(0,0,1,10)
        previous=unary_union([box(0,0,1,1),box(0,9,1,10)])
        result=self.bridge(current,previous,max_span=4)
        self.assertTrue(result["axial"].is_empty)
        self.assertTrue(result["radial"].is_empty)

    def test_periodic_seam_bridge_not_triply_classified(self):
        current=box(6.05,0,6.52,10)
        previous=unary_union([box(6.05,0,6.19,10),box(6.39,0,6.52,10)])
        regions=self.bridge(current,previous,r=15,max_span=8)
        bridge=regions["radial"]
        self.assertFalse(bridge.is_empty)
        self.assertTrue(bridge.covers(Point(6.25,5)))
        self.assertAlmostEqual(bridge.area,
            (6.39-6.19 - (.42/15))*10,delta=.1)
        self.assertLess(bridge.area,2.5)

    def test_mixed_two_anchor_and_one_anchor_lanes_are_separated(self):
        # The left half has two bearings around a 2 mm axial void. The
        # right half has only a lower bearing, so it is a cantilever.
        current = box(0,0,1,10)
        previous = unary_union([
            box(0,0,1,4),
            box(0,6,.45,10),
        ])
        verified = self.bridge(current,previous)
        axial = verified["axial"]
        self.assertFalse(axial.is_empty)
        self.assertTrue(axial.covers(Point(.2,5)))
        self.assertFalse(axial.covers(Point(.7,5)))
        self.assertTrue(verified["radial"].is_empty)
        self.assertLess(axial.area, .5 * current.area)
        self.assertTrue(current.covers(axial))
        for k in (-2,-1,1,4):
            shifted = box(k*2*math.pi,0,1+k*2*math.pi,10)
            moved_bearings = unary_union([
                box(k*2*math.pi,0,1+k*2*math.pi,4),
                box(k*2*math.pi,6,.45+k*2*math.pi,10),
            ])
            result = self.bridge(shifted, moved_bearings)["axial"]
            self.assertAlmostEqual(result.area, axial.area,delta=1e-6)

    def test_sub_bead_width_bridge_patch_remains_supported_not_bridge(self):
        # A bearing pair with just 0.3 mm transverse width cannot fit a
        # .42 mm bridge bead and must not become a micro bridge.
        current=box(0,0,1,10)
        previous=unary_union([box(0,0,1,4),box(0,6,.03,10)])
        regions=self.bridge(current,previous)
        self.assertTrue(all(x.is_empty for x in regions.values()))

    def test_separate_narrow_bridging_strips_do_not_connect_across_cantilever(self):
        current=box(0,0,1,10)
        previous=unary_union([
            box(0,0,1,4),
            box(0,6,.32,10),
            box(.68,6,1,10),
        ])
        regions=self.bridge(current,previous)
        axial=regions["axial"]
        self.assertTrue(axial.covers(Point(.15,5)))
        self.assertTrue(axial.covers(Point(.85,5)))
        self.assertFalse(axial.covers(Point(.5,5)))
        self.assertGreater(len(getattr(axial, "geoms", [axial])),1)

    def test_full_rectangular_bridge_has_exact_bearing_fast_path(self):
        from unittest.mock import patch
        from cyl_slicer_v2.bridges import _anchored_bridge_lanes
        # 180 mm wide rectangle: avoid hundreds of repeated line probes.
        previous = unary_union([
            box(0, 3.5, 180, 4.15),
            box(0, 5.85, 180, 6.5),
        ])
        rectangle = box(0, 4.2, 180, 5.8)
        with patch("cyl_slicer_v2.bridges._lines",
                   side_effect=AssertionError("expected exact fast path")):
            result = _anchored_bridge_lanes(
                rectangle, previous, "axial", .42, 8,
            )
        self.assertTrue(result.equals(rectangle))

    def test_rectangular_gap_with_incomplete_bearings_uses_lane_checks(self):
        from cyl_slicer_v2.bridges import _anchored_bridge_lanes
        candidate = box(0, 4.2, 10, 5.8)
        previous = unary_union([
            box(0, 3.5, 10, 4.15),
            box(0, 5.85, 4, 6.5),
            box(6, 5.85, 10, 6.5),
        ])
        result = _anchored_bridge_lanes(candidate,previous,"axial",.42,8)
        self.assertFalse(result.is_empty)
        self.assertTrue(result.covers(Point(2,5)))
        self.assertTrue(result.covers(Point(8,5)))
        self.assertFalse(result.covers(Point(5,5)))
        self.assertLess(result.area,candidate.area)

    def test_no_previous_layer_has_no_fabricated_support(self):
        regions=self.bridge(box(0,0,1,10),GeometryCollection())
        self.assertTrue(all(x.is_empty for x in regions.values()))


if __name__ == "__main__":
    unittest.main()
