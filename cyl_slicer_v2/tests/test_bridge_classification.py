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

    def test_no_previous_layer_has_no_fabricated_support(self):
        regions=self.bridge(box(0,0,1,10),GeometryCollection())
        self.assertTrue(all(x.is_empty for x in regions.values()))


if __name__ == "__main__":
    unittest.main()
