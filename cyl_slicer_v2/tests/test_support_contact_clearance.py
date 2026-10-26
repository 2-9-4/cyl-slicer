from __future__ import annotations

import math
import unittest
from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.layers import propagate_support
from cyl_slicer_v2.models import LayerSlice


EMPTY=GeometryCollection()


def support_stack(count=4, *, occupied_middle=None):
    geom=box(.2,0,.4,8)
    layers=[
        LayerSlice(i,30+.16*i,
                   occupied_middle if i==count-2 and occupied_middle is not None else EMPTY,
                   EMPTY,EMPTY,EMPTY,EMPTY)
        for i in range(count)
    ]
    layers[-1].needs_support=geom
    layers[-1].support_interface_target=geom
    return layers


class SupportContactClearanceTests(unittest.TestCase):
    def test_default_is_exact_previous_layer_count(self):
        c=SlicerConfig()
        self.assertIsNone(c.support_contact_gap_mm)
        self.assertEqual(c.effective_support_gap_layers,c.support_gap_layers)
        self.assertEqual(c.support_lateral_gap_mm,0)

    def test_explicit_physical_gap_ceil_to_next_whole_layer(self):
        self.assertEqual(SlicerConfig(support_contact_gap_mm=0).effective_support_gap_layers,0)
        self.assertEqual(SlicerConfig(support_contact_gap_mm=.16).effective_support_gap_layers,1)
        self.assertEqual(SlicerConfig(support_contact_gap_mm=.17).effective_support_gap_layers,2)
        self.assertEqual(SlicerConfig(support_contact_gap_mm=.32).effective_support_gap_layers,2)
        self.assertEqual(SlicerConfig(support_contact_gap_mm=.33).effective_support_gap_layers,3)

    def test_support_skips_requested_number_of_physical_radial_layers(self):
        for gap_mm,empty_count in ((0,0),(.16,1),(.30,2)):
            c=SlicerConfig(support_contact_gap_mm=gap_mm,
                support_min_island_area=0,support_min_region_area=0,
                support_bearing_min_area=0)
            layers=support_stack()
            propagate_support(layers,c)
            first=next((i for i in range(2,-1,-1)
                        if not layers[i].support.is_empty),None)
            self.assertIsNotNone(first)
            self.assertEqual(2-first,empty_count)
            for i in range(2,first,-1):
                self.assertTrue(layers[i].support.is_empty)

    def test_lateral_clearance_is_physical_not_theta_radians(self):
        model=box(.2,0,.25,8)
        a=support_stack(3,occupied_middle=model)
        b=support_stack(3,occupied_middle=model)
        settings=dict(support_contact_gap_mm=0,support_min_island_area=0,
                      support_min_region_area=0,support_bearing_min_area=0)
        propagate_support(a,SlicerConfig(**settings))
        propagate_support(b,SlicerConfig(**settings,support_lateral_gap_mm=.3))
        self.assertGreater(a[1].support.area,b[1].support.area)
        self.assertTrue(b[1].support.disjoint(model))
        self.assertGreaterEqual(b[1].support.bounds[0],
                                model.bounds[2]+.3/b[1].radius-2e-5)

    def test_does_not_remove_lower_support_foundation(self):
        layers=support_stack(3,occupied_middle=box(.2,0,.25,8))
        propagate_support(layers,SlicerConfig(
            support_contact_gap_mm=0,support_lateral_gap_mm=.3,
            support_min_island_area=0,support_min_region_area=0,
            support_bearing_min_area=0))
        self.assertFalse(layers[0].support.is_empty)

    def test_cli_resolves_defaults_and_explicit_zero_for_both_profiles(self):
        for profile in ("quality","legacy-cli","library"):
            args=build_parser().parse_args(["--profile",profile])
            self.assertIsNone(_profile_value(args,"support_contact_gap_mm"))
            self.assertEqual(_profile_value(args,"support_lateral_gap_mm"),0)
            args=build_parser().parse_args(["--profile",profile,
                "--support-contact-gap-mm","0.32","--support-lateral-gap-mm",".24"])
            self.assertEqual(_profile_value(args,"support_contact_gap_mm"),.32)
            self.assertEqual(_profile_value(args,"support_lateral_gap_mm"),.24)

    def test_invalid_nonfinite_and_negative_gaps(self):
        for field in ("support_contact_gap_mm","support_lateral_gap_mm"):
            for bad in (-1,math.nan,math.inf):
                with self.subTest(field=field,bad=bad),self.assertRaises(ValueError):
                    SlicerConfig(**{field:bad}).validate()


if __name__=="__main__":
    unittest.main()
