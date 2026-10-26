from __future__ import annotations

import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.layers import propagate_support
from cyl_slicer_v2.models import LayerSlice


EMPTY = GeometryCollection()


def varied_layers():
    """Continuous slabs: 30.00..30.16..30.26..30.36..30.46..30.56."""
    heights = [.16, .10, .10, .10, .10]
    edge = 30.
    layers = []
    for i, thickness in enumerate(heights):
        layers.append(LayerSlice(
            i, edge + thickness/2,
            EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            deposited_thickness=thickness,
        ))
        edge += thickness
    return layers


def config(gap):
    return SlicerConfig(support_contact_gap_mm=gap,
                        support_min_island_area=0,
                        support_min_region_area=0,
                        support_bearing_min_area=0)


class VariableSlabSupportContactTests(unittest.TestCase):
    def test_contact_limit_uses_actual_faces_not_nominal_layer_count(self):
        # Model layer 4 spans [30.46, 30.56].
        geom = box(.2, 0, .4, 8)
        layers = varied_layers()
        layers[4].needs_support = geom
        layers[4].support_interface_target = geom
        propagate_support(layers, config(.15))
        # Candidate layer 2 ends at 30.36: only 0.10 below model.
        # Candidate layer 1 ends at 30.26: 0.20, enough clearance.
        self.assertTrue(layers[3].support.is_empty)
        self.assertTrue(layers[2].support.is_empty)
        self.assertFalse(layers[1].support.is_empty)
        self.assertFalse(layers[0].support.is_empty)
        self.assertGreaterEqual(30.46 - (layers[1].radius + .05), .15)

    def test_exact_boundary_included_and_tiny_gap_excludes_direct_contact(self):
        geom = box(.2, 0, .4, 8)
        for gap, expected in ((0, 3), (.10, 2), (.100001, 1)):
            with self.subTest(gap=gap):
                layers = varied_layers()
                layers[4].needs_support = geom
                layers[4].support_interface_target = geom
                propagate_support(layers, config(gap))
                first = next(i for i in range(3, -1, -1)
                             if not layers[i].support.is_empty)
                self.assertEqual(first, expected)

    def test_separate_overhangs_preserve_independent_contact_planes(self):
        layers = varied_layers()
        early = box(.2, 0, .4, 8)
        late = box(.6, 0, .8, 8)
        layers[3].needs_support = early
        layers[3].support_interface_target = early
        layers[4].needs_support = late
        layers[4].support_interface_target = late
        propagate_support(layers, config(.15))
        # Origin layer 3 underside 30.36 -> support may reach only top
        # 30.21 or less; layer 1 ends 30.26, so layer 0 carries it.
        self.assertTrue(layers[1].support.disjoint(early))
        self.assertFalse(layers[0].support.intersection(early).is_empty)
        # Origin layer 4 underside 30.46 -> layer 1 ends 30.26.
        self.assertFalse(layers[1].support.intersection(late).is_empty)

    def test_explicit_mm_gap_and_implicit_legacy_layers_are_distinct(self):
        geom = box(.2, 0, .4, 8)
        layers = varied_layers()
        layers[4].needs_support = geom
        layers[4].support_interface_target = geom
        legacy = SlicerConfig(support_gap_layers=1,
                              support_min_island_area=0,
                              support_min_region_area=0,
                              support_bearing_min_area=0)
        propagate_support(layers, legacy)
        self.assertTrue(layers[3].support.is_empty)
        self.assertFalse(layers[2].support.is_empty)

    def test_interface_begins_only_on_first_eligible_support_slab(self):
        layers = varied_layers()
        geom = box(.2, 0, .4, 8)
        layers[4].needs_support = geom
        layers[4].support_interface_target = geom
        propagate_support(layers, config(.15))
        self.assertTrue(layers[2].support_interface.is_empty)
        self.assertFalse(layers[1].support_interface.is_empty)

    def test_deposited_thickness_must_be_positive(self):
        layers = varied_layers()
        layers[2].deposited_thickness = -1
        with self.assertRaisesRegex(ValueError, "thickness"):
            propagate_support(layers, config(0))


if __name__ == "__main__":
    unittest.main()
