from __future__ import annotations

import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, gap_support_slices


EMPTY = GeometryCollection()
SHAPE = box(.1, 0, .9, 8)


def model(inner, thickness):
    return LayerSlice(0, inner + thickness/2,
                      SHAPE, EMPTY, EMPTY, EMPTY, EMPTY,
                      deposited_thickness=thickness)


class GapFoundationVariableSlabTests(unittest.TestCase):
    def cfg(self):
        return SlicerConfig(
            bed_radius=30, layer_height=.16, radial_layer_mode="complete-terminal",
            raft_layers=2, generate_gap_support=True
        )

    def test_rebalance_small_remainder_and_touch_model_inner_face(self):
        c = self.cfg()
        # From mandrel R30.0 to underside of first model R30.35.
        first = model(30.35, .095)
        support = gap_support_slices([first], c)
        self.assertEqual(len(support), 3)
        self.assertEqual([round(layer.deposited_thickness, 5) for layer in support],
                         [.16, .095, .095])
        self.assertEqual([round(layer.radius, 5) for layer in support],
                         [30.08, 30.2075, 30.3025])
        edge = c.bed_radius
        for layer in support:
            self.assertAlmostEqual(layer.radius - layer.deposited_thickness/2, edge)
            edge = layer.radius + layer.deposited_thickness/2
            self.assertGreaterEqual(layer.deposited_thickness, .08 - 1e-9)
        self.assertAlmostEqual(edge, first.radius - first.deposited_thickness/2)
        self.assertTrue(support[0].is_raft)
        self.assertTrue(support[1].is_raft)
        self.assertFalse(support[2].is_raft)

    def test_generated_support_toolpath_and_gcode_match_partial_thickness(self):
        c = self.cfg()
        first = model(30.35, .095)
        top = gap_support_slices([first], c)[-1]
        paths = build_layer_toolpaths(top, c)
        self.assertAlmostEqual(paths.nozzle_radius, 30.35)
        self.assertAlmostEqual(paths.deposited_thickness, .095)
        self.assertTrue(paths.paths_by_role.get("SUPPORT"))
        result = emit_gcode([paths], c)
        self.assertIn(";DEPOSITED_THICKNESS:0.0950", result)
        self.assertIn(";Z:30.3500", result)

    def test_exact_full_layers_no_slivers(self):
        c = self.cfg()
        support = gap_support_slices([model(30.32, .16)], c)
        self.assertEqual(len(support), 2)
        self.assertAlmostEqual(support[-1].radius +
                               support[-1].deposited_thickness/2, 30.32)

    def test_gap_shorter_than_half_height_is_rejected(self):
        c = self.cfg()
        with self.assertRaisesRegex(ValueError, "thinner than half"):
            gap_support_slices([model(30.04, .16)], c)

    def test_contact_at_mandrel_requires_no_foundation(self):
        c = self.cfg()
        self.assertEqual(gap_support_slices([model(30., .16)], c), [])

    def test_opt_out_does_not_add_support(self):
        c = SlicerConfig(radial_layer_mode="complete-terminal",
                         generate_gap_support=False)
        self.assertEqual(gap_support_slices([model(30.35, .095)], c), [])

    def test_fixed_mode_keeps_old_center_schedule(self):
        c = SlicerConfig(bed_radius=30, layer_height=.16,
                         radial_layer_mode="fixed")
        # Fixed schedule is intentionally unchanged for existing G-code.
        slabs = gap_support_slices([model(30.35, .095)], c)
        self.assertEqual([round(s.radius, 5) for s in slabs],
                         [30.08, 30.24])
        self.assertTrue(all(s.deposited_thickness is None for s in slabs))


if __name__ == "__main__":
    unittest.main()
