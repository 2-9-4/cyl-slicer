from __future__ import annotations

import math
import re
import unittest

import trimesh
from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.layers import layer_radii, layer_slabs, slice_mesh
from cyl_slicer_v2.models import LayerSlice, LayerToolpaths, Toolpath
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


class TerminalRadialSlabTests(unittest.TestCase):
    @staticmethod
    def mesh():
        mesh = trimesh.creation.box(extents=(20, 20, 10))
        mesh.apply_translation((30, 0, 0))
        return mesh

    def test_fixed_default_preserves_legacy_sampling(self):
        cfg = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.35)
        self.assertEqual(cfg.radial_layer_mode, "fixed")
        self.assertEqual(len(layer_radii(self.mesh(), cfg)), 2)
        self.assertAlmostEqual(layer_radii(self.mesh(), cfg)[-1], 30.24)

    def test_terminal_remainder_rebalanced_no_thin_sliver(self):
        cfg = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.35,
                           radial_layer_mode="complete-terminal")
        slabs = layer_slabs(self.mesh(), cfg)
        self.assertEqual(len(slabs), 3)
        self.assertAlmostEqual(slabs[0][1], .16)
        self.assertAlmostEqual(slabs[1][1], .095)
        self.assertAlmostEqual(slabs[2][1], .095)
        self.assertAlmostEqual(sum(h for _, h in slabs), .35, places=9)
        edge = 30.0
        for radius, height in slabs:
            self.assertAlmostEqual(radius - .5 * height, edge, places=9)
            edge = radius + .5 * height
            self.assertGreaterEqual(height, .08 - 1e-9)
        self.assertAlmostEqual(edge, 30.35, places=9)

    def test_no_remainder_and_single_printable_terminal_slab(self):
        c = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.32,
                         radial_layer_mode="complete-terminal")
        self.assertEqual(len(layer_slabs(self.mesh(), c)), 2)
        c = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.12,
                         radial_layer_mode="complete-terminal")
        self.assertEqual(len(layer_slabs(self.mesh(), c)), 1)
        self.assertAlmostEqual(layer_slabs(self.mesh(), c)[0][1], .12)

    def test_unprintably_short_range_rejected_in_opt_in_mode(self):
        c = SlicerConfig(bed_radius=30, min_radius=30, max_radius=30.04,
                         radial_layer_mode="complete-terminal")
        with self.assertRaisesRegex(ValueError, "half a layer"):
            layer_slabs(self.mesh(), c)

    def test_coupled_bead_area_pitch_and_calibration_multiplier(self):
        cfg = SlicerConfig(extrusion_per_mm=.123, radial_layer_mode="complete-terminal")
        smaller = cfg.for_deposited_thickness(.095)
        self.assertAlmostEqual(smaller.extrusion_per_mm / cfg.extrusion_per_mm,
                               smaller.bead_area / cfg.bead_area, places=12)
        self.assertAlmostEqual(smaller.dense_line_spacing,
                               smaller.bead_area / .095)
        self.assertEqual(cfg.layer_height, .16)
        self.assertIs(cfg.for_deposited_thickness(.16), cfg)

    def test_builder_places_nozzle_on_outer_slab_and_uses_thickness(self):
        geom = box(.1, 0, .7, 10)
        empty = GeometryCollection()
        base = SlicerConfig(bed_radius=30, wall_count=1,
                            radial_layer_mode="complete-terminal")
        layer = LayerSlice(0, 30.3025, geom, geom, empty, empty, empty,
                           deposited_thickness=.095)
        toolpaths = build_layer_toolpaths(layer, base)
        self.assertAlmostEqual(toolpaths.nozzle_radius, 30.35)
        self.assertAlmostEqual(toolpaths.deposited_thickness, .095)
        self.assertTrue(toolpaths.paths_by_role)
        self.assertTrue(all(abs(path.radius - 30.3025) < 1e-9
                            for paths in toolpaths.paths_by_role.values() for path in paths))

    def test_gcode_partial_layer_uses_volume_not_nominal_height(self):
        cfg = SlicerConfig(bed_radius=30, retract_mm=0, cooling_min_layer_time_s=0,
                           cooling_min_island_print_time_s=0,
                           radial_layer_mode="complete-terminal")
        path = Toolpath("SOLID_INFILL", 30.3025, [(0, 0), (0, 10)])
        layer = LayerToolpaths(0, 30.3025, {"SOLID_INFILL": [path]},
                               nozzle_radius=30.35, deposited_thickness=.095)
        code = emit_gcode([layer], cfg)
        self.assertIn(";SAMPLE_RADIUS:30.3025", code)
        self.assertIn(";DEPOSITED_THICKNESS:0.0950", code)
        self.assertIn(";Z:30.3500", code)
        actual = float(re.search(r"(?m)^G1 X[^\\n]*? E(-?[0-9.]+)", code).group(1))
        expected = 10 * cfg.for_deposited_thickness(.095).extrusion_per_mm
        self.assertAlmostEqual(actual, expected, places=4)

    def test_invalid_mode(self):
        with self.assertRaisesRegex(ValueError, "radial_layer_mode"):
            SlicerConfig(radial_layer_mode="broken").validate()


if __name__ == "__main__":
    unittest.main()
