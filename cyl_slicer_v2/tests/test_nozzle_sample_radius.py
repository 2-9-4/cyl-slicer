from __future__ import annotations

import re
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerSlice, LayerToolpaths
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


def layer_at(index, sample_radius, thickness=0.2):
    return LayerToolpaths(index, sample_radius, nozzle_radius=sample_radius + thickness / 2,
                          deposited_thickness=thickness)


class NozzleSampleRadiusTests(unittest.TestCase):
    def test_first_layer_uses_full_thickness_not_midpoint_clearance(self):
        config = SlicerConfig(bed_radius=30, layer_height=0.2)
        code = emit_gcode([layer_at(0, 30.1)], config)
        self.assertIn(";SAMPLE_RADIUS:30.1000", code)
        self.assertIn(";DEPOSITED_THICKNESS:0.2000", code)
        self.assertIn(";Z:30.2000", code)
        self.assertIn("G0 Z0.2000", code)
        self.assertNotIn("G0 Z0.1000", code)

    def test_terminal_layer_envelope_is_sample_plus_half_thickness(self):
        config = SlicerConfig(bed_radius=30, layer_height=0.2)
        code = emit_gcode([layer_at(0, 30.1), layer_at(1, 30.3)], config)
        self.assertIn(";Z:30.4000", code)
        self.assertIn("G0 Z10.4000 F600.0 ; end lift", code)

    def test_metadata_propagates_from_toolpath_builder_without_resampling(self):
        geom = box(0, 0, 0.5, 10)
        layer = LayerSlice(0, 30.1, geom, GeometryCollection(),
                           GeometryCollection(), GeometryCollection(), GeometryCollection())
        config = SlicerConfig(bed_radius=30, layer_height=0.2, wall_count=1)
        toolpaths = build_layer_toolpaths(layer, config)
        self.assertEqual(toolpaths.radius, 30.1)
        self.assertAlmostEqual(toolpaths.nozzle_radius, 30.2)
        self.assertAlmostEqual(toolpaths.deposited_thickness, 0.2)
        self.assertTrue(all(path.radius == 30.1 for paths in toolpaths.paths_by_role.values()
                            for path in paths))

    def test_nozzle_radius_with_explicit_thickness_does_not_depend_on_config_default(self):
        config = SlicerConfig(bed_radius=30, layer_height=0.5)
        code = emit_gcode([layer_at(0, 30.1, thickness=0.2)], config)
        self.assertIn("G0 Z0.2000", code)
        self.assertIn(";DEPOSITED_THICKNESS:0.2000", code)

    def test_manual_toolpath_metadata_defaults_to_configured_half_layer(self):
        code = emit_gcode([LayerToolpaths(0, 30.1)],
                          SlicerConfig(bed_radius=30, layer_height=0.2))
        self.assertIn("G0 Z0.2000", code)

    def test_reject_invalid_nozzle_or_thickness(self):
        config = SlicerConfig(bed_radius=30, layer_height=0.2)
        with self.assertRaises(ValueError):
            emit_gcode([LayerToolpaths(0, 30.1, nozzle_radius=29.99)], config)
        with self.assertRaises(ValueError):
            emit_gcode([LayerToolpaths(0, 30.1, deposited_thickness=-0.1)], config)
