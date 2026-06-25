from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, Point, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import cylindrical_distance
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, build_toolpaths, gap_support_slices


class ToolpathTests(unittest.TestCase):
    def test_toolpath_points_lie_on_or_inside_source(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, wall_count=2)
        layer = next(layer for layer in slice_mesh(make_fixture("solid_box"), config, max_layers=3) if not layer.occupied.is_empty)
        layer_paths = build_layer_toolpaths(layer, config)
        for role, paths in layer_paths.paths_by_role.items():
            if role == "SUPPORT":
                source = layer.support
            else:
                source = layer.occupied.buffer(1e-8)
            for path in paths:
                for theta, axial in path.points:
                    self.assertTrue(source.covers(Point(theta, axial)), role)

    def test_perimeter_count_matches_wall_count_for_wide_box(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, wall_count=2)
        layer = next(layer for layer in slice_mesh(make_fixture("solid_box"), config, max_layers=3) if not layer.occupied.is_empty)
        paths = build_layer_toolpaths(layer, config).paths_by_role
        self.assertGreaterEqual(len(paths.get("EXTERNAL_PERIMETER", [])), 1)
        self.assertGreaterEqual(len(paths.get("INTERNAL_PERIMETER", [])), 1)

    def test_infill_spacing_is_reasonable(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, infill_spacing=2.0)
        layer = next(layer for layer in slice_mesh(make_fixture("tube"), config, max_layers=3) if not layer.occupied.is_empty)
        paths = build_layer_toolpaths(layer, config).paths_by_role
        radial_z = sorted({round(path.points[0][1], 3) for path in paths.get("INFILL_RADIAL", [])})
        if len(radial_z) >= 2:
            gaps = [b - a for a, b in zip(radial_z, radial_z[1:])]
            self.assertTrue(any(abs(gap - config.infill_spacing) <= 0.01 for gap in gaps))

    def test_no_degenerate_paths(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0)
        layer = next(layer for layer in slice_mesh(make_fixture("tube"), config, max_layers=3) if not layer.occupied.is_empty)
        paths = build_layer_toolpaths(layer, config).paths_by_role
        for role_paths in paths.values():
            for path in role_paths:
                self.assertGreaterEqual(len(path.points), 2)
                length = sum(cylindrical_distance(a, b, path.radius) for a, b in zip(path.points, path.points[1:]))
                self.assertGreater(length, config.min_path_length)

    def test_full_wrap_radial_paths_span_one_turn(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, infill_spacing=2.0, wall_count=1)
        occupied = box(0.0, 0.0, 2.0 * math.pi, 10.0)
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=occupied,
            top_surfaces=GeometryCollection(),
            bottom_surfaces=GeometryCollection(),
            needs_support=GeometryCollection(),
            support=GeometryCollection(),
            theta_window=(0.0, 2.0 * math.pi),
            is_full_wrap=True,
        )
        paths = build_layer_toolpaths(layer, config).paths_by_role.get("INFILL_RADIAL", [])
        full = [path for path in paths if abs((path.points[-1][0] - path.points[0][0]) - 2 * math.pi) < 1e-6]
        self.assertTrue(full)

    def test_sparse_infill_is_not_generated_inside_solid_infill(self):
        config = SlicerConfig(wall_count=1, nozzle_width=0.4, infill_spacing=1.0)
        occupied = box(0.0, 0.0, 1.0, 10.0)
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=occupied,
            top_surfaces=occupied,
            bottom_surfaces=GeometryCollection(),
            needs_support=GeometryCollection(),
            support=GeometryCollection(),
            theta_window=(0.0, 1.0),
            is_full_wrap=False,
        )
        paths = build_layer_toolpaths(layer, config).paths_by_role
        self.assertGreater(len(paths.get("SOLID_INFILL", [])), 0)
        self.assertEqual(paths.get("INFILL_RADIAL", []), [])
        self.assertEqual(paths.get("INFILL_AXIAL", []), [])

    def test_gap_support_layers_are_emitted_below_floating_model(self):
        config = SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0, bed_radius=15.0, wall_count=1)
        layers = slice_mesh(make_fixture("solid_box"), config, max_layers=1)
        gap_layers = gap_support_slices(layers, config)
        self.assertTrue(gap_layers)
        self.assertTrue(all(layer.index < 0 for layer in gap_layers))
        self.assertLess(gap_layers[-1].radius, layers[0].radius)

        paths = build_toolpaths(layers, config)
        support_layers = [layer for layer in paths if layer.layer_index < 0]
        self.assertTrue(support_layers)
        self.assertTrue(any(layer.paths_by_role.get("SUPPORT") for layer in support_layers))

    def test_gap_support_can_be_disabled(self):
        config = SlicerConfig(
            theta_step_deg=15,
            axial_step=2.0,
            layer_height=1.0,
            bed_radius=15.0,
            generate_gap_support=False,
        )
        layers = slice_mesh(make_fixture("solid_box"), config, max_layers=1)
        self.assertEqual(gap_support_slices(layers, config), [])

    def test_parallel_toolpaths_match_serial_output(self):
        slice_config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, wall_count=2)
        layers = slice_mesh(make_fixture("tube"), slice_config, max_layers=4)
        serial = build_toolpaths(layers, SlicerConfig(toolpath_workers=1))
        parallel = build_toolpaths(layers, SlicerConfig(toolpath_workers=2))
        self.assertEqual(parallel, serial)


if __name__ == "__main__":
    unittest.main()
