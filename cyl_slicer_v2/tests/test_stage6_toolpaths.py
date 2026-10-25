from __future__ import annotations

import math
from pathlib import Path
import unittest

from shapely.geometry import GeometryCollection, LineString, Point, Polygon, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import cylindrical_distance
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import LayerSlice, Toolpath
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, build_toolpaths, gap_support_slices, _postprocess_role_paths


class ToolpathTests(unittest.TestCase):
    def test_toolpath_points_lie_on_or_inside_source(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, wall_count=2)
        layer = next(layer for layer in slice_mesh(make_fixture("solid_box"), config, max_layers=3) if not layer.occupied.is_empty)
        layer_paths = build_layer_toolpaths(layer, config)
        for role, paths in layer_paths.paths_by_role.items():
            if role in {"SUPPORT", "SUPPORT_INTERFACE"}:
                source = layer.support
            else:
                source = layer.occupied.buffer(1e-8)
            for path in paths:
                for theta, axial in path.points:
                    physical_theta = (
                        theta % (2.0 * math.pi)
                        if role in {"SUPPORT", "SUPPORT_INTERFACE"} and theta > 2.0 * math.pi
                        else theta
                    )
                    self.assertTrue(source.covers(Point(physical_theta, axial)), role)

    def test_perimeter_count_matches_wall_count_for_wide_box(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0, wall_count=2)
        layer = next(layer for layer in slice_mesh(make_fixture("solid_box"), config, max_layers=3) if not layer.occupied.is_empty)
        paths = build_layer_toolpaths(layer, config).paths_by_role
        self.assertGreaterEqual(len(paths.get("EXTERNAL_PERIMETER", [])), 1)
        self.assertGreaterEqual(len(paths.get("INTERNAL_PERIMETER", [])), 1)

    def test_full_wrap_external_perimeter_is_two_continuous_window_runs(self):
        theta_max = 2.0 * math.pi
        occupied = Polygon(
            [
                (0.05, 10.0),
                (1.0, 10.0),
                (theta_max, 10.0),
                (theta_max, -10.0),
                (1.0, -10.0),
                (0.05, -10.0),
                (0.0, -10.0),
                (0.0, 10.0),
                (0.05, 10.0),
            ]
        )
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=occupied,
            top_surfaces=GeometryCollection(),
            bottom_surfaces=GeometryCollection(),
            needs_support=GeometryCollection(),
            support=GeometryCollection(),
            theta_window=(0.0, theta_max),
            is_full_wrap=True,
        )

        paths = build_layer_toolpaths(layer, SlicerConfig(wall_count=1)).paths_by_role.get("EXTERNAL_PERIMETER", [])

        self.assertEqual(len(paths), 2)
        self.assertTrue(all(len(path.points) >= 2 for path in paths))
        self.assertTrue(all(b[0] >= a[0] for path in paths for a, b in zip(path.points, path.points[1:])))
        for path in paths:
            self.assertAlmostEqual(path.points[-1][0] - path.points[0][0], theta_max)
            self.assertTrue(all(point[1] == path.points[0][1] for point in path.points))

    def test_cubeandring_layer_73_external_perimeters_are_continuous_runs(self):
        mesh_path = Path(__file__).resolve().parents[2] / "cubeandring.stl"
        layers = slice_mesh(
            load_mesh(mesh_path),
            SlicerConfig(layer_backend="intersection", intersection_refine_factor=0, layer_height=0.2, slice_workers=4),
            max_layers=74,
        )

        paths = build_layer_toolpaths(layers[73], SlicerConfig(layer_height=0.2)).paths_by_role.get("EXTERNAL_PERIMETER", [])

        self.assertEqual(len(paths), 2)
        boundary = layers[73].occupied.boundary.buffer(1e-8)
        self.assertTrue(all(len(path.points) >= 2 for path in paths))
        self.assertTrue(all(boundary.covers(LineString(path.points)) for path in paths))

    def test_perimeters_keep_tiny_triangle_corners_below_min_path_length(self):
        occupied = Polygon([(0.0, 0.0), (0.01, 0.0), (0.01, 0.01), (0.0, 0.0)])
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=occupied,
            top_surfaces=GeometryCollection(),
            bottom_surfaces=GeometryCollection(),
            needs_support=GeometryCollection(),
            support=GeometryCollection(),
        )

        paths = build_layer_toolpaths(layer, SlicerConfig(wall_count=1, min_path_length=0.6)).paths_by_role.get(
            "EXTERNAL_PERIMETER",
            [],
        )

        self.assertTrue(paths)
        points = {(round(theta, 4), round(axial, 4)) for path in paths for theta, axial in path.points}
        self.assertIn((0.01, 0.0), points)
        self.assertIn((0.01, 0.01), points)

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

    def test_sparse_support_radial_lines_are_meandered_in_local_angle_groups(self):
        empty = GeometryCollection()
        support = box(0.0, 0.0, math.pi / 2.0, 4.0)
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=empty,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=support,
            theta_window=(0.0, 2.0 * math.pi),
            is_full_wrap=True,
        )
        config = SlicerConfig(support_spacing=1.0)

        paths = build_layer_toolpaths(layer, config).paths_by_role.get("SUPPORT", [])

        self.assertTrue(paths)
        meanders = [path for path in paths if len(path.points) > 2]
        self.assertTrue(meanders)
        radial_segments = [
            (a, b)
            for path in meanders
            for a, b in zip(path.points, path.points[1:])
            if abs(a[1] - b[1]) < 1e-8 and abs(a[0] - b[0]) > 1e-8
        ]
        self.assertTrue(radial_segments)
        self.assertTrue(all(abs(a[0] - b[0]) <= math.radians(15.0) + 1e-8 for a, b in radial_segments))
        self.assertTrue(any(b[0] < a[0] for a, b in radial_segments))
        self.assertTrue(
            any(
                abs(a[0] - b[0]) < 1e-8 and abs(a[1] - b[1]) > 0.5
                for path in meanders
                for a, b in zip(path.points, path.points[1:])
            )
        )

    def test_sparse_support_has_coarse_axial_tie_lines(self):
        empty = GeometryCollection()
        support = box(0.0, 0.0, math.pi / 2.0, 4.0)
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=empty,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=support,
        )
        config = SlicerConfig(support_spacing=1.0, support_cross_spacing=2.0)

        paths = build_layer_toolpaths(layer, config).paths_by_role.get("SUPPORT", [])

        self.assertTrue(
            any(
                abs(a[0] - b[0]) < 1e-8 and abs(a[1] - b[1]) > 3.5
                for path in paths
                for a, b in zip(path.points, path.points[1:])
            )
        )

    def test_short_support_lines_are_filtered(self):
        empty = GeometryCollection()
        support = box(0.0, 0.0, 0.01, 2.0)
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=empty,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=support,
        )

        paths = build_layer_toolpaths(layer, SlicerConfig(min_path_length=0.6)).paths_by_role

        self.assertNotIn("SUPPORT", paths)

    def test_straight_lines_contained_in_longer_same_role_lines_are_filtered(self):
        config = SlicerConfig(min_path_length=0.6)
        outer = Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 10.0)])
        inner = Toolpath("SUPPORT", 10.0, [(0.0, 2.0), (0.0, 4.0)])

        self.assertEqual(_postprocess_role_paths([outer, inner], config), [outer])

    def test_support_interface_uses_dedicated_role_and_removes_sparse_overlap(self):
        empty = GeometryCollection()
        interface = box(0.0, 0.0, 1.0, 4.0)
        layer = LayerSlice(
            index=0,
            radius=10.0,
            occupied=empty,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=interface,
            support_interface=interface,
            support_interface_direction="axial",
        )
        config = SlicerConfig(nozzle_width=0.5, support_interface_spacing=0.5)

        paths = build_layer_toolpaths(layer, config).paths_by_role

        self.assertTrue(paths.get("SUPPORT_INTERFACE"))
        self.assertNotIn("SUPPORT", paths)
        self.assertTrue(
            all(abs(path.points[0][0] - path.points[-1][0]) < 1e-8 for path in paths["SUPPORT_INTERFACE"])
        )

    def test_rotatecube_debug_layers_87_to_95_do_not_emit_floating_support(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        config = SlicerConfig(
            layer_backend="intersection",
            intersection_refine_factor=0,
            layer_height=0.2,
            wall_count=3,
            support_min_island_area=10.0,
            support_min_region_area=10.0,
            support_bearing_min_area=25.0,
            support_merge_gap=10.0,
            bed_radius=30.0,
            slice_workers=4,
            toolpath_workers=4,
        )
        layers = slice_mesh(load_mesh(mesh_path), config, max_layers=100)
        paths_by_index = {layer.layer_index: layer.paths_by_role for layer in build_toolpaths(layers, config)}

        for layer_index in range(87, 96):
            with self.subTest(layer=layer_index):
                self.assertNotIn("SUPPORT", paths_by_index[layer_index])
                self.assertNotIn("SUPPORT_INTERFACE", paths_by_index[layer_index])

    def test_cubeandring2_emits_real_support_without_floating_midlayer_support(self):
        mesh_path = Path(__file__).resolve().parents[2] / "cubeandring2.stl"
        config = SlicerConfig(
            layer_backend="intersection",
            intersection_refine_factor=0,
            layer_height=0.2,
            wall_count=3,
            support_min_island_area=10.0,
            support_min_region_area=10.0,
            support_bearing_min_area=25.0,
            support_merge_gap=10.0,
            bed_radius=30.0,
            slice_workers=4,
            toolpath_workers=4,
        )
        layers = slice_mesh(load_mesh(mesh_path), config, max_layers=73)
        paths_by_index = {layer.layer_index: layer.paths_by_role for layer in build_toolpaths(layers, config)}

        self.assertTrue(
            any(paths_by_index[layer_index].get("SUPPORT") for layer_index in range(7, 14))
        )
        for layer_index in range(62, 73):
            with self.subTest(layer=layer_index):
                self.assertNotIn("SUPPORT", paths_by_index[layer_index])
                self.assertNotIn("SUPPORT_INTERFACE", paths_by_index[layer_index])

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

    def test_first_gap_support_layers_are_raft_with_axial_relief_gap(self):
        empty = GeometryCollection()
        occupied = box(0.0, 0.0, 2.0 * math.pi, 2.0).union(box(0.0, 10.0, 2.0 * math.pi, 12.0))
        first_model = LayerSlice(
            index=0,
            radius=12.0,
            occupied=occupied,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=empty,
            theta_window=(0.0, 2.0 * math.pi),
            is_full_wrap=True,
        )
        config = SlicerConfig(layer_height=0.5, bed_radius=10.0, raft_layers=2, nozzle_width=1.0, raft_grid_spacing=2.0)

        gap_layers = gap_support_slices([first_model], config)

        self.assertTrue(gap_layers[0].is_raft)
        self.assertTrue(gap_layers[1].is_raft)
        self.assertFalse(gap_layers[2].is_raft)
        self.assertEqual(gap_layers[0].raft_direction, "radial")
        self.assertEqual(gap_layers[1].raft_direction, "axial")
        self.assertIsNone(gap_layers[2].raft_direction)
        self.assertFalse(gap_layers[0].support.covers(Point(1.0, 6.0)))
        self.assertTrue(gap_layers[0].support.covers(Point(1.0, 1.0)))
        self.assertTrue(gap_layers[0].support.covers(Point(1.0, 11.0)))
        self.assertFalse(gap_layers[2].support.covers(Point(1.0, 6.0)))

        support_paths = build_layer_toolpaths(gap_layers[0], config).paths_by_role.get("SUPPORT", [])
        self.assertTrue(
            any(abs(path.points[-1][0] - path.points[0][0]) > 2.0 * math.pi * 0.95 for path in support_paths)
        )
        self.assertFalse(any(abs(path.points[-1][1] - path.points[0][1]) > 10.0 for path in support_paths))

        support_paths = build_layer_toolpaths(gap_layers[1], config).paths_by_role.get("SUPPORT", [])
        self.assertTrue(any(abs(path.points[-1][1] - path.points[0][1]) > 4.0 for path in support_paths))
        self.assertFalse(
            any(
                min(path.points[0][1], path.points[-1][1]) < 6.0 < max(path.points[0][1], path.points[-1][1])
                for path in support_paths
            )
        )
        self.assertFalse(
            any(abs(path.points[-1][0] - path.points[0][0]) > 2.0 * math.pi * 0.95 for path in support_paths)
        )

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
