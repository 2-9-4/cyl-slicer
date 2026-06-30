from __future__ import annotations

import math
from pathlib import Path
import unittest

import numpy as np
from shapely.geometry import LineString, Point
import trimesh

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import cylindrical_to_cartesian
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.gcode import emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.intersection_backend import _polygon_parts, cylindrical_intersection_segments, sample_intersection_layer
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.mesh_io import axial_bounds, load_mesh, radial_bounds
from cyl_slicer_v2.toolpaths import build_toolpaths


class IntersectionBackendTests(unittest.TestCase):
    def test_cylindrical_edge_crossings_are_exact_for_non_radial_edges(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [3.0, 4.0, 0.0],
                    [8.0, 0.0, 0.0],
                    [3.0, 4.0, 1.0],
                ],
                dtype=float,
            ),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )

        points, edges = cylindrical_intersection_segments(mesh, 6.0)

        self.assertEqual(len(points), 2)
        self.assertEqual(len(edges), 1)
        self.assertTrue(np.allclose(np.hypot(points[:, 0], points[:, 1]), 6.0, atol=1e-9))

    def test_cylindrical_intersection_can_tessellate_curved_face_crossing(self):
        mesh = trimesh.Trimesh(
            vertices=np.array(
                [
                    [3.0, 4.0, 0.0],
                    [8.0, 0.0, 0.0],
                    [3.0, 4.0, 8.0],
                ],
                dtype=float,
            ),
            faces=np.array([[0, 1, 2]]),
            process=False,
        )

        points, edges = cylindrical_intersection_segments(
            mesh,
            6.0,
            max_theta_step=np.radians(1.0),
            max_axial_step=1.0,
        )

        self.assertGreater(len(points), 2)
        self.assertGreater(len(edges), 1)
        self.assertTrue(np.allclose(np.hypot(points[:, 0], points[:, 1]), 6.0, atol=1e-9))

    def test_tube_intersection_layer_is_nonempty_and_full_wrap(self):
        mesh = make_fixture("tube")
        geom, _, full = sample_intersection_layer(mesh, 11.0, SlicerConfig(layer_backend="intersection"))
        self.assertFalse(geom.is_empty)
        self.assertTrue(full)

    def test_solid_box_intersection_layer_is_nonempty_island(self):
        mesh = make_fixture("solid_box")
        geom, _, full = sample_intersection_layer(mesh, 20.0, SlicerConfig(layer_backend="intersection"))
        self.assertFalse(geom.is_empty)
        self.assertFalse(full)

    def test_intersection_slice_mesh_smoke(self):
        layers = slice_mesh(
            make_fixture("solid_box"),
            SlicerConfig(layer_backend="intersection", layer_height=1.0, slice_workers=1),
            max_layers=3,
        )
        self.assertEqual(len(layers), 3)
        self.assertTrue(any(not layer.occupied.is_empty for layer in layers))

    def test_parallel_intersection_matches_serial_area(self):
        mesh = make_fixture("solid_box")
        serial = slice_mesh(
            mesh,
            SlicerConfig(layer_backend="intersection", layer_height=1.0, slice_workers=1),
            max_layers=3,
        )
        parallel = slice_mesh(
            mesh,
            SlicerConfig(layer_backend="intersection", layer_height=1.0, slice_workers=2),
            max_layers=3,
        )
        self.assertEqual([round(layer.occupied.area, 8) for layer in serial], [round(layer.occupied.area, 8) for layer in parallel])

    def test_cubeandring_intersection_layers_64_to_66_no_full_turn_perimeters(self):
        root = Path(__file__).resolve().parents[2]
        config = SlicerConfig(layer_backend="intersection", slice_workers=1)
        layers = slice_mesh(load_mesh(root / "cubeandring.stl"), config, max_layers=67)
        gcode = emit_gcode(build_toolpaths(layers, config), config)
        bad = []
        for layer, role, x0, y0, x1, y1 in parse_extrusion_segments(gcode):
            if layer in {64, 65, 66} and role in {"EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"}:
                if abs(y1 - y0) >= 300.0:
                    bad.append((layer, role, x0, y0, x1, y1))
        self.assertEqual(bad, [])

    def test_cubeandring_intersection_layer_58_at_fine_height_is_not_fully_infilled(self):
        root = Path(__file__).resolve().parents[2]
        mesh = load_mesh(root / "cubeandring.stl")
        config = SlicerConfig(layer_backend="intersection", layer_height=0.2, slice_workers=1, intersection_refine_factor=0)
        layers = slice_mesh(mesh, config, max_layers=60)
        layer = layers[58]
        frame_area = 2.0 * np.pi * (mesh.bounds[1, 2] - mesh.bounds[0, 2])

        sample_count = 160
        theta_values = np.linspace(0.0, 2.0 * np.pi, sample_count, endpoint=False) + np.pi / sample_count
        axial_values = np.linspace(mesh.bounds[0, 2], mesh.bounds[1, 2], sample_count)
        tt, zz = np.meshgrid(theta_values, axial_values, indexing="ij")
        expected = mesh.contains(cylindrical_to_cartesian(tt.ravel(), zz.ravel(), layer.radius))
        actual = (
            np.array([layer.occupied.covers(Point(theta, axial)) for theta, axial in zip(tt.ravel(), zz.ravel())])
            if not layer.occupied.is_empty
            else np.zeros_like(expected)
        )

        self.assertLess(layer.occupied.area, frame_area * 0.25)
        self.assertLess(float(np.mean(~expected & actual)), 0.025)

    def test_cubeandring2_intersection_layers_do_not_drop_or_fill_frames(self):
        root = Path(__file__).resolve().parents[2]
        mesh = load_mesh(root / "cubeandring2.stl")
        config = SlicerConfig(layer_backend="intersection", layer_height=0.25, intersection_refine_factor=0)
        r_min, _ = radial_bounds(mesh)
        frame_area = 2.0 * math.pi * (mesh.bounds[1, 2] - mesh.bounds[0, 2])

        def sample_area(layer_index: int) -> float:
            radius = r_min + config.layer_height * 0.5 + layer_index * config.layer_height
            geom, _, _ = sample_intersection_layer(mesh, radius, config)
            return 0.0 if geom.is_empty else geom.area

        for layer_index in (19, 20, 23, 45):
            with self.subTest(layer=layer_index):
                self.assertGreater(sample_area(layer_index), frame_area * 0.30)

        for layer_index in (21, 22, 24, 25, 44):
            with self.subTest(layer=layer_index):
                self.assertLess(sample_area(layer_index), frame_area * 0.90)

    def test_cubeandring2_layer_20_keeps_intersection_boundary_not_cell_grid(self):
        root = Path(__file__).resolve().parents[2]
        mesh = load_mesh(root / "cubeandring2.stl")
        config = SlicerConfig(layer_backend="intersection", layer_height=0.25, intersection_refine_factor=0)
        r_min, _ = radial_bounds(mesh)
        radius = r_min + config.layer_height * 0.5 + 20 * config.layer_height
        geom, _, _ = sample_intersection_layer(mesh, radius, config)

        self.assertFalse(geom.is_empty)
        self.assertLess(geom.bounds[1], -7.0)
        self.assertGreater(geom.bounds[3], 7.0)

        axis_aligned = 0
        segment_count = 0
        for poly in _polygon_parts(geom):
            for ring in (poly.exterior, *poly.interiors):
                coords = np.asarray(ring.coords)
                deltas = np.diff(coords, axis=0)
                segment_count += len(deltas)
                axis_aligned += int(np.count_nonzero((np.abs(deltas[:, 0]) < 1e-9) | (np.abs(deltas[:, 1]) < 1e-9)))
        self.assertLess(axis_aligned / segment_count, 0.05)

    def test_helix_intersection_does_not_merge_turns_into_full_width_sheet(self):
        root = Path(__file__).resolve().parents[2]
        mesh = load_mesh(root / "helix.stl")
        r_min, _ = radial_bounds(mesh)
        radius = r_min + 2.5
        geom, _, full = sample_intersection_layer(
            mesh,
            radius,
            SlicerConfig(layer_backend="intersection", intersection_refine_factor=0),
        )

        self.assertFalse(geom.is_empty)
        self.assertFalse(full)
        mid_turn = LineString([(0.0, 20.0), (2.0 * np.pi, 20.0)])
        self.assertLess(mid_turn.intersection(geom).length, 2.0 * np.pi * 0.9)

    def test_rotatecube_debug_layers_match_mesh_contains_sampling_without_rasterizing_faces(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(layer_backend="intersection", layer_height=1.0, slice_workers=1, intersection_refine_factor=2)
        layers = slice_mesh(mesh, config, max_layers=13)
        theta_values = np.linspace(0.0, 2.0 * np.pi, 60, endpoint=False) + np.pi / 60.0
        axial_values = np.linspace(mesh.bounds[0, 2], mesh.bounds[1, 2], 60)
        tt, zz = np.meshgrid(theta_values, axial_values, indexing="ij")

        for layer_index in range(7, 13):
            with self.subTest(layer=layer_index):
                layer = layers[layer_index]
                points = cylindrical_to_cartesian(tt.ravel(), zz.ravel(), layer.radius)
                expected = mesh.contains(points)
                actual = (
                    np.array([layer.occupied.covers(Point(theta, axial)) for theta, axial in zip(tt.ravel(), zz.ravel())])
                    if not layer.occupied.is_empty
                    else np.zeros_like(expected)
                )
                mismatch_ratio = float(np.mean(expected != actual))
                self.assertLess(mismatch_ratio, 0.10)

    def test_rotatecube_debug_default_intersection_does_not_fill_entire_layers(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(layer_backend="intersection", layer_height=1.0, slice_workers=1, intersection_refine_factor=0)
        layers = slice_mesh(mesh, config, max_layers=15)
        theta_values = np.linspace(0.0, 2.0 * np.pi, 80, endpoint=False) + np.pi / 80.0
        axial_values = np.linspace(mesh.bounds[0, 2], mesh.bounds[1, 2], 80)
        tt, zz = np.meshgrid(theta_values, axial_values, indexing="ij")
        frame_area = 2.0 * np.pi * (mesh.bounds[1, 2] - mesh.bounds[0, 2])

        for layer_index in (5, 6, 13, 14):
            with self.subTest(layer=layer_index):
                layer = layers[layer_index]
                points = cylindrical_to_cartesian(tt.ravel(), zz.ravel(), layer.radius)
                expected = mesh.contains(points)
                actual = (
                    np.array([layer.occupied.covers(Point(theta, axial)) for theta, axial in zip(tt.ravel(), zz.ravel())])
                    if not layer.occupied.is_empty
                    else np.zeros_like(expected)
                )
                self.assertLess(layer.occupied.area, frame_area * 0.90)
                self.assertLess(float(np.mean(expected != actual)), 0.20)

        layer = layers[14]
        seam_theta_values = np.r_[
            np.linspace(0.005, 0.5, 30),
            np.linspace(2.0 * np.pi - 0.5, 2.0 * np.pi - 0.005, 30),
        ]
        seam_tt, seam_zz = np.meshgrid(seam_theta_values, axial_values, indexing="ij")
        seam_points = cylindrical_to_cartesian(seam_tt.ravel(), seam_zz.ravel(), layer.radius)
        seam_expected = mesh.contains(seam_points)
        seam_actual = (
            np.array([layer.occupied.covers(Point(theta, axial)) for theta, axial in zip(seam_tt.ravel(), seam_zz.ravel())])
            if not layer.occupied.is_empty
            else np.zeros_like(seam_expected)
        )
        self.assertLess(float(np.mean(seam_expected & ~seam_actual)), 0.03)
        self.assertLess(float(np.mean(seam_expected != seam_actual)), 0.05)

    def test_rotatecube_debug_layer_34_at_fine_height_rejects_degenerate_full_span(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(layer_backend="intersection", layer_height=0.2, slice_workers=1, intersection_refine_factor=0)
        r_min, _ = radial_bounds(mesh)
        radius = r_min + config.layer_height * 0.5 + 34 * config.layer_height
        geom, _, full = sample_intersection_layer(mesh, radius, config)
        frame_area = 2.0 * np.pi * (mesh.bounds[1, 2] - mesh.bounds[0, 2])

        theta_values = np.linspace(0.0, 2.0 * np.pi, 100, endpoint=False) + np.pi / 100.0
        axial_values = np.linspace(mesh.bounds[0, 2], mesh.bounds[1, 2], 100)
        tt, zz = np.meshgrid(theta_values, axial_values, indexing="ij")
        expected = mesh.contains(cylindrical_to_cartesian(tt.ravel(), zz.ravel(), radius))
        actual = (
            np.array([geom.covers(Point(theta, axial)) for theta, axial in zip(tt.ravel(), zz.ravel())])
            if not geom.is_empty
            else np.zeros_like(expected)
        )

        self.assertFalse(full)
        self.assertLess(geom.area, frame_area * 0.55)
        self.assertLess(float(np.mean(~expected & actual)), 0.025)
        self.assertLess(float(np.mean(expected != actual)), 0.03)

    def test_rotatecube_export_layer_33_does_not_fill_axial_margin(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2export.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(
            layer_backend="intersection",
            intersection_refine_factor=0,
            layer_height=0.25,
            bed_radius=30.0,
            axial_margin=15.0,
        )
        r_min, _ = radial_bounds(mesh)
        radius = r_min + config.layer_height * 0.5 + 33 * config.layer_height
        geom, _, full = sample_intersection_layer(mesh, radius, config)
        z_min, z_max = axial_bounds(mesh, 0.0)

        self.assertFalse(full)
        self.assertGreaterEqual(geom.bounds[1], z_min - 1e-6)
        self.assertLessEqual(geom.bounds[3], z_max + 1e-6)

    def test_rotatecube_export_layer_44_handles_periodic_normalization(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2export.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(layer_backend="intersection", intersection_refine_factor=0, layer_height=0.2)
        r_min, _ = radial_bounds(mesh)
        radius = r_min + config.layer_height * 0.5 + 44 * config.layer_height
        geom, _, _ = sample_intersection_layer(mesh, radius, config)

        self.assertFalse(geom.is_empty)
        self.assertTrue(geom.is_valid)

    def test_rotatecube_debug_intersection_volume_stays_close_with_smooth_topology_voting(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(layer_backend="intersection", layer_height=1.0, slice_workers=1, intersection_refine_factor=2)
        layers = slice_mesh(mesh, config)
        sliced_volume = sum(
            0.0 if layer.occupied.is_empty else layer.occupied.area * layer.radius * config.layer_height
            for layer in layers
        )
        relative_error = abs(sliced_volume - mesh.volume) / mesh.volume
        self.assertLess(relative_error, 0.05)

    def test_rotatecube_debug_smooth_topology_voting_does_not_emit_raster_boundary(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        mesh = load_mesh(mesh_path)
        config = SlicerConfig(layer_backend="intersection", layer_height=1.0, intersection_refine_factor=2)
        layers = slice_mesh(mesh, config, max_layers=8)
        parts = _polygon_parts(layers[7].occupied)
        self.assertTrue(parts)
        coords = np.asarray(parts[0].exterior.coords)
        deltas = np.diff(coords, axis=0)
        axis_aligned = np.count_nonzero((np.abs(deltas[:, 0]) < 1e-8) | (np.abs(deltas[:, 1]) < 1e-8))
        self.assertLess(axis_aligned / len(deltas), 0.10)

    def test_intersection_refinement_defaults_to_off(self):
        self.assertEqual(SlicerConfig(layer_backend="intersection").intersection_refine_factor, 0)


if __name__ == "__main__":
    unittest.main()
