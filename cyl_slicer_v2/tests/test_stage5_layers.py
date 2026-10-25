from __future__ import annotations

import math
from pathlib import Path
import unittest

import shapely
from shapely.geometry import GeometryCollection, LineString, Point, Polygon, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.intersection_backend import _polygon_parts
from cyl_slicer_v2.layers import (
    _clean,
    _metric_area,
    _needs_support_from_previous,
    _support_regions,
    propagate_solid_surfaces,
    propagate_support,
    slice_mesh,
)
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import LayerSlice


def _corner_count(geom) -> int:
    return sum(max(0, len(poly.exterior.coords) - 1) for poly in _polygon_parts(geom)) + sum(
        max(0, len(ring.coords) - 1) for poly in _polygon_parts(geom) for ring in poly.interiors
    )


class LayerClassificationTests(unittest.TestCase):
    def test_clean_preserves_valid_cubeandring2_top_surface_components(self):
        # Reduced from the layer 48 exposed surface at radius 54.634266 mm.
        # buffer(0) on this valid MultiPolygon drops the large printable patch.
        fixture = Path(__file__).resolve().parent / "fixtures" / "cubeandring2_top_surfaces.wkt"
        surface = shapely.from_wkt(fixture.read_text())
        self.assertTrue(surface.is_valid)

        cleaned = _clean(surface)

        self.assertTrue(cleaned.equals(surface))
        self.assertAlmostEqual(cleaned.area, surface.area, places=12)
        self.assertTrue(cleaned.covers(Point(math.radians(295), -8.0)))

    def test_clean_repairs_invalid_surfaces_and_discards_non_polygonal_parts(self):
        bowtie = Polygon([(0, 0), (2, 2), (0, 2), (2, 0), (0, 0)])
        cleaned = _clean(bowtie)
        self.assertTrue(cleaned.is_valid)
        self.assertAlmostEqual(cleaned.area, 2.0)

        mixed = GeometryCollection([box(0, 0, 1, 1), LineString([(2, 0), (2, 1)]), Point(3, 0)])
        self.assertTrue(_clean(mixed).equals(box(0, 0, 1, 1)))
        self.assertTrue(_clean(LineString([(0, 0), (1, 1)])).is_empty)

    def test_first_occupied_layer_is_bottom_surface(self):
        layers = slice_mesh(make_fixture("solid_box"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0), max_layers=3)
        first = next(layer for layer in layers if not layer.occupied.is_empty)
        self.assertAlmostEqual(first.bottom_surfaces.area, first.occupied.area)

    def test_tube_middle_layers_have_no_false_surface_churn(self):
        config = SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=0.5, top_layers=1, bottom_layers=1)
        layers = slice_mesh(make_fixture("tube"), config)
        occupied = [layer for layer in layers if not layer.occupied.is_empty]
        middle = occupied[len(occupied) // 2]
        self.assertLess(middle.top_surfaces.area, middle.occupied.area * 0.05)
        self.assertLess(middle.bottom_surfaces.area, middle.occupied.area * 0.05)

    def test_solid_surfaces_propagate_in_expected_layer_directions(self):
        empty = GeometryCollection()
        occupied = box(0.0, 0.0, 2.0, 2.0)
        top_skin = box(0.0, 0.0, 1.0, 1.0)
        bottom_skin = box(1.0, 1.0, 2.0, 2.0)
        layers = [
            LayerSlice(
                index=i,
                radius=10.0 + i,
                occupied=occupied,
                top_surfaces=top_skin if i == 3 else empty,
                bottom_surfaces=bottom_skin if i == 0 else empty,
                needs_support=empty,
                support=empty,
            )
            for i in range(5)
        ]

        propagate_solid_surfaces(layers, SlicerConfig(top_layers=3, bottom_layers=2))

        self.assertTrue(layers[1].top_surfaces.covers(Point(0.5, 0.5)))
        self.assertTrue(layers[2].top_surfaces.covers(Point(0.5, 0.5)))
        self.assertTrue(layers[4].top_surfaces.is_empty)
        self.assertTrue(layers[1].bottom_surfaces.covers(Point(1.5, 1.5)))
        self.assertTrue(layers[2].bottom_surfaces.is_empty)

    def test_overhang_step_creates_support(self):
        layers = slice_mesh(make_fixture("overhang_step"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0))
        self.assertTrue(any(not layer.needs_support.is_empty for layer in layers))
        self.assertTrue(any(not layer.support.is_empty for layer in layers))

    def test_support_regions_have_valid_bounds(self):
        layers = slice_mesh(make_fixture("overhang_step"), SlicerConfig(theta_step_deg=10, axial_step=1.0, layer_height=1.0))
        for layer in layers:
            if layer.support.is_empty:
                continue
            t_min, z_min, t_max, z_max = layer.support.bounds
            self.assertGreaterEqual(t_min, 0.0)
            self.assertGreater(t_max, t_min)
            self.assertGreater(z_max, z_min)

    def test_tiny_support_islands_are_removed(self):
        config = SlicerConfig(support_min_island_area=0.05, support_min_region_area=1.0)
        tiny = box(0.0, 0.0, 0.001, 0.5)

        self.assertTrue(_support_regions(tiny, 10.0, config).is_empty)

    def test_surviving_small_support_regions_are_promoted(self):
        config = SlicerConfig(support_min_island_area=0.05, support_min_region_area=1.0)
        small = box(0.0, 0.0, 0.02, 0.5)

        promoted = _support_regions(small, 10.0, config)

        self.assertFalse(promoted.is_empty)
        self.assertGreaterEqual(_metric_area(promoted, 10.0), 0.9)

    def test_tiny_previous_contact_does_not_fully_support_layer_above(self):
        config = SlicerConfig(support_min_island_area=0.0, support_min_region_area=0.0, support_bearing_min_area=1.0)
        current = box(0.0, 0.0, 1.0, 1.0)
        weak_previous = box(0.04, 0.45, 0.06, 0.55)

        needs = _needs_support_from_previous(current, weak_previous, 10.0, config)

        self.assertTrue(needs.covers(Point(0.05, 0.5)))

    def test_tiny_model_contact_does_not_ground_carried_support(self):
        config = SlicerConfig(support_min_island_area=0.0, support_min_region_area=0.0, support_bearing_min_area=1.0)
        empty = GeometryCollection()
        carry_source = box(0.0, 0.0, 1.0, 1.0)
        layers = [
            LayerSlice(
                index=0,
                radius=9.5,
                occupied=empty,
                top_surfaces=empty,
                bottom_surfaces=empty,
                needs_support=empty,
                support=empty,
            ),
            LayerSlice(
                index=1,
                radius=10.0,
                occupied=box(0.04, 0.45, 0.06, 0.55),
                top_surfaces=empty,
                bottom_surfaces=empty,
                needs_support=empty,
                support=empty,
            ),
            LayerSlice(
                index=2,
                radius=10.5,
                occupied=empty,
                top_surfaces=empty,
                bottom_surfaces=empty,
                needs_support=carry_source,
                support=empty,
            ),
        ]

        propagate_support(layers, config)

        self.assertFalse(layers[1].support.covers(Point(0.05, 0.5)))
        self.assertTrue(layers[0].support.covers(Point(0.05, 0.5)))

    def test_support_gap_and_interface_layers_are_propagated_below_overhang(self):
        empty = GeometryCollection()
        needs_support = box(0.0, 0.0, 1.0, 1.0)
        layers = [
            LayerSlice(
                index=i,
                radius=10.0 + i,
                occupied=empty,
                top_surfaces=empty,
                bottom_surfaces=empty,
                needs_support=needs_support if i == 4 else empty,
                support=empty,
            )
            for i in range(5)
        ]
        config = SlicerConfig(
            support_gap_layers=1,
            support_interface_layers=2,
            support_min_island_area=0.0,
            support_min_region_area=0.0,
            support_bearing_min_area=0.0,
        )

        propagate_support(layers, config)

        self.assertTrue(layers[3].support.is_empty)
        self.assertFalse(layers[2].support.is_empty)
        self.assertFalse(layers[2].support_interface.is_empty)
        self.assertEqual(layers[2].support_interface_direction, "axial")
        self.assertFalse(layers[1].support_interface.is_empty)
        self.assertEqual(layers[1].support_interface_direction, "radial")
        self.assertFalse(layers[0].support.is_empty)
        self.assertTrue(layers[0].support_interface.is_empty)

    def test_cubeandring2_layers_do_not_seed_full_frame_geometry(self):
        mesh_path = Path(__file__).resolve().parent / "fixtures" / "cubeandring2.stl"
        mesh = load_mesh(mesh_path)
        layers = slice_mesh(
            mesh,
            SlicerConfig(layer_backend="intersection", layer_height=0.25, slice_workers=4),
        )
        frame_area = 2.0 * math.pi * (mesh.bounds[1, 2] - mesh.bounds[0, 2])

        layer25 = layers[25]
        self.assertGreater(layer25.occupied.area, 20.0)
        self.assertLess(layer25.support.area, frame_area * 0.25)

        layer46 = layers[46]
        self.assertFalse(layer46.is_full_wrap)
        self.assertLess(layer46.occupied.area, frame_area * 0.75)

    def test_cubeandring_early_growth_seeds_support(self):
        mesh_path = Path(__file__).resolve().parents[2] / "cubeandring.stl"
        layers = slice_mesh(
            load_mesh(mesh_path),
            SlicerConfig(layer_backend="intersection", layer_height=0.25, slice_workers=4),
            max_layers=4,
        )
        self.assertTrue(any(not layer.needs_support.is_empty for layer in layers[1:]))
        self.assertFalse(layers[0].support.is_empty)

    def test_cubeandring2_layers_79_and_80_do_not_spike_corner_count(self):
        mesh_path = Path(__file__).resolve().parent / "fixtures" / "cubeandring2.stl"
        layers = slice_mesh(
            load_mesh(mesh_path),
            SlicerConfig(layer_backend="intersection", layer_height=0.25, slice_workers=4),
        )
        counts = {index: _corner_count(layers[index].occupied) for index in range(78, 82)}
        neighbor_limit = max(counts[78], counts[81]) * 2

        self.assertLessEqual(counts[79], neighbor_limit)
        self.assertLessEqual(counts[80], neighbor_limit)

    def test_parallel_slicing_matches_serial_slicing(self):
        mesh = make_fixture("tube")
        serial = slice_mesh(mesh, SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0, slice_workers=1), max_layers=3)
        parallel = slice_mesh(mesh, SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0, slice_workers=2), max_layers=3)
        self.assertEqual(len(serial), len(parallel))
        for a, b in zip(serial, parallel):
            self.assertAlmostEqual(a.occupied.area, b.occupied.area)
            self.assertAlmostEqual(a.top_surfaces.area, b.top_surfaces.area)
            self.assertAlmostEqual(a.bottom_surfaces.area, b.bottom_surfaces.area)
            self.assertEqual(a.is_full_wrap, b.is_full_wrap)


if __name__ == "__main__":
    unittest.main()
