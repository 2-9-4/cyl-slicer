from dataclasses import replace
import math
from pathlib import Path
import unittest

from shapely.geometry import GeometryCollection, Point, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.layers import merge_support_body, propagate_support, slice_mesh
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.periodic import periodic_copies
from cyl_slicer_v2.toolpaths import (
    build_layer_toolpaths, gap_support_slices, meandered_radial_grid_paths,
)


EMPTY = GeometryCollection()
TURN = 2 * math.pi


def stack():
    return [LayerSlice(i, 30.08 + .16*i, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY)
            for i in range(7)]


class SupportBodyMergingTests(unittest.TestCase):
    def cfg(self, **overrides):
        return replace(SlicerConfig(support_min_island_area=0, support_min_region_area=0,
                                    support_bearing_min_area=0), **overrides)

    def test_merge_preserves_original_branches_and_large_gaps(self):
        shape = unary_union([box(.2, 0, .3, 8), box(.34, 0, .44, 8), box(.6, 0, .7, 8)])
        merged = merge_support_body(shape, 30, self.cfg())
        self.assertTrue(merged.covers(shape))
        self.assertEqual(len(merged.geoms), 2)
        self.assertTrue(merged.covers(Point(.32, 4)))
        self.assertTrue(merged.disjoint(box(.45, 0, .59, 8)))

    def test_merge_distance_is_physical_and_zero_disables(self):
        shape = unary_union([box(.2, 0, .3, 8), box(.34, 0, .44, 8)])
        self.assertEqual(merge_support_body(shape, 30, self.cfg()).geom_type, "Polygon")
        self.assertEqual(merge_support_body(shape, 60, self.cfg()).geom_type, "MultiPolygon")
        self.assertTrue(merge_support_body(shape, 30, self.cfg(support_merge_gap=0)).equals(shape))

    def test_seam_merge_keeps_periodic_connections_and_large_cutout(self):
        shape = unary_union([box(0, 0, .1, 8), box(TURN-.1, 0, TURN-.04, 8)])
        merged = merge_support_body(shape, 30, self.cfg())
        self.assertTrue(merged.covers(shape))
        self.assertTrue(merged.covers(Point(TURN-.02, 4)))
        self.assertTrue(merged.disjoint(box(.2, 0, TURN-.2, 8)))
        connected = periodic_copies(merged)
        self.assertTrue(any(poly.covers(Point(-.02, 4)) and poly.covers(Point(.02, 4))
                            for poly in connected.geoms))

    def test_distinct_contact_planes_share_body_but_keep_interface_pads(self):
        left, right = box(.2, 0, .3, 8), box(.34, 0, .44, 8)
        layers = stack()
        layers[4].needs_support = left
        layers[6].needs_support = right
        propagate_support(layers, self.cfg())
        self.assertEqual(layers[2].support.geom_type, "Polygon")
        self.assertLess(layers[2].support_interface.symmetric_difference(left.union(right)).area, 1e-12)
        self.assertTrue(layers[3].support.disjoint(left))
        self.assertTrue(layers[2].support.covers(left.union(right)))

    def test_merging_does_not_fill_another_branches_pending_contact_gap(self):
        sides = box(.2, 0, .3, 8).union(box(.34, 0, .44, 8))
        middle = box(.315, 2, .325, 6)
        for gap in (None, .16):
            with self.subTest(contact_gap=gap):
                layers = stack()
                layers[6].needs_support = sides
                layers[4].needs_support = middle
                propagate_support(layers, self.cfg(support_contact_gap_mm=gap))
                self.assertEqual(layers[3].support.intersection(middle).area, 0)
                self.assertTrue(layers[2].support.covers(middle))
                self.assertTrue(layers[3].support.covers(Point(.32, 1)))

    def test_model_clearance_is_applied_after_merging(self):
        sides = box(.2, 0, .3, 8).union(box(.34, 0, .44, 8))
        model = box(.315, 2, .325, 6)
        layers = stack()
        layers[6].needs_support = sides
        layers[3].occupied = model
        propagate_support(layers, self.cfg(support_lateral_gap_mm=.2))
        self.assertTrue(layers[3].support.disjoint(model))
        self.assertGreaterEqual(layers[3].support.distance(model)*layers[3].radius, .2 - 1e-6)
        self.assertTrue(layers[3].support.covers(Point(.32, 1)))

    def test_offset_cube_foundation_is_one_body_with_fewer_internal_meander_walls(self):
        mesh = load_mesh(Path(__file__).resolve().parents[2] / "offset_cube.stl")
        config = SlicerConfig(axial_min=30)
        separate = replace(config, support_merge_gap=0)
        layers = slice_mesh(mesh, config)
        baseline = slice_mesh(mesh, separate)
        for actual, before in zip(layers, baseline):
            self.assertTrue(actual.occupied.equals(before.occupied))
            self.assertTrue(actual.needs_support.equals(before.needs_support))
            self.assertTrue(actual.support_interface.equals(before.support_interface))
        foundation = gap_support_slices(layers, config)
        unmerged = gap_support_slices(baseline, separate)
        target = foundation[config.raft_layers]
        before = unmerged[config.raft_layers]
        # Contact clearance may alter the number of foundation slabs; support
        # merging must not make the two schedules disagree.
        self.assertEqual(target.index, before.index)
        self.assertLess(target.index, 0)
        self.assertEqual(target.support.geom_type, "Polygon")
        self.assertGreater(len(before.support.geoms), 10)
        self.assertLess(before.support.difference(target.support).area * target.radius, 1e-7)
        planned = build_layer_toolpaths(target, config)
        self.assertEqual(len(planned.island_paths), 1)
        self.assertLess(len(planned.paths_by_role["SUPPORT"]),
                        len(build_layer_toolpaths(before, separate).paths_by_role["SUPPORT"]))
        def meander_wall_length(layer):
            paths = meandered_radial_grid_paths(layer.support, "SUPPORT", layer.radius,
                                                config.support_spacing, config)
            return sum(abs(b[1]-a[1]) for path in paths for a, b in zip(path.points, path.points[1:])
                       if abs(b[0]-a[0]) < 1e-10)
        self.assertLess(meander_wall_length(target), meander_wall_length(before)*.5)
        code = emit_gcode([planned], config)
        old_code = emit_gcode([build_layer_toolpaths(before, separate)], separate)
        self.assertLess(code.count("; retract"), old_code.count("; retract"))

    def test_invalid_merge_distances_rejected(self):
        for distance in (-1, math.inf, math.nan):
            with self.subTest(distance=distance), self.assertRaisesRegex(ValueError, "support_merge_gap"):
                SlicerConfig(support_merge_gap=distance).validate()


if __name__ == "__main__":
    unittest.main()
