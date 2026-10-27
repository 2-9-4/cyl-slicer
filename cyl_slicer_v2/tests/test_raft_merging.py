from dataclasses import replace
import math
from pathlib import Path
import unittest

from shapely.geometry import GeometryCollection, LineString, Point, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, gap_support_slices


EMPTY = GeometryCollection()
TURN = 2 * math.pi


def foundation(shape, config, lower=0):
    model = LayerSlice(0, config.bed_radius + 1, shape, EMPTY, EMPTY, EMPTY, EMPTY,
                       theta_window=(lower, lower + TURN))
    return gap_support_slices([model], config)


class RaftMergingTests(unittest.TestCase):
    def test_small_angular_gaps_merge_only_in_raft(self):
        shape = unary_union([box(1, 0, 1.1, 4), box(1.15, 0, 1.25, 4)])
        config = SlicerConfig(raft_layers=2, raft_axial_gap=0)
        layers = foundation(shape, config)
        for layer in layers[:2]:
            self.assertEqual(layer.support.geom_type, "Polygon")
            self.assertTrue(layer.support.covers(shape))
            self.assertTrue(layer.support.covers(Point(1.125, 2)))
        self.assertTrue(layers[2].support.equals(shape))
        self.assertEqual(shape.geom_type, "MultiPolygon")

    def test_merging_uses_physical_radius_and_preserves_large_cutouts(self):
        shape = unary_union([box(1, 0, 1.1, 4), box(1.15, 0, 1.25, 4), box(1.35, 0, 1.45, 4)])
        for radius, count in ((30, 2), (60, 3)):
            config = SlicerConfig(bed_radius=radius, raft_axial_gap=0)
            raft = foundation(shape, config)[0]
            self.assertEqual(len(raft.support.geoms), count)
            self.assertTrue(raft.support.disjoint(box(1.26, 0, 1.34, 4)))

    def test_zero_disables_raft_merging(self):
        shape = unary_union([box(1, 0, 1.1, 4), box(1.101, 0, 1.2, 4)])
        config = SlicerConfig(raft_merge_gap_mm=0, raft_axial_gap=0)
        self.assertTrue(foundation(shape, config)[0].support.equals(shape))

    def test_merging_across_periodic_cut_preserves_large_gap(self):
        for lower in (0, .3):
            with self.subTest(lower=lower):
                shape = unary_union([box(lower + .02, 0, lower + .4, 4),
                                     box(lower + TURN - .4, 0, lower + TURN - .02, 4)])
                config = SlicerConfig(raft_axial_gap=0)
                raft = foundation(shape, config, lower)[0]
                self.assertTrue(raft.support.covers(Point(lower, 2)))
                self.assertTrue(raft.support.covers(Point(lower + TURN, 2)))
                cutout = box(lower + .5, 0, lower + TURN - .5, 4)
                self.assertTrue(raft.support.disjoint(cutout))
                paths = build_layer_toolpaths(raft, config)
                self.assertEqual(len(paths.island_paths), 1)
                for path in paths.paths_by_role["SUPPORT"]:
                    self.assertTrue(LineString(path.points).disjoint(cutout))
                relieved = foundation(shape, replace(config, raft_axial_gap=1), lower)[0]
                self.assertFalse(relieved.support.covers(Point(lower, 3.5)))
                self.assertFalse(relieved.support.covers(Point(lower + TURN, 3.5)))

    def test_angular_merge_preserves_axial_relief(self):
        shape = unary_union([box(t, x, t + .1, x + 2)
                             for t in (1, 1.15) for x in (0, 8)])
        raft = foundation(shape, SlicerConfig())[0]
        self.assertTrue(raft.support.covers(Point(1.125, 1)))
        self.assertFalse(raft.support.covers(Point(1.125, 5)))

    def test_invalid_merge_distance_rejected(self):
        for distance in (-1, math.inf, math.nan):
            with self.subTest(distance=distance), self.assertRaisesRegex(ValueError, "raft_merge_gap_mm"):
                SlicerConfig(raft_merge_gap_mm=distance).validate()

    def test_zero_radius_fixture_uses_positive_first_sample_radius(self):
        raft = foundation(box(1, 0, 2, 4), SlicerConfig(bed_radius=0))[0]
        self.assertFalse(raft.support.is_empty)

    def test_offset_cube_first_layer_is_one_raft_with_few_retractions(self):
        root = Path(__file__).resolve().parents[2]
        config = SlicerConfig(axial_min=30)
        layers = slice_mesh(load_mesh(root / "offset_cube.stl"), config)
        raft_layers = gap_support_slices(layers, config)
        unmerged = gap_support_slices(layers, replace(config, raft_merge_gap_mm=0))
        fragmented = build_layer_toolpaths(unmerged[0], config)
        merged = build_layer_toolpaths(raft_layers[0], config)
        self.assertGreater(len(fragmented.paths_by_role["SUPPORT"]), 500)
        self.assertEqual(len(merged.island_paths), 1)
        self.assertLess(len(merged.paths_by_role["SUPPORT"]), 100)
        self.assertTrue(raft_layers[config.raft_layers].support.equals(unmerged[config.raft_layers].support))
        code = emit_gcode([build_layer_toolpaths(layer, config) for layer in raft_layers[:config.raft_layers]], config)
        for block in code.split(";LAYER:")[1:]:
            self.assertLessEqual(block.count("; retract"), 2)
        segments = parse_extrusion_segments(code)
        directions = [s[5] - s[3] for s in segments]
        self.assertTrue(any(delta > 0 for delta in directions))
        self.assertTrue(any(delta < 0 for delta in directions))


if __name__ == "__main__":
    unittest.main()
