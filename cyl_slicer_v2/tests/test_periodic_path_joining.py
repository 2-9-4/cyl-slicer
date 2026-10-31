from __future__ import annotations

from dataclasses import replace
import math
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import TWO_PI, cylindrical_distance
from cyl_slicer_v2.gcode import emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerSlice, Toolpath
from cyl_slicer_v2.periodic import fold_to_canonical
from cyl_slicer_v2.toolpaths import _join_periodic_seam_paths, build_layer_toolpaths


class PeriodicPathJoiningTests(unittest.TestCase):
    def test_model_wall_and_fill_are_continuous_across_cut(self):
        empty = GeometryCollection()
        model = fold_to_canonical(box(-0.2, 0, 0.2, 5), 0, 5)
        config = SlicerConfig(wall_count=1, gap_fill_enabled=False,
                              solid_infill_pattern="alternate", sparse_infill_pattern="both")
        for solid in (True, False):
            with self.subTest(solid=solid):
                layer = LayerSlice(1, 30, model, model if solid else empty, empty, empty, empty,
                                   theta_window=(0, TWO_PI), is_full_wrap=False)
                result = build_layer_toolpaths(layer, config)
                walls = result.paths_by_role["EXTERNAL_PERIMETER"]
                self.assertEqual(len(walls), 1)
                self.assertTrue(walls[0].closed)
                role = "SOLID_INFILL" if solid else "INFILL_RADIAL"
                fills = result.paths_by_role[role]
                self.assertTrue(fills)
                for path in fills:
                    for theta, _ in (path.points[0], path.points[-1]):
                        self.assertGreater(abs(math.remainder(theta, TWO_PI)), 0.01)
                    expected_length = (0.4 * 30 - 2 * (
                        config.wall_coverage_depth(1) - config.nozzle_width * config.infill_wall_overlap_fraction
                    ))
                    actual_length = sum(cylindrical_distance(a, b, 30)
                                        for a, b in zip(path.points, path.points[1:]))
                    self.assertAlmostEqual(actual_length, expected_length)
                # Each complete pass has a single entry; no artificial G0 or
                # retract is inserted at the seam inside a deposited path.
                for path in [*walls, *fills]:
                    state = GcodeState(theta_accum=path.points[0][0], axial=path.points[0][1], radial=30)
                    code = emit_toolpath(path, state, config, 0, retract_after=False)
                    self.assertEqual(sum(s.startswith("G0 X") for s in code), 1)
                    self.assertFalse(any("; retract" in s for s in code))

    def test_tiny_clipped_fill_end_is_kept_until_joined(self):
        empty = GeometryCollection()
        config = SlicerConfig(wall_count=1, gap_fill_enabled=False, min_path_length=0.1,
                              solid_infill_pattern="alternate", infill_wall_overlap_fraction=0)
        inset = config.wall_coverage_depth(1) / 30
        model = fold_to_canonical(box(-0.0001 - inset, 0, 0.2 + inset, 5), 0, 5)
        layer = LayerSlice(1, 30, model, model, empty, empty, empty, theta_window=(0, TWO_PI))
        paths = build_layer_toolpaths(layer, config).paths_by_role["SOLID_INFILL"]
        self.assertTrue(paths)
        for path in paths:
            length = sum(cylindrical_distance(a, b, 30) for a, b in zip(path.points, path.points[1:]))
            self.assertAlmostEqual(length, 6.003)

    def test_reversed_variable_width_fragment_keeps_widths_and_e(self):
        config = SlicerConfig(min_path_length=0, cooling_min_layer_time_s=0,
                              cooling_min_island_print_time_s=0)
        for cut in (0.0, 0.7):
            with self.subTest(cut=cut):
                left = Toolpath("GAP_FILL", 30, [(cut + TWO_PI - .2, 0),
                                (cut + TWO_PI - .1, 0), (cut + TWO_PI, 0)],
                                segment_widths=[.25, .30])
                right = Toolpath("GAP_FILL", 30, [(cut + .2, 0), (cut + .1, 0), (cut, 0)],
                                 segment_widths=[.40, .35])
                paths = _join_periodic_seam_paths([left, right], (cut, cut + TWO_PI))
                self.assertEqual(len(paths), 1)
                self.assertEqual(paths[0].segment_widths, [.25, .30, .35, .40])
                self.assertTrue(all(b[0] > a[0] for a, b in zip(paths[0].points, paths[0].points[1:])))
                amounts = []
                for path in [left, right, paths[0]]:
                    state = GcodeState(radial=30)
                    emit_toolpath(path, state, config, 0, retract_after=False)
                    amounts.append(state.extrusion)
                self.assertAlmostEqual(amounts[2], sum(amounts[:2]), places=10)

    def test_voids_branches_corners_and_role_changes_do_not_join(self):
        left = Toolpath("INFILL_RADIAL", 30, [(TWO_PI - .2, 0), (TWO_PI, 0)])
        right = Toolpath("INFILL_RADIAL", 30, [(0, 0), (.2, 0)])
        cases = [
            [left, replace(right, points=[(0, .001), (.2, .001)])],
            [left, right, replace(right, points=[(0, 0), (.2, .1)])],
            [left, replace(right, points=[(0, 0), (.2, .1)])],
            [left, replace(right, role="SOLID_INFILL")],
        ]
        for paths in cases:
            with self.subTest(paths=paths):
                self.assertEqual(len(_join_periodic_seam_paths(paths, (0, TWO_PI))), len(paths))

    def test_full_turns_and_closed_winding_are_unchanged(self):
        for sign in (1, -1):
            path = Toolpath("EXTERNAL_PERIMETER", 30, [(0, 0), (sign * TWO_PI, 0)])
            self.assertEqual(_join_periodic_seam_paths([path], (0, TWO_PI)), [path])
        path = Toolpath("EXTERNAL_PERIMETER", 30,
                        [(0, 0), (.1, 0), (.1, 1), (0, 1), (0, 0)], closed=True)
        self.assertEqual(_join_periodic_seam_paths([path], (0, TWO_PI)), [path])


if __name__ == "__main__":
    unittest.main()
