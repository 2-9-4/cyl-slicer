from __future__ import annotations

from collections import Counter
from dataclasses import replace
import math
import re
import unittest

from shapely.geometry import box

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (
    MAX_MARLIN_Y_DEGREES, _one_direction_layer, _one_direction_path_runs,
    _safe_island_travel_route, emit_gcode, parse_extrusion_segments,
)
from cyl_slicer_v2.models import LayerToolpaths, Toolpath
from cyl_slicer_v2.periodic import periodic_copies
from cyl_slicer_v2.tests.test_marlin_rotary_precision import MarlinRotaryModel


TAU = 2*math.pi


def edges(paths):
    result = Counter()
    for path in paths:
        for i, (a, b) in enumerate(zip(path.points, path.points[1:])):
            endpoints = tuple(sorted((tuple(round(v, 9) for v in (p[0] % TAU, p[1]))
                                      for p in (a, b))))
            width = path.segment_widths[i] if path.segment_widths is not None else None
            result[(path.role, path.radius, endpoints, round(abs(b[0]-a[0]), 9), width)] += 1
    return result


def deposited_e(code):
    current = total = 0.
    for line in code.splitlines():
        match = re.search(r"\bE(-?[\d.]+)", line)
        if match and line.startswith(("G0 ", "G1 ", "G92 ")):
            value = float(match.group(1))
            if line.startswith("G1 X"):
                total += value-current
            current = value
    return total


class OneDirectionFirstLayersTests(unittest.TestCase):
    def config(self, **kwargs):
        return SlicerConfig(one_direction_first_layers=2,
                            cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0, **kwargs)

    def test_split_reverse_keeps_edges_widths_short_fragments_and_full_turns(self):
        config = self.config(min_path_length=.5)
        sources = [
            Toolpath("GAP_FILL", 30, [(6.1, 0), (6.5, 0), (6.4999, .001),
                                      (6.8, 1), (6.8, 1.01)],
                     segment_widths=[.24, .27, .32, .4]),
            Toolpath("SUPPORT", 30, [(0, 2), (-2*TAU, 2)], segment_widths=[.4]),
            Toolpath("EXTERNAL_PERIMETER", 30,
                     [(0, 0), (.5, 0), (.5, 3), (0, 3), (0, 0)], closed=True),
        ]
        for source in sources:
            with self.subTest(role=source.role):
                runs = _one_direction_path_runs(source, config)
                self.assertEqual(edges([source]), edges(runs))
                self.assertTrue(all(b[0] >= a[0] for path in runs
                                    for a, b in zip(path.points, path.points[1:])))
                code = emit_gcode([LayerToolpaths(-3, 30, {source.role: [source]})], config)
                normal = emit_gcode([LayerToolpaths(-3, 30, {source.role: [source]})],
                                    replace(config, one_direction_first_layers=0))
                self.assertAlmostEqual(deposited_e(code), deposited_e(normal), places=4)
                self.assertTrue(all(step >= 0 for step in MarlinRotaryModel().replay(code).move_steps))
        # The short negative corner remains in the output despite the whole-path threshold.
        tiny = _one_direction_path_runs(sources[0], config)[1]
        self.assertLess(math.hypot((tiny.points[-1][0]-tiny.points[0][0])*30,
                                   tiny.points[-1][1]-tiny.points[0][1]), .5)

    def test_first_n_include_raft_support_and_empty_layers_do_not_consume_count(self):
        rows = [Toolpath("SUPPORT", 30, [(.1, x), (6., x)]) for x in (0, .4, .8)]
        layers = [LayerToolpaths(i, 30, {"SUPPORT": rows}, travel_region=box(.1, 0, 6., 1))
                  for i in (-8, -7, -6)]
        layers.insert(1, LayerToolpaths(-100, 30))
        code = emit_gcode(layers, self.config())
        self.assertEqual(code.count(";ONE_DIRECTION_LAYER:"), 2)
        model = MarlinRotaryModel()
        moves = {}
        index = None
        for line in code.splitlines():
            if line.startswith(";LAYER:"):
                index = int(line.split(":")[1])
                moves[index] = []
            before = len(model.move_steps)
            model.feed(line)
            if index is not None:
                moves[index].extend(model.move_steps[before:])
        for index in (-8, -7):
            self.assertTrue(moves[index])
            self.assertTrue(all(step >= 0 for step in moves[index]))
        self.assertTrue(any(step < 0 for step in moves[-6]), "Normal bidirectional output must resume")

    def test_gap_crossings_are_forward_non_extruding_and_retracted_even_for_tiny_slit(self):
        # Gap travel is under two line widths, but its entire connection crosses air.
        rows = [Toolpath("SUPPORT", 30, [(.01, x), (6.27, x)]) for x in (0, .4, .8)]
        config = self.config()
        code = emit_gcode([LayerToolpaths(-20, 30, {"SUPPORT": rows},
                                         travel_region=box(.01, -.3, 6.27, 1.1))], config)
        blocks = code.split(";start of path")[1:]
        self.assertEqual(len(blocks), 3)
        for block in blocks[:-1]:
            self.assertIn("; retract", block)
        for block in blocks[1:]:
            self.assertIn("; unretract", block)
            travel = next(line for line in block.splitlines() if line.startswith("G0 X"))
            self.assertNotIn(" E", travel)
        segments = parse_extrusion_segments(code)
        self.assertEqual(len(segments), 3)
        self.assertTrue(all(abs((s[5]-s[3])-math.degrees(6.26)) < .002 for s in segments))
        self.assertTrue(all(step >= 0 for step in MarlinRotaryModel().replay(code).move_steps))

    def test_native_rotary_resets_keep_direction_and_deposited_volume_bounded(self):
        rows = [Toolpath("SUPPORT", 30, [(6., x), (.1, x)]) for x in range(24)]
        source = LayerToolpaths(-4, 30, {"SUPPORT": rows}, travel_region=box(.1, 0, 6., 24))
        config = self.config(max_volumetric_flow_mm3_s=6, slow_start_layers=10)
        code = emit_gcode([source], config)
        normal = emit_gcode([source], replace(config, one_direction_first_layers=0))
        model = MarlinRotaryModel().replay(code)
        self.assertGreater(model.native_resets, 1)
        self.assertTrue(all(step >= 0 for step in model.move_steps))
        self.assertLessEqual(model.max_native, MAX_MARLIN_Y_DEGREES)
        self.assertLessEqual(model.max_logical, MAX_MARLIN_Y_DEGREES)
        self.assertAlmostEqual(deposited_e(code), deposited_e(normal), places=4)

    def test_layer_preparation_retains_separate_island_views_without_mutating_sources(self):
        a = Toolpath("SOLID_INFILL", 30, [(.5, 1), (.1, 1)], segment_widths=[.3])
        b = Toolpath("SOLID_INFILL", 30, [(2, 10), (1, 10)], segment_widths=[.4])
        groups = [{a.role: [a]}, {b.role: [b]}]
        source = LayerToolpaths(0, 30, {a.role: [a, b]}, island_paths=groups)
        prepared = _one_direction_layer(source, self.config())
        self.assertEqual(len(prepared.island_paths), 2)
        self.assertIs(prepared.island_paths[0][a.role][0], prepared.paths_by_role[a.role][0])
        self.assertEqual(edges([a, b]), edges(prepared.paths_by_role[a.role]))
        self.assertEqual(a.points, [(.5, 1), (.1, 1)])
        self.assertEqual(a.segment_widths, [.3])

    def test_router_rejects_backtracking_and_keeps_retract_hop_fallback(self):
        region = box(0, 0, 1, 10).difference(box(.4, 4, .6, 6))
        periodic = periodic_copies(region)
        start, end = (.5, 2), (.5, 8)
        config = self.config(route_max_unretracted_mm=12)
        normal = replace(config, one_direction_first_layers=0)
        self.assertIsNotNone(_safe_island_travel_route(start, end, 10, periodic, normal))
        self.assertIsNone(_safe_island_travel_route(start, end, 10, periodic, config))
        paths = [Toolpath("SUPPORT", 10, [(.5, 1), start]),
                 Toolpath("SUPPORT", 10, [end, (.5, 9)])]
        # Bed radius must be below this synthetic 10 mm sample radius.
        code = emit_gcode([LayerToolpaths(-1, 10, {"SUPPORT": paths}, travel_region=region)],
                          replace(config, bed_radius=9))
        self.assertIn("; z hop", code)
        self.assertIn("; retract", code)
        self.assertNotIn(";ROUTED_TRAVEL:", code)

    def test_monotonic_top_skin_keeps_widths_matched_when_axial_pass_reverses(self):
        path = Toolpath("TOP_SKIN", 30, [(.2, 3), (.2, 2), (.2, 0)],
                        segment_widths=[.25, .4])
        config = self.config(monotonic_top_skin=True, retract_mm=0)
        code = emit_gcode([LayerToolpaths(0, 30, {path.role: [path]})], config)
        # Compute independently from the physical segment lengths and widths.
        h = config.layer_height
        expected = sum(length*(h*(width-h)+math.pi*h*h/4)
                       for length, width in ((1, .25), (2, .4))) / (math.pi*(1.75/2)**2)
        self.assertAlmostEqual(deposited_e(code), expected, places=5)
        self.assertTrue(all(step >= 0 for step in MarlinRotaryModel().replay(code).move_steps))

    def test_profiles_defaults_override_and_validation(self):
        for profile in ("quality", "legacy-cli", "library"):
            args = build_parser().parse_args(["--profile", profile])
            self.assertEqual(_profile_value(args, "one_direction_first_layers"), 0)
            args = build_parser().parse_args(["--profile", profile, "--one-direction-first-layers", "10"])
            self.assertEqual(_profile_value(args, "one_direction_first_layers"), 10)
        for value in (-1, 1.5, True, math.inf, math.nan):
            with self.subTest(value=value), self.assertRaises(ValueError):
                replace(SlicerConfig(), one_direction_first_layers=value).validate()


if __name__ == "__main__":
    unittest.main()
