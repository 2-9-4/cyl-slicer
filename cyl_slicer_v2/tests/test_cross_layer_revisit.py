from __future__ import annotations

from collections import Counter
from dataclasses import replace
import math
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (_island_revisit_anchor, _revisit_anchor_distance,
    _ordered_layer_paths_with_islands, emit_gcode)
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


def make_layer(index, *, remove_a=False, shift_b=0, radius=31):
    a = Toolpath("SOLID_INFILL", radius, [(0.1, 0), (0.1, 12)])
    b = Toolpath("SOLID_INFILL", radius, [(1.1+shift_b, 0), (1.1+shift_b, 12)])
    groups = ([] if remove_a else [{"SOLID_INFILL": [a]}]) + [{"SOLID_INFILL": [b]}]
    return LayerToolpaths(index, radius, {"SOLID_INFILL": ([] if remove_a else [a]) + [b]},
                          island_paths=groups)


def order_by_layer(gcode):
    return [tuple(int(x) for x in line.split(":", 1)[1].split(","))
            for line in gcode.splitlines() if line.startswith(";COOLING_ISLAND_ORDER:")]


def final_e_per_layer(gcode):
    results = []
    e = 0.0
    for line in gcode.splitlines():
        if line.startswith(";LAYER:"):
            if results:
                results[-1] = e
            results.append(0.0)
            e = 0.0
        elif line.startswith("G1 ") and " E" in (" "+line):
            import re
            m = re.search(r"\bE(-?\d+(?:\.\d+)?)", line)
            if m:
                e = float(m.group(1))
    if results:
        results[-1] = e
    return results


class CrossLayerRevisitTests(unittest.TestCase):
    def cfg(self, **opts):
        return replace(SlicerConfig(bed_radius=30, retract_mm=0,
                cooling_min_layer_time_s=0,
                cooling_min_island_print_time_s=0,
                cooling_revisit_order="avoid-last"), **opts)

    def test_avoids_immediate_revisit_across_layers(self):
        layers = [make_layer(i) for i in range(3)]
        code = emit_gcode(layers, self.cfg())
        self.assertEqual(order_by_layer(code), [(0,1), (0,1), (0,1)])
        nearest = emit_gcode(layers, self.cfg(cooling_revisit_order="nearest"))
        self.assertEqual(order_by_layer(nearest), [])
        # Nearest-first on layer 1 would restart on previously last island B.
        previous_last = _island_revisit_anchor(layers[0].island_paths[1])
        selected = _ordered_layer_paths_with_islands(
            layers[1], GcodeState(theta_accum=1.1, axial=12),
            self.cfg(), avoid_first_anchor=previous_last)
        self.assertEqual([i for _,i in selected], [0,1])
        baseline = _ordered_layer_paths_with_islands(
            layers[1], GcodeState(theta_accum=1.1, axial=12),
            self.cfg(cooling_revisit_order="nearest"))
        self.assertEqual([i for _,i in baseline], [1,0])
        self.assertEqual(final_e_per_layer(code), final_e_per_layer(nearest))
        self.assertNotIn("\nG4",code)

    def test_island_identity_survives_permutation_and_periodic_rotary_shift(self):
        a = _island_revisit_anchor(make_layer(0).island_paths[1])
        b = _island_revisit_anchor(make_layer(1,shift_b=2*math.pi).island_paths[1])
        self.assertAlmostEqual(_revisit_anchor_distance(a,b,31), 0, places=7)
        # Re-ordered group list: choose its current index, not stale index.
        original = make_layer(0)
        next_layer = make_layer(1)
        next_layer.island_paths = list(reversed(next_layer.island_paths))
        planned = _ordered_layer_paths_with_islands(
            next_layer, GcodeState(theta_accum=1.1, axial=12),
            self.cfg(), avoid_first_anchor=a)
        self.assertEqual([i for _,i in planned], [1,0])

    def test_one_island_or_no_matching_island_falls_back(self):
        one = make_layer(1,remove_a=True)
        last = _island_revisit_anchor(make_layer(0).island_paths[1])
        paths = _ordered_layer_paths_with_islands(one,GcodeState(theta_accum=1.1),
                 self.cfg(),avoid_first_anchor=last)
        self.assertEqual([i for _,i in paths],[0])
        far = make_layer(1,shift_b=.65)
        alt = _ordered_layer_paths_with_islands(far,
                 GcodeState(theta_accum=1.75, axial=0),
                 self.cfg(cooling_revisit_match_mm=5),avoid_first_anchor=last)
        self.assertEqual([i for _,i in alt],[1,0])

    def test_separate_model_and_support_never_match(self):
        layer = make_layer(0)
        support = {"SUPPORT":[Toolpath("SUPPORT",31,[(1.1,0),(1.1,12)])]}
        self.assertTrue(math.isinf(_revisit_anchor_distance(
            _island_revisit_anchor(layer.island_paths[1]),
            _island_revisit_anchor(support),31)))

    def test_invalid_inputs_and_cli_profiles(self):
        for profile in ("quality","legacy-cli","library"):
            parser = build_parser()
            args=parser.parse_args(["--profile",profile])
            self.assertEqual(_profile_value(args,"cooling_revisit_order"),"nearest")
            args=parser.parse_args(["--profile",profile,
                "--cooling-revisit-order","avoid-last",
                "--cooling-revisit-match-mm","8"])
            self.assertEqual(_profile_value(args,"cooling_revisit_order"),"avoid-last")
            self.assertEqual(_profile_value(args,"cooling_revisit_match_mm"),8)
        for bad in ({"cooling_revisit_order":"spread"},
                    {"cooling_revisit_match_mm":0},
                    {"cooling_revisit_match_mm":-1},
                    {"cooling_revisit_match_mm":math.nan}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                replace(self.cfg(), **bad).validate()


if __name__=="__main__":
    unittest.main()
