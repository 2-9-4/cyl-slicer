from __future__ import annotations

from collections import Counter
from dataclasses import replace
from itertools import groupby
import re
import unittest

from shapely.geometry import box

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _island_role_sequence, _ordered_layer_paths_with_islands, emit_gcode
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath
from cyl_slicer_v2.toolpaths import ROLE_ORDER


def line(role, theta, x, width=None):
    return Toolpath(role, 31, [(theta, x), (theta, x + 4)],
                    segment_widths=[width] if width is not None else None)


def fixture():
    a = {
        "EXTERNAL_PERIMETER": [line("EXTERNAL_PERIMETER", .1, 0)],
        "INTERNAL_PERIMETER": [line("INTERNAL_PERIMETER", .12, 0)],
        "TOP_SKIN": [line("TOP_SKIN", .15, 0)],
        "SOLID_INFILL": [line("SOLID_INFILL", .17, 0)],
        "INFILL_RADIAL": [line("INFILL_RADIAL", .18, 0)],
        "INFILL_AXIAL": [line("INFILL_AXIAL", .19, 0)],
        "GAP_FILL": [line("GAP_FILL", .21, 0, .3)],
    }
    b = {
        "EXTERNAL_PERIMETER": [line("EXTERNAL_PERIMETER", 1.1, 15)],
        "INTERNAL_PERIMETER": [line("INTERNAL_PERIMETER", 1.12, 15)],
        "INFILL_RADIAL": [line("INFILL_RADIAL", 1.2, 15)],
    }
    support = {"SUPPORT": [line("SUPPORT", 2, 30)],
               "SUPPORT_INTERFACE": [line("SUPPORT_INTERFACE", 2.1, 30)]}
    by_role = {}
    for island in (a, b, support):
        for role, paths in island.items():
            by_role.setdefault(role, []).extend(paths)
    return LayerToolpaths(0, 31, by_role, island_paths=[a, b, support],
                          travel_region=box(.01, -1, 2.2, 37))


def deposition(code):
    current = total = 0.0
    for line in code.splitlines():
        if line.startswith("G92 E0"):
            current = 0.0
        if line.startswith("G1 "):
            match = re.search(r"\bE(-?[0-9.]+)", line)
            if match:
                e = float(match.group(1))
                if " X" in line:
                    total += e - current
                current = e
    return total


class IslandStageOrderTests(unittest.TestCase):
    def config(self, **kwargs):
        return SlicerConfig(bed_radius=30, retract_mm=0, z_hop_height=0,
                            cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0, **kwargs)

    def test_default_and_legacy_wall_order_preserved(self):
        self.assertIsNone(self.config().island_print_order)
        self.assertEqual(_island_role_sequence(None), ROLE_ORDER)
        self.assertEqual(_island_role_sequence(self.config()), ROLE_ORDER)
        legacy = _island_role_sequence(self.config(wall_order="inner-first"))
        self.assertEqual(legacy[:2], ["INTERNAL_PERIMETER", "EXTERNAL_PERIMETER"])
        self.assertEqual(legacy[2:], ROLE_ORDER[2:])
        self.assertNotIn(";ISLAND_PRINT_ORDER", emit_gcode([fixture()], self.config()))

    def test_every_model_island_has_its_own_three_stage_order(self):
        layer = fixture()
        c = self.config(island_print_order=("infill", "outer", "inner"))
        seq = _ordered_layer_paths_with_islands(layer, GcodeState(radial=31), c)
        by_island = {}
        for path, ident in seq:
            by_island.setdefault(ident, []).append(path.role)
        self.assertEqual(by_island[0], [
            "TOP_SKIN", "SOLID_INFILL", "INFILL_RADIAL", "INFILL_AXIAL",
            "GAP_FILL", "EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"])
        self.assertEqual(by_island[1], [
            "INFILL_RADIAL", "EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"])
        self.assertEqual(by_island[2], ["SUPPORT", "SUPPORT_INTERFACE"])
        self.assertEqual([key for key, _ in groupby([i for _, i in seq])],
                         list(dict.fromkeys(i for _, i in seq)))
        # Entry planning can reverse a line or rotate a closed seam, which
        # creates a new Toolpath object. Compare deposited roles, not ids.
        self.assertEqual(Counter(p.role for p, _ in seq),
                         Counter(p.role for paths in layer.paths_by_role.values() for p in paths))

    def test_arbitrary_stage_permutations_override_legacy_wall_order(self):
        for order, expected in (
            (("infill", "outer", "inner"), ["TOP_SKIN", "EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"]),
            (("outer", "infill", "inner"), ["EXTERNAL_PERIMETER", "TOP_SKIN", "INTERNAL_PERIMETER"]),
            (("inner", "outer", "infill"), ["INTERNAL_PERIMETER", "EXTERNAL_PERIMETER", "TOP_SKIN"]),
        ):
            with self.subTest(order=order):
                roles = _island_role_sequence(self.config(wall_order="inner-first", island_print_order=order))
                self.assertEqual([r for r in roles if r in expected], expected)

    def test_full_gcode_same_roles_and_filament(self):
        layer = fixture()
        a = emit_gcode([layer], self.config())
        b = emit_gcode([layer], self.config(island_print_order=("infill","outer","inner")))
        self.assertIn(";ISLAND_PRINT_ORDER:infill,outer,inner", b)
        tags = lambda code: Counter(line for line in code.splitlines() if line.startswith(";TYPE:"))
        self.assertEqual(tags(a), tags(b))
        self.assertAlmostEqual(deposition(a), deposition(b), places=4)

    def test_missing_stages_and_one_wall_are_safe(self):
        one = line("EXTERNAL_PERIMETER", .1, 0)
        layer = LayerToolpaths(0, 31, {one.role: [one]})
        a = _ordered_layer_paths_with_islands(layer, GcodeState(), self.config())
        b = _ordered_layer_paths_with_islands(layer, GcodeState(),
            self.config(island_print_order=("infill","outer","inner")))
        self.assertEqual([path.points for path, _ in a], [path.points for path, _ in b])

    def test_cli_profiles_and_invalid_orders(self):
        for profile in ("quality", "legacy-cli", "library"):
            args = build_parser().parse_args(["--profile", profile])
            self.assertIsNone(_profile_value(args, "island_print_order"))
            specified = build_parser().parse_args(
                ["--profile", profile, "--island-print-order", "infill, outer, inner"])
            self.assertEqual(_profile_value(specified, "island_print_order"),
                             ("infill","outer","inner"))
        for invalid in (("infill","outer"), ("infill","infill","inner"),
                        ("fill","outer","inner"), ("outer","inner","infill","bridge"),
                        "infill,outer,inner", []):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "island_print_order"):
                self.config(island_print_order=invalid).validate()


if __name__ == "__main__":
    unittest.main()
