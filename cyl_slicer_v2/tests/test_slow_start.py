from __future__ import annotations

from dataclasses import replace
import math
import re
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _surface_feed_for_marlin, emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def layer(index, radius, role="SOLID_INFILL", *, widths=None):
    path = Toolpath(role, radius, [(.1, 1), (.2, 6), (.3, 11)], segment_widths=widths)
    return LayerToolpaths(index, radius, {role: [path]})


def layer_blocks(code):
    return re.split(r"(?m)^;LAYER:[^\n]*\n", code)[1:]


def print_feed(block):
    line = next(line for line in block.splitlines() if line.startswith("G1 X"))
    return float(re.search(r"\bF([\d.]+)", line).group(1))


class SlowStartTests(unittest.TestCase):
    def config(self, **kwargs):
        return SlicerConfig(feed_mm_per_min=4500, travel_mm_per_min=12000,
                            retract_mm=0, cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0, **kwargs)

    def test_first_ten_nonempty_layers_include_negative_support_and_restore_speeds(self):
        layers = [layer(-20+i, 30.1+i*.16, "SUPPORT" if i < 4 else "SOLID_INFILL")
                  for i in range(12)]
        # Empty output layers must not consume the ten printed-layer allowance.
        layers.insert(3, LayerToolpaths(-18, 30.5, {}))
        config = self.config(slow_start_layers=10)
        code = emit_gcode(layers, config)
        printed = 0
        for item, block in zip(layers, layer_blocks(code)):
            if not item.paths_by_role:
                self.assertNotIn(";SLOW_START_LAYER:", block)
                continue
            printed += 1
            slow = printed <= 10
            expected = _surface_feed_for_marlin(5, .1, item.radius, 1200 if slow else 4500)
            self.assertAlmostEqual(print_feed(block), expected, delta=.051)
            radial_move = next(line for line in block.splitlines() if line.startswith("G0 Z"))
            self.assertIn(f" F{4800 if slow else 12000}.0", radial_move)
            if slow:
                self.assertIn(f";SLOW_START_LAYER:{printed}/10", block)
            else:
                self.assertNotIn(";SLOW_START_LAYER:", block)
        self.assertEqual(code.count(";SLOW_START_LAYER:"), 10)

    def test_slow_start_preserves_all_coordinates_e_widths_retraction_and_firmware_commands(self):
        layers = [layer(-3, 31, "GAP_FILL", widths=[.24, .38]), layer(0, 31.16)]
        base = replace(self.config(), retract_mm=1.4, print_accel=1200, travel_accel=2000)
        a = emit_gcode(layers, base)
        b = emit_gcode(layers, replace(base, slow_start_layers=1))

        def geometry_and_commands(code):
            return [re.sub(r" F[\d.]+", "", line) for line in code.splitlines()
                    if line and not line.startswith(";")]

        self.assertEqual(geometry_and_commands(a), geometry_and_commands(b))
        self.assertEqual(a.count("; retract"), b.count("; retract"))
        self.assertEqual(a.count("; unretract"), b.count("; unretract"))
        self.assertRegex(b, r"M204 P[\d.]+ T[\d.]+ ; layer surface acceleration")

    def test_cooling_uses_slow_start_ceiling_and_role_and_flow_limits_still_apply(self):
        layers = [layer(0, 31, "SUPPORT_INTERFACE"), layer(1, 31.16)]
        base = self.config(slow_start_layers=1, support_interface_print_feed_mm_per_min=600,
                           max_volumetric_flow_mm3_s=1)
        config = replace(base, cooling_min_layer_time_s=.1, cooling_min_island_print_time_s=.1)
        code = emit_gcode(layers, config)
        blocks = layer_blocks(code)
        self.assertIn(";COOLING_PRINT_SPEED_MM_MIN:1200.0000", blocks[0])
        self.assertIn(";COOLING_PRINT_SPEED_MM_MIN:4500.0000", blocks[1])
        for item, block, role_cap in zip(layers, blocks, (600, 4500)):
            flow_cap = 60/config.bead_area
            expected = _surface_feed_for_marlin(5, .1, item.radius, min(role_cap, flow_cap))
            self.assertLessEqual(print_feed(block), expected+.051)
        self.assertNotIn("\nG4", code)

    def test_above_normal_ceiling_never_speeds_up_and_disabled_is_identical(self):
        layers = [layer(0, 31)]
        base = self.config()
        code = emit_gcode(layers, base)
        high = replace(base, slow_start_layers=10,
                       slow_start_print_feed_mm_per_min=9000,
                       slow_start_travel_feed_mm_per_min=24000)
        enabled = emit_gcode(layers, high)
        self.assertEqual(print_feed(layer_blocks(code)[0]), print_feed(layer_blocks(enabled)[0]))
        self.assertIn(" F12000.0", next(line for line in layer_blocks(enabled)[0].splitlines()
                                       if line.startswith("G0 Z")))
        self.assertEqual(emit_gcode(layers, replace(high, slow_start_layers=0)), code)

    def test_profiles_overrides_and_invalid_values(self):
        for profile in ("quality", "legacy-cli", "library"):
            plain = build_parser().parse_args(["--profile", profile])
            explicit = build_parser().parse_args(["--profile", profile,
                "--slow-start-layers", "10", "--slow-start-print-feed", "900",
                "--slow-start-travel-feed", "3600"])
            for key, expected in (("slow_start_layers", 10),
                                  ("slow_start_print_feed_mm_per_min", 900),
                                  ("slow_start_travel_feed_mm_per_min", 3600)):
                self.assertEqual(_profile_value(plain, key), getattr(SlicerConfig(), key))
                self.assertEqual(_profile_value(explicit, key), expected)
        for key, invalid in (("slow_start_layers", (-1, 1.5, True, math.nan, math.inf)),
                             ("slow_start_print_feed_mm_per_min", (-1, 0, math.nan, math.inf)),
                             ("slow_start_travel_feed_mm_per_min", (-1, 0, math.nan, math.inf))):
            for value in invalid:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    replace(self.config(), **{key: value}).validate()


if __name__ == "__main__":
    unittest.main()
