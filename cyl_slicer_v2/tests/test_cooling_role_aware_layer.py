from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _cooling_print_speed, emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def layer(roles: dict[str, float], *, radius=30.5):
    return LayerToolpaths(0, radius, {
        role: [Toolpath(role, radius, [(0, 0), (0, length)])]
        for role, length in roles.items()
    })


class RoleAwareLayerCoolingTests(unittest.TestCase):
    def test_role_cap_already_satisfies_layer_time(self):
        cfg = SlicerConfig(cooling_min_layer_time_s=8,
                           cooling_min_island_print_time_s=0,
                           feed_mm_per_min=2700,
                           bridge_print_feed_mm_per_min=900)
        # 120 mm of bridge takes 8 seconds at its own 900 mm/min role cap.
        chosen, estimate = _cooling_print_speed(layer({"BRIDGE": 120}), cfg)
        self.assertAlmostEqual(chosen, 2700)
        self.assertAlmostEqual(estimate, 8)

    def test_combined_ordinary_and_bridge_slow_only_as_needed(self):
        cfg = SlicerConfig(cooling_min_layer_time_s=8,
                           cooling_min_island_print_time_s=0,
                           feed_mm_per_min=2700,
                           bridge_print_feed_mm_per_min=900,
                           cooling_min_print_speed_mm_per_min=300)
        paths = layer({"BRIDGE": 45, "TOP_SKIN": 90})
        feed, seconds = _cooling_print_speed(paths, cfg)
        # Bridge already accounts for 3 s, so ordinary 90 mm needs 5 s.
        self.assertAlmostEqual(feed, 1080, places=4)
        self.assertAlmostEqual(seconds, 8, places=5)
        code = emit_gcode([paths], cfg)
        self.assertIn(";COOLING_PRINT_SPEED_MM_MIN:1080.0000", code)
        self.assertIn(";COOLING_PRINT_DURATION_ESTIMATE_S:8.0000", code)
        # A lower role-specific limit remains intact after global slowdown.
        extrude = [line for line in code.splitlines()
                   if line.startswith("G1 X") and " E" in line]
        f_values = [float(re.search(r" F([0-9.]+)", s).group(1)) for s in extrude]
        self.assertTrue(any(abs(f-900) < .1 for f in f_values))
        self.assertTrue(any(abs(f-1080) < .1 for f in f_values))

    def test_support_interface_cap_counts_even_when_cooling_disabled(self):
        cfg = SlicerConfig(cooling_min_layer_time_s=0,
                           support_interface_print_feed_mm_per_min=600)
        speed, seconds = _cooling_print_speed(
            layer({"SUPPORT_INTERFACE": 30}), cfg
        )
        self.assertEqual(speed, cfg.feed_mm_per_min)
        self.assertAlmostEqual(seconds, 3)

    def test_unreachable_target_still_respects_minimum_speed_and_no_dwell(self):
        cfg = SlicerConfig(cooling_min_layer_time_s=100,
                           cooling_min_island_print_time_s=0,
                           cooling_min_print_speed_mm_per_min=300)
        p = layer({"TOP_SKIN": 2})
        feed, duration = _cooling_print_speed(p, cfg)
        self.assertEqual(feed, 300)
        self.assertAlmostEqual(duration, .4)
        self.assertNotIn("\nG4", emit_gcode([p], cfg))

    def test_physical_arc_and_axial_roles_have_equal_time(self):
        cfg = SlicerConfig(cooling_min_layer_time_s=3,
                           bridge_print_feed_mm_per_min=900)
        r = 60
        a = layer({"BRIDGE": 45}, radius=r)
        b = LayerToolpaths(0, r, {"BRIDGE": [Toolpath(
            "BRIDGE", r, [(0, 0), (45/r, 0)]
        )]})
        self.assertEqual(_cooling_print_speed(a, cfg),
                         _cooling_print_speed(b, cfg))


if __name__ == "__main__":
    unittest.main()
