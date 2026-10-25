from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _cooling_print_speed, _layer_fan_pwm, emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def layer(index, r, mm, *, role="EXTERNAL_PERIMETER"):
    return LayerToolpaths(index, r, {
        role: [Toolpath(role, r, [(0, 0), (0, mm)])]
    })


class LayerCoolingTests(unittest.TestCase):
    def test_disabled_matches_default_print_feed(self):
        c = SlicerConfig(feed_mm_per_min=1800)
        chosen, seconds = _cooling_print_speed(layer(0, 30, 3), c)
        self.assertEqual(chosen, 1800)
        self.assertAlmostEqual(seconds, 0.1)
        code = emit_gcode([layer(0, 30, 3)], c)
        self.assertNotIn(";COOLING_PRINT_SPEED_MM_MIN:", code)

    def test_short_print_slowed_to_minimum_time_in_physical_units(self):
        c = SlicerConfig(feed_mm_per_min=1800,
                         cooling_min_layer_time_s=2,
                         cooling_min_print_speed_mm_per_min=30,
                         bed_radius=29)
        speed, estimated = _cooling_print_speed(layer(0, 30, 3), c)
        self.assertAlmostEqual(speed, 90)
        self.assertAlmostEqual(estimated, 2)
        code = emit_gcode([layer(0, 30, 3)], c)
        self.assertIn(";COOLING_PRINT_SPEED_MM_MIN:90.0000", code)
        self.assertIn(";COOLING_PRINT_DURATION_ESTIMATE_S:2.0000", code)
        self.assertIn("G1 X18.000 Y0.000", code)
        self.assertIn(" F90.0", code)

    def test_minimum_print_speed_prevents_extreme_slowdown(self):
        c = SlicerConfig(feed_mm_per_min=1800,
                         cooling_min_layer_time_s=20,
                         cooling_min_print_speed_mm_per_min=300)
        speed, duration = _cooling_print_speed(layer(0, 30, 3), c)
        self.assertEqual(speed, 300)
        self.assertAlmostEqual(duration, .6)

    def test_long_layer_keeps_configured_speed(self):
        c = SlicerConfig(cooling_min_layer_time_s=2, feed_mm_per_min=1800)
        speed, _ = _cooling_print_speed(layer(0, 30, 150), c)
        self.assertEqual(speed, 1800)

    def test_arc_mm_not_rotary_degrees_drives_cooling(self):
        r = 60
        arc = Toolpath("INFILL_RADIAL", r, [(0, 0), (3/r, 0)])
        config = SlicerConfig(cooling_min_layer_time_s=1, feed_mm_per_min=1800,
                              cooling_min_print_speed_mm_per_min=30)
        linear = layer(0, r, 3)
        rotary = LayerToolpaths(0, r, {"INFILL_RADIAL": [arc]})
        for target in (linear, rotary):
            speed, duration = _cooling_print_speed(target, config)
            self.assertAlmostEqual(speed, 180)
            self.assertAlmostEqual(duration, 1)

    def test_only_layer_feedrate_changes_not_extrusion(self):
        p = layer(0, 30, 3)
        a = emit_gcode([p], SlicerConfig(bed_radius=29, retract_mm=0))
        b = emit_gcode([p], SlicerConfig(
            bed_radius=29, retract_mm=0, cooling_min_layer_time_s=2,
            cooling_min_print_speed_mm_per_min=30))
        e = lambda code: re.findall(r"(?m)^G1 X[^\n]*? E(-?[0-9.]+)", code)
        self.assertEqual(e(a), e(b))

    def test_fan_ramp_over_radial_nozzle_height(self):
        cfg = SlicerConfig(fan_speed=200, fan_start_height=1,
                           fan_ramp_height_mm=2, bed_radius=30, layer_height=.2)
        self.assertEqual(_layer_fan_pwm(.5, cfg), 0)
        self.assertEqual(_layer_fan_pwm(1, cfg), 0)
        self.assertEqual(_layer_fan_pwm(1.5, cfg), 50)
        self.assertEqual(_layer_fan_pwm(3, cfg), 200)
        shells = [layer(i, r, 3) for i, r in
                  enumerate([30.9, 31.4, 31.9, 32.9])]
        code = emit_gcode(shells, cfg)
        self.assertIn("M106 S50 ; part cooling ramp", code)
        self.assertIn("M106 S100 ; part cooling ramp", code)
        self.assertIn("M106 S200 ; part cooling on", code)
        self.assertEqual(code.count("M106 S200 ; part cooling on"), 1)

    def test_invalid_cooling_configuration_fails_fast(self):
        for kwargs in ({"fan_ramp_height_mm": -1},
                       {"fan_ramp_height_mm": math.nan},
                       {"cooling_min_layer_time_s": -1},
                       {"cooling_min_layer_time_s": math.inf},
                       {"cooling_min_print_speed_mm_per_min": 0},
                       {"cooling_min_print_speed_mm_per_min": math.nan},
                       {"cooling_min_layer_time_s": 1,
                        "cooling_min_print_speed_mm_per_min": 2000}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                SlicerConfig(**kwargs).validate()
