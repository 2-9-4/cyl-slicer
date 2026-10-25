from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.cli import build_parser
from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import _role_print_settings, _surface_feed_for_marlin, emit_gcode, emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


def single(role, *, rotary=False):
    points = [(0, 0), (0.2, 0)] if rotary else [(0, 0), (0, 10)]
    return Toolpath(role, 30.1, points)


def output(path, config):
    return emit_toolpath(path, GcodeState(radial=30.3), config, 15, retract_after=False)


def last_extrusion(lines):
    commands = [line for line in lines if line.startswith("G1 X")]
    assert len(commands) == 1, commands
    e = float(re.search(r"\bE(-?\d+(?:\.\d+)?)", commands[0]).group(1))
    f = float(re.search(r"\bF(-?\d+(?:\.\d+)?)", commands[0]).group(1))
    return e, f


class SupportRoleFlowSpeedTests(unittest.TestCase):
    def test_defaults_leave_support_and_model_deposition_identical(self):
        c = SlicerConfig(layer_height=.2, feed_mm_per_min=1200)
        for role in ("SUPPORT", "SUPPORT_INTERFACE", "EXTERNAL_PERIMETER"):
            self.assertEqual(_role_print_settings(role, c), (1200, 1.0))
            e, f = last_extrusion(output(single(role), c))
            self.assertAlmostEqual(e, 10 * c.extrusion_per_mm, delta=.00001)
            self.assertAlmostEqual(f, 1200, delta=.01)

    def test_support_interface_independent_speed_and_flow(self):
        c = SlicerConfig(
            layer_height=.2, feed_mm_per_min=1200,
            support_print_feed_mm_per_min=900,
            support_interface_print_feed_mm_per_min=450,
            support_flow_ratio=.85,
            support_interface_flow_ratio=1.1,
        )
        for role, speed, flow in (
            ("SUPPORT",900,.85),
            ("SUPPORT_INTERFACE",450,1.1),
            ("SOLID_INFILL",1200,1),
            ("TOP_SKIN",1200,1),
        ):
            with self.subTest(role=role):
                e, f = last_extrusion(output(single(role), c))
                self.assertAlmostEqual(e, 10 * c.extrusion_per_mm * flow, delta=.00001)
                self.assertAlmostEqual(f, speed, delta=.01)

    def test_mixed_rotary_moves_keep_physical_speed_conversion(self):
        c = SlicerConfig(
            layer_height=.2, feed_mm_per_min=1800,
            support_interface_print_feed_mm_per_min=600,
            support_interface_flow_ratio=.75,
        )
        path = single("SUPPORT_INTERFACE", rotary=True)
        e, f = last_extrusion(output(path, c))
        actual_length = .2 * path.radius
        self.assertAlmostEqual(e, actual_length*c.extrusion_per_mm*.75, delta=.00001)
        expected = _surface_feed_for_marlin(0, .2, path.radius, 600)
        self.assertAlmostEqual(f, expected, delta=.06)

    def test_cooling_layer_speed_caps_support_overrides(self):
        c = SlicerConfig(
            feed_mm_per_min=300, support_print_feed_mm_per_min=900,
            support_interface_print_feed_mm_per_min=1000,
        )
        self.assertEqual(_role_print_settings("SUPPORT", c), (300, 1))
        self.assertEqual(_role_print_settings("SUPPORT_INTERFACE", c), (300, 1))
        self.assertEqual(_role_print_settings("TOP_SKIN", c), (300, 1))

    def test_gcode_metadata_and_model_extrusion_preserved(self):
        p = Toolpath("EXTERNAL_PERIMETER", 30.1, [(0, 0), (0, 5)])
        layers = [LayerToolpaths(0, 30.1, {"EXTERNAL_PERIMETER": [p]})]
        base = SlicerConfig(bed_radius=30, retract_mm=0)
        extra = SlicerConfig(
            bed_radius=30, retract_mm=0,
            support_print_feed_mm_per_min=600,
            support_interface_print_feed_mm_per_min=300,
            support_flow_ratio=.9,
            support_interface_flow_ratio=1.1,
        )
        a = emit_gcode(layers, base)
        b = emit_gcode(layers, extra)
        self.assertEqual(
            [line for line in a.splitlines() if line.startswith("G1 X")],
            [line for line in b.splitlines() if line.startswith("G1 X")],
        )
        self.assertIn(";SUPPORT_INTERFACE_FLOW_RATIO:1.1", b)
        self.assertIn(";SUPPORT_PRINT_FEED_MM_MIN:600", b)

    def test_cli_settings_and_invalid_parameters(self):
        args = build_parser().parse_args([
            "--support-print-feed", "700",
            "--support-interface-print-feed", "400",
            "--support-flow-ratio", ".9",
            "--support-interface-flow-ratio", "1.1",
        ])
        self.assertEqual(args.support_print_feed, 700)
        self.assertEqual(args.support_interface_print_feed, 400)
        self.assertEqual(args.support_flow_ratio, .9)
        self.assertEqual(args.support_interface_flow_ratio, 1.1)
        for bad in (
            {"support_print_feed_mm_per_min":0},
            {"support_interface_print_feed_mm_per_min":math.inf},
            {"support_interface_print_feed_mm_per_min":math.nan},
            {"support_flow_ratio":0},
            {"support_flow_ratio":math.nan},
            {"support_interface_flow_ratio":2.1},
        ):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                SlicerConfig(**bad).validate()


if __name__ == "__main__":
    unittest.main()
