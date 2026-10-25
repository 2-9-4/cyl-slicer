from __future__ import annotations

import math
import re
import subprocess
import sys
import unittest

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def sample():
    return [LayerToolpaths(0, 30.1, {
        "EXTERNAL_PERIMETER": [Toolpath(
            "EXTERNAL_PERIMETER", 30.1, [(0, 0), (.5, 0), (.5, 1)]
        )]
    })]


class MarlinLinearAdvanceTests(unittest.TestCase):
    def test_default_has_no_m900_commands(self):
        result = emit_gcode(sample(), SlicerConfig(bed_radius=30))
        self.assertNotIn("M900 ", result)
        self.assertIn(";LINEAR_ADVANCE:disabled", result)

    def test_explicit_firmware_and_k_are_scoped_to_print(self):
        cfg = SlicerConfig(bed_radius=30, linear_advance_k=0.17,
                           linear_advance_firmware="marlin-m900", material_profile="petg")
        g = emit_gcode(sample(), cfg)
        self.assertIn(";MATERIAL_PROFILE:petg", g)
        self.assertIn(";LINEAR_ADVANCE:Marlin M900 K0.17", g)
        self.assertEqual(len(re.findall(r"(?m)^M900 K", g)), 2)
        self.assertLess(g.index("M109 "), g.index("\nM900 K0.17"))
        self.assertLess(g.index("\nM900 K0.17"), g.index(";LAYER:0"))
        self.assertLess(g.index("; end lift"), g.index("M900 K0 ; clear"))
        self.assertLess(g.index("M900 K0 ; clear"), g.index("M84 ;"))

    def test_zero_is_an_explicit_k_not_an_unconfigured_value(self):
        cfg = SlicerConfig(linear_advance_k=0, linear_advance_firmware="marlin-m900")
        g = emit_gcode([], cfg)
        self.assertEqual(len(re.findall(r"(?m)^M900 K0(?:\s|$)", g)), 2)

    def test_k_does_not_alter_physical_e_or_surface_speed(self):
        plain = emit_gcode(sample(), SlicerConfig(bed_radius=30, retract_mm=0))
        active = emit_gcode(sample(), SlicerConfig(
            bed_radius=30, retract_mm=0, linear_advance_firmware="marlin-m900",
            linear_advance_k=.05))
        def motion(code):
            return [(x[0],x[1],x[2],x[3],x[4],x[5]) for x in parse_extrusion_segments(code)]
        self.assertEqual(motion(plain), motion(active))
        self.assertEqual([s for s in plain.splitlines() if s.startswith("G1 X")],
                         [s for s in active.splitlines() if s.startswith("G1 X")])

    def test_rejects_missing_firmware_capability_and_invalid_material(self):
        examples = [
            {"linear_advance_k": .2},
            {"linear_advance_firmware": "marlin-m900"},
            {"linear_advance_k": -0.1, "linear_advance_firmware": "marlin-m900"},
            {"linear_advance_k": math.inf, "linear_advance_firmware": "marlin-m900"},
            {"linear_advance_k": math.nan, "linear_advance_firmware": "marlin-m900"},
            {"linear_advance_k": 11, "linear_advance_firmware": "marlin-m900"},
            {"linear_advance_k": .2, "linear_advance_firmware": "unknown"},
            {"material_profile": "fantasy"},
        ]
        for args in examples:
            with self.subTest(args=args), self.assertRaises(ValueError):
                SlicerConfig(**args).validate()

    def test_cli_does_not_allow_unsupported_k(self):
        args = [sys.executable, "-m", "cyl_slicer_v2.cli", "--fixture", "tube",
                "--linear-advance-k", "0.15"]
        result = subprocess.run(args, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires explicit Marlin M900", result.stderr)
