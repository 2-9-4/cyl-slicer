from __future__ import annotations

import math
import subprocess
import sys
import unittest
from pathlib import Path

from cyl_slicer_v2.cli import CLI_PROFILE_OVERRIDES, _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode


class HardwareProfileTests(unittest.TestCase):
    def test_original_cli_profile_is_explicit_not_the_default(self):
        args = build_parser().parse_args(["--profile", "legacy-cli"])
        self.assertEqual(args.profile, "legacy-cli")
        for key, old in CLI_PROFILE_OVERRIDES.items():
            self.assertEqual(_profile_value(args, key), old)

    def test_library_profile_matches_python_config_and_flags_override(self):
        plain = SlicerConfig()
        args = build_parser().parse_args(["--profile", "library"])
        for key in CLI_PROFILE_OVERRIDES:
            self.assertEqual(_profile_value(args, key), getattr(plain, key))
        override = build_parser().parse_args(["--profile", "library",
            "--wall-count", "5", "--build-cylinder-radius", "42",
            "--support-merge-gap", "7"])
        self.assertEqual(_profile_value(override, "wall_count"), 5)
        self.assertEqual(_profile_value(override, "bed_radius"), 42.0)
        self.assertEqual(_profile_value(override, "support_merge_gap"), 7.0)

    def test_physical_nozzle_is_separate_from_nominal_line_width(self):
        config = SlicerConfig(nozzle_diameter=0.6, nozzle_width=0.4, layer_height=0.2)
        config.validate()
        self.assertAlmostEqual(config.extrusion_per_mm,
            config.bead_area / (math.pi * (1.75 / 2) ** 2))
        gcode = emit_gcode([], config)
        self.assertIn(";LINE_WIDTH:0.4", gcode)
        self.assertIn(";NOZZLE_DIAMETER:0.6", gcode)
        self.assertIn(";NOZZLE_DIAMETER:unknown",
            emit_gcode([], SlicerConfig(nozzle_diameter=None)))

    def test_explicit_nozzle_rejects_thick_or_invalid_layer(self):
        for kw in ({"nozzle_diameter": 0.4, "layer_height": 0.5},
                   {"nozzle_diameter": 0.0}, {"nozzle_diameter": -1.0},
                   {"nozzle_diameter": math.nan}, {"nozzle_diameter": math.inf}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                SlicerConfig(**kw).validate()

    def test_round_bead_profile_and_explicit_nozzle_coexist(self):
        config = SlicerConfig(bead_model="rounded", layer_height=0.2,
                              nozzle_width=0.4, nozzle_diameter=0.4)
        config.validate()
        self.assertAlmostEqual(config.dense_line_spacing, 0.3570796326795)

    def test_cli_rejects_bad_hardware_before_slicing(self):
        root = Path(__file__).resolve().parents[2]
        cmd = [sys.executable, "-m", "cyl_slicer_v2.cli", "--fixture",
               "solid_box", "--layer-height", "0.5", "--nozzle-diameter", "0.4"]
        run = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("layer_height cannot exceed", run.stderr)
        self.assertNotIn("Timings:", run.stdout)
