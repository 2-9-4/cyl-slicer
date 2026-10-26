from __future__ import annotations

import re
import subprocess
import sys
import unittest

from cyl_slicer_v2.calibration import RetractionSweep, sweep_from_args
from cyl_slicer_v2.cli import build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def layers():
    return [LayerToolpaths(i, radius, {
        "EXTERNAL_PERIMETER": [
            Toolpath("EXTERNAL_PERIMETER", radius,
                     [(0, 0), (.2, 0), (.2, 4)])
        ]})
        for i, radius in enumerate((30.42, 31.42, 32.42, 33.42))]


def motion(gcode):
    return [line for line in gcode.splitlines()
            if line.startswith(("G0 ", "G1 "))]


class LinearAdvanceBandCalibrationTests(unittest.TestCase):
    def test_band_k_and_base_reset_order(self):
        sweep = RetractionSweep("linear-advance-k", (.04, .08, .12),
                                base_height=1., band_height=1.)
        config = SlicerConfig(
            linear_advance_firmware="marlin-m900", linear_advance_k=0.,
            retraction_sweep=sweep, bed_radius=30,
            retract_mm=0, cooling_min_layer_time_s=0,
            cooling_min_island_print_time_s=0,
        )
        g = emit_gcode(layers(), config)
        m900 = [line for line in g.splitlines() if line.startswith("M900 K")]
        self.assertEqual([re.match(r"M900 K([^ ]+)", s).group(1) for s in m900],
                         ["0", "0.04", "0.08", "0.12", "0"])
        self.assertLess(g.index("M109 "), g.index("\nM900 K0 ; enable"))
        for i, k in ((1, "0.04"), (2, "0.08"), (3, "0.12")):
            section = g.split(f";LAYER:{i}\n", 1)[1].split(";LAYER:", 1)[0]
            self.assertIn(f"M900 K{k} ; calibrated linear advance radial band", section)
        self.assertLess(g.index("M900 K0 ; clear"), g.index("M84 ;"))

    def test_command_only_sweep_does_not_affect_e_or_feed_or_rotary(self):
        active = SlicerConfig(
            linear_advance_firmware="marlin-m900", linear_advance_k=0.,
            retraction_sweep=RetractionSweep(
                "linear-advance-k", (0., .05, .12), 1, 1),
            bed_radius=30, retract_mm=0,
            cooling_min_layer_time_s=0, cooling_min_island_print_time_s=0,
        )
        baseline = SlicerConfig(
            linear_advance_firmware="marlin-m900", linear_advance_k=0.,
            bed_radius=30, retract_mm=0,
            cooling_min_layer_time_s=0, cooling_min_island_print_time_s=0,
        )
        self.assertEqual(motion(emit_gcode(layers(), active)),
                         motion(emit_gcode(layers(), baseline)))

    def test_sweep_requires_declared_firmware_and_base_k(self):
        sweep = RetractionSweep("linear-advance-k", (.02, .06))
        for extra in ({}, {"linear_advance_k": .0},
                      {"linear_advance_firmware": "marlin-m900"}):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                SlicerConfig(retraction_sweep=sweep, **extra).validate()
        argv = [sys.executable, "-m", "cyl_slicer_v2.cli",
                "--calibration", "linear-advance-k",
                "--calibration-values", "0,.05,.1",
                "--fixture", "tube"]
        result = subprocess.run(argv, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires --marlin-linear-advance", result.stderr)

    def test_k_values_and_height_bands_validated(self):
        self.assertEqual(RetractionSweep("linear-advance-k", (0., .05, .2))
                         .overrides(5), {"linear_advance_k": .05})
        for bad in ((-1,), (10.01,), (float("inf"),),
                    (float("nan"),)):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                RetractionSweep("linear-advance-k", bad)

    def test_cli_accepts_sweep_parameters(self):
        parser = build_parser()
        args = parser.parse_args(["--marlin-linear-advance",
                                  "--calibration", "linear-advance-k",
                                  "--calibration-values", "0,.1,.2"])
        sweep = sweep_from_args(args, parser)
        self.assertEqual(sweep.parameter, "linear-advance-k")
        self.assertEqual(sweep.values, (0., .1, .2))


if __name__ == "__main__":
    unittest.main()
