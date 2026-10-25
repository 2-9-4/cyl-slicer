from dataclasses import replace
from contextlib import redirect_stderr
import io
import re
import unittest

from cyl_slicer_v2.calibration import RetractionSweep, sweep_from_args
from cyl_slicer_v2.cli import build_parser
from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


class CalibrationTests(unittest.TestCase):
    def test_feed_settings_must_be_finite_and_positive(self):
        for field in ("feed_mm_per_min", "travel_mm_per_min"):
            for value in (0, -1, float("nan"), float("inf")):
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    SlicerConfig(**{field: value}).validate()

    def test_band_boundaries_base_and_final_hold(self):
        sweep = RetractionSweep("retract-mm", (0, 2, 4))
        for height, expected in [(0.2, None), (1.8, None), (2, 0),
                                 (4.8, 0), (4.999999999999, 1), (5, 1),
                                 (8, 2), (100, 2)]:
            with self.subTest(height=height):
                self.assertEqual(sweep.band_at(height), expected)
        self.assertEqual(sweep.overrides(1), {})
        self.assertEqual(sweep.overrides(2), {"retract_mm": 0})

    def test_speed_units_and_invalid_inputs(self):
        self.assertEqual(RetractionSweep("retract-speed", (25,)).overrides(2),
                         {"retract_speed_mm_per_min": 1500})
        for kwargs in [dict(values=()), dict(values=(-1,)), dict(values=(float("nan"),)),
                       dict(values=(float("inf"),)), dict(band_height=0),
                       dict(base_height=-1), dict(parameter="temperature"),
                       dict(parameter="retract-speed", values=(0,))]:
            params = dict(parameter="retract-mm", values=(1, 2))
            params.update(kwargs)
            with self.subTest(params=params), self.assertRaises(ValueError):
                RetractionSweep(**params)

    def test_cli_requires_explicit_values_and_parameter(self):
        parser = build_parser()
        args = parser.parse_args(["--calibration", "retract-mm", "--calibration-values", "0,1,2"])
        self.assertEqual(sweep_from_args(args, parser).values, (0, 1, 2))
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            sweep_from_args(parser.parse_args(["--calibration", "retract-mm"]), parser)

    def test_gcode_balances_retraction_across_resets_and_zero_band(self):
        sweep = RetractionSweep("retract-mm", (3, 0, 1), base_height=0, band_height=1)
        layers = [LayerToolpaths(i, 30.2+i, {"SUPPORT": [
            Toolpath("SUPPORT", 30.2+i, [(0, 0), (0, 6)]),
            Toolpath("SUPPORT", 30.2+i, [(0, 20), (0, 26)]),
        ]}) for i in range(3)]
        config = SlicerConfig(bed_radius=30, retraction_sweep=sweep, z_hop_height=0,
                              max_axial_segment=20, travel_mm_per_min=4200,
                              retract_speed_mm_per_min=1500)
        gcode = emit_gcode(layers, config)
        e, pending, layer = 0.0, 0.0, -1
        counts = [0, 0, 0]
        for line in gcode.splitlines():
            if line.startswith(";LAYER:"):
                layer = int(line.split(":")[1])
            if line.startswith("G92 E0"):
                e = 0
                continue
            match = re.search(r"\bE(-?[\d.]+)", line)
            if not line.startswith("G1 ") or not match:
                continue
            new_e = float(match[1])
            if line.endswith("; retract"):
                self.assertEqual(pending, 0)
                pending = e - new_e
                self.assertAlmostEqual(pending, (3, 0, 1)[layer], places=4)
                self.assertIn("F1500.0", line)
                counts[layer] += 1
            elif line.endswith("; unretract"):
                self.assertGreater(pending, 0)
                self.assertAlmostEqual(new_e - e, pending, places=4)
                self.assertIn("F1500.0", line)
                pending = 0
            else:
                self.assertEqual(pending, 0, "printing while still retracted")
            e = new_e
        self.assertEqual(counts, [2, 0, 2])
        # Sweeping retraction must not alter total deposited filament.
        baseline = emit_gcode(layers, replace(config, retraction_sweep=None, retract_mm=0))
        def deposited(code):
            previous = total = 0
            for line in code.splitlines():
                if line.startswith("G92 E0"):
                    previous = 0
                m = re.search(r"\bE(-?[\d.]+)", line)
                if line.startswith("G1 ") and m:
                    current = float(m[1])
                    if " X" in line or " Y" in line:
                        total += current - previous
                    previous = current
            return total
        self.assertAlmostEqual(deposited(gcode), deposited(baseline), places=4)

    def test_speed_sweep_uses_old_speed_to_repay_pending_retraction(self):
        layers = [LayerToolpaths(i, 30.2+i, {"SUPPORT": [
            Toolpath("SUPPORT", 30.2+i, [(0, 0), (0, 6)])]}) for i in range(2)]
        code = emit_gcode(layers, SlicerConfig(
            bed_radius=30, z_hop_height=0,
            retraction_sweep=RetractionSweep("retract-speed", (20, 40), 0, 1)))
        moves = [line for line in code.splitlines() if line.endswith(("; retract", "; unretract"))]
        self.assertIn("F1200.0 ; retract", moves[0])
        self.assertIn("F1200.0 ; unretract", moves[1])
        self.assertIn("F2400.0 ; retract", moves[2])
