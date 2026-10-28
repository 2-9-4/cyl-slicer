from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _island_cooling_print_speed, emit_gcode, emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


ROOT = Path(__file__).resolve().parents[2]


def print_moves(lines):
    return [line for line in lines if line.startswith("G1 X")]


def fields(line):
    return {key: float(value) for key, value in
            re.findall(r"\b([XYEF])(-?[\d.]+)", line)}


class FlowLimitTests(unittest.TestCase):
    def config(self, **kwargs):
        return SlicerConfig(layer_height=.28, nozzle_width=.48,
                            gap_fill_min_width_mm=.28,
                            feed_mm_per_min=4500, min_path_length=0,
                            cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0,
                            retract_mm=0, z_hop_height=0, **kwargs)

    def test_actual_flow_capped_with_widths_role_flow_calibration_and_ramps(self):
        for bead_model in ("rounded", "rectangular"):
            for role in ("GAP_FILL", "SUPPORT_INTERFACE", "SOLID_INFILL"):
                for radius in (30., 70.):
                    with self.subTest(bead_model=bead_model, role=role, radius=radius):
                        base = self.config(bead_model=bead_model,
                            support_interface_flow_ratio=1.3,
                            extrusion_feed_ramp_mm=1.5)
                        base = replace(base, extrusion_per_mm=base.extrusion_per_mm*1.2)
                        limited = replace(base, max_volumetric_flow_mm3_s=6)
                        path = Toolpath(role, radius,
                            [(0, 0), (0, 10), (10/radius, 10), (20/radius, 20)],
                            segment_widths=[.30, .40, .48])
                        a = print_moves(emit_toolpath(path, GcodeState(radial=radius), base, 0))
                        b = print_moves(emit_toolpath(path, GcodeState(radial=radius), limited, 0))
                        self.assertEqual([s.split(" F")[0] for s in a],
                                         [s.split(" F")[0] for s in b])
                        self.assertNotEqual(a, b)
                        previous = dict(X=0., Y=0., E=0.)
                        for line in b:
                            move = fields(line)
                            planner_length = math.hypot(move["X"]-previous["X"],
                                                        move["Y"]-previous["Y"])
                            seconds = planner_length*60/move["F"]
                            flow = (move["E"]-previous["E"])*math.pi*(1.75/2)**2/seconds
                            # Account for output coordinate/E/F decimal rounding.
                            self.assertLessEqual(flow, 6.003)
                            previous = move

    def test_partial_layer_uses_its_actual_bead_and_calibrated_e(self):
        base = self.config()
        full = replace(base, max_volumetric_flow_mm3_s=6)
        partial = full.for_deposited_thickness(.14)
        path = Toolpath("SOLID_INFILL", 31, [(0, 0), (0, 10)])
        for config in (full, partial):
            line, = print_moves(emit_toolpath(path, GcodeState(radial=31), config, 0))
            move = fields(line)
            flow = move["E"]*math.pi*(1.75/2)**2/(10*60/move["F"])
            self.assertLessEqual(flow, 6.001)
            self.assertAlmostEqual(move["E"], 10*config.extrusion_per_mm, delta=.00001)
            if config is partial:
                self.assertAlmostEqual(move["F"], 4500)

    def test_output_rounding_cannot_cause_flow_spikes_on_tiny_moves(self):
        base = self.config()
        limited = replace(base, max_volumetric_flow_mm3_s=6)
        # Includes a corner that rounds to a positive E-only move.
        path = Toolpath("SOLID_INFILL", 31,
            [(0, 0), (0, 10.0011), (.000005, 10.0012), (.1, 10.5)],
            segment_widths=[.48, .48, .48])
        a = print_moves(emit_toolpath(path, GcodeState(radial=31), base, 0))
        b = print_moves(emit_toolpath(path, GcodeState(radial=31), limited, 0))
        self.assertEqual([s.split(" F")[0] for s in a],
                         [s.split(" F")[0] for s in b])
        previous = dict(X=0., Y=0., E=0.)
        saw_rounded_e_only = False
        for line in b:
            move = fields(line)
            delta_e = move["E"]-previous["E"]
            distance = math.hypot(move["X"]-previous["X"], move["Y"]-previous["Y"])
            if delta_e > 0:
                saw_rounded_e_only |= distance == 0
                seconds = (distance or delta_e)*60/move["F"]
                self.assertLessEqual(delta_e*math.pi*(1.75/2)**2/seconds, 6+1e-10)
            previous = move
        self.assertTrue(saw_rounded_e_only)

    def test_unrepresentably_slow_flow_is_rejected_before_writing_zero_feed(self):
        config = self.config(max_volumetric_flow_mm3_s=1e-8)
        path = Toolpath("SOLID_INFILL", 31, [(0, 0), (0, 10)])
        with self.assertRaisesRegex(ValueError, "below G-code precision"):
            emit_toolpath(path, GcodeState(radial=31), config, 0)

    def test_cooling_counts_actual_width_role_flow_and_volume_ceiling(self):
        config = self.config(max_volumetric_flow_mm3_s=1,
                             support_interface_flow_ratio=1.3)
        config = replace(config, cooling_min_island_print_time_s=2,
                         extrusion_per_mm=config.extrusion_per_mm*1.2)
        path = Toolpath("SUPPORT_INTERFACE", 31,
                        [(0, 0), (0, 50), (0, 100)], segment_widths=[.30, .48])
        feed, duration = _island_cooling_print_speed([path], config)
        area = sum((w-.28)*.28+math.pi*.28**2/4 for w in (.30, .48))
        expected = area*50*1.3*1.2/1
        self.assertEqual(feed, config.feed_mm_per_min)
        self.assertAlmostEqual(duration, expected)
        code = emit_gcode([LayerToolpaths(0, 31, {path.role: [path]})], config)
        self.assertNotIn("\nG4", code)

    def test_new_options_resolve_from_config_and_validate_in_all_profiles(self):
        for profile in ("quality", "legacy-cli", "library"):
            plain = build_parser().parse_args(["--profile", profile])
            explicit = build_parser().parse_args(["--profile", profile,
                "--max-volumetric-flow", "6", "--print-accel", "1200", "--travel-accel", "2000"])
            for key, value in (("max_volumetric_flow_mm3_s", 6),
                               ("print_accel", 1200), ("travel_accel", 2000)):
                self.assertEqual(_profile_value(plain, key), getattr(SlicerConfig(), key))
                self.assertEqual(_profile_value(explicit, key), value)
        for key in ("max_volumetric_flow_mm3_s", "print_accel", "travel_accel"):
            for value in (-1, math.nan, math.inf):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    replace(self.config(), **{key: value}).validate()

    def test_acceleration_is_explicit_and_does_not_change_motion_or_filament(self):
        path = Toolpath("EXTERNAL_PERIMETER", 31, [(0, 0), (.1, 10)])
        layer = LayerToolpaths(0, 31, {path.role: [path]})
        config = self.config()
        a = emit_gcode([layer], config)
        b = emit_gcode([layer], replace(config, print_accel=1200, travel_accel=2000))
        command = "M204 P1200 T2000 ; native planner acceleration, firmware axis limits apply"
        self.assertNotIn("\nM204", a)
        self.assertEqual(a, b.replace(command+"\n", ""))
        self.assertLess(b.index(command), b.index(";LAYER:"))
        self.assertNotRegex(b, r"(?m)^M(?:201|203|205|500)\b")
        for key, letter in (("print_accel", "P"), ("travel_accel", "T")):
            code = emit_gcode([], replace(config, **{key: 1234}))
            self.assertIn(f"\nM204 {letter}1234 ;", code)


class FastScriptTests(unittest.TestCase):
    def test_script_from_another_directory_accepts_relative_paths_and_overrides(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = ROOT / "offset_cube.stl"
            (Path(tmp)/"input with spaces.stl").write_bytes(source.read_bytes())
            command = ["bash", str(ROOT/"slice_fast.sh"), "input with spaces.stl",
                       "--out", "fast output.gcode", "--max-layers", "1",
                       "--slice-workers", "1", "--toolpath-workers", "1"]
            subprocess.run(command, cwd=tmp, check=True, capture_output=True, text=True)
            code = (Path(tmp)/"fast output.gcode").read_text()
            for expected in (";LAYER_HEIGHT:0.28", ";LINE_WIDTH:0.48",
                             ";MAX_VOLUMETRIC_FLOW_MM3_S:6", "M204 P1200 T2000",
                             ";PRINT_FEED:4500", ";TRAVEL_FEED:12000"):
                self.assertIn(expected, code)
            self.assertTrue(print_moves(code.splitlines()))
            self.assertNotIn("\nM900", code)
            subprocess.run([*command, "--layer-height", "0.24", "--print-feed", "2400",
                            "--max-volumetric-flow", "4", "--print-accel", "0",
                            "--travel-accel", "0"], cwd=tmp,
                           check=True, capture_output=True, text=True)
            overridden = (Path(tmp)/"fast output.gcode").read_text()
            for expected in (";LAYER_HEIGHT:0.24", ";PRINT_FEED:2400",
                             ";MAX_VOLUMETRIC_FLOW_MM3_S:4"):
                self.assertIn(expected, overridden)
            self.assertNotIn("\nM204", overridden)


if __name__ == "__main__":
    unittest.main()
