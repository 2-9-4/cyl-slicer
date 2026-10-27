"""Replay Y commands with Marlin 2.0.7's float32 workspace/planner semantics.

References: Marlin/src/gcode/geometry/G92.cpp, module/motion.h and
module/planner.cpp in ../printer/Mizar-S/Marlin2.0.7. The configured rotary
scale is 52.7777778 steps/degree. This is a numeric model, not a firmware build.
The reduced model has no leveling, home offset, or Y software endstops, matching
the rotary configuration. G92.9 requires POWER_LOSS_RECOVERY in that firmware.
"""
from __future__ import annotations

import math
import re
import unittest
from io import StringIO

import numpy as np

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import MAX_MARLIN_Y_DEGREES, emit_gcode, emit_toolpath, parse_extrusion_segments
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


class MarlinRotaryModel:
    def __init__(self):
        self.native = np.float32(0)
        self.shift = np.float32(0)
        self.steps_per_degree = np.float32(52.7777778)
        self.steps = 0
        self.physical_steps = 0
        self.max_native = 0.0
        self.max_logical = 0.0
        self.max_steps = 0
        self.move_steps: list[int] = []
        self.last_command = ""
        self.native_resets = 0
        self.relative = False

    def _planner_steps(self, native):
        # Marlin uses LROUND(float32(position * axis_steps_per_mm)).
        value = float(np.float32(native * self.steps_per_degree))
        steps = math.floor(value + 0.5) if value >= 0 else math.ceil(value - 0.5)
        if not -(2**31) <= steps < 2**31:
            raise OverflowError("Marlin signed 32-bit planner position overflow")
        return steps

    def feed(self, line):
        tokens = line.split(";", 1)[0].split()
        if not tokens:
            return
        command = tokens[0]
        if command in {"G90", "G91"}:
            self.relative = command == "G91"
        y = next((np.float32(token[1:]) for token in tokens[1:] if token.startswith("Y")), None)
        if y is not None:
            if command == "G92":
                # G92.cpp: v=LOGICAL_TO_NATIVE(l), d=v-current_position;
                # position_shift+=d. It does not change native Y or steps.
                native_target = np.float32(y - self.shift)
                delta = np.float32(native_target - self.native)
                self.shift = np.float32(self.shift + delta)
            elif command == "G92.9":
                if self.last_command != "M400":
                    raise AssertionError("Native rotary reset must wait for queued moves")
                self.native = y
                self.steps = self._planner_steps(y)
                self.native_resets += 1
            elif command in {"G0", "G1"}:
                self.max_logical = max(self.max_logical, abs(float(y)))
                self.native = np.float32(self.native + y if self.relative else y - self.shift)
                new_steps = self._planner_steps(self.native)
                delta = new_steps - self.steps
                self.physical_steps += delta
                self.move_steps.append(delta)
                self.steps = new_steps
        self.max_native = max(self.max_native, abs(float(self.native)))
        self.max_steps = max(self.max_steps, abs(self.steps))
        self.last_command = command

    def replay(self, code):
        for line in StringIO(code):
            self.feed(line)
        return self


class MarlinRotaryPrecisionTests(unittest.TestCase):
    def test_compact_multilayer_print_keeps_one_rotary_coordinate_frame(self):
        # The path revisits the same physical arc on every layer. Repeated
        # G92.9/G92 Y rebases add a firmware dependency without bounding
        # anything further, and can introduce a layer-to-layer phase jump.
        layers = [LayerToolpaths(i, 30.0, {"EXTERNAL_PERIMETER": [
            Toolpath("EXTERNAL_PERIMETER", 30.0, [
                (math.radians(60), 0.0), (math.radians(120), 0.0),
            ])
        ]}) for i in range(100)]
        code = emit_gcode(layers, SlicerConfig(retract_mm=0, z_hop_height=0))
        marlin = MarlinRotaryModel().replay(code)

        self.assertEqual(code.count("G92.9 Y"), 1)
        self.assertEqual(code.count("G92 E0 Y0"), 1)
        self.assertEqual(code.count("G92 E0 ; reset extrusion"), len(layers))
        self.assertLess(marlin.max_native, 181)
        self.assertAlmostEqual(marlin.physical_steps / float(marlin.steps_per_degree), 120, delta=0.02)

    def test_generated_program_uses_absolute_rotary_moves(self):
        layer = LayerToolpaths(0, 30.0, {"SUPPORT": [
            Toolpath("SUPPORT", 30.0, [(0.0, 0.0), (math.radians(10), 0.0)]),
        ]})
        marlin = MarlinRotaryModel().replay("G91\n")
        marlin.replay(emit_gcode([layer], SlicerConfig(retract_mm=0, z_hop_height=0)))
        self.assertFalse(marlin.relative)
        self.assertAlmostEqual(float(marlin.native), 10.0)

    def test_plain_g92_leaves_native_precision_problem_in_place(self):
        marlin = MarlinRotaryModel().replay("G0 Y1000000\nG92 Y0\nG1 Y0.026\n")
        self.assertEqual(float(marlin.native), 1_000_000)
        self.assertEqual(marlin.move_steps[-1], 0, "Small move is lost despite G92 Y0")
        self.assertEqual(np.float32(80_000_000) + np.float32(0.026), np.float32(80_000_000))
        with self.assertRaises(OverflowError):
            MarlinRotaryModel().feed("G0 Y80000000")

    def test_long_print_bounds_native_and_logical_y_and_preserves_turns(self):
        # Both many layers and many full turns in a single layer. Neither layer
        # G92s alone nor a layer-only native reset suffices for this program.
        turns_per_layer, count = 32, 7000  # 80,640,000 degrees of physical rotation.
        layers = [LayerToolpaths(i, 30.0, {"EXTERNAL_PERIMETER": [
            Toolpath("EXTERNAL_PERIMETER", 30.0, [(0.0, 0.0), (2 * math.pi, 0.0)])
        ] * turns_per_layer}) for i in range(count)]
        code = emit_gcode(layers, SlicerConfig(retract_mm=0, z_hop_height=0))
        marlin = MarlinRotaryModel().replay(code)

        self.assertLessEqual(marlin.max_native, MAX_MARLIN_Y_DEGREES)
        self.assertLessEqual(marlin.max_logical, MAX_MARLIN_Y_DEGREES)
        self.assertLess(marlin.max_steps, 2**31)
        self.assertGreater(marlin.native_resets, count + 1, "Must reset inside large layers too")
        self.assertEqual(marlin.physical_steps, count * turns_per_layer * 19000)
        self.assertGreater(marlin.physical_steps / float(marlin.steps_per_degree), 80_000_000)
        # Resetting counters is not a physical unwind or a lost full turn.
        self.assertEqual(len([steps for steps in marlin.move_steps if steps > 0]), count * turns_per_layer)

    def test_native_and_workspace_resets_do_not_rotate_an_already_offset_axis(self):
        marlin = MarlinRotaryModel().replay("G0 Y123.25\nG92 Y-45\n")
        physical_before = marlin.physical_steps
        self.assertNotEqual(float(marlin.shift), 0)

        marlin.replay(emit_gcode([LayerToolpaths(0, 30.0)], SlicerConfig()))

        self.assertEqual(float(marlin.native), 0)
        self.assertEqual(float(marlin.shift), 0)
        self.assertEqual(marlin.steps, 0)
        self.assertEqual(marlin.physical_steps, physical_before)

    def test_travel_crossing_the_precision_bound_resets_before_moving(self):
        for sign in (1, -1):
            with self.subTest(sign=sign):
                current = math.radians(sign * 3590)
                state = GcodeState(theta_accum=current, radial=30.0)
                path = Toolpath("SUPPORT", 30.0, [
                    (math.radians(sign * 10), 0.0), (math.radians(sign * 11), 0.0),
                ])
                lines = emit_toolpath(path, state, SlicerConfig(retract_mm=0, min_path_length=0), 0.0)
                marlin = MarlinRotaryModel().replay(f"G0 Y{sign * 3590}\n")
                physical_before = marlin.physical_steps
                marlin.replay("\n".join(lines))

                self.assertEqual(marlin.native_resets, 1)
                self.assertLessEqual(marlin.max_native, MAX_MARLIN_Y_DEGREES)
                degrees_moved = (marlin.physical_steps - physical_before) / float(marlin.steps_per_degree)
                self.assertAlmostEqual(degrees_moved, sign * 21, delta=1 / float(marlin.steps_per_degree))

    def test_mid_path_reset_preserves_small_reverse_moves_and_extrusion(self):
        path = Toolpath("SUPPORT", 30.0, [
            (0.0, 0.0), (math.radians(3599.999), 0.0),
            (math.radians(3599.999), 1.0), (math.radians(3600.025), 1.0),
            (math.radians(3599.999), 1.0),
        ])
        code = emit_gcode([LayerToolpaths(0, 30.0, {"SUPPORT": [path]})],
                          SlicerConfig(retract_mm=0, z_hop_height=0, extrusion_per_mm=0.1, min_path_length=0))
        marlin = MarlinRotaryModel().replay(code)
        self.assertLessEqual(marlin.max_native, MAX_MARLIN_Y_DEGREES)
        self.assertEqual(marlin.move_steps[-2:], [1, -1])
        expected_e = (math.radians(3599.999 + 2 * 0.026) * 30 + 1) * 0.1
        self.assertAlmostEqual(float(re.findall(r"\bE([0-9.]+)", code)[-1]), expected_e, places=5)
        segments = parse_extrusion_segments(code)
        self.assertAlmostEqual(segments[-2][5] - segments[-2][3], 0.026, places=3)
        self.assertAlmostEqual(segments[-1][5] - segments[-1][3], -0.026, places=3)

    def test_huge_internal_winding_never_leaks_into_marlin_coordinates(self):
        for sign in (1, -1):
            with self.subTest(sign=sign):
                current = math.radians(sign * 80_000_000)
                state = GcodeState(theta_accum=current, theta_origin=0, radial=30.0)
                theta = current % (2 * math.pi)
                path = Toolpath("SUPPORT", 30.0, [(theta, 0.0), (theta + math.radians(sign * 0.026), 0.0)])
                lines = emit_toolpath(path, state, SlicerConfig(retract_mm=0, min_path_length=0), 0.0)
                marlin = MarlinRotaryModel().replay("\n".join(lines))
                self.assertEqual(marlin.native_resets, 1)
                self.assertLess(marlin.max_native, 181.0)
                self.assertLess(marlin.max_logical, 1.0)
                self.assertGreater(marlin.move_steps[-1] * sign, 0)

    def test_fractional_layer_endpoints_do_not_accumulate_step_rounding_drift(self):
        for start, step in ((0.0, 0.026), (170.0, 0.026), (-170.0, -0.026)):
            with self.subTest(start=start, step=step):
                layers = [LayerToolpaths(i, 30.0, {"EXTERNAL_PERIMETER": [
                    Toolpath("EXTERNAL_PERIMETER", 30.0,
                             [(math.radians(start + i * step), 0.0),
                              (math.radians(start + (i + 1) * step), 0.0)])
                ]}) for i in range(700)]
                marlin = MarlinRotaryModel().replay(emit_gcode(layers, SlicerConfig(retract_mm=0, z_hop_height=0)))
                expected_steps = marlin._planner_steps(np.float32(start + 700 * step))
                self.assertAlmostEqual(marlin.physical_steps, expected_steps, delta=1)

    def test_unrepresentable_rotary_move_is_rejected_before_writing_gcode(self):
        for sign in (1, -1):
            path = Toolpath("SUPPORT", 30.0, [(0.0, 0.0), (math.radians(sign * 80_000_000), 0.0)])
            with self.subTest(sign=sign), self.assertRaisesRegex(ValueError, "Marlin precision budget"):
                emit_gcode([LayerToolpaths(0, 30.0, {"SUPPORT": [path]})], SlicerConfig())
