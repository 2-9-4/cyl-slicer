from contextlib import redirect_stderr
from dataclasses import replace
import io
import math
from pathlib import Path
import re
import tempfile
import unittest

from cyl_slicer_v2.calibration import CalibrationSweep, RetractionSweep, sweep_from_args
from cyl_slicer_v2.cli import build_parser, main
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode, emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath
from cyl_slicer_v2.tests.test_marlin_rotary_precision import MarlinRotaryModel


def fields(line):
    return {key: float(value) for key, value in
            re.findall(r"\b([XYZEPFT])(-?\d+(?:\.\d+)?(?:e[+-]?\d+)?)", line.split(";")[0])}


def measured_moves(code, bed_radius=30, axis_limits=None):
    """Reconstruct physical acceleration from commands, independent of emitter math."""
    position = dict(X=0., Y=0., Z=0., E=0.)
    acceleration = {}
    limits = dict(axis_limits) if axis_limits is not None else None
    layer = None
    for line in code.splitlines():
        values = fields(line)
        if line.startswith(";LAYER:"):
            layer = int(line.split(":")[1])
        elif line.startswith("M204 "):
            acceleration.update(values)
        elif line.startswith("M201 ") and limits is not None:
            limits.update(values)
        elif line.startswith("G92 "):
            position.update(values)
        elif line.startswith(("G0 ", "G1 ")):
            endpoint = {axis: values.get(axis, current) for axis, current in position.items()}
            dx, dy, dz = (endpoint[k] - position[k] for k in "XYZ")
            native = math.sqrt(dx*dx + dy*dy + dz*dz)
            radius = bed_radius + position["Z"]
            physical = math.sqrt(dx*dx + (dy * math.pi / 180 * radius)**2 + dz*dz)
            printing = endpoint["E"] > position["E"] and native > 0
            key = "P" if printing else "T"
            if native > 1e-12:
                requested = acceleration.get(key, math.nan)
                if limits is not None:
                    for axis, maximum in limits.items():
                        delta = abs(endpoint[axis] - position[axis])
                        if delta > 0:
                            requested = min(requested, maximum*native/delta)
                yield layer, printing, requested * physical/native, line, values
            position = endpoint


def calibration_layers():
    result = []
    # Exact nozzle heights exercise base, boundary, next band and final hold.
    for i, height in enumerate((.5, 1., 2., 4.)):
        radius = 30 + height - .08
        path = Toolpath("GAP_FILL", radius, [(0, 0), (0, 5), (.15, 5), (.3, 0)],
                        segment_widths=[.24, .32, .4])
        result.append(LayerToolpaths(i, radius, {path.role: [path]}, nozzle_radius=30+height))
    return result


class PhysicalAccelerationTests(unittest.TestCase):
    def config(self, **overrides):
        return replace(SlicerConfig(axial_min=0, retract_mm=0, z_hop_height=0,
                                   cooling_min_layer_time_s=0, cooling_min_island_print_time_s=0),
                       **overrides)

    def test_constant_physical_acceleration_for_axial_rotary_diagonal_and_radial_moves(self):
        layers = calibration_layers()
        # Different sample/slab radii ensure conversion uses actual nozzle radius.
        layers[-1] = replace(layers[-1], nozzle_radius=65.)
        config = self.config(print_accel=700, travel_accel=1300,
                             retract_mm=1.4, z_hop_height=.4, z_hop_min_travel=0)
        code = emit_gcode(layers, config)
        counts = {False: 0, True: 0}
        for layer, printing, actual, line, _ in measured_moves(code):
            with self.subTest(layer=layer, line=line):
                self.assertAlmostEqual(actual, 700 if printing else 1300, delta=.00002)
            counts[printing] += 1
        self.assertGreater(counts[True], 5)
        self.assertGreater(counts[False], 5)
        baseline = emit_gcode(layers, replace(config, print_accel=0, travel_accel=0))
        self.assertEqual(baseline, "\n".join(s for s in code.splitlines() if not s.startswith(("M201 ", "M204 "))) + "\n")
        self.assertNotRegex(code, r"(?m)^M(?:203|205|500)\b|^M204 .*\bR|^M201 .*\bE")

    def test_reasserts_each_layer_including_empty_and_partial_layers(self):
        layers = [LayerToolpaths(-1, 30.03, deposited_thickness=.06),
                  LayerToolpaths(0, 31., nozzle_radius=31.2),
                  LayerToolpaths(1, 31.1, nozzle_radius=31.2)]
        code = emit_gcode(layers, self.config(print_accel=600, travel_accel=900))
        for layer, block in zip(layers, code.split(";LAYER:")[1:]):
            first_command = next(s for s in block.splitlines()[1:] if not s.startswith(";"))
            self.assertTrue(first_command.startswith("M204 P"))
            radius = layer.deposition_radius(.16)
            for key, target in (("P", 600), ("T", 900)):
                self.assertAlmostEqual(fields(first_command)[key]*math.pi/180*radius, target, places=5)
        self.assertNotIn("M204", emit_gcode([], self.config(print_accel=600)))

    def test_no_redundant_acceleration_blocks_on_straight_variable_width_moves(self):
        path = Toolpath("GAP_FILL", 31, [(0, 0), (0, 5), (0, 10), (0, 15)],
                        segment_widths=[.24, .32, .4])
        lines = emit_toolpath(path, GcodeState(radial=31), self.config(print_accel=600), 0)
        self.assertEqual(sum(s.startswith("M204 ") for s in lines), 1)
        self.assertEqual(sum(s.startswith("G1 X") for s in lines), 3)

    def test_rotary_resets_and_bidirectional_motion_preserve_acceleration_and_e(self):
        # Vary widths to retain individual segments through simplification.
        paths = [Toolpath("GAP_FILL", 31, [(i*.1, x) for i in range(700)],
                          segment_widths=[.24 if i % 2 else .32 for i in range(699)]) for x in (0, 5)]
        layer = LayerToolpaths(0, 31, {"GAP_FILL": paths})
        for one_direction in (0, 1):
            with self.subTest(one_direction=one_direction):
                base = self.config(one_direction_first_layers=one_direction)
                code = emit_gcode([layer], replace(base, print_accel=600, travel_accel=900))
                self.assertGreater(code.count("G92.9 Y"), 1)
                for _, printing, actual, _, _ in measured_moves(code):
                    self.assertAlmostEqual(actual, 600 if printing else 900, delta=.00002)
                baseline = emit_gcode([layer], base)
                self.assertEqual(baseline, "\n".join(s for s in code.splitlines() if not s.startswith(("M201 ", "M204 "))) + "\n")
                self.assertEqual(MarlinRotaryModel().replay(code).move_steps,
                                 MarlinRotaryModel().replay(baseline).move_steps)

    def test_acceleration_commands_use_decimal_not_scientific_notation(self):
        path = Toolpath("EXTERNAL_PERIMETER", 31, [(0, 0), (0, 5)])
        for acceleration in (1e-5, 1e10):
            lines = emit_toolpath(path, GcodeState(radial=31), self.config(print_accel=acceleration), 0)
            command = next(s.split(";")[0].strip() for s in lines if s.startswith("M204 "))
            self.assertRegex(command, r"^M204 P\d+(?:\.\d+)?$")

    def test_single_startup_m201_covers_print_travel_and_all_declared_bands(self):
        layers = calibration_layers()[:1]  # Highest band isn't reached in this slice.
        for mode, values in (("print-accel", (400., 1600.)),
                             ("speed-accel", ((20., 400.), (80., 1600.)))):
            config = self.config(print_accel=200, travel_accel=900,
                                 retraction_sweep=CalibrationSweep(mode, values, 2, 3))
            code = emit_gcode(layers, config)
            commands = [s for s in code.splitlines() if s.startswith("M201 ")]
            self.assertEqual(len(commands), 1)
            self.assertEqual(fields(commands[0]), dict(X=1601, Y=math.ceil(1600*180/math.pi/30.5)+1, Z=901))
            self.assertLess(code.index(commands[0]), code.index(";LAYER:"))
            self.assertLess(code.index(commands[0]), code.index("\nM204 "))
            self.assertNotRegex(code, r"(?m)^M500\b|^M201 .*\bE|^M204 .*\bR")

    def test_startup_ceilings_prevent_xyz_firmware_limits_clipping_targets(self):
        sweep = CalibrationSweep("speed-accel", ((20., 400.), (40., 800.), (60., 1600.)), 1, 1)
        config = self.config(print_accel=200, travel_accel=900, retraction_sweep=sweep,
                             retract_mm=1.4, z_hop_height=.4, z_hop_min_travel=0)
        layers = calibration_layers()
        code = emit_gcode(layers, config)
        defaults = dict(X=500., Y=1000., Z=100., E=1000.)
        clipped = False
        without_m201 = "\n".join(s for s in code.splitlines() if not s.startswith("M201 "))
        for index, printing, actual, _, _ in measured_moves(without_m201, axis_limits=defaults):
            expected = (200, 400, 800, 1600)[index] if printing else 900
            clipped |= actual < expected-.01
        self.assertTrue(clipped, "Fixture must exercise the original firmware axis caps")
        for index, printing, actual, line, _ in measured_moves(code, axis_limits=defaults):
            expected = (200, 400, 800, 1600)[index] if printing else 900
            with self.subTest(line=line):
                self.assertAlmostEqual(actual, expected, delta=.00002)

    def test_startup_m201_handles_partial_radius_unsorted_layers_and_hops(self):
        layers = [LayerToolpaths(0, 60., nozzle_radius=60.16),
                  LayerToolpaths(1, 30.03, deposited_thickness=.06)]
        code = emit_gcode(layers, self.config(print_accel=1200, travel_accel=2000, z_hop_height=2))
        command = next(s for s in code.splitlines() if s.startswith("M201 "))
        self.assertEqual(fields(command), dict(X=2001, Y=math.ceil(2000*180/math.pi/30.06)+1, Z=2001))

    def test_startup_m201_requires_an_acceleration_target_and_preserves_e_only_settings(self):
        layers = calibration_layers()
        for extra in ({}, {"retraction_sweep": CalibrationSweep("print-speed", (20., 40.))}):
            self.assertNotIn("\nM201 ", emit_gcode(layers, self.config(**extra)))
        self.assertNotIn("\nM201 ", emit_gcode([], self.config(print_accel=1200, travel_accel=2000)))
        print_only = emit_gcode(layers, self.config(print_accel=1200))
        self.assertEqual(set(fields(next(s for s in print_only.splitlines() if s.startswith("M201 ")))), {"X", "Y"})
        travel_only = emit_gcode(layers, self.config(travel_accel=2000))
        self.assertEqual(set(fields(next(s for s in travel_only.splitlines() if s.startswith("M201 ")))), {"X", "Y", "Z"})

    def test_startup_m201_rejects_native_integer_overflow(self):
        with self.assertRaisesRegex(ValueError, "32-bit M201"):
            emit_gcode(calibration_layers(), self.config(print_accel=2**32))


class SpeedAccelerationCalibrationTests(unittest.TestCase):
    def config(self, sweep, **overrides):
        return replace(SlicerConfig(axial_min=0, retract_mm=0, z_hop_height=0,
                                   feed_mm_per_min=900, print_accel=200, travel_accel=900,
                                   cooling_min_layer_time_s=0, cooling_min_island_print_time_s=0,
                                   retraction_sweep=sweep), **overrides)

    def test_modes_units_boundaries_and_final_hold(self):
        for mode, values, expected in (
            ("print-speed", (20., 40.), {"feed_mm_per_min": 2400}),
            ("print-accel", (400., 800.), {"print_accel": 800}),
            ("speed-accel", ((20., 400.), (40., 800.)), {"feed_mm_per_min": 2400, "print_accel": 800}),
        ):
            sweep = CalibrationSweep(mode, values, 1, 1)
            self.assertEqual(sweep.overrides(.9), {})
            self.assertEqual(sweep.band_at(1), 0)
            self.assertEqual(sweep.band_at(2-1e-12), 1)
            self.assertEqual(sweep.overrides(2), expected)
            self.assertEqual(sweep.overrides(100), expected)
        self.assertIs(RetractionSweep, CalibrationSweep)

    def test_paired_sweep_changes_only_speed_and_acceleration_preserving_width_e_and_coordinates(self):
        sweep = CalibrationSweep("speed-accel", ((20., 400.), (40., 800.)), 1, 1)
        config = self.config(sweep, z_hop_height=2, z_hop_min_travel=0, retract_mm=1.4)
        layers = calibration_layers()
        code = emit_gcode(layers, config)
        baseline = emit_gcode(layers, replace(config, retraction_sweep=None))
        def geometry(g):
            return [s.split(" F")[0] for s in g.splitlines() if s.startswith(("G0 ", "G1 ", "G92"))]
        self.assertEqual(geometry(code), geometry(baseline))
        for index, printing, actual, line, _ in measured_moves(code):
            self.assertAlmostEqual(actual, (200, 400, 800, 800)[index] if printing else 900, delta=.00002)
            if printing:
                # First extrusion is axial, so F is directly the physical target.
                first = next(s for s in code.split(f";LAYER:{index}\n")[1].splitlines() if s.startswith("G1 X"))
                self.assertEqual(fields(first)["F"], (900, 1200, 2400, 2400)[index])
        self.assertIn(";CALIBRATION:base height=0.5000", code)
        self.assertIn("value=20:400", code)
        self.assertIn("values=20:400,40:800", code)

    def test_cooling_slow_start_role_and_flow_caps_remain_active(self):
        sweep = CalibrationSweep("speed-accel", ((100., 800.),), 0, 1)
        for options, ceiling in ((dict(slow_start_layers=5, slow_start_print_feed_mm_per_min=600), 600),
                                 (dict(cooling_min_layer_time_s=100), 301),
                                 (dict(max_volumetric_flow_mm3_s=.2), 601),
                                 (dict(support_print_feed_mm_per_min=450), 450)):
            radius = 31
            path = Toolpath("SUPPORT", radius, [(0, 0), (0, 10)])
            code = emit_gcode([LayerToolpaths(0, radius, {path.role: [path]})], self.config(sweep, **options))
            move = next(s for s in code.splitlines() if s.startswith("G1 X"))
            self.assertLessEqual(fields(move)["F"], ceiling)
            self.assertNotIn("\nG4", code)

    def test_invalid_and_cli_values(self):
        parser = build_parser()
        for mode, text, expected in (("print-speed", "20,40", (20, 40)),
                                     ("print-accel", "400,800", (400, 800)),
                                     ("speed-accel", "20:400,40:800", ((20, 400), (40, 800)))):
            args = parser.parse_args(["--calibration", mode, "--calibration-values", text])
            self.assertEqual(sweep_from_args(args, parser).values, expected)
        for mode, text in (("print-speed", "0"), ("print-speed", "1e308"),
                           ("print-accel", "-1"), ("print-accel", "nan"),
                           ("speed-accel", "20"), ("speed-accel", "20:400:5"),
                           ("speed-accel", "20:0"), ("speed-accel", "20:inf")):
            with self.subTest(mode=mode, text=text), redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                sweep_from_args(parser.parse_args(["--calibration", mode, "--calibration-values", text]), parser)

    def test_cli_generates_paired_calibration_gcode(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()):
            output = Path(tmp)/"speed-accel.gcode"
            self.assertEqual(main(["--fixture", "tube", "--max-layers", "3", "--out", str(output),
                                   "--build-cylinder-radius", "8",
                                   "--slice-workers", "1", "--toolpath-workers", "1",
                                   "--calibration", "speed-accel", "--calibration-values", "20:400,40:800",
                                   "--calibration-base-height", "0", "--calibration-band-height", ".16"]), 0)
            code = output.read_text()
            self.assertIn(";CALIBRATION_SCHEDULE:speed-accel values=20:400,40:800", code)
            self.assertIn("; layer surface acceleration", code)
            self.assertEqual(code.count("\nM201 "), 1)


if __name__ == "__main__":
    unittest.main()
