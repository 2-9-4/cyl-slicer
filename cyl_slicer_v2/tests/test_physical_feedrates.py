from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import _surface_feed_for_marlin, emit_toolpath
from cyl_slicer_v2.models import GcodeState, Toolpath


def _f_value(line: str) -> float:
    match = re.search(r"\bF([0-9.]+)", line)
    if match is None:
        raise AssertionError(f"Missing F: {line}")
    return float(match.group(1))


class PhysicalSurfaceFeedTests(unittest.TestCase):
    def test_theory_x_only_rotary_only_and_mixed(self):
        target = 1800.0  # physical mm/min
        for radius in (30.0, 60.0):
            for dx, surface_arc in ((10.0, 0.0), (0.0, 10.0), (10.0, 10.0)):
                dtheta = surface_arc / radius
                with self.subTest(radius=radius, dx=dx, arc=surface_arc):
                    commanded_f = _surface_feed_for_marlin(dx, dtheta, radius, target)
                    physical_length = math.hypot(dx, surface_arc)
                    marlin_length = math.hypot(dx, math.degrees(dtheta))
                    self.assertAlmostEqual(marlin_length / (commanded_f / 60),
                                           physical_length / (target / 60))
                    if surface_arc == 0:
                        self.assertAlmostEqual(commanded_f, target)

    def test_emitted_print_moves_preserve_requested_surface_speed_and_e(self):
        config = SlicerConfig(axial_min=0.0, min_path_length=0.0, retract_mm=0.0,
                              feed_mm_per_min=1800.0)
        for radius in (30.0, 60.0):
            for dx, arc in ((10.0, 0.0), (0.0, 10.0), (10.0, 10.0)):
                with self.subTest(radius=radius, dx=dx, arc=arc):
                    state = GcodeState(radial=radius)
                    theta = arc / radius
                    path = Toolpath("EXTERNAL_PERIMETER", radius, [(0.0, 0.0), (theta, dx)])
                    lines = emit_toolpath(path, state, config, 0.0, retract_after=False)
                    printed, = [s for s in lines if s.startswith("G1 X")]
                    emitted_f = _f_value(printed)
                    marlin_length = math.hypot(dx, math.degrees(theta))
                    physical_length = math.hypot(dx, arc)
                    self.assertAlmostEqual(marlin_length / (emitted_f / 60),
                                           physical_length / (config.feed_mm_per_min / 60),
                                           delta=0.0002)
                    self.assertAlmostEqual(state.extrusion,
                                           physical_length * config.extrusion_per_mm,
                                           places=8)

    def test_xy_travel_corrected_without_changing_radial_hop_speed(self):
        config = SlicerConfig(axial_min=0.0, min_path_length=0.0,
                              z_hop_height=0.0, retract_mm=0.0,
                              travel_mm_per_min=3000)
        radius = 30.0
        path = Toolpath("EXTERNAL_PERIMETER", radius, [(0.2, 10.0), (0.4, 10.0)])
        lines = emit_toolpath(path, GcodeState(radial=radius), config, 0.0,
                              retract_after=False)
        travel, = [s for s in lines if s.startswith("G0 X")]
        expected = _surface_feed_for_marlin(10.0, 0.2, radius, 3000)
        self.assertAlmostEqual(_f_value(travel), expected, delta=0.05)

    def test_zero_length_surface_move_does_not_divide_by_zero(self):
        self.assertEqual(_surface_feed_for_marlin(0, 0, 30, 1800), 1800)
