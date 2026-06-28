from __future__ import annotations

import math
import unittest

import numpy as np

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.coordinates import (
    TWO_PI,
    cartesian_to_cylindrical,
    cylindrical_distance,
    cylindrical_to_cartesian,
    unwrap_path_from_current,
)


class CoordinatesAndConfigTests(unittest.TestCase):
    def test_cylindrical_cartesian_round_trip(self):
        theta = np.array([0.0, math.pi / 2.0, math.pi, 1.75 * math.pi])
        axial = np.array([0.0, 2.0, 4.0, 6.0])
        cart = cylindrical_to_cartesian(theta, axial, 12.0)
        cyl = cartesian_to_cylindrical(cart)
        self.assertTrue(np.allclose(cyl[:, 0], 12.0))
        self.assertTrue(np.allclose(cyl[:, 1], theta % TWO_PI))
        self.assertTrue(np.allclose(cyl[:, 2], axial))

    def test_arc_length_for_pure_rotation(self):
        self.assertAlmostEqual(
            cylindrical_distance((0.0, 4.0), (math.pi / 2.0, 4.0), 10.0),
            5.0 * math.pi,
        )

    def test_unwrapped_thetas_never_move_backward(self):
        values = unwrap_path_from_current([math.radians(350), math.radians(10), math.radians(30)], math.radians(340))
        self.assertGreaterEqual(values[0], math.radians(340))
        self.assertTrue(all(b >= a for a, b in zip(values, values[1:])))

    def test_unwrap_preserves_intentional_full_turn(self):
        values = unwrap_path_from_current([0.0, TWO_PI], 0.0)
        self.assertAlmostEqual(values[-1] - values[0], TWO_PI)

    def test_unwrap_does_not_turn_same_theta_drift_into_full_turn(self):
        theta = math.radians(355.0)
        values = unwrap_path_from_current([theta, theta - 1e-12], math.radians(350.0))
        self.assertAlmostEqual(values[-1] - values[0], 0.0)

    def test_unwrap_starts_at_nearest_equivalent_to_current(self):
        current = math.radians(358.038)
        values = unwrap_path_from_current([math.radians(357.634), math.radians(357.926)], current)
        self.assertLess(values[0], current)
        self.assertLess(abs(values[0] - current), math.radians(1.0))
        self.assertAlmostEqual(values[1] - values[0], math.radians(0.292), places=9)

    def test_config_validation(self):
        config = SlicerConfig()
        config.validate()
        self.assertEqual(config.top_layers, 3)
        self.assertEqual(config.bottom_layers, 3)
        self.assertEqual(config.raft_layers, 4)
        self.assertIsNone(config.raft_grid_spacing)
        self.assertEqual(config.support_gap_layers, 1)
        self.assertEqual(config.support_interface_layers, 2)
        self.assertIsNone(config.support_interface_spacing)
        self.assertEqual(config.support_meander_max_angle_deg, 15.0)
        self.assertIsNone(config.support_cross_spacing)
        self.assertEqual(config.axial_min, 15.0)
        self.assertEqual(config.raft_axial_gap, 1.0)
        self.assertEqual(config.min_path_length, 0.6)
        self.assertEqual(config.retract_mm, 1.0)
        expected_extrusion = config.nozzle_width * config.layer_height / (
            math.pi * (config.filament_diameter / 2.0) ** 2
        )
        self.assertAlmostEqual(config.extrusion_per_mm, expected_extrusion)

        override = SlicerConfig(extrusion_per_mm=0.035)
        self.assertEqual(override.extrusion_per_mm, 0.035)
        with self.assertRaises(ValueError):
            SlicerConfig(layer_height=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(nozzle_width=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(filament_diameter=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(extrusion_per_mm=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(top_layers=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(bottom_layers=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(raft_layers=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(raft_grid_spacing=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(support_gap_layers=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(support_interface_layers=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(support_interface_spacing=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(support_meander_max_angle_deg=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(support_cross_spacing=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(raft_axial_gap=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(min_path_length=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(retract_mm=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(axial_min=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(axial_step=0).validate()


if __name__ == "__main__":
    unittest.main()
