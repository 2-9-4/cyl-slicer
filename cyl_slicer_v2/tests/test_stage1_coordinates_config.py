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

    def test_config_validation(self):
        SlicerConfig().validate()
        with self.assertRaises(ValueError):
            SlicerConfig(layer_height=0).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(nozzle_width=-1).validate()
        with self.assertRaises(ValueError):
            SlicerConfig(axial_step=0).validate()


if __name__ == "__main__":
    unittest.main()
