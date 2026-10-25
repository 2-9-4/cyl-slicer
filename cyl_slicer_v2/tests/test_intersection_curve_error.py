from __future__ import annotations

import math
import unittest

import numpy as np
import trimesh

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.intersection_backend import cylindrical_intersection_segments


def sloping_triangle():
    # z = x + y plane: the intersection with radius 10 is a genuine sinusoid
    # in unwrapped physical (R*theta, axial) surface coordinates.
    return trimesh.Trimesh(vertices=np.asarray([(5., 0., 5.), (20., 0., 20.), (0., 20., 20.)]),
                           faces=np.asarray([[0, 1, 2]]), process=False)


class IntersectionCurveErrorTests(unittest.TestCase):
    def test_physical_chord_error_is_below_requested_tolerance(self):
        radius, tolerance = 10.0, 0.02
        points, edges = cylindrical_intersection_segments(
            sloping_triangle(), radius, max_curve_error_mm=tolerance
        )
        self.assertGreater(len(edges), 2)
        for first, second in edges:
            a, b = points[first], points[second]
            theta0 = math.atan2(a[1], a[0])
            theta1 = math.atan2(b[1], b[0])
            dt = (theta1 - theta0 + math.pi) % (2 * math.pi) - math.pi
            tm = theta0 + dt * 0.5
            true_z = radius * (math.cos(tm) + math.sin(tm))
            chord_z = (float(a[2]) + float(b[2])) * 0.5
            self.assertLessEqual(abs(true_z - chord_z), tolerance + 1e-8)
            self.assertAlmostEqual(math.hypot(a[0], a[1]), radius, places=8)
            self.assertAlmostEqual(math.hypot(b[0], b[1]), radius, places=8)

    def test_tighter_tolerance_produces_more_segments(self):
        mesh = sloping_triangle()
        coarse = cylindrical_intersection_segments(mesh, 10, max_curve_error_mm=0.4)
        fine = cylindrical_intersection_segments(mesh, 10, max_curve_error_mm=0.002)
        self.assertGreater(len(fine[1]), len(coarse[1]))

    def test_flat_face_yields_one_straight_block(self):
        mesh = trimesh.Trimesh(vertices=np.asarray([
            (5, 0, 4), (20, 0, 4), (0, 20, 4)]), faces=[[0, 1, 2]], process=False)
        points, edges = cylindrical_intersection_segments(
            mesh, 10, max_axial_step=0.01, max_curve_error_mm=0.025
        )
        self.assertEqual(len(edges), 1)

    def test_invalid_curve_error_rejected(self):
        for tolerance in (0, -0.001, math.inf, math.nan):
            with self.subTest(tolerance=tolerance):
                with self.assertRaises(ValueError):
                    SlicerConfig(intersection_curve_error_mm=tolerance).validate()
                with self.assertRaises(ValueError):
                    cylindrical_intersection_segments(
                        sloping_triangle(), 10, max_curve_error_mm=tolerance
                    )
