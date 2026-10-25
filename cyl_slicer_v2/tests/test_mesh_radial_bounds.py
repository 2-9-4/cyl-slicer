from __future__ import annotations

import math
import unittest

import numpy as np
import trimesh

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.layers import layer_radii
from cyl_slicer_v2.mesh_io import radial_bounds


def _triangle(points):
    return trimesh.Trimesh(vertices=np.asarray(points, dtype=float),
                           faces=[[0, 1, 2]], process=False)


class TrueRadialBoundsTests(unittest.TestCase):
    def test_offset_box_minimum_occurs_on_face_not_a_vertex(self):
        mesh = trimesh.creation.box(extents=(10, 10, 10))
        mesh.apply_translation((30, 0, 0))
        minimum, maximum = radial_bounds(mesh)
        self.assertAlmostEqual(minimum, 25.0, places=9)
        self.assertAlmostEqual(maximum, math.hypot(35, 5), places=9)
        # The missed first slab previously started at ~25.495 mm instead of 25.
        radii = layer_radii(mesh, SlicerConfig(layer_height=0.2))
        self.assertAlmostEqual(radii[0], 25.1, places=8)

    def test_origin_inside_triangle_projection(self):
        mesh = _triangle([(-2, -1, 3), (2, -1, 5), (0, 3, 8)])
        self.assertAlmostEqual(radial_bounds(mesh)[0], 0.0)

    def test_origin_on_edge_in_projection(self):
        mesh = _triangle([(-2, 0, 0), (2, 0, 3), (2, 4, 8)])
        self.assertAlmostEqual(radial_bounds(mesh)[0], 0.0)

    def test_projected_degenerate_triangle_away_from_origin_is_not_inside(self):
        mesh = _triangle([(2, 2, 0), (3, 3, 1), (4, 4, 2)])
        self.assertAlmostEqual(radial_bounds(mesh)[0], math.sqrt(8), places=9)

    def test_nearest_edge_point_between_its_vertices(self):
        mesh = _triangle([(25, -5, 0), (25, 5, 3), (35, 5, 5)])
        self.assertAlmostEqual(radial_bounds(mesh)[0], 25.0, places=9)

    def test_negative_quadrant_and_axis_translation(self):
        mesh = trimesh.creation.box(extents=(10, 10, 10))
        mesh.apply_translation((-30, -30, 0))
        self.assertAlmostEqual(radial_bounds(mesh)[0], math.hypot(25, 25), places=9)

    def test_mesh_without_faces_retains_vertex_fallback(self):
        mesh = trimesh.Trimesh(vertices=[[4, 0, 0], [0, 7, 0]],
                               faces=[], process=False)
        self.assertAlmostEqual(radial_bounds(mesh)[0], 4.0)
        self.assertAlmostEqual(radial_bounds(mesh)[1], 7.0)
