from __future__ import annotations

import math
import unittest

import numpy as np
import trimesh

from cyl_slicer_v2.nonplanar_surface import (
    RadialSurfaceSampler, sample_surface_ribbon, segment_surface_length,
)


TAU = 2 * math.pi


def cylindrical_shell(radius, *, theta_segments=96, axial_segments=8,
                      theta_min=0.0, theta_max=TAU, skip=None,
                      inward=False):
    """Make outward-wound, possibly open patches in CAD XYZ coordinates."""
    vertices = []
    for j in range(axial_segments + 1):
        x = j * 1.0
        for i in range(theta_segments + 1):
            theta = theta_min + (theta_max-theta_min)*i/theta_segments
            r = radius(theta, x)
            vertices.append((r*math.cos(theta),r*math.sin(theta),x))
    triangles = []
    for j in range(axial_segments):
        for i in range(theta_segments):
            theta = theta_min+(theta_max-theta_min)*(i+.5)/theta_segments
            x = j+.5
            if skip is not None and skip(theta,x):
                continue
            a = j*(theta_segments+1)+i
            b = a+1
            c = a+(theta_segments+1)
            d = c+1
            triangles.extend(((a,b,c),(b,d,c)))
    if inward:
        triangles = [(a,c,b) for a,b,c in triangles]
    return trimesh.Trimesh(vertices=vertices, faces=triangles, process=False)


def together(*meshes):
    vertices, faces = [], []
    for mesh in meshes:
        faces.extend((mesh.faces + len(vertices)).tolist())
        vertices.extend(mesh.vertices.tolist())
    return trimesh.Trimesh(vertices=vertices,faces=faces,process=False)


class RadialSurfaceTests(unittest.TestCase):
    def test_cylinder_unique_surface_and_periodic_seam(self):
        mesh = cylindrical_shell(lambda t,x: 30)
        sampler = RadialSurfaceSampler(mesh)
        for theta in (0, .02, 1.9, TAU-.00001, TAU):
            one=sampler.sample(theta, 3.5)
            two=sampler.sample(theta + 2*TAU, 3.5)
            self.assertIsNotNone(one)
            self.assertIsNotNone(two)
            self.assertAlmostEqual(one.radius,two.radius,places=6)
            self.assertLess(abs(one.radius-30),.017)
            self.assertGreater(one.outward_normal_radial,.95)
        self.assertIsNone(sampler.sample(1.,10.5))
        with self.assertRaisesRegex(ValueError,"finite"):
            sampler.sample(math.nan,3)

    def test_taper_ribbon_follows_true_surface_and_has_no_default_gcode(self):
        sampler = RadialSurfaceSampler(cylindrical_shell(lambda t,x: 30+.08*x))
        ribbon = sample_surface_ribbon(sampler,[(.2,1),(.3,5)],max_segment_mm=.5)
        self.assertGreater(len(ribbon),8)
        self.assertAlmostEqual(ribbon[0].radius,30.08,delta=.018)
        self.assertAlmostEqual(ribbon[-1].radius,30.4,delta=.018)
        self.assertTrue(all(segment_surface_length(a,b) <= .50001
                            for a,b in zip(ribbon,ribbon[1:])))
        self.assertTrue(any(abs(b.radius-a.radius)>.001 for a,b in zip(ribbon,ribbon[1:])))

    def test_shallow_crown_midpoint_error_is_bounded(self):
        mesh=cylindrical_shell(lambda t,x: 30+.18*math.cos(2*t)+.08*math.sin(.4*x),
                               theta_segments=128,axial_segments=12)
        sampler=RadialSurfaceSampler(mesh)
        ribbon=sample_surface_ribbon(sampler,[(.02,2),(.42,4)],
                       max_segment_mm=.65,max_radial_chord_error_mm=.002,
                       max_radial_slope=.3)
        self.assertGreater(len(ribbon),5)
        for a,b in zip(ribbon,ribbon[1:]):
            mid=sampler.sample((a.theta+b.theta)/2,(a.axial+b.axial)/2)
            self.assertIsNotNone(mid)
            self.assertLessEqual(abs(mid.radius-(a.radius+b.radius)/2),.002001)

    def test_holes_do_not_get_bridge_ribbons(self):
        mesh=cylindrical_shell(lambda t,x: 30,theta_segments=120,
               skip=lambda t,x: .25<t<.5 and 2<x<5)
        sampler=RadialSurfaceSampler(mesh)
        self.assertIsNone(sampler.sample(.35,3))
        with self.assertRaisesRegex(ValueError,"hole|surface"):
            sample_surface_ribbon(sampler,[(.1,3),(.7,3)],
                                  max_segment_mm=.15)

    def test_two_outer_sheets_and_reversed_surface_are_rejected(self):
        outward = cylindrical_shell(lambda t,x: 31)
        hidden = cylindrical_shell(lambda t,x: 29)
        duplicate=together(outward,hidden)
        sampler=RadialSurfaceSampler(duplicate)
        self.assertIsNone(sampler.sample(.4,3))
        self.assertIsNone(RadialSurfaceSampler(cylindrical_shell(
            lambda t,x: 30,inward=True)).sample(.4,3))

    def test_slope_excursion_budget_and_invalid_bounds_fail(self):
        sampler=RadialSurfaceSampler(cylindrical_shell(lambda t,x: 30+.5*x))
        with self.assertRaisesRegex(ValueError,"slope"):
            sample_surface_ribbon(sampler,[(.2,1),(.2,4)],
                                  max_segment_mm=.5,max_radial_slope=.2,
                                  max_total_radial_excursion_mm=3)
        with self.assertRaisesRegex(ValueError,"excursion"):
            sample_surface_ribbon(sampler,[(.2,1),(.2,4)],
                                  max_segment_mm=.5,max_radial_slope=.8,
                                  max_total_radial_excursion_mm=.5)
        with self.assertRaisesRegex(ValueError,"finite and positive"):
            sample_surface_ribbon(sampler,[(.2,1),(.2,2)],max_segment_mm=0)
        with self.assertRaisesRegex(ValueError,"subdivision budget"):
            sample_surface_ribbon(sampler,[(.2,1),(.2,7)],max_segment_mm=.1,
                                  max_points=6,max_radial_slope=.9,
                                  max_total_radial_excursion_mm=5)

    def test_nonplanar_candidate_is_seam_translation_invariant(self):
        sampler=RadialSurfaceSampler(cylindrical_shell(lambda t,x: 30+.1*x))
        baseline=sample_surface_ribbon(sampler,[(TAU-.06,2),(TAU+.05,3)],
                                       max_segment_mm=.35)
        for turns in (-2,-1,1,3):
            shifted=sample_surface_ribbon(sampler,
                 [(TAU-.06+turns*TAU,2),(TAU+.05+turns*TAU,3)],
                 max_segment_mm=.35)
            self.assertEqual(len(baseline),len(shifted))
            for a,b in zip(baseline,shifted):
                self.assertAlmostEqual(a.radius,b.radius,places=6)
                self.assertAlmostEqual(b.theta-a.theta,turns*TAU,places=7)
                self.assertAlmostEqual(a.axial,b.axial,places=8)


if __name__=="__main__":
    unittest.main()
