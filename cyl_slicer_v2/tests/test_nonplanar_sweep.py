from __future__ import annotations

import math
import unittest

import numpy as np
import trimesh

from cyl_slicer_v2.nonplanar_sweep import (
    RadialHeadSection, audit_continuous_head_sweep, cone_head_sections,
)
from cyl_slicer_v2.tests.test_nonplanar_surface import cylindrical_shell


def audit(mesh, poses, sections=None, **kw):
    if sections is None:
        sections=cone_head_sections(flat_tip_radius_mm=.2,
                                   cone_half_angle_deg=35,
                                   cone_height_mm=2,
                                   min_height_mm=.15,
                                   section_count=8)
    return audit_continuous_head_sweep(mesh,poses,sections,
                                      **kw)


class NonplanarContinuousSweepTests(unittest.TestCase):
    def test_smooth_cylinder_has_provable_clearance_for_declared_raised_sections(self):
        mesh=cylindrical_shell(lambda theta,x:30,theta_segments=60)
        result=audit(mesh,[(.25,2,30),(.25,5,30)])
        self.assertEqual(result["status"],"modeled_sections_clear")
        self.assertTrue(result["modeled_sections_clear"])
        self.assertFalse(result["full_head_collision_certified"])
        self.assertGreater(result["checked_triangle_roots"],100)
        self.assertEqual(result["unresolved_regions"],0)
        self.assertFalse(result["mesh_watertight"])  # open surface fixture

    def test_moving_rotary_head_on_smooth_cylinder_can_clear_body(self):
        mesh=cylindrical_shell(lambda theta,x:30,theta_segments=72)
        result=audit(mesh,[(.2,2,30),(.35,3.2,30)],
                     max_nodes=200000)
        self.assertEqual(result["status"],"modeled_sections_clear")
        self.assertTrue(result["modeled_sections_clear"])

    def test_mid_motion_boss_is_witnessed_despite_safe_endpoints(self):
        # Only the central axial mesh row is raised: true narrow boss.
        mesh=cylindrical_shell(
            lambda theta,x: 31 if abs(x-4) < .1 else 30,
            theta_segments=48,axial_segments=9)
        section=[RadialHeadSection(.2,1.5,.65)]
        left=audit(mesh,[(.4,2,30)],section)
        right=audit(mesh,[(.4,6,30)],section)
        combined=audit(mesh,[(.4,2,30),(.4,6,30)],section)
        self.assertTrue(left["modeled_sections_clear"])
        self.assertTrue(right["modeled_sections_clear"])
        self.assertEqual(combined["status"],"witnessed_interference")
        self.assertFalse(combined["modeled_sections_clear"])
        self.assertEqual(combined["witness"]["segment_index"],0)
        self.assertGreater(combined["witness"]["head_local_radial_mm"],.2)
        self.assertLess(combined["witness"]["head_local_radial_mm"],1.5)

    def test_wide_heater_section_detects_obstruction_tip_could_avoid(self):
        mesh=cylindrical_shell(
            lambda theta,x: 31.5 if abs(x-4)<.1 else 30,
            theta_segments=48,axial_segments=9)
        result=audit(mesh,[(.4,2.5,30)],
                     [RadialHeadSection(.8,2.1,2.3)])
        self.assertEqual(result["status"],"witnessed_interference")
        self.assertEqual(result["witness"]["section_index"],0)

    def test_conservative_unresolved_is_not_promoted_to_clear(self):
        mesh=trimesh.Trimesh(vertices=np.asarray([
             [31,0,0],[31,0,8],[31,1,0]],dtype=float),
             faces=[[0,1,2]],process=False)
        section=[RadialHeadSection(.5,1.5,.15)]
        result=audit(mesh,[(0,2.3,30)],section,
                     max_time_depth=0,max_face_depth=0)
        self.assertEqual(result["status"],"unresolved_geometry")
        self.assertFalse(result["modeled_sections_clear"])
        self.assertFalse(result["full_head_collision_certified"])
        self.assertGreater(result["unresolved_regions"],0)
        self.assertIsNotNone(result["first_uncertain"])

    def test_triangle_spatial_refinement_rejects_false_box_overlap(self):
        mesh=trimesh.Trimesh(vertices=np.asarray([
            [31,0,0],[31,1,0],[31,1,8]],dtype=float),
            faces=[[0,1,2]],process=False)
        section=[RadialHeadSection(.5,1.5,.12)]
        # Axial 4.5 is outside this triangle at theta 0; coarse triangle
        # AABB overlaps section but conservative spatial splitting separates.
        coarse=audit(mesh,[(0,4.5,30)],section,
                     max_face_depth=0,max_time_depth=0)
        finer=audit(mesh,[(0,4.5,30)],section,
                    max_face_depth=8,max_time_depth=0,max_nodes=30000)
        self.assertEqual(coarse["status"],"unresolved_geometry")
        self.assertEqual(finer["status"],"modeled_sections_clear")

    def test_budget_exhaustion_is_explicit_failure(self):
        mesh=trimesh.Trimesh(vertices=np.asarray([
            [31,0,0],[31,0,8],[31,1,0]],dtype=float),
            faces=[[0,1,2]],process=False)
        result=audit(mesh,[(0,2.3,30)],
                     [RadialHeadSection(.5,1.5,.15)],
                     max_nodes=1,max_face_depth=4)
        self.assertEqual(result["status"],"budget_exceeded")
        self.assertFalse(result["modeled_sections_clear"])

    def test_winding_translation_preserves_collision_and_clearance(self):
        mesh=cylindrical_shell(lambda theta,x:30,theta_segments=60)
        sections=[RadialHeadSection(.3,1.5,.7)]
        base=audit(mesh,[(.1,2,30),(.1,5,30)],sections)
        for turns in (-3,-1,1,4):
            shifted=audit(mesh,[(.1+turns*2*math.pi,2,30),
                                (.1+turns*2*math.pi,5,30)],sections)
            self.assertEqual(base["status"],shifted["status"])
            self.assertEqual(base["modeled_sections_clear"],
                             shifted["modeled_sections_clear"])

    def test_cone_envelope_is_outer_bound_at_band_top(self):
        sections=cone_head_sections(flat_tip_radius_mm=.2,
                    cone_half_angle_deg=30,cone_height_mm=1.6,
                    min_height_mm=.2,section_count=7)
        self.assertEqual(len(sections),7)
        for sec in sections:
            self.assertAlmostEqual(
                sec.radius_mm,.2+sec.max_height_mm*math.tan(math.pi/6))
            self.assertGreater(sec.min_height_mm,0)
            self.assertGreater(sec.max_height_mm,sec.min_height_mm)

    def test_bad_geometry_settings_and_models_fail_closed(self):
        mesh=cylindrical_shell(lambda theta,x:30)
        pose=[(.2,2,30)]
        for kwargs in (
            {"radial_margin_mm":-1},
            {"lateral_margin_mm":float("nan")},
            {"max_time_depth":-1},
            {"max_face_depth":11},
            {"max_nodes":0},
        ):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                audit(mesh,pose,**kwargs)
        for sec in (
            RadialHeadSection(-.1,1,1),
            RadialHeadSection(.2,.1,1),
            RadialHeadSection(.2,1,0),
            RadialHeadSection(.2,1,math.inf),
        ):
            with self.subTest(sec=sec),self.assertRaises(ValueError):
                audit(mesh,pose,[sec])
        for coord in ((.2,2,-1),(float("nan"),2,30),(.2,2,float("inf"))):
            with self.subTest(coord=coord),self.assertRaises(ValueError):
                audit(mesh,[coord])
        with self.assertRaises(ValueError):
            audit(mesh,[],[RadialHeadSection(.1,1,.3)])
        with self.assertRaises(ValueError):
            audit(mesh,pose,[])
        for opts in (
            {"flat_tip_radius_mm":0},
            {"cone_half_angle_deg":90},
            {"cone_height_mm":.1},
            {"section_count":0},
        ):
            with self.subTest(opts=opts),self.assertRaises(ValueError):
                args={"flat_tip_radius_mm":.2,"cone_half_angle_deg":30,
                      "cone_height_mm":1,"min_height_mm":.2}
                args.update(opts)
                cone_head_sections(**args)


if __name__ == "__main__":
    unittest.main()
