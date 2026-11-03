from __future__ import annotations

from dataclasses import replace
import math
import unittest

from cyl_slicer_v2.nonplanar_clearance import NozzleTipProfile, probe_nozzle_ribbon
from cyl_slicer_v2.nonplanar_surface import RadialSurfaceSampler, TAU
from cyl_slicer_v2.tests.test_nonplanar_surface import cylindrical_shell


def profile(**overrides):
    return replace(NozzleTipProfile(scan_radius_mm=1.0,
                   neighbor_step_mm=.25,path_step_mm=.4,
                   contact_allowance_mm=.03,safety_margin_mm=0),**overrides)


class NozzleInterferenceTests(unittest.TestCase):
    def test_cylinder_has_no_interference_in_sampled_conical_tip_model(self):
        sampler=RadialSurfaceSampler(cylindrical_shell(lambda t,x:30,
                                      theta_segments=120))
        result=probe_nozzle_ribbon(sampler,[(.5,4)],profile())
        self.assertEqual(result["status"],"no_sampled_interference")
        self.assertFalse(result["candidate_is_collision_certified"])
        self.assertEqual(result["unknown_samples"],0)
        self.assertGreater(result["probes"],30)
        self.assertGreaterEqual(result["min_sampled_clearance_mm"],-1e-8)

    def test_neighboring_radial_boss_protrudes_into_tip(self):
        mesh=cylindrical_shell(lambda t,x:30 if x<4 else 31,
                               theta_segments=120)
        sampler=RadialSurfaceSampler(mesh)
        result=probe_nozzle_ribbon(sampler,[(.4,3.0)],profile())
        self.assertEqual(result["status"],"sampled_interference")
        self.assertGreater(result["first_interference"]["radial_penetration_mm"],.01)
        self.assertFalse(result["candidate_is_collision_certified"])

    def test_intermediate_path_pose_is_checked_not_only_endpoints(self):
        mesh=cylindrical_shell(lambda t,x:30 if x<4 or x>5 else 31,
                               theta_segments=120)
        sampler=RadialSurfaceSampler(mesh)
        result=probe_nozzle_ribbon(sampler,[(.4,2.4),(.4,6.6)],profile())
        self.assertEqual(result["status"],"sampled_interference")
        self.assertGreater(result["path_samples"],1)

    def test_unresolved_external_coverage_is_not_assumed_safe(self):
        mesh=cylindrical_shell(lambda t,x:30,
                               theta_segments=32,
                               theta_min=.36,theta_max=.44)
        sampler=RadialSurfaceSampler(mesh)
        result=probe_nozzle_ribbon(sampler,[(.4,4)],profile())
        self.assertEqual(result["status"],"unresolved_neighbors")
        self.assertGreater(result["unknown_samples"],0)
        self.assertFalse(result["candidate_is_collision_certified"])

    def test_periodic_rotary_offsets_are_equivalent(self):
        sampler=RadialSurfaceSampler(cylindrical_shell(
                  lambda t,x:30,theta_segments=120))
        base=probe_nozzle_ribbon(sampler,[(.03,4)],profile())
        for turns in (-2,-1,1,3):
            shifted=probe_nozzle_ribbon(sampler,[(.03+turns*TAU,4)],
                                         profile())
            self.assertEqual(base["status"],shifted["status"])
            self.assertEqual(base["unknown_samples"],shifted["unknown_samples"])
            self.assertAlmostEqual(base["min_sampled_clearance_mm"],
                                   shifted["min_sampled_clearance_mm"],delta=1e-8)

    def test_budget_exhaustion_cannot_claim_clearance(self):
        sampler=RadialSurfaceSampler(cylindrical_shell(lambda t,x:30))
        result=probe_nozzle_ribbon(sampler,[(.5,4)],
                                    profile(max_total_probes=10))
        self.assertEqual(result["status"],"budget_exceeded")
        self.assertFalse(result["candidate_is_collision_certified"])

    def test_invalid_profiles_and_centers_fail_closed(self):
        sampler=RadialSurfaceSampler(cylindrical_shell(lambda t,x:30))
        for changed in (
            {"flat_tip_radius_mm":0},
            {"flat_tip_radius_mm":2},
            {"cone_half_angle_deg":90},
            {"cone_half_angle_deg":0},
            {"scan_radius_mm":0},
            {"neighbor_step_mm":float("nan")},
            {"contact_allowance_mm":-1},
            {"max_total_probes":0},
        ):
            with self.subTest(profile=changed),self.assertRaises(ValueError):
                probe_nozzle_ribbon(sampler,[(.5,4)],profile(**changed))
        for positions in ([],[(float("nan"),4)],[(.1,)],[(.1,4,30)]):
            with self.subTest(positions=positions),self.assertRaises(ValueError):
                probe_nozzle_ribbon(sampler,positions,profile())


if __name__=="__main__":
    unittest.main()
