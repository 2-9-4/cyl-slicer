from __future__ import annotations

import math
import unittest
from dataclasses import replace

from cyl_slicer_v2.nonplanar_surface import RadialSurfaceSampler, TAU
from cyl_slicer_v2.nonplanar_patches import scan_surface_patches
from cyl_slicer_v2.tests.test_nonplanar_surface import cylindrical_shell, together


def scan(mesh, **options):
    return scan_surface_patches(mesh, theta_min=0, theta_max=TAU,
                                axial_min=1, axial_max=7,
                                theta_cells=24, axial_cells=6,
                                max_chord_error_mm=.045, **options)


class PatchScannerTests(unittest.TestCase):
    def test_full_turn_cylinder_is_one_periodic_patch_without_duplicate_seam(self):
        mesh=cylindrical_shell(lambda t,x: 30, theta_segments=96)
        result=scan(mesh)
        self.assertEqual(result["status"],"geometry_candidates_only")
        self.assertFalse(result["nonplanar_extrusion_enabled"])
        self.assertFalse(result["collision_certified"])
        self.assertEqual(result["accepted_cells"],144)
        self.assertEqual(result["rejected_cells"],0)
        self.assertEqual(len(result["patches"]),1)
        self.assertTrue(result["patches"][0]["crosses_periodic_seam"])
        self.assertEqual(result["patches"][0]["cell_count"],144)
        self.assertGreater(result["patches"][0]["approx_unwrapped_area_mm2"],1100)
        self.assertEqual(len(result["candidate_mask"]),6)
        self.assertTrue(all(row == "1"*24 for row in result["candidate_mask"]))

    def test_gradual_topography_is_eligible_without_flattening(self):
        mesh=cylindrical_shell(lambda t,x: 30+.08*x + .04*math.sin(t),
                               theta_segments=96)
        result=scan(mesh, max_cell_excursion_mm=.5)
        self.assertGreater(result["accepted_cells"],100)
        self.assertEqual(len(result["patches"]),1)
        self.assertGreater(result["patches"][0]["max_radius_mm"]
                           -result["patches"][0]["min_radius_mm"],.35)

    def test_tiny_hole_between_sample_nodes_is_rejected_by_mesh_edges(self):
        # Theta point probes occur at 0.4 and 0.45, but the missing
        # facet band occupies theta ~0.419..0.445 within the *same* cell.
        mesh=cylindrical_shell(
            lambda t,x:30, theta_segments=240,
            skip=lambda t,x: .422<t<.444 and 3<x<4)
        result=scan_surface_patches(
            mesh,theta_min=.2,theta_max=.8,axial_min=2.2,axial_max=5.8,
            theta_cells=6,axial_cells=4,max_chord_error_mm=.03)
        self.assertIn("mesh_boundary_or_radial_crease",result["rejection_counts"])
        self.assertGreater(result["rejected_cells"],0)
        self.assertGreater(result["accepted_cells"],0)
        self.assertNotEqual(result["candidate_mask"][1][2],"1")

    def test_ambiguous_double_sheet_rejected_and_no_live_print(self):
        mesh=together(cylindrical_shell(lambda t,x:31),
                      cylindrical_shell(lambda t,x:29))
        result=scan(mesh)
        self.assertEqual(result["accepted_cells"],0)
        self.assertEqual(result["patches"],[])
        self.assertFalse(result["nonplanar_extrusion_enabled"])
        self.assertGreater(result["rejection_counts"].get(
            "missing_or_ambiguous_surface",0),0)

    def test_excessive_slope_and_excursion_are_independent_guards(self):
        mesh=cylindrical_shell(lambda t,x:30+.65*x)
        default=scan(mesh)
        self.assertEqual(default["accepted_cells"],0)
        self.assertGreater(default["rejection_counts"].get("radial_excursion",0),0)
        # Allow the larger within-cell height range but retain the same
        # steep-surface slope ceiling. Now the slope guard should reject.
        result=scan(mesh,max_cell_excursion_mm=2.0)
        self.assertEqual(result["accepted_cells"],0)
        self.assertGreater(result["rejection_counts"].get("radial_slope",0),0)

    def test_full_turn_translation_is_invariant(self):
        mesh=cylindrical_shell(lambda t,x:30+.08*x)
        baseline=scan(mesh)
        shifted=scan_surface_patches(mesh,
                      theta_min=-2*TAU,theta_max=-TAU,
                      axial_min=1,axial_max=7,
                      theta_cells=24,axial_cells=6,max_chord_error_mm=.045)
        self.assertEqual(baseline["candidate_mask"],shifted["candidate_mask"])
        self.assertEqual(baseline["accepted_cells"],shifted["accepted_cells"])
        self.assertEqual(len(baseline["patches"]),len(shifted["patches"]))
        self.assertEqual(baseline["patches"][0]["cell_count"],
                         shifted["patches"][0]["cell_count"])

    def test_invalid_budget_and_scan_ranges_are_rejected(self):
        mesh=cylindrical_shell(lambda t,x:30)
        cases=[
            {"theta_min":0,"theta_max":TAU+0.1},
            {"axial_min":3,"axial_max":2},
            {"max_cells":50},
            {"max_slope":0},
            {"theta_cells":0},
            {"axial_cells":1.5},
            {"sharp_dihedral_deg":181},
            {"max_chord_error_mm":float("nan")},
        ]
        for change in cases:
            opts={"theta_min":0,"theta_max":TAU,"axial_min":1,"axial_max":7,
                  "theta_cells":24,"axial_cells":6}
            opts.update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                scan_surface_patches(mesh,**opts)

    def test_small_patches_can_be_removed_without_claiming_collision_safety(self):
        mesh=cylindrical_shell(lambda t,x:30,theta_segments=96,
               skip=lambda t,x: .3<t<.6 and 1<x<7)
        result=scan_surface_patches(mesh,theta_min=.1,theta_max=.8,
                      axial_min=1.2,axial_max=6.8,
                      theta_cells=14,axial_cells=8,
                      max_chord_error_mm=.05,min_patch_cells=40)
        self.assertFalse(result["collision_certified"])
        self.assertEqual(result["accepted_cells"],
                         sum(p["cell_count"] for p in result["patches"]))


class AxialIndexTests(unittest.TestCase):
    def test_index_matches_exact_unindexed_triangle_queries(self):
        # Lots of axial strips ensure the conservative bin path is active.
        mesh=cylindrical_shell(lambda t,x:30+.1*x,
                  theta_segments=28,axial_segments=160)
        sampler=RadialSurfaceSampler(mesh)
        self.assertIsNotNone(sampler._axial_bins)
        index=sampler._axial_bins
        sampler._axial_bins=None
        positions=[(a,b) for a in (-TAU-.04,0,.1,3.,TAU-.01,TAU+.05)
                        for b in (0, .025, 1.2, 4.0, 7.98, 8.0, 10)]
        baseline=[sampler.sample(a,b) for a,b in positions]
        sampler._axial_bins=index
        actual=[sampler.sample(a,b) for a,b in positions]
        self.assertEqual([a is None for a in baseline],
                         [a is None for a in actual])
        for a,b in zip(baseline,actual):
            if a is not None:
                self.assertAlmostEqual(a.radius,b.radius,places=7)
                self.assertAlmostEqual(a.outward_normal_radial,
                                       b.outward_normal_radial,places=7)


if __name__=="__main__":
    unittest.main()
