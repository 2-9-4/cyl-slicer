from __future__ import annotations

import math
import unittest

from cyl_slicer_v2.nonplanar_patches import scan_surface_patches
from cyl_slicer_v2.nonplanar_raster import propose_axial_conformal_rasters
from cyl_slicer_v2.tests.test_nonplanar_surface import cylindrical_shell


TAU=2*math.pi


def scan(mesh, *, theta_min=0.,theta_max=TAU,theta_cells=12,
         axial_cells=6,axial_min=1.,axial_max=7.):
    return scan_surface_patches(mesh,theta_min=theta_min,theta_max=theta_max,
             axial_min=axial_min,axial_max=axial_max,
             theta_cells=theta_cells,axial_cells=axial_cells,
             max_chord_error_mm=.05)


def rasters(mesh, report, **kwargs):
    settings=dict(nominal_surface_pitch_mm=4.,
                  edge_inset_mm=.1,min_path_length_mm=.5,
                  max_segment_mm=.5,max_radial_chord_error_mm=.02)
    settings.update(kwargs)
    return propose_axial_conformal_rasters(mesh,report,**settings)


class CandidateAxialRasterTests(unittest.TestCase):
    def test_full_turn_has_registered_tracks_without_duplicate_seam(self):
        mesh=cylindrical_shell(lambda t,x:30,theta_segments=96)
        report=scan(mesh)
        result=rasters(mesh,report)
        self.assertEqual(result["status"],"candidate_axial_rasters_only")
        self.assertFalse(result["nonplanar_extrusion_enabled"])
        self.assertFalse(result["actual_surface_coverage_certified"])
        self.assertFalse(result["full_head_collision_certified"])
        self.assertEqual(result["candidate_count"],
                         result["global_angular_track_count"])
        self.assertGreater(result["candidate_count"],40)
        paths=result["paths"]
        self.assertEqual(len({p["theta_deg"] for p in paths}),len(paths))
        self.assertGreater(paths[0]["theta_deg"],0)
        self.assertLess(paths[-1]["theta_deg"],360)
        self.assertTrue(all(p["axial_row_span"]==[0,6] for p in paths))
        self.assertTrue(all(p["physical_length_mm"]>5 for p in paths))
        differences=[b["theta_deg"]-a["theta_deg"]
                     for a,b in zip(paths,paths[1:])]
        self.assertAlmostEqual(min(differences),max(differences),delta=1e-10)
        self.assertLess(result["projected_angular_pitch_at_outer_radius_mm"],
                        4 / math.sqrt(1+.3**2)+1e-9)

    def test_curved_cylindrical_mesh_produces_true_varying_radius_ribbons(self):
        mesh=cylindrical_shell(lambda t,x:30+.1*x+.05*math.sin(t),
                               theta_segments=96)
        report=scan(mesh)
        result=rasters(mesh,report,include_points=True)
        self.assertGreater(result["candidate_count"],40)
        self.assertGreater(result["total_centerline_length_mm"],200)
        for path in result["paths"]:
            points=path["points"]
            self.assertEqual(len(points),path["point_count"])
            self.assertGreater(points[-1]["radius_mm"]-points[0]["radius_mm"],.3)
            self.assertTrue(all(abs(p["theta_deg"]-path["theta_deg"])<1e-8
                                for p in points))
            self.assertAlmostEqual(points[0]["axial_mm"],
                                   path["axial_range_mm"][0],places=8)

    def test_mask_holes_split_runs_and_no_candidate_crosses_rejected_cells(self):
        mesh=cylindrical_shell(
            lambda t,x:30,theta_segments=240,
            skip=lambda t,x: .4<t<.6 and 2<x<5)
        report=scan(mesh,theta_cells=24,axial_cells=12)
        self.assertGreater(report["rejected_cells"],0)
        result=rasters(mesh,report,nominal_surface_pitch_mm=2.0)
        self.assertGreater(result["candidate_count"],0)
        spans_by_col={}
        for path in result["paths"]:
            col=path["scan_column"]
            row0,row1=path["axial_row_span"]
            self.assertTrue(all(report["candidate_mask"][j][col]=="1"
                                for j in range(row0,row1)))
            spans_by_col.setdefault(col,[]).append((row0,row1))
        # A thin hole cannot be traversed by a single drawn axial raster.
        self.assertTrue(any(len(spans)>1 for spans in spans_by_col.values()))

    def test_shift_by_full_turn_only_changes_angles_not_paths(self):
        mesh=cylindrical_shell(lambda t,x:30+.07*x,theta_segments=96)
        initial=rasters(mesh,scan(mesh))
        shifted=rasters(mesh,scan(mesh,theta_min=TAU,theta_max=2*TAU))
        self.assertEqual(initial["candidate_count"],shifted["candidate_count"])
        self.assertEqual(initial["global_angular_track_count"],
                         shifted["global_angular_track_count"])
        for a,b in zip(initial["paths"],shifted["paths"]):
            self.assertAlmostEqual(b["theta_deg"]-a["theta_deg"],360,delta=1e-9)
            self.assertEqual(a["axial_row_span"],b["axial_row_span"])
            self.assertAlmostEqual(a["physical_length_mm"],
                                   b["physical_length_mm"],delta=1e-8)

    def test_real_path_and_point_budgets_are_enforced(self):
        mesh=cylindrical_shell(lambda t,x:30,theta_segments=96)
        report=scan(mesh)
        with self.assertRaisesRegex(ValueError,"path budget"):
            rasters(mesh,report,max_paths=1)
        with self.assertRaisesRegex(ValueError,"point budget"):
            rasters(mesh,report,max_total_points=10)
        with self.assertRaisesRegex(ValueError,"track budget"):
            rasters(mesh,report,nominal_surface_pitch_mm=.0005)
        for opts in (
            {"nominal_surface_pitch_mm":0},
            {"allowed_surface_slope":-1},
            {"edge_inset_mm":-1},
            {"max_paths":0},
            {"max_total_points":0},
        ):
            with self.subTest(opts=opts),self.assertRaises(ValueError):
                propose_axial_conformal_rasters(mesh,report,**opts)

    def test_scan_metadata_must_represent_a_valid_analysis_only_mask(self):
        mesh=cylindrical_shell(lambda t,x:30)
        report=scan(mesh)
        for change in (
            {"status":"geometry_candidate_only"},
            {"nonplanar_extrusion_enabled":True},
            {"candidate_mask":["1"]},
            {"theta_range_deg":[0,720]},
            {"axial_range_mm":[7,1]},
            {"theta_cells":-1},
        ):
            with self.subTest(change=change),self.assertRaises(ValueError):
                propose_axial_conformal_rasters(mesh,{**report,**change})

    def test_empty_safe_region_returns_no_paths_not_a_false_surface(self):
        mesh=cylindrical_shell(lambda t,x:30+.8*x,theta_segments=96)
        report=scan(mesh)
        self.assertEqual(report["accepted_cells"],0)
        result=rasters(mesh,report)
        self.assertEqual(result["candidate_count"],0)
        self.assertEqual(result["paths"],[])
        self.assertFalse(result["nonplanar_extrusion_enabled"])


if __name__=="__main__":
    unittest.main()
