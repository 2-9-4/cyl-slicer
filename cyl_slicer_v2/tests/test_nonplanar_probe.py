from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import io
import json
import math
from unittest.mock import patch
import unittest

from cyl_slicer_v2.nonplanar_probe import analyze_ribbon, degree_point, main
from cyl_slicer_v2.tests.test_nonplanar_surface import cylindrical_shell


class NonplanarProbeTests(unittest.TestCase):
    def test_degrees_and_axial_input(self):
        theta, axial = degree_point("370,19.5")
        self.assertAlmostEqual(theta, math.radians(370))
        self.assertEqual(axial,19.5)
        for value in ("1", "1,2,3", "nan,20", "30,inf", "not_a_number,5"):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                degree_point(value)

    def test_candidate_json_does_not_claim_printable_or_change_gcode(self):
        mesh=cylindrical_shell(lambda theta,axial: 30+axial*.07)
        result=analyze_ribbon(mesh, [degree_point("14,2"),degree_point("20,6")])
        self.assertEqual(result["status"],"geometry_candidate_only")
        self.assertFalse(result["nonplanar_extrusion_enabled"])
        self.assertGreater(result["path_length_3d_mm"],4)
        self.assertGreater(result["radial_excursion_mm"],.2)
        self.assertTrue(all("radius_mm" in p for p in result["samples"]))
        self.assertLessEqual(result["max_segment_3d_mm"],.50001)

    def test_automatic_grid_cli_is_analysis_only_and_does_not_need_points(self):
        mesh=cylindrical_shell(lambda t,x:30+.05*x,theta_segments=96)
        output=io.StringIO()
        with patch("cyl_slicer_v2.nonplanar_probe.load_mesh",return_value=mesh):
            with redirect_stdout(output):
                code=main(["unused.stl","--scan",
                    "--axial-range-mm","1","7",
                    "--theta-cells","16","--axial-cells","4",
                    "--max-chord-error-mm",".06"])
        self.assertEqual(code,0)
        result=json.loads(output.getvalue())
        self.assertEqual(result["status"],"geometry_candidates_only")
        self.assertFalse(result["nonplanar_extrusion_enabled"])
        self.assertFalse(result["collision_certified"])
        self.assertGreater(result["accepted_cells"],0)
        self.assertTrue(result["patches"])

    def test_scan_and_explicit_ribbon_are_mutually_exclusive(self):
        with self.assertRaises(SystemExit):
            main(["unused.stl","--scan","--point","20,1","--point","30,2"])

    def test_motion_plan_json_has_bounded_e_and_native_y_degrees(self):
        mesh=cylindrical_shell(lambda theta,axial:30+.25*axial,
                               theta_segments=120)
        output=io.StringIO()
        with patch("cyl_slicer_v2.nonplanar_probe.load_mesh",return_value=mesh):
            with redirect_stdout(output):
                code=main(["unused.stl",
                           "--point","12,2","--point","48,5",
                           "--max-segment-mm","50",
                           "--max-excursion-mm","2",
                           "--motion-plan","--bead-area-mm2",".0672",
                           "--motion-max-e-error-fraction",".0001"])
        self.assertEqual(code,0)
        result=json.loads(output.getvalue())
        self.assertEqual(result["status"],"geometry_candidate_only")
        self.assertFalse(result["nonplanar_extrusion_enabled"])
        motion=result["motion_plan"]
        self.assertEqual(motion["status"],"analysis_only_no_gcode")
        self.assertEqual(motion["native_marlin_y_units"],"degrees")
        self.assertGreater(motion["split_count"],0)
        self.assertGreater(motion["filament_length_mm"],0)
        self.assertLessEqual(motion["max_local_e_error_fraction"],.00010000001)
        self.assertAlmostEqual(motion["segments"][-1]["cumulative_filament_mm"],
                               motion["filament_length_mm"],places=10)
        self.assertTrue(all(seg["marlin_feed_mm_per_min"]>0
                            for seg in motion["segments"]))

    def test_motion_plan_requires_explicit_bead_geometry_and_ribbon(self):
        for args in (
            ["--point","12,2","--point","18,3","--motion-plan"],
            ["--scan","--motion-plan","--bead-area-mm2",".07"],
            ["--point","12,2","--point","18,3","--bead-area-mm2",".07"],
        ):
            with self.subTest(args=args), self.assertRaises(SystemExit):
                main(["unused.stl",*args])

    def test_rejected_candidates_have_explicit_reason(self):
        mesh=cylindrical_shell(lambda theta,axial: 30+.4*axial)
        result=analyze_ribbon(mesh,[degree_point("10,1"),degree_point("10,6")],
                              max_slope=.1)
        self.assertEqual(result["status"],"rejected")
        self.assertIn("slope",result["reason"])
        self.assertFalse(result["nonplanar_extrusion_enabled"])


if __name__=="__main__":
    unittest.main()
