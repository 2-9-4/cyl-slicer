from __future__ import annotations

import math
import unittest

from cyl_slicer_v2.nonplanar_probe import analyze_ribbon, degree_point
from cyl_slicer_v2.tests.test_nonplanar_surface import cylindrical_shell


class NonplanarProbeTests(unittest.TestCase):
    def test_degrees_and_axial_input(self):
        theta, axial = degree_point("370,19.5")
        self.assertAlmostEqual(theta, math.radians(370))
        self.assertEqual(axial,19.5)
        for value in ("1", "1,2,3", "nan,20", "30,inf", "not_a_number,5"):
            with self.subTest(value=value), self.assertRaises(ValueError):
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

    def test_rejected_candidates_have_explicit_reason(self):
        mesh=cylindrical_shell(lambda theta,axial: 30+.4*axial)
        result=analyze_ribbon(mesh,[degree_point("10,1"),degree_point("10,6")],
                              max_slope=.1)
        self.assertEqual(result["status"],"rejected")
        self.assertIn("slope",result["reason"])
        self.assertFalse(result["nonplanar_extrusion_enabled"])


if __name__=="__main__":
    unittest.main()
