from __future__ import annotations

import math
import unittest

from scipy.integrate import quad

from cyl_slicer_v2.nonplanar_motion import evaluate_conformal_move


def evaluate(a,b,**kwargs):
    return evaluate_conformal_move(
        a,b,physical_print_feed_mm_per_min=2700,
        bead_area_mm2=.056,filament_diameter_mm=1.75,**kwargs)


class NonplanarMotionTests(unittest.TestCase):
    def test_constant_radius_matches_existing_cylindrical_length_and_feed(self):
        start=(.2,5,30)
        end=(.42,8,30)
        result=evaluate(start,end)
        true=math.hypot(3,30*.22)
        native=math.hypot(3,math.degrees(.22))
        self.assertAlmostEqual(result.length_mm,true,places=10)
        self.assertAlmostEqual(result.marlin_native_length,native,places=10)
        self.assertAlmostEqual(result.marlin_feed_mm_per_min,2700*native/true,places=9)
        self.assertAlmostEqual(result.max_physical_speed_mm_per_min,2700,places=8)
        self.assertAlmostEqual(result.min_physical_speed_mm_per_min,2700,places=8)
        self.assertFalse(result.needs_shorter_segments)
        self.assertAlmostEqual(result.max_local_volume_error_fraction,0,places=12)

    def test_radial_only_uses_exact_physical_length(self):
        result=evaluate((0,0,30),(0,0,30.5))
        self.assertAlmostEqual(result.length_mm,.5,places=12)
        self.assertAlmostEqual(result.marlin_native_length,.5,places=12)
        self.assertAlmostEqual(result.marlin_feed_mm_per_min,2700,places=12)
        self.assertAlmostEqual(result.filament_mm,
            .056*.5/(math.pi*(1.75*.5)**2),places=12)

    def test_true_polar_xyz_arc_matches_independent_integration(self):
        start=(1.2,3,30)
        end=(1.65,6,32)
        result=evaluate(start,end,max_local_volume_error_fraction=.8)
        dtheta=end[0]-start[0]
        dx=end[1]-start[1]
        dr=end[2]-start[2]
        expected=quad(lambda t: math.sqrt(dx*dx+dr*dr
                            +(dtheta*(start[2]+dr*t))**2),0,1,
                      epsabs=1e-11,epsrel=1e-11)[0]
        self.assertAlmostEqual(result.length_mm,expected,delta=1e-9)
        self.assertGreater(result.length_mm,
                           math.hypot(dx,(start[2]+end[2])*.5*dtheta))
        self.assertAlmostEqual(result.filament_mm,
                           .056*result.length_mm/(math.pi*(1.75*.5)**2),
                           places=10)
        self.assertLessEqual(result.max_physical_speed_mm_per_min,2700+1e-7)

    def test_high_polar_rate_variation_requires_subdivision(self):
        a=(0,0,30)
        b=(.7,1,42)
        result=evaluate(a,b)
        self.assertTrue(result.needs_shorter_segments)
        self.assertGreater(result.max_local_volume_error_fraction,.01)
        self.assertLessEqual(result.max_physical_speed_mm_per_min,2700+1e-7)
        self.assertGreater(result.max_physical_speed_mm_per_min,
                           result.min_physical_speed_mm_per_min)
        short=evaluate((0,0,30),(.002,.2,30.02))
        self.assertFalse(short.needs_shorter_segments)

    def test_positive_negative_winding_and_turn_offsets_keep_motion_identical(self):
        a=(.1,2,31)
        b=(.3,3,31.2)
        base=evaluate(a,b)
        for turns in (-5,-2,1,7):
            shifted=evaluate((a[0]+turns*2*math.pi,a[1],a[2]),
                             (b[0]+turns*2*math.pi,b[1],b[2]))
            self.assertAlmostEqual(base.length_mm,shifted.length_mm,places=10)
            self.assertAlmostEqual(base.filament_mm,shifted.filament_mm,places=10)
            self.assertAlmostEqual(base.marlin_feed_mm_per_min,
                                   shifted.marlin_feed_mm_per_min,places=9)
        reverse=evaluate(b,a)
        self.assertAlmostEqual(base.length_mm,reverse.length_mm,places=10)
        self.assertAlmostEqual(base.filament_mm,reverse.filament_mm,places=10)

    def test_flow_limit_reduces_native_feed_without_altering_e(self):
        a=(.1,5,30)
        b=(.2,7,30.6)
        unlimited=evaluate(a,b)
        limited=evaluate(a,b,max_volumetric_flow_mm3_s=.1)
        self.assertEqual(unlimited.filament_mm,limited.filament_mm)
        self.assertLess(limited.marlin_feed_mm_per_min,
                        unlimited.marlin_feed_mm_per_min)
        self.assertLessEqual(limited.volumetric_rate_mm3_s,.100000001)
        self.assertLess(limited.max_physical_speed_mm_per_min,2700)
        self.assertEqual(limited.needs_shorter_segments,
                         unlimited.needs_shorter_segments)

    def test_machine_angle_degrees_not_physical_radian_native(self):
        result=evaluate((0,0,40),(math.pi,0,40))
        self.assertAlmostEqual(result.marlin_native_length,180,places=10)
        self.assertAlmostEqual(result.length_mm,40*math.pi,places=10)
        self.assertAlmostEqual(result.marlin_feed_mm_per_min,
                               2700*180/(40*math.pi),places=8)

    def test_bad_input_fails_closed(self):
        a=(.1,1,30)
        b=(.2,2,31)
        cases=[
            ((.1,1,0),b,{}),
            ((.1,1,-30),b,{}),
            ((.1,float("nan"),30),b,{}),
            (a,a,{}),
            (a,b,{"physical_print_feed_mm_per_min":0}),
            (a,b,{"bead_area_mm2":0}),
            (a,b,{"filament_diameter_mm":0}),
            (a,b,{"max_local_volume_error_fraction":1.1}),
            (a,b,{"max_volumetric_flow_mm3_s":-1}),
            (a,b,{"max_volumetric_flow_mm3_s":float("inf")}),
        ]
        for p,q,values in cases:
            kw={"physical_print_feed_mm_per_min":2700,
                "bead_area_mm2":.056,
                "filament_diameter_mm":1.75}
            kw.update(values)
            with self.subTest(args=(p,q,values)),self.assertRaises(ValueError):
                evaluate_conformal_move(p,q,**kw)


if __name__=="__main__":
    unittest.main()
