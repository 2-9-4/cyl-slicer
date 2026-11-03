from __future__ import annotations

import math
import unittest

from cyl_slicer_v2.nonplanar_motion import evaluate_conformal_move
from cyl_slicer_v2.nonplanar_plan import plan_conformal_motion
from cyl_slicer_v2.nonplanar_surface import RadialSurfaceSample


def point(theta, axial, radius):
    return RadialSurfaceSample(theta, axial, radius, 0, 1.0)


def plan(samples, **kw):
    args=dict(physical_feed_mm_per_min=2700,bead_area_mm2=.0672,
              filament_diameter_mm=1.75,max_local_volume_error_fraction=.01)
    args.update(kw)
    return plan_conformal_motion(samples,**args)


def unsplit_length_e(a,b):
    return evaluate_conformal_move(
        (a.theta,a.axial,a.radius),(b.theta,b.axial,b.radius),
        physical_print_feed_mm_per_min=2700,bead_area_mm2=.0672,
        filament_diameter_mm=1.75)


class ConformalSegmentationTests(unittest.TestCase):
    def test_constant_radius_needs_no_e_split_even_for_rotary(self):
        a,b=point(.1,2,30),point(.6,4,30)
        result=plan([a,b])
        self.assertEqual(result.split_count,0)
        self.assertEqual(len(result.segments),1)
        self.assertAlmostEqual(result.filament_length_mm,
                               unsplit_length_e(a,b).filament_mm,places=11)
        self.assertAlmostEqual(result.physical_length_mm,
                               unsplit_length_e(a,b).length_mm,places=11)

    def test_varying_radius_subdivides_without_changing_exact_3d_geometry_or_e(self):
        a,b=point(.1,2,30),point(.8,5,43)
        original=unsplit_length_e(a,b)
        result=plan([a,b])
        self.assertTrue(original.needs_shorter_segments)
        self.assertGreater(result.split_count,4)
        self.assertEqual(len(result.segments),result.split_count+1)
        self.assertLessEqual(result.largest_volume_error_fraction,.0100000000001)
        self.assertAlmostEqual(result.filament_length_mm,
                               original.filament_mm,delta=1e-10)
        self.assertAlmostEqual(result.physical_length_mm,
                               original.length_mm,delta=1e-9)
        self.assertEqual(result.segments[0].start,(a.theta,a.axial,a.radius))
        self.assertEqual(result.segments[-1].end,(b.theta,b.axial,b.radius))
        for i,seg in enumerate(result.segments):
            self.assertFalse(seg.metrics.needs_shorter_segments)
            self.assertLessEqual(seg.metrics.max_physical_speed_mm_per_min,
                                 2700+1e-9)
            self.assertGreater(seg.metrics.filament_mm,0)
            if i:
                prev=result.segments[i-1]
                self.assertEqual(prev.end,seg.start)
                self.assertAlmostEqual(prev.cumulative_e_end_mm,
                                       seg.cumulative_e_start_mm,places=12)
            # Every inserted point stays on the ORIGINAL machine-linear
            # motion, not on a resampled/altered radial curve.
            fraction=(seg.end[0]-a.theta)/(b.theta-a.theta)
            self.assertAlmostEqual(seg.end[1],a.axial+fraction*(b.axial-a.axial),
                                   delta=1e-10)
            self.assertAlmostEqual(seg.end[2],a.radius+fraction*(b.radius-a.radius),
                                   delta=1e-10)

    def test_cumulative_e_across_multiple_original_segments(self):
        points=[point(.1,2,30),point(.3,3,36),point(.4,5,37)]
        expected=sum(unsplit_length_e(a,b).filament_mm
                     for a,b in zip(points,points[1:]))
        result=plan(points)
        self.assertEqual(result.input_segment_count,2)
        self.assertGreater(result.split_count,0)
        self.assertAlmostEqual(result.filament_length_mm,expected,delta=1e-10)
        self.assertAlmostEqual(result.segments[0].cumulative_e_start_mm,0)
        self.assertAlmostEqual(result.segments[-1].cumulative_e_end_mm,
                               result.filament_length_mm,places=12)

    def test_extreme_flow_cap_preserves_e_and_reduces_feed(self):
        a,b=point(.05,2,32),point(.4,3,36)
        baseline=plan([a,b])
        cap=plan([a,b],max_volumetric_flow_mm3_s=.15)
        self.assertAlmostEqual(baseline.filament_length_mm,
                               cap.filament_length_mm,delta=1e-10)
        self.assertEqual(baseline.split_count,cap.split_count)
        self.assertTrue(all(seg.metrics.volumetric_rate_mm3_s<=.15000001
                            for seg in cap.segments))
        self.assertLess(cap.segments[0].metrics.marlin_feed_mm_per_min,
                        baseline.segments[0].metrics.marlin_feed_mm_per_min)

    def test_winding_and_reverse_are_physically_equivalent(self):
        a,b=point(.1,2,30),point(.5,5,35)
        base=plan([a,b])
        for turns in (-4,-1,2,5):
            shifted=plan([point(a.theta+2*math.pi*turns,a.axial,a.radius),
                          point(b.theta+2*math.pi*turns,b.axial,b.radius)])
            self.assertEqual(len(base.segments),len(shifted.segments))
            self.assertAlmostEqual(base.filament_length_mm,
                                   shifted.filament_length_mm,delta=1e-10)
        reverse=plan([b,a])
        self.assertAlmostEqual(base.filament_length_mm,
                               reverse.filament_length_mm,delta=1e-10)

    def test_recursion_budgets_and_invalid_input_fail_closed(self):
        a,b=point(0,1,30),point(.8,4,43)
        for args in (
            {"max_segments":1},
            {"max_segments":2},
            {"max_split_depth":0},
            {"max_segments":0},
            {"max_split_depth":-1},
            {"max_local_volume_error_fraction":0},
        ):
            # zero local E error is mathematically valid but can be
            # unachievable within bounded depth; must reject explicitly.
            with self.subTest(args=args),self.assertRaises(ValueError):
                plan([a,b],**args)
        with self.assertRaises(ValueError):
            plan([a])
        with self.assertRaises(ValueError):
            plan([a,a])
        with self.assertRaises(ValueError):
            plan([a,point(.1,2,-1)])
        with self.assertRaises(ValueError):
            plan([a,point(float("nan"),2,30)])


if __name__=="__main__":
    unittest.main()
