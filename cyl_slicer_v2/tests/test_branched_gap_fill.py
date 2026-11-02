from __future__ import annotations

from dataclasses import replace
import math
import unittest

from shapely import affinity
from shapely.geometry import GeometryCollection, LineString, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths
from cyl_slicer_v2.branch_gap_fill import branched_gap_paths


R=30
EMPTY=GeometryCollection()


def shape():
    return unary_union([box(5.85,0,6.15,8),box(4,7.85,8,8.15)])


def paths(physical, config=None):
    c=config or SlicerConfig(min_path_length=.1)
    unwrapped=affinity.scale(physical,xfact=1/R,yfact=1,origin=(0,0))
    layer=LayerSlice(0,R,unwrapped,EMPTY,EMPTY,EMPTY,EMPTY)
    return build_layer_toolpaths(layer,c).paths_by_role


def deposited_footprints(items):
    results=[]
    for path in items:
        parts=[
            LineString([(a[0]*R,a[1]),(b[0]*R,b[1])]).buffer(w*.5,cap_style="flat")
            for a,b,w in zip(path.points,path.points[1:],path.segment_widths)
        ]
        results.append(unary_union(parts))
    return results


class BranchGapFillTests(unittest.TestCase):
    def test_t_junction_recovers_all_arms_without_overlap_or_overshoot(self):
        original=shape()
        roles=paths(original)
        self.assertFalse(roles.get("EXTERNAL_PERIMETER"))
        self.assertFalse(roles.get("INTERNAL_PERIMETER"))
        result=roles.get("GAP_FILL",[])
        self.assertEqual(len(result),2)
        sweeps=deposited_footprints(result)
        for path,footprint in zip(result,sweeps):
            self.assertEqual(len(path.segment_widths),len(path.points)-1)
            self.assertTrue(all(.22<=width<=.42 for width in path.segment_widths))
            self.assertLessEqual(footprint.difference(original).area,.012)
        self.assertLessEqual(sweeps[0].intersection(sweeps[1]).area,.0011)
        self.assertGreater(sum(footprint.area for footprint in sweeps),2.7)
        self.assertGreater(max(LineString([(theta*R,axial) for theta,axial in p.points]).length
                               for p in result),9)

    def test_variable_width_e_is_preserved_through_emission(self):
        c=SlicerConfig(retract_mm=0,z_hop_height=0,cooling_min_layer_time_s=0)
        for path in paths(shape())["GAP_FILL"]:
            state=GcodeState(theta_accum=path.points[0][0],
                             axial=path.points[0][1],radial=R+c.layer_height/2)
            commands=emit_toolpath(path,state,c,0,retract_after=False)
            moves=[line for line in commands if line.startswith("G1 X") and " E" in line]
            self.assertTrue(moves)
            actual=float(moves[-1].split(" E",1)[1].split()[0])
            expected=0
            for a,b,w in zip(path.points,path.points[1:],path.segment_widths):
                area=c.layer_height*(w-c.layer_height)+math.pi*c.layer_height**2/4
                expected+=math.hypot((b[0]-a[0])*R,b[1]-a[1])*area/c.bead_area*c.extrusion_per_mm
            self.assertAlmostEqual(actual,expected,delta=.0002)

    def test_periodic_theta_copies_have_same_safe_two_arm_layout(self):
        baseline=paths(shape())["GAP_FILL"]
        length=sum(LineString([(t*R,x) for t,x in p.points]).length for p in baseline)
        for shift in (-2*math.pi,2*math.pi):
            translated=affinity.translate(shape(),xoff=shift*R)
            result=paths(translated)["GAP_FILL"]
            self.assertEqual(len(result),2)
            self.assertAlmostEqual(sum(
                LineString([(t*R,x) for t,x in p.points]).length for p in result
            ),length,delta=.02)
            sweeps=deposited_footprints(result)
            self.assertLessEqual(sweeps[0].intersection(sweeps[1]).area,.0011)

    def test_four_arm_cross_uses_one_trunk_and_two_nonoverlapping_arms(self):
        cross=unary_union([box(5.85,0,6.15,8),box(3,3.85,9,4.15)])
        roles=paths(cross)
        self.assertFalse(roles.get("EXTERNAL_PERIMETER"),
                         "a microscopic eroded junction is not a printable wall")
        self.assertFalse(roles.get("INTERNAL_PERIMETER"))
        result=roles.get("GAP_FILL",[])
        self.assertEqual(len(result),3)
        sweeps=deposited_footprints(result)
        self.assertGreater(sum(footprint.area for footprint in sweeps),3.2)
        for p,footprint in zip(result,sweeps):
            self.assertEqual(len(p.points)-1,len(p.segment_widths))
            self.assertTrue(all(.22<=w<=.42 for w in p.segment_widths))
            self.assertLess(footprint.difference(cross).area,.015)
        for i in range(3):
            for j in range(i+1,3):
                self.assertLessEqual(sweeps[i].intersection(sweeps[j]).area,.0011)
        for shift in (-2*math.pi,2*math.pi):
            other=paths(affinity.translate(cross,xoff=shift*R)).get("GAP_FILL",[])
            self.assertEqual(len(other),3)
            for a,b in zip(result,other):
                self.assertAlmostEqual(
                    LineString([(t*R,x) for t,x in a.points]).length,
                    LineString([(t*R,x) for t,x in b.points]).length,delta=.02)

    def test_four_arm_extrusion_matches_individual_bead_area(self):
        cross=unary_union([box(5.85,0,6.15,8),box(3,3.85,9,4.15)])
        c=SlicerConfig(retract_mm=0,z_hop_height=0,cooling_min_layer_time_s=0)
        for path in paths(cross,c).get("GAP_FILL",[]):
            state=GcodeState(theta_accum=path.points[0][0],axial=path.points[0][1],
                             radial=R+c.layer_height/2)
            moves=[line for line in emit_toolpath(path,state,c,0,retract_after=False)
                   if line.startswith("G1 X") and " E" in line]
            self.assertTrue(moves)
            actual=float(moves[-1].split(" E",1)[1].split()[0])
            expected=sum(math.hypot((b[0]-a[0])*R,b[1]-a[1]) *
                         (c.layer_height*(w-c.layer_height)+math.pi*c.layer_height**2/4)
                         /c.bead_area*c.extrusion_per_mm
                         for a,b,w in zip(path.points,path.points[1:],path.segment_widths))
            self.assertAlmostEqual(actual,expected,delta=.0002)

    def test_multi_junction_and_holed_regions_are_not_fabricated(self):
        c=SlicerConfig()
        multiple=unary_union([shape(),box(5.85,3,8,3.3)])
        def directly(geom):
            return branched_gap_paths(geom,R,c.nozzle_width,c.gap_fill_min_width_mm,
                                      c.gap_fill_sample_step_mm,c.gap_fill_min_length_mm)
        self.assertEqual(directly(multiple),[])
        holed=box(3,0,9,9).difference(box(4,2,8,7))
        self.assertEqual(directly(holed),[])

    def test_opt_out_disables_recovery(self):
        c=replace(SlicerConfig(),gap_fill_enabled=False)
        self.assertFalse(paths(shape(),c).get("GAP_FILL"))


if __name__=="__main__":
    unittest.main()
