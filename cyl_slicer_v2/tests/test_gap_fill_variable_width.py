from __future__ import annotations

import math
import unittest

from shapely import affinity
from shapely.geometry import GeometryCollection, LineString, Polygon, box
from shapely.ops import unary_union

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _orient_path_for_current_position, emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerSlice, Toolpath
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


R = 30.0
EMPTY = GeometryCollection()


def layer(region):
    return LayerSlice(0,R,region,EMPTY,EMPTY,EMPTY,EMPTY)


def metric(line):
    return LineString([(p[0]*R,p[1]) for p in line.points])


def gcode_e(line):
    return float(line.split(" E")[1].split()[0])


class GapFillVariableWidthTests(unittest.TestCase):
    def setUp(self):
        self.c = SlicerConfig(wall_count=4, min_path_length=.1)

    def test_default_quality_enabled_and_legacy_explicitly_disabled(self):
        self.assertTrue(self.c.gap_fill_enabled)
        self.assertTrue(_profile_value(build_parser().parse_args([]),"gap_fill_enabled"))
        self.assertFalse(_profile_value(build_parser().parse_args(
            ["--profile","legacy-cli"]),"gap_fill_enabled"))
        self.assertFalse(_profile_value(build_parser().parse_args(
            ["--no-gap-fill-enabled"]),"gap_fill_enabled"))

    def test_isolated_sub_nominal_strip_recovers_material_without_wall_overlap(self):
        region=box(.1,0,.1+.30/R,10)
        roles=build_layer_toolpaths(layer(region),self.c).paths_by_role
        self.assertFalse(roles.get("EXTERNAL_PERIMETER"))
        self.assertFalse(roles.get("INTERNAL_PERIMETER"))
        strips=roles.get("GAP_FILL",[])
        self.assertTrue(strips)
        for path in strips:
            self.assertEqual(len(path.segment_widths),len(path.points)-1)
            self.assertTrue(all(.22 <= width <= .31 for width in path.segment_widths))
            for a,b,w in zip(path.points,path.points[1:],path.segment_widths):
                nominal=LineString([(a[0]*R,a[1]),(b[0]*R,b[1])]).buffer(
                    w/2,cap_style="flat")
                real=region
                from shapely import affinity
                actual=affinity.scale(real,xfact=R,yfact=1,origin=(0,0))
                self.assertLess(nominal.difference(actual).area, .01)

    def test_tapered_tip_receives_width_changes_without_printing_subminimum(self):
        geom=Polygon([(.1,0),(.1+.32/R,0),(.1+.32/R,7),(.1,10)])
        paths=build_layer_toolpaths(layer(geom),self.c).paths_by_role.get("GAP_FILL",[])
        self.assertTrue(paths)
        widths=[w for p in paths for w in p.segment_widths]
        self.assertGreater(max(widths)-min(widths),.025)
        self.assertTrue(all(self.c.gap_fill_min_width_mm <= w <= .42 for w in widths))

    def test_gently_curved_narrow_strip_recovers_local_width_without_overshoot(self):
        # The bounding rectangle is wider than a bead due to the bend,
        # while the actual strand stays about 0.30 mm across everywhere.
        spine = LineString([(3, 0), (3.1, 2), (3.45, 4), (3.55, 8)])
        physical = spine.buffer(.15, cap_style="flat")
        self.assertGreater(
            min(math.dist(a, b) for a, b in zip(
                list(physical.minimum_rotated_rectangle.exterior.coords)[:4],
                list(physical.minimum_rotated_rectangle.exterior.coords)[1:5]
            )), self.c.nozzle_width)
        region = affinity.scale(physical, xfact=1/R, yfact=1, origin=(0, 0))
        roles = build_layer_toolpaths(layer(region), self.c).paths_by_role
        self.assertFalse(roles.get("EXTERNAL_PERIMETER"))
        fills = roles.get("GAP_FILL", [])
        self.assertTrue(fills)
        self.assertGreater(sum(metric(p).length for p in fills), 3)
        for path in fills:
            self.assertEqual(len(path.segment_widths), len(path.points)-1)
            for a,b,w in zip(path.points, path.points[1:], path.segment_widths):
                self.assertGreaterEqual(w, self.c.gap_fill_min_width_mm)
                self.assertLessEqual(w, self.c.nozzle_width)
                swept = LineString([(a[0]*R,a[1]),(b[0]*R,b[1])]).buffer(
                    w*.5, cap_style="flat")
                self.assertLessEqual(swept.difference(physical).area,
                                     max(1e-8, swept.length*.006))

    def test_wide_bent_region_remains_ineligible(self):
        from cyl_slicer_v2.gap_fill import _narrow_polygon_paths
        spine = LineString([(3, 0), (3.1, 2), (3.45, 4), (3.55, 8)])
        wide = spine.buffer(.5, cap_style="flat")
        self.assertEqual(_narrow_polygon_paths(
            wide, R, self.c.nozzle_width,
            self.c.gap_fill_min_width_mm, self.c.gap_fill_sample_step_mm,
            self.c.gap_fill_min_length_mm,
        ), [])

    def test_regular_walls_and_sparse_interior_do_not_become_gap_fill(self):
        region=box(.1,0,1.1,20)
        roles=build_layer_toolpaths(layer(region),self.c).paths_by_role
        self.assertTrue(roles.get("EXTERNAL_PERIMETER"))
        self.assertFalse(roles.get("GAP_FILL"),"do not solidify intentionally sparse infill")

    def test_empty_hole_is_not_crossed_even_in_narrow_parts(self):
        region=unary_union([box(.1,0,.1+.3/R,10),
                            box(.1+1/R,0,.1+1.3/R,10)])
        roles=build_layer_toolpaths(layer(region),self.c).paths_by_role
        fill=roles.get("GAP_FILL",[])
        self.assertGreaterEqual(len(fill),2)
        hole=box(.1+.35/R,0,.1+.95/R,10)
        for p in fill:
            self.assertTrue(metric(p).disjoint(
                LineString([((.1+.5/R)*R,0),((.1+.5/R)*R,10)])))

    def test_periodic_theta_window_kept_without_spurious_second_traversal(self):
        region=box(2*math.pi-.08,0,2*math.pi+.08, .30)
        roles=build_layer_toolpaths(layer(region),self.c).paths_by_role
        fill=roles.get("GAP_FILL",[])
        # A narrow axial strip spanning theta=2pi must not be folded
        # into a disconnected fragment or duplicated by periodic copies.
        self.assertFalse(any(any(p[0]<0 for p in q.points) for q in fill))
        self.assertLessEqual(len(fill),1)

    def test_each_segment_e_is_scaled_by_its_own_rounded_bead_area(self):
        widths=[.24,.36,.30]
        path=Toolpath("GAP_FILL",R,[(.2,0),(.2,2),(.2,4),(.2,6)],
                      segment_widths=widths)
        config=SlicerConfig(retract_mm=0,z_hop_height=0,cooling_min_layer_time_s=0)
        state=GcodeState(theta_accum=.2,axial=0,radial=R+.08)
        commands=emit_toolpath(path,state,config,0,retract_after=False)
        moves=[line for line in commands if line.startswith("G1 X")]
        self.assertEqual(len(moves),3)
        area=lambda w:config.layer_height*(w-config.layer_height)+math.pi*config.layer_height**2/4
        expected=sum(2*area(w)/config.bead_area*config.extrusion_per_mm for w in widths)
        self.assertAlmostEqual(gcode_e(moves[-1]),expected,delta=1e-5)

    def test_no_implicit_path_reversal_can_detach_segment_widths(self):
        path=Toolpath("GAP_FILL",R,[(.1,0),(.1,2),(.1,4)],
                      segment_widths=[.22,.41])
        state=GcodeState(theta_accum=.1,axial=4.1)
        self.assertIs(_orient_path_for_current_position(path,state,self.c),path)

    def test_invalid_gap_parameters_or_width_metadata_rejected(self):
        for vals in ({"gap_fill_min_width_mm":0},
                     {"gap_fill_min_width_mm":.1},
                     {"gap_fill_min_width_mm":.6},
                     {"gap_fill_sample_step_mm":0},
                     {"gap_fill_min_length_mm":-1}):
            with self.subTest(vals=vals),self.assertRaises(ValueError):
                SlicerConfig(**vals).validate()
        path=Toolpath("GAP_FILL",R,[(.2,0),(.2,1),(.2,2)],segment_widths=[.3])
        with self.assertRaisesRegex(ValueError,"Variable-width"):
            emit_toolpath(path,GcodeState(),self.c,0,retract_after=False)


if __name__=="__main__":
    unittest.main()
