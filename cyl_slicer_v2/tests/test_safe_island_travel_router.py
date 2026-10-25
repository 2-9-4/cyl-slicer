from __future__ import annotations

import math
import unittest

from shapely.geometry import LineString, Point, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (
    _safe_island_travel_route, _should_retract_after_path, emit_toolpath,
)
from cyl_slicer_v2.models import GcodeState, Toolpath
from cyl_slicer_v2.periodic import periodic_copies


def with_hole():
    return box(0, 0, 1, 10).difference(box(.4, 4, .6, 6))


def physical(route, start, radius):
    points = [start, *route]
    return sum(math.hypot((b[0]-a[0])*radius,b[1]-a[1])
               for a,b in zip(points, points[1:]))


class SafeIslandRoutingTests(unittest.TestCase):
    def setUp(self):
        self.cfg = SlicerConfig(route_max_unretracted_mm=12)
        self.region = periodic_copies(with_hole())
        self.start = (.2, 5)
        self.end = (.8, 5)

    def test_hole_detour_stays_inside_physical_island(self):
        route = _safe_island_travel_route(self.start, self.end, 10,
                                          self.region, self.cfg)
        self.assertIsNotNone(route)
        self.assertGreater(len(route), 1)
        self.assertEqual(route[-1], self.end)
        for a,b in zip([self.start,*route],route):
            self.assertTrue(self.region.buffer(1e-8).covers(LineString([a,b])),
                            msg=f"travel not contained: {a} -> {b}")
        direct = math.hypot((self.end[0]-self.start[0])*10, 0)
        self.assertGreater(physical(route,self.start,10), direct)
        self.assertLessEqual(physical(route,self.start,10), direct*1.65+1e-9)

    def test_concave_outer_edge_detour_is_contained(self):
        region = box(0,0,1,10).difference(box(.4,5,.6,10.1))
        start,end = (.2,7),(.8,7)
        route = _safe_island_travel_route(start,end,10,periodic_copies(region),
                                          SlicerConfig(route_max_extra_mm=10,
                                                       route_max_factor=3))
        self.assertIsNotNone(route)
        for a,b in zip([start,*route],route):
            self.assertTrue(periodic_copies(region).buffer(1e-8).covers(LineString([a,b])))

    def test_periodic_seam_and_large_accumulated_phase(self):
        region = box(6.0,0,6.6,10).difference(box(6.2,4,6.4,6))
        periodic = periodic_copies(region)
        for revolutions in (0, 1, 7, -3):
            shift = revolutions * 2*math.pi
            start,end=(6.05+shift,5),(6.55+shift,5)
            route = _safe_island_travel_route(start,end,10,periodic,self.cfg)
            self.assertIsNotNone(route)
            self.assertAlmostEqual(route[-1][0],end[0],places=8)
            for a,b in zip([start,*route],route):
                original = LineString([(a[0]-shift,a[1]),(b[0]-shift,b[1])])
                self.assertTrue(periodic.buffer(1e-8).covers(original))

    def test_disconnected_islands_never_get_routed(self):
        region = unary_union([box(0,0,.2,10), box(.8,0,1,10)])
        self.assertIsNone(_safe_island_travel_route((.1,5),(.9,5),10,
                           periodic_copies(region),self.cfg))

    def test_bounded_length_and_vertex_count_fall_back(self):
        c = SlicerConfig(route_max_factor=1.01,route_max_extra_mm=.05)
        self.assertIsNone(_safe_island_travel_route(self.start,self.end,10,
                           self.region,c))
        c = SlicerConfig(route_max_vertices=4)
        self.assertIsNone(_safe_island_travel_route(self.start,self.end,10,
                           self.region,c))
        c = SlicerConfig(route_travel_inside_island=False)
        self.assertIsNone(_safe_island_travel_route(self.start,self.end,10,
                           self.region,c))

    def test_existing_safe_straight_travel_is_not_replaced(self):
        self.assertIsNone(_safe_island_travel_route((.1,1),(.9,1),10,
                          self.region,self.cfg))

    def test_short_routed_travel_can_avoid_retraction_but_longer_still_retracts(self):
        current=Toolpath("SOLID_INFILL",10,[(.1,5),self.start])
        following=Toolpath("SOLID_INFILL",10,[self.end,(.9,5)])
        self.assertFalse(_should_retract_after_path(
            current, following, GcodeState(theta_accum=.1,axial=5),self.cfg,
            self.region))
        self.assertTrue(_should_retract_after_path(
            current, following, GcodeState(theta_accum=.1,axial=5),
            SlicerConfig(route_max_unretracted_mm=2),self.region))
        self.assertTrue(_should_retract_after_path(
            current, following, GcodeState(theta_accum=.1,axial=5),
            SlicerConfig(route_travel_inside_island=False),self.region))

    def test_emitted_g0_waypoints_do_not_extrude_and_e_is_conserved(self):
        following=Toolpath("SOLID_INFILL",10,[self.end,(.8,8)])
        before=GcodeState(theta_accum=.2,axial=5,radial=10.08)
        routed=emit_toolpath(following,before,self.cfg,0,
                             retract_after=False,travel_region=self.region)
        straight=emit_toolpath(following,
            GcodeState(theta_accum=.2,axial=5,radial=10.08),
            SlicerConfig(route_travel_inside_island=False),0,retract_after=False,
            travel_region=self.region)
        self.assertTrue(any(s.startswith(";ROUTED_TRAVEL:") for s in routed))
        self.assertGreater(len([s for s in routed if s.startswith("G0 X")]),1)
        self.assertEqual([s for s in straight if s.startswith("G1 X")],
                         [s for s in routed if s.startswith("G1 X")])
        self.assertAlmostEqual(before.extrusion,
            3*self.cfg.extrusion_per_mm,places=6)

    def test_validation(self):
        for fields in (
            {"route_max_vertices":3},
            {"route_max_factor":.9},
            {"route_max_extra_mm":-1},
            {"route_max_unretracted_mm":math.nan},
            {"route_travel_inside_island":"sure"},
        ):
            with self.subTest(fields=fields),self.assertRaises(ValueError):
                SlicerConfig(**fields).validate()


if __name__ == "__main__":
    unittest.main()
