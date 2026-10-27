from __future__ import annotations

import math
import unittest

from shapely.geometry import LineString, Point, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (
    _safe_island_travel_route, _should_retract_after_path, emit_gcode, emit_toolpath,
)
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath
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

    def test_boundary_travel_is_primed_in_every_rotary_frame_and_role(self):
        lo, hi = 2.844154778070864, 3.4390475507779064
        region = periodic_copies(box(lo, 0, hi, 2))
        config = SlicerConfig()
        for role in ("SUPPORT", "SUPPORT_INTERFACE", "SOLID_INFILL", "INFILL_RADIAL",
                     "EXTERNAL_PERIMETER"):
            for turns in (-250000, -20, -3, -1, 0, 1, 7, 20, 250000):
                with self.subTest(role=role, turns=turns):
                    current = Toolpath(role, 30.4, [(hi, 1), (lo, 1)])
                    following = Toolpath("SOLID_INFILL", 30.4, [(lo, 1.4), (hi, 1.4)])
                    state = GcodeState(theta_accum=hi + turns * 2*math.pi, axial=1)
                    self.assertFalse(_should_retract_after_path(
                        current, following, state, config, region))

    def test_short_boundary_travel_still_retracts_across_a_real_gap(self):
        lo, hi = 2.844154778070864, 3.4390475507779064
        # The 0.002 mm notch is much smaller than output travel spacing, but
        # is a real opening, not numerical drift from a rotary winding.
        region = periodic_copies(box(lo, 0, hi, 2).difference(box(lo-.1, 1.199, lo+.1, 1.201)))
        for turns in (-7, 0, 7):
            current = Toolpath("SUPPORT", 30.4, [(hi, 1), (lo, 1)])
            following = Toolpath("SOLID_INFILL", 30.4, [(lo, 1.4), (hi, 1.4)])
            state = GcodeState(theta_accum=hi + turns * 2*math.pi, axial=1)
            self.assertTrue(_should_retract_after_path(
                current, following, state, SlicerConfig(), region))

    def test_numerical_boundary_travel_uses_no_visibility_graph_or_waypoints(self):
        lo, hi = 2.844154778070864, 3.4390475507779064
        region = periodic_copies(box(lo, 0, hi, 2))
        for radius in (1, 30.4, 100):
            for turns in (-20, -1, 0, 7):
                shift = turns * 2*math.pi
                # A physical 1e-9 mm error at the boundary is below the
                # tolerance, independent of radius or accumulated turns.
                start = (lo + shift - 1e-9/radius, 1)
                end = (lo + shift, 1.4)
                self.assertIsNone(_safe_island_travel_route(
                    start, end, radius, region, SlicerConfig(route_max_vertices=4)))
                path = Toolpath("SUPPORT", radius, [(lo, 1.4), (hi, 1.4)])
                lines = emit_toolpath(path, GcodeState(theta_accum=start[0], axial=1),
                                      SlicerConfig(), 0, travel_region=region, retract_after=False)
                self.assertFalse(any(line.startswith(";ROUTED_TRAVEL:") for line in lines))

    def test_repeated_bidirectional_layers_do_not_retract_between_boundary_passes(self):
        lo, hi = 2.844154778070864, 3.4390475507779064
        for role in ("SUPPORT", "SUPPORT_INTERFACE", "SOLID_INFILL", "INFILL_RADIAL"):
            with self.subTest(role=role):
                paths = [Toolpath(role, 30.4, [(lo, x), (hi, x)]) for x in (0, .4, .8, 1.2)]
                layers = [LayerToolpaths(i, 30.4, {role: paths}, travel_region=box(lo, 0, hi, 1.2))
                          for i in range(4)]
                code = emit_gcode(layers, SlicerConfig())
                for layer in code.split(";LAYER:")[1:]:
                    self.assertEqual(layer.count("; retract"), 1, "Layer end only")

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

    def test_retraction_and_emitted_detour_agree_at_accumulated_rotary_positions(self):
        current = Toolpath("SOLID_INFILL", 10, [(.1, 5), self.start])
        following = Toolpath("SUPPORT", 10, [self.end, (.9, 5)])
        for turns in (-250000, -7, 0, 7, 250000):
            with self.subTest(turns=turns):
                shift = turns * 2*math.pi
                state = GcodeState(theta_accum=.1 + shift, axial=5)
                self.assertFalse(_should_retract_after_path(
                    current, following, state, self.cfg, self.region))
                emitted = emit_toolpath(following,
                    GcodeState(theta_accum=self.start[0] + shift, axial=5),
                    self.cfg, 0, travel_region=self.region, retract_after=False)
                self.assertTrue(any(line.startswith(";ROUTED_TRAVEL:") for line in emitted))

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
