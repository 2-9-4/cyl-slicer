from collections import Counter
from dataclasses import replace
import math
from pathlib import Path
import re
import unittest

from shapely import affinity
from shapely.geometry import LineString, Point, Polygon, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (
    MAX_MARLIN_Y_DEGREES,
    _ordered_paths_for_current_position, _orient_path_for_current_position,
    _safe_island_travel_route, _should_retract_after_path,
    emit_gcode, emit_toolpath, parse_extrusion_segments,
)
from cyl_slicer_v2.layers import slice_mesh
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import GcodeState, Toolpath
from cyl_slicer_v2.periodic import periodic_copies
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, gap_support_slices
from cyl_slicer_v2.tests.test_marlin_rotary_precision import MarlinRotaryModel


def cylindrical(shape, radius):
    return affinity.scale(shape, xfact=1/radius, yfact=1, origin=(0, 0))


def route_length(start, route, radius):
    return sum(math.hypot((b[0]-a[0])*radius, b[1]-a[1])
               for a, b in zip([start, *route], route))


class LocalSupportTravelTests(unittest.TestCase):
    def assert_inside(self, shape, start, route, radius):
        self.assertIsNotNone(route)
        for a, b in zip([start, *route], route):
            segment = LineString([(a[0]*radius, a[1]), (b[0]*radius, b[1])])
            self.assertTrue(shape.buffer(1e-7).covers(segment), (a, b))

    def test_remote_outline_detail_does_not_consume_local_graph_budget(self):
        # Hundreds of real corners far from the move must not prevent routing
        # around a small notch. The useful graph needs just its two corners.
        bottom = [(20*i/300, .05*(i % 2)) for i in range(301)]
        shape = Polygon([*bottom, (20, 10), (11, 10), (11, 5),
                         (9, 5), (9, 10), (0, 10)])
        radius = 30
        start, end = (8/radius, 7), (12/radius, 7)
        config = SlicerConfig(route_max_vertices=4)
        route = _safe_island_travel_route(start, end, radius,
                                          periodic_copies(cylindrical(shape, radius)), config)
        self.assert_inside(shape, start, route, radius)
        self.assertEqual(route[-1], end)
        self.assertLessEqual(len(route)+1, 4)
        self.assertLessEqual(route_length(start, route, radius), 4*config.route_max_factor + 1e-8)

    def test_dense_curved_boundary_is_followed_without_crossing_its_notch(self):
        shape = box(0, 0, 20, 10).difference(Point(10, 10).buffer(2, quad_segs=128))
        radius = 30
        points = sorted((x, y) for x, y in shape.exterior.coords if 8.5 < x < 10 and y < 9)
        a, b = points[0], points[-1]
        start, end = (a[0]/radius, a[1]), (b[0]/radius, b[1])
        self.assertFalse(shape.covers(LineString([a, b])))
        region = periodic_copies(cylindrical(shape, radius))
        route = _safe_island_travel_route(start, end, radius, region, SlicerConfig())
        self.assert_inside(shape, start, route, radius)
        direct = math.dist(a, b)
        self.assertLessEqual(route_length(start, route, radius), direct*1.65 + 1e-8)
        self.assertLessEqual(len(route)+1, 96)
        # The same local obstacle can exceed a deliberately small budget;
        # fallback must remain available rather than cutting through the notch.
        small = SlicerConfig(route_max_vertices=4)
        self.assertIsNone(_safe_island_travel_route(start, end, radius, region, small))
        current = Toolpath("SUPPORT", radius, [(start[0], start[1]-.2), start])
        following = Toolpath("SUPPORT", radius, [end, (end[0], end[1]-.2)])
        state = GcodeState(theta_accum=start[0], axial=start[1]-.2)
        self.assertTrue(_should_retract_after_path(current, following, state, small, region))
        lines = emit_toolpath(following, GcodeState(theta_accum=start[0], axial=start[1], radial=radius+.08),
                              replace(small, z_hop_min_travel=.1), 0, travel_region=region)
        self.assertTrue(any("; z hop" in line for line in lines))

    def test_local_routing_respects_distance_budget_and_periodic_frame(self):
        shape = box(0, 0, 20, 10).difference(Point(10, 10).buffer(2, quad_segs=64))
        radius = 30
        boundary = sorted((x, y) for x, y in shape.exterior.coords if 8.5 < x < 10 and y < 9)
        a, b = boundary[0], boundary[-1]
        base = cylindrical(shape, radius)
        region = periodic_copies(affinity.translate(base, xoff=2*math.pi-.3))
        for turns in (-7, 0, 9):
            shift = 2*math.pi-.3 + turns*2*math.pi
            start, end = (a[0]/radius+shift, a[1]), (b[0]/radius+shift, b[1])
            route = _safe_island_travel_route(start, end, radius, region, SlicerConfig())
            self.assertIsNotNone(route)
            local = [(p[0]-shift, p[1]) for p in route]
            self.assert_inside(shape, (start[0]-shift, start[1]), local, radius)
            self.assertEqual(route[-1], end)
            self.assertIsNone(_safe_island_travel_route(start, end, radius, region,
                              SlicerConfig(route_max_factor=1.001, route_max_extra_mm=.001)))

    def test_support_order_uses_verified_travel_length_instead_of_air_distance(self):
        region = periodic_copies(box(0, 0, 1, 10).difference(box(.4, 4, .6, 6)))
        near = Toolpath("SUPPORT", 10, [(.61, 5), (.9, 5)])
        reachable = Toolpath("SUPPORT", 10, [(.3, 8.2), (.3, 9.2)])
        state = GcodeState(theta_accum=.3, axial=5)
        # The nearer endpoint needs a longer detour around the hole.
        old = _ordered_paths_for_current_position([near, reachable], state)
        new = _ordered_paths_for_current_position([near, reachable], state, SlicerConfig(), region)
        self.assertIs(old[0], near)
        self.assertIs(new[0], reachable)

    def test_routed_waypoints_rebase_marlin_safely_in_both_directions(self):
        region = periodic_copies(box(-.1, 0, .1, 10).difference(box(-.005, 4.9, .005, 5.1)))
        for sign in (-1, 1):
            with self.subTest(sign=sign):
                state = GcodeState(theta_accum=math.radians(sign*3599), axial=5, radial=30.08)
                following = Toolpath("SUPPORT", 30,
                                     [(math.radians(sign), 5), (math.radians(sign), 5.3)])
                lines = emit_toolpath(following, state, SlicerConfig(), 0,
                                      travel_region=region, retract_after=False)
                self.assertTrue(any(line.startswith(";ROUTED_TRAVEL:") for line in lines))
                marlin = MarlinRotaryModel().replay(f"G0 Y{sign*3599}\n")
                before = marlin.physical_steps
                marlin.replay("\n".join(lines))
                self.assertEqual(marlin.native_resets, 1)
                self.assertLessEqual(marlin.max_native, MAX_MARLIN_Y_DEGREES)
                self.assertAlmostEqual((marlin.physical_steps-before)/float(marlin.steps_per_degree),
                                       sign*2, delta=1/float(marlin.steps_per_degree))

    def test_short_contained_support_travel_stays_primed_but_real_gaps_do_not(self):
        current = Toolpath("SUPPORT", 30, [(.1, 0), (.2, 0)])
        following = Toolpath("SUPPORT_INTERFACE", 30, [(.2, 3), (.3, 3)])
        state = GcodeState(theta_accum=.1, axial=0)
        region = periodic_copies(box(.1, 0, .3, 3))
        self.assertFalse(_should_retract_after_path(current, following, state, SlicerConfig(), region))
        self.assertTrue(_should_retract_after_path(current, following, state,
                                                 SlicerConfig(route_max_unretracted_mm=2), region))
        broken = region.difference(box(.09, 1, .31, 2))
        self.assertTrue(_should_retract_after_path(current, following, state, SlicerConfig(), broken))

    def test_reversing_support_meanders_preserves_widths_endpoints_and_filament(self):
        for role in ("SUPPORT", "SUPPORT_INTERFACE"):
            points = [(math.radians(t), x) for t, x in ((350, 0), (370, 0), (370, 1), (350, 1))]
            path = Toolpath(role, 30, points, segment_widths=[.25, .3, .4])
            state = GcodeState(theta_accum=points[-1][0], axial=points[-1][1])
            oriented = _orient_path_for_current_position(path, state)
            self.assertEqual(oriented.points, list(reversed(points)))
            self.assertEqual(oriented.segment_widths, [.4, .3, .25])
            amounts = []
            for candidate in (path, oriented):
                before = replace(state)
                emit_toolpath(candidate, before, SlicerConfig(), 0, retract_after=False)
                amounts.append(before.extrusion)
            self.assertAlmostEqual(*amounts, places=9)
            self.assertEqual(path.points, points)
            self.assertEqual(path.segment_widths, [.25, .3, .4])

    def test_offset_cube_support_layer_reduces_retractions_without_changing_deposition(self):
        config = SlicerConfig(axial_min=30)
        mesh = load_mesh(Path(__file__).resolve().parents[2] / "offset_cube.stl")
        layer = gap_support_slices(slice_mesh(mesh, config), config)[4]
        planned = build_layer_toolpaths(layer, config)
        routed = emit_gcode([planned], config)
        fallback = emit_gcode([planned], replace(config, route_travel_inside_island=False))
        self.assertLessEqual(routed.count("; retract"), 3)
        self.assertGreater(fallback.count("; retract"), 10)
        def deposition(code):
            def point(x, y):
                return (round(x, 3), round(y % 360, 3) % 360)
            return Counter((role, *sorted((point(x0, y0), point(x1, y1))))
                           for _, role, x0, y0, x1, y1 in parse_extrusion_segments(code))
        self.assertEqual(deposition(routed), deposition(fallback))
        def filament(code):
            previous = total = 0
            for line in code.splitlines():
                match = re.search(r"\bE(-?\d+(?:\.\d+)?)", line)
                if match:
                    current = float(match[1])
                    if line.startswith("G1 X"):
                        total += current - previous
                    previous = current
            return total
        self.assertAlmostEqual(filament(routed), filament(fallback), places=4)


if __name__ == "__main__":
    unittest.main()
