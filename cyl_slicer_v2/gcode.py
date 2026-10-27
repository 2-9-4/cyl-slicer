from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import heapq
import math
import re

from shapely import affinity

from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from shapely.geometry.base import BaseGeometry
from shapely.strtree import STRtree

from .config import SlicerConfig
from .coordinates import coalesce_collinear_points, cylindrical_distance, unwrap_path_from_current
from .models import GcodeState, LayerToolpaths, Toolpath
from .periodic import periodic_copies
from .toolpaths import MIN_GEOMETRY_LENGTH, MIN_OUTPUT_LENGTH_GUARD, PERIMETER_ROLES, ROLE_ORDER
from .wall_support import measure_wall_support


TWO_PI = 2.0 * math.pi
PATH_ORDER_OPTIMIZE_LIMIT = 1024
BIDIRECTIONAL_OUTPUT_ROLES = PERIMETER_ROLES | frozenset({"SUPPORT", "SUPPORT_INTERFACE", "INFILL_RADIAL", "SOLID_INFILL", "GAP_FILL"})
# Precision budget, not a Marlin travel limit: float32 spacing at 3600 degrees
# is 0.000244 degrees, below our 0.001-degree output resolution.
MAX_MARLIN_Y_DEGREES = 3600.0
# Numerical boundary tolerance in physical mm, far below output resolution.
TRAVEL_CONTAINMENT_TOLERANCE_MM = 1e-7


@dataclass(frozen=True)
class PerimeterIndex:
    segments: list[LineString]
    tree: STRtree | None
    period: float


def _format_move(cmd: str, x: float | None = None, y_deg: float | None = None, z: float | None = None, e: float | None = None, f: float | None = None) -> str:
    parts = [cmd]
    if x is not None:
        parts.append(f"X{x:.3f}")
    if y_deg is not None:
        parts.append(f"Y{y_deg:.3f}")
    if z is not None:
        parts.append(f"Z{z:.4f}")
    if e is not None:
        parts.append(f"E{e:.5f}")
    if f is not None:
        parts.append(f"F{f:.1f}")
    return " ".join(parts)


def _surface_feed_for_marlin(
    dx_mm: float, dtheta_radians: float, radius_mm: float, physical_feed_mm_per_min: float
) -> float:
    """Convert physical surface speed to Marlin's mixed X(mm)/Y(degrees) F.

    Marlin computes the planar planner distance from its X and Y coordinate
    numbers, which are not both millimetres for this cylindrical machine.
    The motion duration is L_planner / F, so scale F by the ratio of planner
    to physical length. Z-only motions and E-only retracts remain in mm/min.
    """
    physical_length = math.hypot(dx_mm, radius_mm * dtheta_radians)
    if physical_length <= 1e-12:
        return physical_feed_mm_per_min
    planner_length = math.hypot(dx_mm, math.degrees(dtheta_radians))
    return physical_feed_mm_per_min * planner_length / physical_length


def _emitted_y_degrees(theta: float, state: GcodeState) -> float:
    return math.degrees(theta - state.theta_origin)


def _reset_rotary_axis(
    state: GcodeState,
    *,
    reset_extrusion: bool = False,
    comment: str = "",
) -> list[str]:
    # The target Marlin firmware enables POWER_LOSS_RECOVERY (G92.9). Ordinary
    # G92 only shifts the workspace, leaving native/planner coordinates growing.
    # Rebase by whole turns, retaining angular phase so layer resets don't
    # repeatedly round a fractional motor step into a new zero. The following
    # G92 keeps logical Y layer-local while native Y starts at the wrapped phase.
    native_phase = float(f"{math.degrees(math.remainder(state.theta_accum, TWO_PI)):.3f}")
    state.theta_origin = state.theta_accum
    state.native_y_origin_degrees = native_phase
    if reset_extrusion:
        state.extrusion = 0.0
    reset = "G92 E0 Y0" if reset_extrusion else "G92 Y0"
    return [
        "M400 ; finish moves before rotary coordinate reset",
        f"G92.9 Y{native_phase:.3f} ; rebase Marlin native rotary position",
        reset + (" ; " + comment if comment else ""),
    ]


def _ensure_rotary_range(theta: float, state: GcodeState) -> list[str]:
    target = _emitted_y_degrees(theta, state)
    if not math.isfinite(target):
        raise ValueError("Rotary coordinate must be finite")
    # Leave room for formatting to 0.001 degrees and float32 addition of the
    # workspace offset. Bound the native planner target as well as logical Y.
    limit = MAX_MARLIN_Y_DEGREES - 0.001
    if max(abs(target), abs(target + state.native_y_origin_degrees)) <= limit:
        return []
    delta = math.degrees(theta - state.theta_accum)
    phase = float(f"{math.degrees(math.remainder(state.theta_accum, TWO_PI)):.3f}")
    if max(abs(delta), abs(delta + phase)) > limit:
        raise ValueError("A single rotary move exceeds the Marlin precision budget (3600 degrees)")
    return _reset_rotary_axis(state, comment="reset rotary accumulator within layer")


def _emitted_x(axial: float, output_x_offset: float) -> float:
    return axial + output_x_offset


def _emitted_z(radial: float, config: SlicerConfig) -> float:
    return radial - config.bed_radius


def _output_x_offset(layers: list[LayerToolpaths], config: SlicerConfig) -> float:
    axial_values = [
        point[1]
        for layer in layers
        for paths in layer.paths_by_role.values()
        for path in paths
        for point in path.points
    ]
    if not axial_values:
        return config.axial_min
    return config.axial_min - min(axial_values)


def _nearest_theta(theta: float, reference: float) -> float:
    delta = (theta - reference + math.pi) % TWO_PI - math.pi
    return reference + delta


def _travel_distance_to_point(path: Toolpath, point: tuple[float, float], state: GcodeState) -> float:
    theta = _nearest_theta(point[0], state.theta_accum)
    return math.hypot((theta - state.theta_accum) * path.radius, point[1] - state.axial)


def _metric_point(theta: float, axial: float, radius: float) -> tuple[float, float]:
    return theta * radius, axial


def _path_metric_segments(path: Toolpath) -> list[LineString]:
    if len(path.points) < 2:
        return []
    thetas = _unwrap_toolpath_thetas(path, path.points[0][0])
    pts = [_metric_point(theta, axial, path.radius) for theta, (_, axial) in zip(thetas, path.points)]
    return [
        LineString([a, b])
        for a, b in zip(pts[:-1], pts[1:])
        if math.hypot(b[0] - a[0], b[1] - a[1]) > 1e-8
    ]


def _layer_perimeter_index(layer: LayerToolpaths) -> PerimeterIndex | None:
    segments = []
    for role in ("EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"):
        for path in layer.paths_by_role.get(role, []):
            segments.extend(_path_metric_segments(path))
    if not segments:
        return None
    radius = next(
        path.radius
        for role in ("EXTERNAL_PERIMETER", "INTERNAL_PERIMETER")
        for path in layer.paths_by_role.get(role, [])
        if len(path.points) >= 2
    )
    return PerimeterIndex(segments, STRtree(segments), TWO_PI * radius)


def _query_perimeter_index(index: PerimeterIndex, geom: LineString):
    if index.tree is None:
        yield from index.segments
        return
    for item in index.tree.query(geom):
        try:
            yield index.segments[int(item)]
        except (TypeError, ValueError, IndexError):
            yield item


def _intersection_is_away_from_endpoints(
    geom: BaseGeometry,
    start: Point,
    end: Point,
    tol: float = 1e-6,
) -> bool:
    if geom.is_empty:
        return False
    if isinstance(geom, Point):
        return geom.distance(start) > tol and geom.distance(end) > tol
    if hasattr(geom, "geoms"):
        return any(_intersection_is_away_from_endpoints(part, start, end, tol) for part in geom.geoms)
    # Following a perimeter is not crossing it. In particular, adjacent solid
    # infill lines end on the same inner wall and travel along it to the next row.
    return False


def _travel_crosses_perimeter(
    path: Toolpath,
    point: tuple[float, float],
    state: GcodeState,
    perimeter_index: PerimeterIndex,
) -> bool:
    if not perimeter_index.segments:
        return False
    target_theta = _nearest_theta(point[0], state.theta_accum)
    start_xy = _metric_point(state.theta_accum, state.axial, path.radius)
    end_xy = _metric_point(target_theta, point[1], path.radius)
    tx, ty = end_xy[0] - start_xy[0], end_xy[1] - start_xy[1]
    travel_length = math.hypot(tx, ty)
    if travel_length <= 1e-8:
        return False

    travel = LineString([start_xy, end_xy])
    travel_min_x, _, travel_max_x, _ = travel.bounds
    period = perimeter_index.period
    shift_min = math.floor((travel_min_x - period) / period) - 1
    shift_max = math.ceil(travel_max_x / period) + 1
    for shift in range(shift_min, shift_max + 1):
        shifted_travel = LineString([(x - shift * period, y) for x, y in travel.coords])
        shifted_start = Point(shifted_travel.coords[0])
        shifted_end = Point(shifted_travel.coords[-1])
        for segment in _query_perimeter_index(perimeter_index, shifted_travel):
            # Collinear fragments can also intersect at a shared vertex. Those
            # point contacts are still travel along the wall, not across it.
            (sx0, sy0), (sx1, sy1) = segment.coords
            sx, sy = sx1 - sx0, sy1 - sy0
            if abs(tx * sy - ty * sx) <= 1e-12 * travel_length * math.hypot(sx, sy):
                continue
            inter = shifted_travel.intersection(segment)
            if _intersection_is_away_from_endpoints(inter, shifted_start, shifted_end):
                return True
    return False



def _travel_polygon_parts(geom: BaseGeometry):
    """Polygonal components of the periodic island footprint, with holes."""
    if isinstance(geom, Polygon):
        yield geom
    elif isinstance(geom, MultiPolygon) or hasattr(geom, "geoms"):
        for part in geom.geoms:
            yield from _travel_polygon_parts(part)


@lru_cache(maxsize=8)
def _metric_travel_islands(travel_region: BaseGeometry, radius: float):
    """Reuse physical footprints for containment checks and visibility routing.

    Buffer each island separately so numerical boundary tolerance cannot join
    separate islands. Shapely geometries are immutable; bound the layer cache.
    """
    result = []
    for polygon in _travel_polygon_parts(travel_region):
        surface = affinity.scale(polygon, xfact=radius, yfact=1.0, origin=(0, 0))
        result.append((surface, surface.buffer(TRAVEL_CONTAINMENT_TOLERANCE_MM)))
    return tuple(result)


def _local_travel_points(start, destination, radius):
    # Shift BOTH endpoints into the footprint's periodic neighborhood. Never
    # combine a wrapped destination with an accumulated rotary start position.
    shift = math.floor(start[0] / TWO_PI) * TWO_PI
    return (shift, ((start[0] - shift) * radius, start[1]),
            ((destination[0] - shift) * radius, destination[1]))


def _travel_is_contained(start, destination, radius, travel_region):
    _, start_mm, end_mm = _local_travel_points(start, destination, radius)
    travel = LineString([start_mm, end_mm])
    return any(clearance.covers(travel)
               for _, clearance in _metric_travel_islands(travel_region, radius))


def _safe_island_travel_route(
    start: tuple[float, float],
    destination: tuple[float, float],
    radius: float,
    travel_region: BaseGeometry | None,
    config: SlicerConfig,
) -> list[tuple[float, float]] | None:
    """Shortest bounded inside-island polyline, excluding the starting point.

    The theta coordinate is unwrapped and measured in radians. Search/limit
    lengths in physical (theta*radius, axial) millimeters; never connect two
    disjoint polygonal islands. An unsupported or too-expensive graph returns
    None, preserving existing retraction/hop safety.
    """
    if (not config.route_travel_inside_island or travel_region is None
            or travel_region.is_empty or radius <= 0):
        return None
    direct = math.hypot((destination[0]-start[0])*radius,
                        destination[1]-start[1])
    if direct < 1e-7:
        return None
    shift, start_mm, end_mm = _local_travel_points(start, destination, radius)
    s = Point(start_mm)
    e = Point(end_mm)
    straight = LineString([start_mm, end_mm])
    for surface, clearance in _metric_travel_islands(travel_region, radius):
        if not clearance.covers(s) or not clearance.covers(e):
            continue
        if clearance.covers(straight):
            return None  # Existing straight travel is already safe.
        # Route around holes or concave cutouts by visible boundary vertices.
        coords = list(surface.exterior.coords)[:-1]
        for ring in surface.interiors:
            coords.extend(list(ring.coords)[:-1])
        if len(coords) + 2 > config.route_max_vertices:
            return None
        nodes = [start_mm, end_mm, *coords]
        n = len(nodes)
        max_length = min(direct * config.route_max_factor,
                         direct + config.route_max_extra_mm)
        edges = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                dist = math.dist(nodes[i], nodes[j])
                if dist > max_length + 1e-8:
                    continue
                if not clearance.covers(LineString([nodes[i], nodes[j]])):
                    continue
                edges[i].append((j, dist))
                edges[j].append((i, dist))
        distances = [math.inf] * n
        distances[0] = 0.0
        parents = [-1] * n
        heap = [(0.0, 0)]
        while heap:
            distance, i = heapq.heappop(heap)
            if distance > distances[i] + 1e-10 or distance > max_length + 1e-8:
                continue
            if i == 1:
                order = []
                while i != 0:
                    order.append(i)
                    i = parents[i]
                order.reverse()
                route = [(nodes[j][0]/radius + shift, nodes[j][1]) for j in order]
                route[-1] = destination
                # A tolerance buffer admits endpoints on the actual contour;
                # never accept a segment crossing the actual empty space.
                unwrapped = [start_mm, *[nodes[j] for j in order]]
                if all(clearance.covers(LineString([a,b]))
                       for a,b in zip(unwrapped, unwrapped[1:])):
                    return route
                return None
            for j, step in edges[i]:
                candidate = distance + step
                if candidate + 1e-10 < distances[j] and candidate <= max_length + 1e-8:
                    distances[j] = candidate
                    parents[j] = i
                    heapq.heappush(heap, (candidate, j))
        return None
    return None


def _should_z_hop_to_path(
    path: Toolpath,
    state: GcodeState,
    config: SlicerConfig,
    perimeter_index: PerimeterIndex | None,
) -> bool:
    if config.z_hop_height <= 0 or len(path.points) < 1:
        return False
    target = path.points[0]
    travel_distance = _travel_distance_to_point(path, target, state)
    if travel_distance <= 1e-8:
        return False
    if travel_distance >= config.z_hop_min_travel:
        return True
    return (
        config.z_hop_on_perimeter_crossing
        and perimeter_index is not None
        and _travel_crosses_perimeter(path, target, state, perimeter_index)
    )


def _can_reverse_for_output(path: Toolpath, tol: float = 1e-8) -> bool:
    if path.role in PERIMETER_ROLES or path.role in {"TOP_SKIN", "GAP_FILL"}:
        return False
    if path.closed or len(path.points) < 2:
        return False
    theta0, axial0 = path.points[0]
    return all(abs(point[0] - theta0) <= tol for point in path.points[1:]) or (
        path.role in {"INFILL_RADIAL", "SOLID_INFILL", "SUPPORT", "SUPPORT_INTERFACE"}
        and all(abs(point[1] - axial0) <= tol for point in path.points[1:])
    )


def _orient_path_for_current_position(
    path: Toolpath, state: GcodeState, config: SlicerConfig | None = None
) -> Toolpath:
    if path.role in PERIMETER_ROLES and len(path.points) > 3 and cylindrical_distance(path.points[0], path.points[-1], path.radius) <= 1e-8:
        return _rotate_closed_perimeter(path, state, config)
    if not _can_reverse_for_output(path):
        return path
    start_distance = _travel_distance_to_point(path, path.points[0], state)
    end_distance = _travel_distance_to_point(path, path.points[-1], state)
    if end_distance + 1e-9 >= start_distance:
        return path
    return replace(
        path,
        points=list(reversed(path.points)),
        segment_widths=list(reversed(path.segment_widths)) if path.segment_widths is not None else None,
    )


def _rotate_closed_perimeter(
    path: Toolpath, state: GcodeState, config: SlicerConfig | None = None
) -> Toolpath:
    if config is not None and (
        config.seam_policy != "nearest" or config.seam_exclude_deg is not None
    ):
        return _rotate_preferred_seam(path, state, config)
    loop = path.points[:-1]
    center = (min(p[0] for p in loop) + max(p[0] for p in loop)) * 0.5
    current_theta = _nearest_theta(state.theta_accum, center)
    best = None
    for shift in (0.0, -TWO_PI, TWO_PI):
        current = ((current_theta + shift) * path.radius, state.axial)
        for i, (a, b) in enumerate(zip(path.points[:-1], path.points[1:])):
            ax, ay = a[0] * path.radius, a[1]
            dx, dy = (b[0] - a[0]) * path.radius, b[1] - a[1]
            length2 = dx * dx + dy * dy
            if length2 <= 1e-16:
                continue
            fraction = max(0.0, min(1.0, ((current[0] - ax) * dx + (current[1] - ay) * dy) / length2))
            distance2 = (current[0] - ax - fraction * dx) ** 2 + (current[1] - ay - fraction * dy) ** 2
            if best is None or distance2 < best[0]:
                best = (distance2, i, fraction)
    if best is None:
        return path
    _, index, fraction = best
    if fraction <= 1e-10 or fraction >= 1.0 - 1e-10:
        index = (index + (fraction >= 1.0 - 1e-10)) % len(loop)
        points = loop[index:] + loop[:index] + [loop[index]]
    else:
        a, b = path.points[index:index + 2]
        point = (a[0] + (b[0] - a[0]) * fraction, a[1] + (b[1] - a[1]) * fraction)
        points = [point, *loop[index + 1:], *loop[:index + 1], point]
    return replace(path, points=points)



def _seam_angle_distance(theta: float, angle_rad: float) -> float:
    """Shortest periodic theta distance, independent of rotary winding."""
    return abs((theta - angle_rad + math.pi) % TWO_PI - math.pi)


def _seam_angle_excluded(theta: float, window: tuple[float, float] | None) -> bool:
    if window is None:
        return False
    angle = math.degrees(theta) % 360.0
    begin, end = window
    if begin < end:
        return begin - 1e-9 <= angle < end - 1e-9
    return angle >= begin - 1e-9 or angle < end - 1e-9


def _rotate_preferred_seam(path: Toolpath, state: GcodeState, config: SlicerConfig) -> Toolpath:
    """Choose a phase-aligned or high-curvature seam without changing loop geometry.

    Candidate points include every original vertex, the closest point on
    each edge in the physical cylindrical metric, and intersections with the
    desired fixed-angle ray (and exclusion-window boundaries). Only choose
    points on real edges. A cut simply rotates/splits the original polyline,
    retaining exact oriented geometry, including non-zero rotary winding.
    """
    loop = path.points[:-1]
    n = len(loop)
    radius = path.radius
    target = math.radians(config.seam_angle_deg)
    center = 0.5 * (min(p[0] for p in loop) + max(p[0] for p in loop))
    current_theta = _nearest_theta(state.theta_accum, center)
    current = (current_theta * radius, state.axial)
    best = None

    for index, (a, b) in enumerate(zip(path.points[:-1], path.points[1:])):
        ax, ay = a[0] * radius, a[1]
        dx, dy = (b[0] - a[0]) * radius, b[1] - a[1]
        length2 = dx * dx + dy * dy
        if length2 < 1e-16:
            continue
        fractions = [0.0, 1.0]
        for shift in (0, -1, 1):
            here = ((current_theta + shift * TWO_PI) * radius, state.axial)
            fractions.append(max(0.0, min(1.0,
                ((here[0] - ax) * dx + (here[1] - ay) * dy) / length2)))

        ray_degrees = [config.seam_angle_deg]
        if config.seam_exclude_deg is not None:
            # The upper exclusion boundary is printable. The lower boundary
            # is blocked by design, but its neighboring vertices are tested.
            ray_degrees.append(config.seam_exclude_deg[1])
        for degrees in ray_degrees:
            phase = math.radians(degrees)
            lower = math.floor((min(a[0], b[0]) - phase) / TWO_PI) - 1
            upper = math.ceil((max(a[0], b[0]) - phase) / TWO_PI) + 1
            if abs(b[0] - a[0]) <= 1e-14:
                continue
            for k in range(lower, upper + 1):
                fraction = (phase + TWO_PI * k - a[0]) / (b[0] - a[0])
                if 0 <= fraction <= 1:
                    fractions.append(fraction)
        for fraction in fractions:
            theta = a[0] + fraction * (b[0] - a[0])
            axial = a[1] + fraction * (b[1] - a[1])
            if _seam_angle_excluded(theta, config.seam_exclude_deg):
                continue
            angular_cost = radius * _seam_angle_distance(theta, target)
            travel = min(
                math.hypot((theta + shift * TWO_PI) * radius - current[0], axial - current[1])
                for shift in (-1, 0, 1)
            )
            turn = 0.0
            if fraction < 1e-10 or fraction > 1.0 - 1e-10:
                vertex = index if fraction < 1e-10 else (index + 1) % n
                prev_pt = loop[(vertex - 1) % n]
                vertex_pt = loop[vertex]
                next_pt = loop[(vertex + 1) % n]
                u = ((vertex_pt[0] - prev_pt[0]) * radius, vertex_pt[1] - prev_pt[1])
                v = ((next_pt[0] - vertex_pt[0]) * radius, next_pt[1] - vertex_pt[1])
                denom = math.hypot(*u) * math.hypot(*v)
                if denom > 1e-12:
                    turn = 1.0 - max(-1.0, min(1.0,
                        (u[0] * v[0] + u[1] * v[1]) / denom))
            if config.seam_policy == "corner":
                cost = (-round(turn, 9), angular_cost, travel)
            elif config.seam_policy == "aligned":
                cost = (angular_cost, travel, -turn)
            else:
                cost = (travel, angular_cost, -turn)
            if best is None or cost < best[0]:
                best = (cost, index, fraction)

    if best is None:
        raise ValueError("No eligible perimeter seam point outside excluded angular window")
    _, index, fraction = best
    if fraction <= 1e-10 or fraction >= 1.0 - 1e-10:
        index = (index + (fraction >= 1.0 - 1e-10)) % n
        points = loop[index:] + loop[:index] + [loop[index]]
    else:
        a, b = path.points[index:index + 2]
        point = (a[0] + (b[0] - a[0]) * fraction, a[1] + (b[1] - a[1]) * fraction)
        points = [point, *loop[index + 1:], *loop[:index + 1], point]
    return replace(path, points=points)


def _shadow_after_path(path: Toolpath, state: GcodeState) -> GcodeState:
    accum = _unwrap_toolpath_thetas(path, state.theta_accum)
    if not accum:
        return state
    return GcodeState(
        theta_accum=accum[-1],
        theta_origin=state.theta_origin,
        native_y_origin_degrees=state.native_y_origin_degrees,
        axial=path.points[-1][1],
        radial=state.radial,
        extrusion=state.extrusion,
        retracted=state.retracted,
        retracted_mm=state.retracted_mm,
        retracted_speed_mm_per_min=state.retracted_speed_mm_per_min,
    )


def _monotonic_skin_paths(paths: list[Toolpath]) -> list[Toolpath]:
    """Order independent parallel stripes; every deposited pass has one sign.

    Do not connect across openings or reorder islands. The existing G-code
    travel/retraction guards still handle the intervening non-extruding moves.
    """
    if not paths:
        return []
    # Paths of a skin layer share the same primary axis; alternate solid
    # direction controls whether the constant coordinate is theta or axial.
    first = paths[0].points
    axial_direction = abs(first[-1][0] - first[0][0]) <= 1e-8
    output = []
    for path in paths:
        points = list(path.points)
        start = points[0][1 if axial_direction else 0]
        end = points[-1][1 if axial_direction else 0]
        if start > end:
            points.reverse()
        output.append(replace(path, points=points))
    output.sort(key=lambda path: (
        round(path.points[0][0 if axial_direction else 1], 9),
        path.points[0][1 if axial_direction else 0],
    ))
    return output


def _ordered_paths_for_current_position(
    paths: list[Toolpath], state: GcodeState, config: SlicerConfig | None = None
) -> list[Toolpath]:
    if len(paths) > PATH_ORDER_OPTIMIZE_LIMIT:
        return paths

    remaining = list(paths)
    ordered = []
    shadow = GcodeState(
        theta_accum=state.theta_accum,
        theta_origin=state.theta_origin,
        native_y_origin_degrees=state.native_y_origin_degrees,
        axial=state.axial,
        radial=state.radial,
        extrusion=state.extrusion,
        retracted=state.retracted,
        retracted_mm=state.retracted_mm,
        retracted_speed_mm_per_min=state.retracted_speed_mm_per_min,
    )
    while remaining:
        best_index = 0
        best_path = _orient_path_for_current_position(remaining[0], shadow, config)
        best_distance = _travel_distance_to_point(best_path, best_path.points[0], shadow)
        for index, path in enumerate(remaining[1:], start=1):
            oriented = _orient_path_for_current_position(path, shadow, config)
            distance = _travel_distance_to_point(oriented, oriented.points[0], shadow)
            if distance < best_distance:
                best_index = index
                best_path = oriented
                best_distance = distance
        remaining.pop(best_index)
        ordered.append(best_path)
        shadow = _shadow_after_path(best_path, shadow)
    return ordered


def _should_retract_after_path(
    path: Toolpath,
    next_path: Toolpath | None,
    state: GcodeState,
    config: SlicerConfig,
    travel_region: BaseGeometry | None = None,
) -> bool:
    if config.retract_mm <= 0:
        return False
    if next_path is None:
        return True

    end_state = _shadow_after_path(path, state)
    next_path = _orient_path_for_current_position(next_path, end_state, config)
    travel_distance = _travel_distance_to_point(next_path, next_path.points[0], end_state)
    inside = False
    if travel_region is not None and travel_distance > 1e-8:
        start = (end_state.theta_accum, end_state.axial)
        destination = (_nearest_theta(next_path.points[0][0], end_state.theta_accum),
                       next_path.points[0][1])
        inside = _travel_is_contained(start, destination, next_path.radius, travel_region)
        if not inside:
            if travel_distance > config.route_max_unretracted_mm:
                return True  # A detour cannot be shorter than direct travel.
            route = _safe_island_travel_route(
                start, destination,
                next_path.radius, travel_region, config
            )
            if route is None:
                return True
            route_length = sum(
                cylindrical_distance(a,b,next_path.radius)
                for a,b in zip(
                    [(end_state.theta_accum,end_state.axial),*route],
                    route
                )
            )
            # A bounded in-island detour can remain primed for short moves,
            # but longer travel still retracts to avoid oozing.
            return route_length > config.route_max_unretracted_mm
    if travel_distance < config.nozzle_width * 2.0:
        return False
    fill_spacing = {
        "SOLID_INFILL": config.nozzle_width,
        "TOP_SKIN": config.dense_line_spacing,
        "INFILL_RADIAL": config.infill_spacing,
        "INFILL_AXIAL": config.infill_spacing,
        "SUPPORT": config.support_spacing,
        "SUPPORT_INTERFACE": config.support_interface_spacing or config.nozzle_width,
        "EXTERNAL_PERIMETER": config.nozzle_width,
        "INTERNAL_PERIMETER": config.nozzle_width,
    }
    if path.role in fill_spacing and next_path.role in fill_spacing and travel_region is not None:
        local_distance = max(fill_spacing[path.role], fill_spacing[next_path.role]) + config.nozzle_width * 0.5
        if travel_distance <= local_distance + 1e-8:
            # Only omit retraction when the entire connection stays in material.
            # The tiny residual allowance handles floating-point boundary drift,
            # including paths touching the periodic seam, without bridging holes.
            if inside:
                return False
    return True


def _ordered_layer_paths_with_islands(
    layer: LayerToolpaths, state: GcodeState, config: SlicerConfig | None = None
) -> list[tuple[Toolpath, int]]:
    # Keep the existing default order for callers of this geometry helper;
    # the outer/inner swap happens *within each island*, not across islands.
    roles = list(ROLE_ORDER)
    if config is not None and config.wall_order == "inner-first":
        external = roles.index("EXTERNAL_PERIMETER")
        internal = roles.index("INTERNAL_PERIMETER")
        roles[external], roles[internal] = roles[internal], roles[external]
    ordered = []
    shadow = state
    groups = list(layer.island_paths) if layer.island_paths else [layer.paths_by_role]
    remaining = list(enumerate(groups))
    while remaining:
        def entry_distance(group):
            first_role = next((role for role in roles if group.get(role)), None)
            if first_role is None:
                return math.inf
            candidates = [_orient_path_for_current_position(path, shadow, config) for path in group[first_role]]
            return min(_travel_distance_to_point(path, path.points[0], shadow) for path in candidates)

        nearest = (min(range(len(remaining)), key=lambda i: entry_distance(remaining[i][1]))
                   if len(remaining) <= PATH_ORDER_OPTIMIZE_LIMIT else 0)
        group_index, group = remaining.pop(nearest)
        for role in roles:
            role_paths = (
                _monotonic_skin_paths(group.get(role, []))
                if role == "TOP_SKIN" and config is not None and config.monotonic_top_skin
                else _ordered_paths_for_current_position(group.get(role, []), shadow, config)
            )
            for path in role_paths:
                path = _orient_path_for_current_position(path, shadow, config)
                ordered.append((path, group_index))
                shadow = _shadow_after_path(path, shadow)
    return ordered


def _ordered_layer_paths(
    layer: LayerToolpaths, state: GcodeState, config: SlicerConfig | None = None
) -> list[Toolpath]:
    """Legacy public order-only helper; no change to path orientation/order."""
    return [path for path, _ in _ordered_layer_paths_with_islands(layer, state, config)]


def _retract(state: GcodeState, config: SlicerConfig) -> list[str]:
    if config.retract_mm <= 0 or state.retracted:
        return []
    state.extrusion -= config.retract_mm
    state.retracted = True
    state.retracted_mm = config.retract_mm
    state.retracted_speed_mm_per_min = config.retract_speed_mm_per_min
    return [_format_move("G1", e=state.extrusion, f=config.retract_speed_mm_per_min) + " ; retract"]


def _unretract(state: GcodeState, config: SlicerConfig) -> list[str]:
    if not state.retracted:
        return []
    state.extrusion += state.retracted_mm
    state.retracted = False
    state.retracted_mm = 0.0
    return [_format_move("G1", e=state.extrusion, f=state.retracted_speed_mm_per_min) + " ; unretract"]


def _unwrap_bidirectional_path_from_current(thetas: list[float], current_accum: float) -> list[float]:
    if not thetas:
        return []
    out = [_nearest_theta(thetas[0], current_accum)]
    previous_raw = thetas[0]
    for theta in thetas[1:]:
        raw_delta = theta - previous_raw
        if abs(raw_delta) <= 1e-9:
            raw_delta = 0.0
        out.append(out[-1] + raw_delta)
        previous_raw = theta
    return out


def _unwrap_toolpath_thetas(path: Toolpath, current_accum: float) -> list[float]:
    raw_thetas = [point[0] for point in path.points]
    allows_reverse = path.role in BIDIRECTIONAL_OUTPUT_ROLES
    if allows_reverse and any(b < a - 1e-9 for a, b in zip(raw_thetas, raw_thetas[1:])):
        return _unwrap_bidirectional_path_from_current(raw_thetas, current_accum)
    return unwrap_path_from_current(raw_thetas, current_accum)


def _min_emitted_path_length(path: Toolpath, config: SlicerConfig) -> float:
    """A length threshold selects whole features, never individual corners."""
    if path.role in PERIMETER_ROLES:
        return MIN_GEOMETRY_LENGTH
    if config.min_path_length <= 0:
        return MIN_GEOMETRY_LENGTH
    return config.min_path_length + MIN_OUTPUT_LENGTH_GUARD


def _role_print_settings(role: str, config: SlicerConfig) -> tuple[float, float]:
    """Physical print feed (mm/min), extrusion multiplier, for a path role.

    Per-layer cooling may already have reduced config.feed_mm_per_min.
    The role's speed is an upper bound, never a way to override cooling.
    """
    if role == "BRIDGE":
        return min(config.feed_mm_per_min, config.bridge_print_feed_mm_per_min), config.bridge_flow_ratio
    if role == "SUPPORT":
        maximum = config.support_print_feed_mm_per_min
        return min(config.feed_mm_per_min, maximum) if maximum is not None else config.feed_mm_per_min, config.support_flow_ratio
    if role == "SUPPORT_INTERFACE":
        maximum = config.support_interface_print_feed_mm_per_min
        return min(config.feed_mm_per_min, maximum) if maximum is not None else config.feed_mm_per_min, config.support_interface_flow_ratio
    return config.feed_mm_per_min, 1.0


def _ramped_extrusion_steps(
    pts: list[tuple[float, float]], radius: float, config: SlicerConfig,
    *, seam_end: bool = False,
) -> list[tuple[tuple[float, float], int, float]]:
    """Split only at bounded physical-mm feed ramps; preserve endpoints and E.

    Start ramp applies to all extrusions when configured. End ramp applies
    only to closed wall paths (caller sets seam_end). Feed changes do NOT
    reduce deposited volume, path length or source segment-width metadata.
    """
    lengths = [
        cylindrical_distance(a, b, radius)
        for a, b in zip(pts, pts[1:])
    ]
    total = sum(lengths)
    start_span = min(total, config.extrusion_feed_ramp_mm)
    end_span = min(total, config.seam_end_ramp_mm) if seam_end else 0.0
    if total <= MIN_GEOMETRY_LENGTH or (start_span <= 0 and end_span <= 0):
        return [(pt, i, 1.0) for i, pt in enumerate(pts[1:])]

    if start_span + end_span > total:
        # Both ramps remain monotonic; never overlay conflicting start and
        # finish tuning on a too-short closed loop.
        scale = total / (start_span + end_span)
        start_span *= scale
        end_span *= scale
    edges = []
    if start_span > 0:
        stages = max(1, min(4, int(start_span / 0.8)))
        edges.extend(start_span * j / stages for j in range(1, stages + 1))
    if end_span > 0:
        stages = max(1, min(4, int(end_span / 0.8)))
        edges.extend(total - end_span + end_span * j / stages
                     for j in range(stages))
    edges = sorted(set(x for x in edges if 0 < x < total))
    result = []
    traversed = 0.0
    eps = 1e-9
    for i, ((theta0, x0), (theta1, x1)) in enumerate(zip(pts, pts[1:])):
        length = lengths[i]
        if length <= MIN_GEOMETRY_LENGTH:
            result.append(((theta1, x1), i, 1.0))
            traversed += length
            continue
        splits = [
            boundary for boundary in edges
            if traversed + eps < boundary < traversed + length - eps
        ]
        cumulative = [0.0, *[p - traversed for p in splits], length]
        for begin, finish in zip(cumulative, cumulative[1:]):
            fraction = finish / length
            point = (
                theta0 + (theta1-theta0)*fraction,
                x0 + (x1-x0)*fraction,
            )
            midpoint = traversed + (begin+finish)*0.5
            speed_fraction = 1.0
            if start_span > 0 and midpoint < start_span:
                progress = midpoint / start_span
                speed_fraction = min(speed_fraction,
                    config.extrusion_feed_ramp_start_fraction +
                    (1.0-config.extrusion_feed_ramp_start_fraction)*progress)
            if end_span > 0 and midpoint > total - end_span:
                progress = (total - midpoint) / end_span
                speed_fraction = min(speed_fraction,
                    config.seam_end_ramp_fraction +
                    (1.0-config.seam_end_ramp_fraction)*progress)
            result.append((point, i, speed_fraction))
        traversed += length
    return result


def emit_toolpath(
    path: Toolpath,
    state: GcodeState,
    config: SlicerConfig,
    output_x_offset: float,
    *,
    retract_after: bool = True,
    perimeter_index: PerimeterIndex | None = None,
    travel_region: BaseGeometry | None = None,
) -> list[str]:
    if len(path.points) < 2:
        return []
    accum = _unwrap_toolpath_thetas(path, state.theta_accum)
    pts = list(zip(accum, [p[1] for p in path.points]))
    if path.segment_widths is None:
        pts = coalesce_collinear_points(pts, path.radius)
    elif (len(path.segment_widths) != len(pts)-1
          or any(not math.isfinite(w) or w < config.layer_height or w > config.nozzle_width
                 for w in path.segment_widths)):
        raise ValueError("Variable-width segments must be finite, printable, and match the path")
    if len(pts) < 2:
        return []
    total_length = sum(cylindrical_distance(a, b, path.radius) for a, b in zip(pts, pts[1:]))
    if total_length < _min_emitted_path_length(path, config):
        return []

    lines = [";start of path"]
    t0, x0 = pts[0]
    route = _safe_island_travel_route(
        (state.theta_accum, state.axial), (t0,x0),
        path.radius, travel_region, config
    )
    if route is None:
        lines.extend(_ensure_rotary_range(t0, state))
    z_hop = _should_z_hop_to_path(path, state, config, perimeter_index)
    if route is not None:
        length = sum(cylindrical_distance(a,b,path.radius)
                     for a,b in zip([(state.theta_accum,state.axial),*route],route))
        if length < config.z_hop_min_travel:
            z_hop = False
    if z_hop:
        # state.radial tracks physical nozzle radius (outer face of the slab),
        # while path.radius is the geometric slice sample/metric radius.
        state.radial += config.z_hop_height
        lines.append(_format_move("G0", z=_emitted_z(state.radial, config), f=config.travel_mm_per_min) + " ; z hop")
    # XY travel uses Marlin's mixed-unit planner metric too. Radial hop
    # motions are separate Z-only commands and retain their native mm/min.
    if route is None:
        travel_feed = _surface_feed_for_marlin(
            x0 - state.axial, t0 - state.theta_accum, path.radius, config.travel_mm_per_min
        )
        lines.append(_format_move("G0", x=_emitted_x(x0, output_x_offset), y_deg=_emitted_y_degrees(t0, state), f=travel_feed))
        state.theta_accum = t0
        state.axial = x0
    else:
        lines.append(f";ROUTED_TRAVEL:{len(route)}")
        for theta, axial in route:
            lines.extend(_ensure_rotary_range(theta, state))
            feed = _surface_feed_for_marlin(
                axial-state.axial, theta-state.theta_accum,
                path.radius, config.travel_mm_per_min
            )
            lines.append(_format_move("G0", x=_emitted_x(axial, output_x_offset),
                                      y_deg=_emitted_y_degrees(theta, state), f=feed))
            state.theta_accum = theta
            state.axial = axial
    if z_hop:
        state.radial -= config.z_hop_height
        lines.append(_format_move("G0", z=_emitted_z(state.radial, config), f=config.travel_mm_per_min) + " ; z unhop")
    lines.extend(_unretract(state, config))

    extruded = False
    role_feed, role_flow = _role_print_settings(path.role, config)
    for (t, x), index, ramp_factor in _ramped_extrusion_steps(
        pts, path.radius, config,
        seam_end=path.closed and path.role in PERIMETER_ROLES,
    ):
        # E must describe the move that will ACTUALLY be emitted. If a tiny
        # numerical segment is omitted, the next move starts at the last
        # emitted endpoint, not at the skipped point.
        dist = cylindrical_distance((state.theta_accum, state.axial), (t, x), path.radius)
        if dist <= MIN_GEOMETRY_LENGTH:
            continue
        lines.extend(_ensure_rotary_range(t, state))
        width_ratio = 1.0
        if path.segment_widths is not None:
            width = path.segment_widths[index]
            if config.bead_model == "rounded":
                bead_area = (width - config.layer_height)*config.layer_height + math.pi*config.layer_height**2/4
            else:
                bead_area = width*config.layer_height
            # Use the ratio with nominal bead area so any custom calibrated
            # extrusion_per_mm multiplier is preserved exactly.
            width_ratio = bead_area/config.bead_area
        state.extrusion += dist * config.extrusion_per_mm * role_flow * width_ratio
        extruded = True
        lines.append(
            _format_move(
                "G1",
                x=_emitted_x(x, output_x_offset),
                y_deg=_emitted_y_degrees(t, state),
                e=state.extrusion,
                f=_surface_feed_for_marlin(
                    x - state.axial, t - state.theta_accum, path.radius,
                    role_feed * ramp_factor
                ),
            )
        )
        state.theta_accum = t
        state.axial = x
    if extruded and retract_after:
        lines.extend(_retract(state, config))
    lines.append(";end of path")
    return lines


def _role_lengths(paths: list[Toolpath]) -> dict[str, float]:
    """Count only deposited physical-mm lengths, grouped by actual speed cap."""
    lengths: dict[str, float] = {}
    for path in paths:
        length = sum(
            cylindrical_distance(a, b, path.radius)
            for a, b in zip(path.points, path.points[1:])
        )
        if length > MIN_GEOMETRY_LENGTH:
            lengths[path.role] = lengths.get(path.role, 0.0) + length
    return lengths


def _cooling_feed_for_target(
    lengths: dict[str, float], config: SlicerConfig, target_seconds: float,
) -> tuple[float, float]:
    """Fastest base F meeting a conservative deposited-motion time target.

    Actual role caps (bridge/support/interface) must already count toward
    layer or island time. Ignoring them slows other features unnecessarily.
    Travel, retract, acceleration and junction stalls are intentionally
    excluded to maintain a conservative *lower bound* on elapsed time.
    """
    ceiling = config.feed_mm_per_min
    if not lengths:
        return ceiling, 0.0

    def duration(feed: float) -> float:
        trial = replace(config, feed_mm_per_min=feed)
        return sum(
            length * 60.0 / _role_print_settings(role, trial)[0]
            for role, length in lengths.items()
        )

    natural = duration(ceiling)
    if target_seconds <= 0 or natural >= target_seconds:
        return ceiling, natural
    floor = min(ceiling, config.cooling_min_print_speed_mm_per_min)
    maximum = duration(floor)
    if maximum < target_seconds or floor >= ceiling:
        return floor, maximum
    low, high = floor, ceiling
    for _ in range(35):
        middle = (low + high) * 0.5
        if duration(middle) >= target_seconds:
            low = middle
        else:
            high = middle
    return low, duration(low)


def _cooling_print_speed(layer: LayerToolpaths, config: SlicerConfig) -> tuple[float, float]:
    """Conservative layer print duration in physical cylindrical millimeters.

    Role-limited bridge and support moves count at their actual ceiling, so
    the unrestricted paths slow only when the whole layer needs it. No
    stationary dwell or guessed firmware acceleration is included.
    """
    lengths = _role_lengths([
        path for paths in layer.paths_by_role.values() for path in paths
    ])
    return _cooling_feed_for_target(lengths, config, config.cooling_min_layer_time_s)


def _island_cooling_print_speed(
    paths: list[Toolpath], config: SlicerConfig
) -> tuple[float, float]:
    """Apply the same role-aware conservative timing to each connected island."""
    return _cooling_feed_for_target(
        _role_lengths(paths), config, config.cooling_min_island_print_time_s
    )


_MATERIAL_FAN_PWM_CEILING = {
    "pla": 255, "petg": 180, "abs": 0, "tpu": 180, "unprofiled": 255,
}


def _effective_material_fan_pwm(config: SlicerConfig) -> int:
    """Opt-in fan ceiling; never increase the configured PWM."""
    if config.material_fan_policy == "manual":
        return config.fan_speed
    return min(config.fan_speed, _MATERIAL_FAN_PWM_CEILING[config.material_profile])


def _layer_fan_pwm(height_mm: float, config: SlicerConfig) -> int:
    maximum = _effective_material_fan_pwm(config)
    if maximum <= 0 or height_mm < config.fan_start_height - 1e-8:
        return 0
    if config.fan_ramp_height_mm <= 0:
        return maximum
    fraction = min(1.0, max(0.0,
                          (height_mm - config.fan_start_height) / config.fan_ramp_height_mm))
    return round(maximum * fraction)


def emit_gcode(layers: list[LayerToolpaths], config: SlicerConfig) -> str:
    config.validate()
    state = GcodeState(radial=config.bed_radius)
    output_x_offset = _output_x_offset(layers, config)
    out = [
        ";start of slice",
        ";MARLIN_REQUIRES:POWER_LOSS_RECOVERY (G92.9 native Y reset)",
        f";LAYER_HEIGHT:{config.layer_height}",
        f";LINE_WIDTH:{config.nozzle_width}",
        f";MATERIAL_PROFILE:{config.material_profile}",
        f";LINEAR_ADVANCE:{('Marlin M900 K' + f'{config.linear_advance_k:g}') if config.linear_advance_k is not None else 'disabled'}",
        f";EXTRUSION_FEED_RAMP_MM:{config.extrusion_feed_ramp_mm:g}",
        f";EXTRUSION_FEED_RAMP_START_FRACTION:{config.extrusion_feed_ramp_start_fraction:g}",
        f";NOZZLE_DIAMETER:{config.nozzle_diameter if config.nozzle_diameter is not None else 'unknown'}",
        f";BEAD_MODEL:{config.bead_model}",
        f";WALL_ORDER:{config.wall_order}",
        f";ROUTE_TRAVEL_INSIDE_ISLAND:{str(config.route_travel_inside_island).lower()}",
        f";ROUTE_MAX_EXTRA_MM:{config.route_max_extra_mm:g}",
        f";ROUTE_MAX_FACTOR:{config.route_max_factor:g}",
        f";MONOTONIC_TOP_SKIN:{str(config.monotonic_top_skin).lower()}",
        f";GAP_FILL_ENABLED:{str(config.gap_fill_enabled).lower()}",
        f";GAP_FILL_MIN_WIDTH_MM:{config.gap_fill_min_width_mm:g}",
        f";SEAM_POLICY:{config.seam_policy}",
        f";SEAM_ANGLE_DEG:{config.seam_angle_deg:g}",
        f";SEAM_EXCLUDE_DEG:{config.seam_exclude_deg}",
        f";SEAM_END_RAMP_MM:{config.seam_end_ramp_mm:g}",
        f";SEAM_END_RAMP_FRACTION:{config.seam_end_ramp_fraction:g}",
        f";BEAD_AREA_MM2:{config.bead_area:.6f}",
        f";WALL_SUPPORT_REPORT:{str(config.wall_support_report).lower()}",
        f";BRIDGE_MODE:{config.bridge_mode}",
        f";BRIDGE_MAX_SPAN_MM:{config.bridge_max_span_mm:g}",
        f";BRIDGE_FEED_MM_MIN:{config.bridge_print_feed_mm_per_min:g}",
        f";BRIDGE_FLOW_RATIO:{config.bridge_flow_ratio:g}",
        f";BRIDGE_FAN_PWM:{config.bridge_fan_pwm}",
        f";SUPPORT_PRINT_FEED_MM_MIN:{config.support_print_feed_mm_per_min}",
        f";SUPPORT_INTERFACE_PRINT_FEED_MM_MIN:{config.support_interface_print_feed_mm_per_min}",
        f";SUPPORT_FLOW_RATIO:{config.support_flow_ratio:g}",
        f";SUPPORT_INTERFACE_FLOW_RATIO:{config.support_interface_flow_ratio:g}",
        f";DENSE_SPACING_MM:{config.dense_line_spacing:.6f}",
        f";BUILD_CYLINDER_RADIUS:{config.bed_radius:g}",
        f";RETRACT_MM:{config.retract_mm:g}",
        f";RETRACT_SPEED_MM_PER_S:{config.retract_speed_mm_per_min / 60:g}",
        f";HOTEND_TEMP:{config.hotend_temp}",
        f";PRINT_FEED:{config.feed_mm_per_min:g}",  # requested physical surface mm/min
        f";TRAVEL_FEED:{config.travel_mm_per_min:g}",  # requested physical surface mm/min for XY travel
        f";FAN_SPEED:{config.fan_speed}",
        f";MATERIAL_FAN_POLICY:{config.material_fan_policy}",
        f";EFFECTIVE_FAN_PWM:{_effective_material_fan_pwm(config)}",
        f";FAN_START_HEIGHT:{config.fan_start_height:g}",
        f";FAN_RAMP_HEIGHT_MM:{config.fan_ramp_height_mm:g}",
        f";COOLING_MIN_LAYER_TIME_S:{config.cooling_min_layer_time_s:g}",
        f";COOLING_MIN_ISLAND_PRINT_TIME_S:{config.cooling_min_island_print_time_s:g}",
        f";COOLING_MIN_PRINT_SPEED_MM_MIN:{config.cooling_min_print_speed_mm_per_min:g}",
        "G90 ; absolute positioning",
        "M82 ; absolute extrusion mode",
        *_reset_rotary_axis(state, reset_extrusion=True),
        "M106 S0 ; part cooling off for adhesion",
        f"M104 S{config.hotend_temp}",
        f"M109 S{config.hotend_temp}",
    ]
    if config.linear_advance_k is not None:
        # Only an explicitly verified Marlin LIN_ADVANCE build may see M900.
        # This is a firmware-level pressure offset; absolute E accounting
        # and per-move speed-corrected motion remain untouched.
        out.append(f"M900 K{config.linear_advance_k:g} ; enable calibrated Marlin linear advance")
    if config.retraction_sweep:
        out.append(config.retraction_sweep.header().rstrip())

    fan_pwm = 0
    current_linear_advance_k = config.linear_advance_k
    for layer in layers:
        nozzle_radius = layer.deposition_radius(config.layer_height)
        if not math.isfinite(nozzle_radius) or nozzle_radius < config.bed_radius:
            raise ValueError("Nozzle radius must be finite and at/above bed radius")
        deposited_thickness = (layer.deposited_thickness if layer.deposited_thickness is not None
                               else config.layer_height)
        if not math.isfinite(deposited_thickness) or deposited_thickness <= 0:
            raise ValueError("Deposited layer thickness must be finite and positive")
        height = nozzle_radius - config.bed_radius
        overrides = config.retraction_sweep.overrides(height) if config.retraction_sweep else {}
        layer_config = config.for_deposited_thickness(deposited_thickness)
        if overrides:
            layer_config = replace(layer_config, **overrides)
        cooling_speed, estimated_print_time = _cooling_print_speed(layer, layer_config)
        if abs(cooling_speed - layer_config.feed_mm_per_min) > 1e-10:
            layer_config = replace(layer_config, feed_mm_per_min=cooling_speed)
        out.append(f";LAYER:{layer.layer_index}")
        if config.cooling_min_layer_time_s > 0:
            out.append(f";COOLING_PRINT_DURATION_ESTIMATE_S:{estimated_print_time:.4f}")
            out.append(f";COOLING_PRINT_SPEED_MM_MIN:{cooling_speed:.4f}")
        out.append(f";SAMPLE_RADIUS:{layer.radius:.4f}")
        out.append(f";DEPOSITED_THICKNESS:{deposited_thickness:.4f}")
        out.append(f";Z:{nozzle_radius:.4f}")
        if (config.retraction_sweep is not None
                and config.retraction_sweep.parameter == "linear-advance-k"):
            requested_k = layer_config.linear_advance_k
            if requested_k != current_linear_advance_k:
                # The swept K changes between layers only, after the base
                # M900 has already been emitted on a declared-capable build.
                # Filament E and feed are never modified by calibration.
                out.append(f"M900 K{requested_k:g} ; calibrated linear advance radial band")
                current_linear_advance_k = requested_k
        if config.retraction_sweep:
            out.append(config.retraction_sweep.comment(height).rstrip())
        desired_pwm = _layer_fan_pwm(height, layer_config)
        if desired_pwm != fan_pwm:
            note = "part cooling on" if desired_pwm == _effective_material_fan_pwm(layer_config) else "part cooling ramp"
            out.append(f"M106 S{desired_pwm} ; {note}")
            fan_pwm = desired_pwm
        out.append(_format_move("G0", z=_emitted_z(nozzle_radius, layer_config), f=layer_config.travel_mm_per_min))
        # Reset E every layer, but leave Y in its current coordinate frame.
        # Repeated native/workspace Y rebases are unnecessary for a compact
        # rotary print and expose it to firmware-specific G92 behavior. The
        # rotary frame is rebased by _ensure_rotary_range only when a move
        # would otherwise exceed the Marlin precision budget.
        state.extrusion = 0.0
        out.append("G92 E0 ; reset extrusion accumulator for layer")
        state.radial = nozzle_radius
        perimeter_index = (
            _layer_perimeter_index(layer)
            if layer_config.z_hop_height > 0 and layer_config.z_hop_on_perimeter_crossing
            else None
        )
        travel_region = periodic_copies(layer.travel_region) if layer.travel_region is not None else None
        tagged_paths = _ordered_layer_paths_with_islands(layer, state, layer_config)
        groups_by_id: dict[int, list[Toolpath]] = {}
        for path, group_index in tagged_paths:
            groups_by_id.setdefault(group_index, []).append(path)
        island_configs: dict[int, SlicerConfig] = {}
        if layer_config.cooling_min_island_print_time_s > 0:
            for group_index, group_paths in groups_by_id.items():
                feed, duration = _island_cooling_print_speed(group_paths, layer_config)
                if feed + 1e-8 < layer_config.feed_mm_per_min:
                    island_configs[group_index] = replace(layer_config, feed_mm_per_min=feed)
                out.append(f";COOLING_ISLAND:{group_index} PRINT_DURATION_ESTIMATE_S:{duration:.4f} SPEED_MM_MIN:{feed:.4f}")
        for index, (path, group_index) in enumerate(tagged_paths):
            path_config = island_configs.get(group_index, layer_config)
            path = _orient_path_for_current_position(path, state, path_config)
            next_path = tagged_paths[index + 1][0] if index + 1 < len(tagged_paths) else None
            retract_after = _should_retract_after_path(path, next_path, state, path_config, travel_region)
            out.append(f";TYPE:{path.role}")
            if (path_config.wall_support_report and path.role in PERIMETER_ROLES
                    and layer.previous_occupied is not None):
                # Diagnostics only: independent of travel ordering, print
                # feed, E, seam policy, gravity and support generation.
                wall = measure_wall_support(
                    path, layer.previous_occupied, path_config.nozzle_width
                )
                out.append(
                    f";WALL_SUPPORT:{wall.category} COVERED_FRACTION:{wall.covered_fraction:.4f}"
                    f" MAX_UNSUPPORTED_MM:{wall.max_unsupported_run_mm:.4f}"
                    f" TOTAL_MM:{wall.total_length_mm:.4f}"
                )
            bridge_fan = (
                path.role == "BRIDGE" and path_config.bridge_mode == "print"
                and path_config.bridge_fan_pwm is not None
                and path_config.bridge_fan_pwm != fan_pwm
            )
            if bridge_fan:
                out.append(f"M106 S{path_config.bridge_fan_pwm} ; bridge fan override")
            out.extend(
                emit_toolpath(
                    path,
                    state,
                    path_config,
                    output_x_offset,
                    retract_after=retract_after,
                    perimeter_index=perimeter_index,
                    travel_region=travel_region,
                )
            )
            if bridge_fan:
                out.append(f"M106 S{fan_pwm} ; restore layer fan after bridge")

    if layers:
        # Z is radial lift above the mandrel. Park relative to the last emitted
        # radius, not at a fixed Z that could drive the nozzle into a tall print.
        out.extend(
            [
                "G90 ; absolute positioning for end lift",
                _format_move("G0", z=_emitted_z(state.radial, config) + 10.0, f=600.0)
                + " ; end lift 10 mm radially away from print",
                "M400 ; wait for end lift before disabling steppers",
            ]
        )

    if config.linear_advance_k is not None:
        # Return to the known no-advance state after the print rather than
        # surprising subsequent jobs with a persistent modified K factor.
        out.append("M900 K0 ; clear job-scoped linear advance")
    out.extend(
        [
            "M104 S0 ; hotend off",
            "M106 S0 ; fan off",
            "M84 ; disable steppers",
            ";end of slice",
        ]
    )
    return "\n".join(out) + "\n"


X_RE = re.compile(r"\bX(-?\d+(?:\.\d+)?)")
Y_RE = re.compile(r"\bY(-?\d+(?:\.\d+)?)")
E_RE = re.compile(r"\bE(-?\d+(?:\.\d+)?)")


def parse_extrusion_segments(gcode: str) -> list[tuple[int, str, float, float, float, float]]:
    layer = -9999
    role = "UNKNOWN"
    x = None
    y = None
    x_offset = 0.0
    y_offset = 0.0
    segments = []
    for line in gcode.splitlines():
        if line.startswith(";LAYER:"):
            try:
                layer = int(line.split(":", 1)[1])
            except ValueError:
                layer = -9999
            continue
        if line.startswith(";TYPE:"):
            role = line.split(":", 1)[1]
            continue
        if line.startswith("G92 "):
            mx = X_RE.search(line)
            my = Y_RE.search(line)
            if mx:
                raw_x = float(mx.group(1))
                if x is None:
                    x = raw_x
                    x_offset = 0.0
                else:
                    x_offset = x - raw_x
            if my:
                raw_y = float(my.group(1))
                if y is None:
                    y = raw_y
                    y_offset = 0.0
                else:
                    y_offset = y - raw_y
            continue
        if not line.startswith(("G0 ", "G1 ")):
            continue
        mx = X_RE.search(line)
        my = Y_RE.search(line)
        me = E_RE.search(line)
        nx = float(mx.group(1)) + x_offset if mx else x
        ny = float(my.group(1)) + y_offset if my else y
        moved_xy = (
            x is not None
            and y is not None
            and nx is not None
            and ny is not None
            and (abs(nx - x) > 1e-9 or abs(ny - y) > 1e-9)
        )
        if line.startswith("G1 ") and me and moved_xy:
            segments.append((layer, role, x, y, nx, ny))
        x, y = nx, ny
    return segments
