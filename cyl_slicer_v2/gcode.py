from __future__ import annotations

from dataclasses import dataclass, replace
import math
import re

from shapely.geometry import LineString, Point
from shapely.geometry.base import BaseGeometry
from shapely.strtree import STRtree

from .config import SlicerConfig
from .coordinates import coalesce_collinear_points, cylindrical_distance, unwrap_path_from_current
from .models import GcodeState, LayerToolpaths, Toolpath
from .toolpaths import MIN_GEOMETRY_LENGTH, MIN_OUTPUT_LENGTH_GUARD, PERIMETER_ROLES, ROLE_ORDER


TWO_PI = 2.0 * math.pi
PATH_ORDER_OPTIMIZE_LIMIT = 1024
BIDIRECTIONAL_OUTPUT_ROLES = PERIMETER_ROLES | frozenset({"SUPPORT", "SUPPORT_INTERFACE"})


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


def _emitted_y_degrees(theta: float, state: GcodeState) -> float:
    return math.degrees(theta - state.theta_origin)


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
    if path.role in PERIMETER_ROLES:
        return False
    if path.closed or len(path.points) < 2:
        return False
    theta0 = path.points[0][0]
    return all(abs(point[0] - theta0) <= tol for point in path.points[1:])


def _orient_path_for_current_position(path: Toolpath, state: GcodeState) -> Toolpath:
    if not _can_reverse_for_output(path):
        return path
    start_distance = _travel_distance_to_point(path, path.points[0], state)
    end_distance = _travel_distance_to_point(path, path.points[-1], state)
    if end_distance + 1e-9 >= start_distance:
        return path
    return Toolpath(role=path.role, radius=path.radius, points=list(reversed(path.points)), closed=path.closed)


def _shadow_after_path(path: Toolpath, state: GcodeState) -> GcodeState:
    accum = _unwrap_toolpath_thetas(path, state.theta_accum)
    if not accum:
        return state
    return GcodeState(
        theta_accum=accum[-1],
        theta_origin=state.theta_origin,
        axial=path.points[-1][1],
        radial=state.radial,
        extrusion=state.extrusion,
        retracted=state.retracted,
        retracted_mm=state.retracted_mm,
        retracted_speed_mm_per_min=state.retracted_speed_mm_per_min,
    )


def _ordered_paths_for_current_position(paths: list[Toolpath], state: GcodeState) -> list[Toolpath]:
    if len(paths) > PATH_ORDER_OPTIMIZE_LIMIT:
        return paths

    remaining = list(paths)
    ordered = []
    shadow = GcodeState(
        theta_accum=state.theta_accum,
        theta_origin=state.theta_origin,
        axial=state.axial,
        radial=state.radial,
        extrusion=state.extrusion,
        retracted=state.retracted,
        retracted_mm=state.retracted_mm,
        retracted_speed_mm_per_min=state.retracted_speed_mm_per_min,
    )
    while remaining:
        best_index = 0
        best_path = _orient_path_for_current_position(remaining[0], shadow)
        best_distance = _travel_distance_to_point(best_path, best_path.points[0], shadow)
        for index, path in enumerate(remaining[1:], start=1):
            oriented = _orient_path_for_current_position(path, shadow)
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
) -> bool:
    if config.retract_mm <= 0:
        return False
    if next_path is None:
        return True

    end_state = _shadow_after_path(path, state)
    travel_distance = _travel_distance_to_point(next_path, next_path.points[0], end_state)
    return travel_distance >= config.nozzle_width * 2.0


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


def _min_emitted_segment_length(path: Toolpath, config: SlicerConfig) -> float:
    if path.role in PERIMETER_ROLES:
        return MIN_GEOMETRY_LENGTH
    if config.min_path_length <= 0:
        return 0.0
    return config.min_path_length + MIN_OUTPUT_LENGTH_GUARD


def emit_toolpath(
    path: Toolpath,
    state: GcodeState,
    config: SlicerConfig,
    output_x_offset: float,
    *,
    retract_after: bool = True,
    perimeter_index: PerimeterIndex | None = None,
) -> list[str]:
    if len(path.points) < 2:
        return []
    accum = _unwrap_toolpath_thetas(path, state.theta_accum)
    pts = list(zip(accum, [p[1] for p in path.points]))
    pts = coalesce_collinear_points(pts, path.radius)

    lines = [";start of path"]
    t0, x0 = pts[0]
    z_hop = _should_z_hop_to_path(path, state, config, perimeter_index)
    if z_hop:
        state.radial = path.radius + config.z_hop_height
        lines.append(_format_move("G0", z=_emitted_z(state.radial, config), f=config.travel_mm_per_min) + " ; z hop")
    lines.append(_format_move("G0", x=_emitted_x(x0, output_x_offset), y_deg=_emitted_y_degrees(t0, state), f=config.travel_mm_per_min))
    state.theta_accum = t0
    state.axial = x0
    if z_hop:
        state.radial = path.radius
        lines.append(_format_move("G0", z=_emitted_z(state.radial, config), f=config.travel_mm_per_min) + " ; z unhop")
    lines.extend(_unretract(state, config))

    extruded = False
    min_length = _min_emitted_segment_length(path, config)
    for (t_prev, x_prev), (t, x) in zip(pts[:-1], pts[1:]):
        dist = cylindrical_distance((t_prev, x_prev), (t, x), path.radius)
        if dist <= min_length:
            continue
        state.extrusion += dist * config.extrusion_per_mm
        extruded = True
        lines.append(
            _format_move(
                "G1",
                x=_emitted_x(x, output_x_offset),
                y_deg=_emitted_y_degrees(t, state),
                e=state.extrusion,
                f=config.feed_mm_per_min,
            )
        )
        state.theta_accum = t
        state.axial = x
    if extruded and retract_after:
        lines.extend(_retract(state, config))
    lines.append(";end of path")
    return lines


def emit_gcode(layers: list[LayerToolpaths], config: SlicerConfig) -> str:
    config.validate()
    state = GcodeState(radial=config.bed_radius)
    output_x_offset = _output_x_offset(layers, config)
    out = [
        ";start of slice",
        f";LAYER_HEIGHT:{config.layer_height}",
        f";LINE_WIDTH:{config.nozzle_width}",
        f";BUILD_CYLINDER_RADIUS:{config.bed_radius:g}",
        f";RETRACT_MM:{config.retract_mm:g}",
        f";RETRACT_SPEED_MM_PER_S:{config.retract_speed_mm_per_min / 60:g}",
        f";HOTEND_TEMP:{config.hotend_temp}",
        f";FAN_SPEED:{config.fan_speed}",
        f";FAN_START_HEIGHT:{config.fan_start_height:g}",
        "M82 ; absolute extrusion mode",
        "G92 E0 Y0",
        "M106 S0 ; part cooling off for adhesion",
        f"M104 S{config.hotend_temp}",
        f"M109 S{config.hotend_temp}",
    ]
    if config.retraction_sweep:
        out.append(config.retraction_sweep.header().rstrip())

    fan_started = False
    for layer in layers:
        height = layer.radius - config.bed_radius
        overrides = config.retraction_sweep.overrides(height) if config.retraction_sweep else {}
        layer_config = replace(config, **overrides) if overrides else config
        out.append(f";LAYER:{layer.layer_index}")
        out.append(f";Z:{layer.radius:.4f}")
        if config.retraction_sweep:
            out.append(config.retraction_sweep.comment(height).rstrip())
        if not fan_started and config.fan_speed > 0 and height >= config.fan_start_height - 1e-8:
            out.append(f"M106 S{config.fan_speed} ; part cooling on")
            fan_started = True
        out.append(_format_move("G0", z=_emitted_z(layer.radius, layer_config), f=layer_config.travel_mm_per_min))
        out.append("G92 E0 Y0 ; reset E/Y accumulators for layer")
        state.extrusion = 0.0
        state.theta_origin = state.theta_accum
        state.radial = layer.radius
        perimeter_index = (
            _layer_perimeter_index(layer)
            if layer_config.z_hop_height > 0 and layer_config.z_hop_on_perimeter_crossing
            else None
        )
        for role in ROLE_ORDER:
            role_paths = _ordered_paths_for_current_position(layer.paths_by_role.get(role, []), state)
            for index, path in enumerate(role_paths):
                path = _orient_path_for_current_position(path, state)
                next_path = role_paths[index + 1] if index + 1 < len(role_paths) else None
                retract_after = _should_retract_after_path(path, next_path, state, layer_config)
                out.append(f";TYPE:{role}")
                out.extend(
                    emit_toolpath(
                        path,
                        state,
                        layer_config,
                        output_x_offset,
                        retract_after=retract_after,
                        perimeter_index=perimeter_index,
                    )
                )

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
