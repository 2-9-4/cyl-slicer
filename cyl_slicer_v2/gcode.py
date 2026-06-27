from __future__ import annotations

import math
import re

from .config import SlicerConfig
from .coordinates import cylindrical_distance, unwrap_path_from_current
from .models import GcodeState, LayerToolpaths, Toolpath
from .toolpaths import MIN_OUTPUT_LENGTH_GUARD, ROLE_ORDER


TWO_PI = 2.0 * math.pi


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


def _can_reverse_for_output(path: Toolpath, tol: float = 1e-8) -> bool:
    if path.role in {"EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"}:
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


def _retract(state: GcodeState, config: SlicerConfig) -> list[str]:
    if config.retract_mm <= 0 or state.retracted:
        return []
    state.extrusion -= config.retract_mm
    state.retracted = True
    return [_format_move("G1", e=state.extrusion, f=config.travel_mm_per_min) + " ; retract"]


def _unretract(state: GcodeState, config: SlicerConfig) -> list[str]:
    if config.retract_mm <= 0 or not state.retracted:
        return []
    state.extrusion += config.retract_mm
    state.retracted = False
    return [_format_move("G1", e=state.extrusion, f=config.travel_mm_per_min) + " ; unretract"]


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
    allows_reverse = path.role in {"SUPPORT", "SUPPORT_INTERFACE"}
    if allows_reverse and any(b < a - 1e-9 for a, b in zip(raw_thetas, raw_thetas[1:])):
        return _unwrap_bidirectional_path_from_current(raw_thetas, current_accum)
    return unwrap_path_from_current(raw_thetas, current_accum)


def emit_toolpath(path: Toolpath, state: GcodeState, config: SlicerConfig, output_x_offset: float) -> list[str]:
    if len(path.points) < 2:
        return []
    accum = _unwrap_toolpath_thetas(path, state.theta_accum)
    pts = list(zip(accum, [p[1] for p in path.points]))
    pts = _split_long_axial_segments(pts, config.max_axial_segment)

    lines = [";start of path"]
    t0, x0 = pts[0]
    lines.append(_format_move("G0", x=_emitted_x(x0, output_x_offset), y_deg=_emitted_y_degrees(t0, state), f=config.travel_mm_per_min))
    state.theta_accum = t0
    state.axial = x0
    lines.extend(_unretract(state, config))

    extruded = False
    min_length = config.min_path_length + MIN_OUTPUT_LENGTH_GUARD if config.min_path_length > 0 else 0.0
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
    if extruded:
        lines.extend(_retract(state, config))
    lines.append(";end of path")
    return lines


def _split_long_axial_segments(
    pts: list[tuple[float, float]],
    max_axial_segment: float,
) -> list[tuple[float, float]]:
    if len(pts) < 2:
        return pts
    out = [pts[0]]
    for (t0, x0), (t1, x1) in zip(pts[:-1], pts[1:]):
        dx = x1 - x0
        if abs(dx) > max_axial_segment:
            steps = int(math.ceil(abs(dx) / max_axial_segment))
            for i in range(1, steps + 1):
                frac = i / steps
                out.append((t0 + (t1 - t0) * frac, x0 + dx * frac))
        else:
            out.append((t1, x1))
    return out


def emit_gcode(layers: list[LayerToolpaths], config: SlicerConfig) -> str:
    config.validate()
    state = GcodeState(radial=config.bed_radius)
    output_x_offset = _output_x_offset(layers, config)
    out = [
        ";start of slice",
        f";LAYER_HEIGHT:{config.layer_height}",
        f";LINE_WIDTH:{config.nozzle_width}",
        "M82 ; absolute extrusion mode",
        "G92 E0 Y0",
        f"M104 S{config.hotend_temp}",
        f"M109 S{config.hotend_temp}",
    ]

    for layer in layers:
        out.append(f";LAYER:{layer.layer_index}")
        out.append(f";Z:{layer.radius:.4f}")
        radial_lift = layer.radius - config.bed_radius
        out.append(_format_move("G0", z=radial_lift, f=config.travel_mm_per_min))
        out.append("G92 E0 Y0 ; reset E/Y accumulators for layer")
        state.extrusion = 0.0
        state.theta_origin = state.theta_accum
        state.radial = layer.radius
        for role in ROLE_ORDER:
            for path in layer.paths_by_role.get(role, []):
                path = _orient_path_for_current_position(path, state)
                out.append(f";TYPE:{role}")
                out.extend(emit_toolpath(path, state, config, output_x_offset))

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
