from __future__ import annotations

import math
import re

from .config import SlicerConfig
from .coordinates import cylindrical_distance, unwrap_path_from_current
from .models import GcodeState, LayerToolpaths, Toolpath
from .toolpaths import ROLE_ORDER


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


def emit_toolpath(path: Toolpath, state: GcodeState, config: SlicerConfig) -> list[str]:
    if len(path.points) < 2:
        return []
    raw_thetas = [p[0] for p in path.points]
    accum = unwrap_path_from_current(raw_thetas, state.theta_accum)
    pts = list(zip(accum, [p[1] for p in path.points]))
    pts = _split_long_axial_segments(pts, config.max_axial_segment)

    lines = [";start of path"]
    t0, x0 = pts[0]
    lines.append(_format_move("G0", x=x0, y_deg=_emitted_y_degrees(t0, state), f=config.travel_mm_per_min))
    state.theta_accum = t0
    state.axial = x0

    for (t_prev, x_prev), (t, x) in zip(pts[:-1], pts[1:]):
        dist = cylindrical_distance((t_prev, x_prev), (t, x), path.radius)
        if dist <= config.min_path_length:
            continue
        state.extrusion += dist * config.extrusion_per_mm
        lines.append(
            _format_move(
                "G1",
                x=x,
                y_deg=_emitted_y_degrees(t, state),
                e=state.extrusion,
                f=config.feed_mm_per_min,
            )
        )
        state.theta_accum = t
        state.axial = x
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
                out.append(f";TYPE:{role}")
                out.extend(emit_toolpath(path, state, config))

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
        if not line.startswith(("G0 ", "G1 ")):
            continue
        mx = X_RE.search(line)
        my = Y_RE.search(line)
        me = E_RE.search(line)
        nx = float(mx.group(1)) if mx else x
        ny = float(my.group(1)) if my else y
        if line.startswith("G1 ") and me and x is not None and y is not None and nx is not None and ny is not None:
            segments.append((layer, role, x, y, nx, ny))
        x, y = nx, ny
    return segments
