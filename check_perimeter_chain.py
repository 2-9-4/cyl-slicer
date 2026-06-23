#!/usr/bin/env python3
"""Validate emitted multi-turn perimeters plus solid/support infill spans."""

from __future__ import annotations

import math
import re
import runpy
import sys
from pathlib import Path

import shapely

X_RE = re.compile(r"X(-?\d+(?:\.\d+)?)")
Y_RE = re.compile(r"Y(-?\d+(?:\.\d+)?)")


def parse_gcode_paths(gcode_path: Path, section_name: str) -> dict[int, list[list[tuple[float, float]]]]:
    text = gcode_path.read_text().splitlines()
    layer = None
    section = None
    in_path = False
    points: list[tuple[float, float]] = []
    by_layer: dict[int, list[list[tuple[float, float]]]] = {}

    def flush() -> None:
        nonlocal points
        if in_path and layer is not None and section == section_name and len(points) >= 2:
            by_layer.setdefault(layer, []).append(points)
        points = []

    for line in text:
        if line.startswith(";LAYER:"):
            flush()
            try:
                layer = int(line.split(":", 1)[1])
            except ValueError:
                layer = None
            section = None
            in_path = False
        elif line.startswith(";TYPE:"):
            flush()
            section = line.split(":", 1)[1]
            in_path = False
        elif line.startswith(";start of path"):
            points = []
            in_path = True
        elif line.startswith(";end of path"):
            flush()
            in_path = False
        elif in_path and line.startswith(("G0 ", "G1 ")):
            mx = X_RE.search(line)
            my = Y_RE.search(line)
            if mx and my:
                points.append((float(mx.group(1)), float(my.group(1))))

    flush()
    return by_layer


def _same_path(a: list[tuple[float, float]], b: list[tuple[float, float]], tol: float = 1e-6) -> bool:
    if len(a) != len(b):
        return False
    return all(abs(ax - bx) <= tol and abs(ay - by) <= tol for (ax, ay), (bx, by) in zip(a, b))


def collapse_repeated_path_block(paths: list[list[tuple[float, float]]]) -> list[list[tuple[float, float]]]:
    n = len(paths)
    for block_size in range(1, n // 2 + 1):
        if n % block_size != 0:
            continue
        repeats = n // block_size
        if repeats < 2:
            continue
        first = paths[:block_size]
        if all(
            _same_path(first[idx % block_size], paths[idx])
            for idx in range(n)
        ):
            return first
    return paths


def load_env(stl_path: str) -> dict:
    old_argv = sys.argv[:]
    try:
        sys.argv = ["slicer.py", stl_path]
        return runpy.run_path("slicer.py", run_name="__main__")
    finally:
        sys.argv = old_argv


def raw_external_spans(env: dict) -> dict[int, list[float]]:
    layer_data = env["layer_data"]
    polygon_to_wrapped_lines = env["polygon_to_wrapped_lines"]
    shapely_to_paths = env["shapely_to_paths"]

    spans: dict[int, list[float]] = {}
    for layer, data in sorted(layer_data.items()):
        radius = data["radius"]
        raw_lines = polygon_to_wrapped_lines(data["external_peris"])
        paths = [p for p in shapely_to_paths(raw_lines, radius) if p is not None]
        spans[layer] = []
        for path in paths:
            if not path.points or len(path.points) < 2:
                continue
            theta = [pt[0] for pt in path.points]
            spans[layer].append(max(theta) - min(theta))
    return spans


def expected_solid_spans(env: dict) -> dict[int, float]:
    layer_data = env["layer_data"]
    slicing_params = env["slicing_params"]
    keep_canonical_polygons = env["keep_canonical_polygons"]

    expected: dict[int, float] = {}
    for layer, data in sorted(layer_data.items()):
        current_r = data["radius"]

        def do_transform(x):
            return shapely.transform(x, lambda p: p * [current_r, 1])

        def do_transform_inverse(x):
            return shapely.transform(x, lambda p: p * [1 / current_r, 1])

        current_ext_print = keep_canonical_polygons(data["external_peris"])
        current_top_print = keep_canonical_polygons(data.get("top_surfaces", shapely.GeometryCollection()))
        current_bottom_print = keep_canonical_polygons(data.get("bottom_surfaces", shapely.GeometryCollection()))

        final_perimeter = do_transform(current_ext_print)
        for wall_number in range(1, slicing_params["wall_count"]):
            final_perimeter = do_transform(current_ext_print).buffer(
                -slicing_params["nozzle_size"] * wall_number,
                join_style="mitre",
            )

        final_top_surface = do_transform(current_top_print).buffer(
            slicing_params["nozzle_size"] * slicing_params["wall_count"],
            join_style="mitre",
        ).intersection(final_perimeter)
        final_bottom_surface = do_transform(current_bottom_print).buffer(
            slicing_params["nozzle_size"] * slicing_params["wall_count"],
            join_style="mitre",
        ).intersection(final_perimeter)

        solid_target = do_transform_inverse(final_top_surface.union(final_bottom_surface))
        if solid_target.is_empty:
            continue
        bounds = shapely.bounds(solid_target)
        expected[layer] = bounds[2] - bounds[0]

    return expected


def actual_solid_spans(env: dict) -> dict[int, float]:
    spans: dict[int, float] = {}
    for layer, data in sorted(env["layer_data"].items()):
        solid = data.get("solid_infill", [])
        if not solid:
            continue
        theta_min = min(shapely.bounds(seg)[0] for seg in solid if not seg.is_empty)
        theta_max = max(shapely.bounds(seg)[2] for seg in solid if not seg.is_empty)
        spans[layer] = theta_max - theta_min
    return spans


def expected_support_spans(env: dict) -> dict[int, float]:
    layer_data = env["layer_data"]
    slicing_params = env["slicing_params"]
    keep_canonical_polygons = env["keep_canonical_polygons"]
    gap_support_layers = env["gap_support_layers"]

    expected: dict[int, float] = {}

    n_gap = len(gap_support_layers)
    for gap_layer_idx, gap_layer in enumerate(gap_support_layers):
        layer_num = gap_layer_idx - n_gap
        support_print = keep_canonical_polygons(gap_layer["support"])
        if support_print.is_empty:
            continue
        bounds = shapely.bounds(support_print)
        expected[layer_num] = bounds[2] - bounds[0]

    gap_n = slicing_params.get("support_gap_layers", 0)

    for layer, data in sorted(layer_data.items()):
        support_geom = data.get("support", shapely.GeometryCollection())
        if support_geom.is_empty:
            continue

        is_gap_layer = gap_n > 0 and layer_data.get(layer + gap_n, {}).get("support", shapely.GeometryCollection()).is_empty
        if is_gap_layer:
            continue

        support_print = keep_canonical_polygons(support_geom)
        if support_print.is_empty:
            continue
        bounds = shapely.bounds(support_print)
        expected[layer] = bounds[2] - bounds[0]

    return expected


def actual_support_spans(env: dict) -> dict[int, float]:
    spans: dict[int, float] = {}

    for gap_layer_idx, gap_layer in enumerate(env["gap_support_layers"]):
        layer_num = gap_layer_idx - len(env["gap_support_layers"])
        support_infill = gap_layer.get("support_infill", [])
        if not support_infill:
            continue
        theta_min = min(shapely.bounds(seg)[0] for seg in support_infill if not seg.is_empty)
        theta_max = max(shapely.bounds(seg)[2] for seg in support_infill if not seg.is_empty)
        spans[layer_num] = theta_max - theta_min

    for layer, data in sorted(env["layer_data"].items()):
        support_infill = data.get("support_infill", [])
        if not support_infill:
            continue
        theta_min = min(shapely.bounds(seg)[0] for seg in support_infill if not seg.is_empty)
        theta_max = max(shapely.bounds(seg)[2] for seg in support_infill if not seg.is_empty)
        spans[layer] = theta_max - theta_min

    return spans


def main() -> int:
    stl_path = sys.argv[1] if len(sys.argv) > 1 else "helix.stl"
    gcode_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("output.gcode")

    env = load_env(stl_path)
    expected_spans = raw_external_spans(env)
    expected_solid = expected_solid_spans(env)
    actual_solid = actual_solid_spans(env)
    expected_support = expected_support_spans(env)
    actual_support = actual_support_spans(env)
    gcode_external = parse_gcode_paths(gcode_path, "EXTERNAL_PERIMETER")

    failures: list[str] = []
    perimeter_span_tol_deg = 1.0
    solid_span_tol_deg = 8.0
    support_span_tol_deg = 8.0

    for layer, expected in sorted(expected_spans.items()):
        actual_paths = collapse_repeated_path_block(gcode_external.get(layer, []))
        if len(expected) != len(actual_paths):
            failures.append(
                f"layer {layer}: expected {len(expected)} external paths, found {len(actual_paths)} in G-code"
            )
            continue

        for idx, (expected_span_rad, path) in enumerate(zip(expected, actual_paths)):
            actual_span_deg = path[-1][1] - path[0][1]
            expected_span_deg = math.degrees(expected_span_rad)
            if abs(actual_span_deg - expected_span_deg) > perimeter_span_tol_deg:
                failures.append(
                    f"layer {layer}: perimeter path {idx} span mismatch "
                    f"(expected {expected_span_deg:.3f} deg, got {actual_span_deg:.3f} deg)"
                )
            min_dy = min((y1 - y0 for (_, y0), (_, y1) in zip(path, path[1:])), default=0.0)
            if min_dy < -1e-6:
                failures.append(f"layer {layer}: perimeter path {idx} moves backward in Y")

    for layer, expected_span_rad in sorted(expected_solid.items()):
        expected_span_deg = math.degrees(expected_span_rad)
        if expected_span_deg <= 360.0 + solid_span_tol_deg:
            continue
        actual_span_rad = actual_solid.get(layer)
        if actual_span_rad is None:
            failures.append(
                f"layer {layer}: expected multi-turn solid infill span {expected_span_deg:.3f} deg, found none"
            )
            continue
        actual_span_deg = math.degrees(actual_span_rad)
        if abs(actual_span_deg - expected_span_deg) > solid_span_tol_deg:
            failures.append(
                f"layer {layer}: solid infill span mismatch "
                f"(expected {expected_span_deg:.3f} deg, got {actual_span_deg:.3f} deg)"
            )

    for layer, expected_span_rad in sorted(expected_support.items()):
        expected_span_deg = math.degrees(expected_span_rad)
        actual_span_rad = actual_support.get(layer)
        if actual_span_rad is None:
            failures.append(
                f"layer {layer}: expected support infill span {expected_span_deg:.3f} deg, found none"
            )
            continue
        actual_span_deg = math.degrees(actual_span_rad)
        if expected_span_deg > 360.0 + support_span_tol_deg and abs(actual_span_deg - expected_span_deg) > support_span_tol_deg:
            failures.append(
                f"layer {layer}: support infill span mismatch "
                f"(expected {expected_span_deg:.3f} deg, got {actual_span_deg:.3f} deg)"
            )

    if failures:
        print("FAIL: toolpath emission regressed")
        for line in failures:
            print(line)
        return 1

    print("PASS: external perimeters plus solid/support multi-turn infill match source geometry")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
