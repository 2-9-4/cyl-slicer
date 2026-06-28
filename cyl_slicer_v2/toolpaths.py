from __future__ import annotations

import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import shapely
from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPolygon, Polygon, box
from shapely.ops import unary_union

from .config import SlicerConfig
from .coordinates import TWO_PI, cylindrical_distance
from .models import LayerSlice, LayerToolpaths, Toolpath
from .periodic import periodic_copies

ROLE_ORDER = [
    "EXTERNAL_PERIMETER",
    "INTERNAL_PERIMETER",
    "SOLID_INFILL",
    "SUPPORT",
    "SUPPORT_INTERFACE",
    "INFILL_RADIAL",
    "INFILL_AXIAL",
]

_WORKER_CONFIG: SlicerConfig | None = None
MIN_OUTPUT_LENGTH_GUARD = 1e-3


def _min_printable_length(config: SlicerConfig) -> float:
    if config.min_path_length <= 0:
        return 0.0
    return config.min_path_length + MIN_OUTPUT_LENGTH_GUARD


def _to_metric(geom, radius: float):
    return shapely.transform(geom, lambda p: p * [radius, 1.0])


def _from_metric(geom, radius: float):
    return shapely.transform(geom, lambda p: p * [1.0 / radius, 1.0])


def _clean_polygonal(geom):
    if geom.is_empty:
        return GeometryCollection()
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    geom = geom.buffer(0)
    if geom.is_empty or geom.area <= 1e-12:
        return GeometryCollection()
    return geom


def _safe_polygon_difference(a, b):
    if a is None or a.is_empty:
        return GeometryCollection()
    if b is None or b.is_empty:
        return _clean_polygonal(a)
    try:
        return _clean_polygonal(a.difference(b))
    except shapely.errors.GEOSException:
        left = _clean_polygonal(a)
        right = _clean_polygonal(b)
        if left.is_empty:
            return GeometryCollection()
        if right.is_empty:
            return left
        return _clean_polygonal(shapely.difference(left, right, grid_size=1e-8))


def _polygon_parts(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return list(geom.geoms)
    if isinstance(geom, GeometryCollection):
        result = []
        for part in geom.geoms:
            result.extend(_polygon_parts(part))
        return result
    return []


def _line_parts(geom) -> list[LineString]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, LineString):
        return [geom]
    if isinstance(geom, MultiLineString):
        return list(geom.geoms)
    if isinstance(geom, GeometryCollection):
        result = []
        for part in geom.geoms:
            result.extend(_line_parts(part))
        return result
    return []


def _merged_axial_intervals(geom, tol: float = 1e-8) -> list[tuple[float, float]]:
    intervals = sorted((poly.bounds[1], poly.bounds[3]) for poly in _polygon_parts(geom))
    merged: list[tuple[float, float]] = []
    for lo, hi in intervals:
        if hi - lo <= tol:
            continue
        if not merged or lo > merged[-1][1] + tol:
            merged.append((lo, hi))
        else:
            prev_lo, prev_hi = merged[-1]
            merged[-1] = (prev_lo, max(prev_hi, hi))
    return merged


def _raft_footprint_with_axial_gap(footprint, config: SlicerConfig):
    theta_min, z_min, theta_max, z_max = footprint.bounds
    raft = _clean_polygonal(box(theta_min, z_min, theta_max, z_max))
    gap = config.raft_axial_gap
    if gap <= 0 or z_max - z_min <= gap:
        return raft

    source_gaps = []
    intervals = _merged_axial_intervals(footprint)
    for (_, prev_hi), (next_lo, _) in zip(intervals, intervals[1:]):
        if next_lo - prev_hi >= gap:
            source_gaps.append((prev_hi, next_lo))

    if source_gaps:
        gap_lo, gap_hi = max(source_gaps, key=lambda item: item[1] - item[0])
        center = (gap_lo + gap_hi) * 0.5
        cut_lo = center - gap * 0.5
        cut_hi = center + gap * 0.5
    else:
        cut_lo = z_max - gap
        cut_hi = z_max

    cut = box(theta_min, cut_lo, theta_max, cut_hi)
    relieved = _safe_polygon_difference(raft, cut)
    return relieved if not relieved.is_empty else raft


def _is_window_boundary_edge(
    a: tuple[float, float],
    b: tuple[float, float],
    theta_window: tuple[float, float] | None,
    tol: float = 1e-8,
) -> bool:
    if theta_window is None:
        return False
    lower, upper = theta_window
    if upper - lower <= tol:
        return False
    return (abs(a[0] - lower) <= tol and abs(b[0] - lower) <= tol) or (
        abs(a[0] - upper) <= tol and abs(b[0] - upper) <= tol
    )


def _ring_to_paths(
    ring,
    role: str,
    radius: float,
    config: SlicerConfig,
    theta_window: tuple[float, float] | None = None,
) -> list[Toolpath]:
    pts = [(float(t), float(z)) for t, z in ring.coords]
    if _is_full_theta_window(theta_window):
        return _ring_to_window_runs(pts, role, radius, config, theta_window)
    return _ring_to_edge_paths(pts, role, radius, config, theta_window)


def _is_full_theta_window(theta_window: tuple[float, float] | None, tol: float = 1e-8) -> bool:
    if theta_window is None:
        return False
    lower, upper = theta_window
    return upper - lower >= TWO_PI - tol


def _theta_is_nondecreasing(pts: list[tuple[float, float]], tol: float = 1e-9) -> bool:
    return all(b[0] + tol >= a[0] for a, b in zip(pts, pts[1:]))


def _ring_to_window_runs(
    pts: list[tuple[float, float]],
    role: str,
    radius: float,
    config: SlicerConfig,
    theta_window: tuple[float, float] | None,
) -> list[Toolpath]:
    runs: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    for a, b in zip(pts[:-1], pts[1:]):
        if _is_window_boundary_edge(a, b, theta_window):
            if len(current) >= 2:
                runs.append(current)
            current = []
            continue
        if not current:
            current = [a, b]
        elif abs(current[-1][0] - a[0]) <= 1e-9 and abs(current[-1][1] - a[1]) <= 1e-9:
            current.append(b)
        else:
            if len(current) >= 2:
                runs.append(current)
            current = [a, b]
    if len(current) >= 2:
        runs.append(current)
    if len(runs) > 1 and _points_close(runs[-1][-1], runs[0][0]):
        runs[0] = [*runs[-1], *runs[0][1:]]
        runs.pop()

    paths = []
    for run in runs:
        if not _theta_is_nondecreasing(run):
            run = list(reversed(run))
        path = _make_path(role, radius, run, False, config)
        if path is not None:
            paths.append(path)
    return sorted(paths, key=lambda p: (p.points[0][0], p.points[0][1], p.points[-1][0], p.points[-1][1]))


def _ring_to_edge_paths(
    pts: list[tuple[float, float]],
    role: str,
    radius: float,
    config: SlicerConfig,
    theta_window: tuple[float, float] | None = None,
) -> list[Toolpath]:
    paths = []
    for a, b in zip(pts[:-1], pts[1:]):
        if _is_window_boundary_edge(a, b, theta_window):
            continue
        edge = [a, b] if b[0] >= a[0] else [b, a]
        path = _make_path(role, radius, edge, False, config)
        if path is not None:
            paths.append(path)
    return sorted(paths, key=lambda p: (p.points[0][0], p.points[0][1], p.points[-1][0], p.points[-1][1]))


def _make_path(role: str, radius: float, pts: list[tuple[float, float]], closed: bool, config: SlicerConfig) -> Toolpath | None:
    if len(pts) < 2:
        return None
    min_length = _min_printable_length(config)
    cleaned = [pts[0]]
    for pt in pts[1:]:
        if cylindrical_distance(cleaned[-1], pt, radius) >= min_length:
            cleaned.append(pt)
    if len(cleaned) < 2:
        return None
    return Toolpath(role=role, radius=radius, points=cleaned, closed=closed)


def perimeter_paths(
    geom,
    role: str,
    radius: float,
    config: SlicerConfig,
    theta_window: tuple[float, float] | None = None,
) -> list[Toolpath]:
    paths = []
    for poly in _polygon_parts(geom):
        paths.extend(_ring_to_paths(poly.exterior, role, radius, config, theta_window))
        for interior in poly.interiors:
            paths.extend(_ring_to_paths(interior, role, radius, config, theta_window))
    return paths


def inset_geometry(
    geom,
    radius: float,
    distance: float,
    *,
    periodic: bool = False,
    theta_window: tuple[float, float] | None = None,
):
    if geom.is_empty:
        return geom
    source = periodic_copies(geom, turns=1) if periodic else geom
    metric = _to_metric(source, radius)
    inset = metric.buffer(-distance, join_style="mitre")
    if inset.is_empty:
        return GeometryCollection()
    raw = _from_metric(shapely.make_valid(inset), radius)
    if periodic:
        lower, upper = theta_window if theta_window is not None else (geom.bounds[0], geom.bounds[2])
        raw = raw.intersection(box(lower, geom.bounds[1], upper, geom.bounds[3]))
    return raw


def line_grid_lines(geom, radius: float, spacing: float, *, direction: str) -> list[LineString]:
    if geom.is_empty:
        return []
    theta_min, z_min, theta_max, z_max = geom.bounds
    lines = []
    if direction == "axial":
        step = spacing / radius
        start = math.floor(theta_min / step) * step
        theta = start
        while theta <= theta_max + step * 0.5:
            lines.append(LineString([(theta, z_min - spacing), (theta, z_max + spacing)]))
            theta += step
    elif direction == "radial":
        start = math.floor(z_min / spacing) * spacing
        z = start
        while z <= z_max + spacing * 0.5:
            lines.append(LineString([(theta_min, z), (theta_max, z)]))
            z += spacing
    else:
        raise ValueError(f"unknown line direction {direction!r}")

    result = []
    for source in lines:
        inter = source.intersection(geom)
        for line in _line_parts(inter):
            if line.length > 0:
                result.append(line)
    if direction == "axial":
        return sorted(result, key=lambda line: (line.coords[0][0], line.bounds[1], line.bounds[3]))
    return sorted(result, key=lambda line: (line.coords[0][1], line.bounds[0], line.bounds[2]))


def line_grid_paths(geom, role: str, radius: float, spacing: float, config: SlicerConfig, *, direction: str) -> list[Toolpath]:
    result = []
    for line in line_grid_lines(geom, radius, spacing, direction=direction):
        pts = [(float(t), float(z)) for t, z in line.coords]
        path = _make_path(role, radius, pts, False, config)
        if path is not None:
            result.append(path)
    return result


def _connector_is_covered(geom, theta: float, z0: float, z1: float, tol: float = 1e-8) -> bool:
    if abs(z1 - z0) <= tol:
        return True
    lo, hi = sorted((z0, z1))
    candidates = []
    for candidate in (theta, theta % TWO_PI):
        if not any(abs(candidate - existing) <= tol for existing in candidates):
            candidates.append(candidate)
    for candidate in candidates:
        connector = LineString([(candidate, lo), (candidate, hi)])
        if geom.covers(connector):
            return True
    return False


def _path_length(path: Toolpath) -> float:
    return sum(cylindrical_distance(a, b, path.radius) for a, b in zip(path.points, path.points[1:]))


def _straight_path_interval(path: Toolpath, tol: float = 1e-8):
    if path.closed or len(path.points) != 2:
        return None
    (t0, z0), (t1, z1) = path.points
    if abs(t1 - t0) <= tol:
        lo, hi = sorted((z0, z1))
        return ("axial", round(t0, 8)), lo, hi
    if abs(z1 - z0) <= tol:
        lo, hi = sorted((t0, t1))
        return ("radial", round(z0, 8)), lo, hi
    return None


def _remove_contained_straight_paths(paths: list[Toolpath], tol: float = 1e-8) -> list[Toolpath]:
    groups: dict[tuple[str, float], list[tuple[int, float, float]]] = defaultdict(list)
    for index, path in enumerate(paths):
        interval = _straight_path_interval(path, tol)
        if interval is None:
            continue
        key, lo, hi = interval
        groups[key].append((index, lo, hi))

    remove = set()
    for intervals in groups.values():
        max_hi = -math.inf
        for index, lo, hi in sorted(intervals, key=lambda item: (item[1], -item[2], item[0])):
            if hi <= max_hi + tol:
                remove.add(index)
            else:
                max_hi = hi

    return [path for index, path in enumerate(paths) if index not in remove]


def _points_close(a: tuple[float, float], b: tuple[float, float], tol: float = 1e-8) -> bool:
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def _can_reverse_join(path: Toolpath) -> bool:
    return path.role not in {"EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"}


def _join_touching_paths(paths: list[Toolpath]) -> list[Toolpath]:
    remaining = [path for path in paths if not path.closed]
    joined = [path for path in paths if path.closed]
    while remaining:
        chain = remaining.pop(0)
        changed = True
        while changed:
            changed = False
            for index, candidate in enumerate(remaining):
                if chain.role != candidate.role or abs(chain.radius - candidate.radius) > 1e-9:
                    continue
                if _points_close(chain.points[-1], candidate.points[0]):
                    chain = Toolpath(chain.role, chain.radius, [*chain.points, *candidate.points[1:]], False)
                elif _points_close(candidate.points[-1], chain.points[0]):
                    chain = Toolpath(chain.role, chain.radius, [*candidate.points[:-1], *chain.points], False)
                elif _can_reverse_join(chain) and _points_close(chain.points[-1], candidate.points[-1]):
                    chain = Toolpath(chain.role, chain.radius, [*chain.points, *reversed(candidate.points[:-1])], False)
                elif _can_reverse_join(chain) and _points_close(candidate.points[0], chain.points[0]):
                    chain = Toolpath(chain.role, chain.radius, [*reversed(candidate.points[1:]), *chain.points], False)
                else:
                    continue
                remaining.pop(index)
                changed = True
                break
        joined.append(chain)
    return joined


def _postprocess_role_paths(paths: list[Toolpath], config: SlicerConfig) -> list[Toolpath]:
    kept = [path for path in paths if _path_length(path) >= _min_printable_length(config)]
    return _remove_contained_straight_paths(_join_touching_paths(kept))


def _radial_line_specs(geom, radius: float, spacing: float, config: SlicerConfig) -> list[tuple[float, float, float]]:
    specs = []
    line_tol = max(1e-6, 1e-6 / max(radius, 1e-9))
    max_theta_span = math.radians(config.support_meander_max_angle_deg)
    for line in line_grid_lines(geom, radius, spacing, direction="radial"):
        coords = [(float(t), float(z)) for t, z in line.coords]
        if len(coords) < 2:
            continue
        start = coords[0]
        end = coords[-1]
        if end[0] < start[0]:
            start, end = end, start
        if abs(end[1] - start[1]) <= line_tol:
            span = end[0] - start[0]
            split_count = max(1, int(math.ceil(span / max_theta_span)))
            for split in range(split_count):
                split_start = start[0] + span * split / split_count
                split_end = start[0] + span * (split + 1) / split_count
                if (split_end - split_start) * radius >= _min_printable_length(config):
                    specs.append((split_start, split_end, start[1]))
    return sorted(specs, key=lambda item: (item[0], item[1], item[2]))


def _can_connect_radial_specs(
    geom,
    previous: tuple[float, float, float],
    current: tuple[float, float, float],
    *,
    previous_forward: bool,
    current_forward: bool,
    spacing: float,
    radius: float,
) -> bool:
    prev_start, prev_end, prev_z = previous
    cur_start, cur_end, cur_z = current
    if cur_z - prev_z > spacing * 1.5 + 1e-8:
        return False

    previous_exit = prev_end if previous_forward else prev_start
    current_entry = cur_start if current_forward else cur_end
    endpoint_tol = max(1e-5, 1e-5 / max(radius, 1e-9))
    if abs(previous_exit - current_entry) > endpoint_tol:
        return False
    return _connector_is_covered(geom, previous_exit, prev_z, cur_z)


def meandered_radial_grid_paths(
    geom,
    role: str,
    radius: float,
    spacing: float,
    config: SlicerConfig,
) -> list[Toolpath]:
    if geom.is_empty:
        return []

    specs = _radial_line_specs(geom, radius, spacing, config)
    result = []
    groups: list[list[tuple[float, float, float]]] = []
    current: list[tuple[float, float, float]] = []
    for spec in specs:
        if not current:
            current = [spec]
            continue
        index = len(current)
        can_connect = _can_connect_radial_specs(
            geom,
            current[-1],
            spec,
            previous_forward=(index - 1) % 2 == 0,
            current_forward=index % 2 == 0,
            spacing=spacing,
            radius=radius,
        )
        if can_connect:
            current.append(spec)
        else:
            groups.append(current)
            current = [spec]
    if current:
        groups.append(current)

    for group in groups:
        if not group:
            continue
        pts = []
        for index, (start, end, z) in enumerate(group):
            line_pts = [(start, z), (end, z)] if index % 2 == 0 else [(end, z), (start, z)]
            pts.append(line_pts[0])
            pts.append(line_pts[1])
        path = _make_path(role, radius, pts, False, config)
        if path is not None:
            result.append(path)
    return result


def build_layer_toolpaths(layer: LayerSlice, config: SlicerConfig) -> LayerToolpaths:
    paths: dict[str, list[Toolpath]] = defaultdict(list)
    model = layer.occupied
    radius = layer.radius

    if not model.is_empty:
        periodic = layer.is_full_wrap
        theta_window = layer.theta_window
        paths["EXTERNAL_PERIMETER"].extend(perimeter_paths(model, "EXTERNAL_PERIMETER", radius, config, theta_window))
        last_inset = model
        for wall in range(1, config.wall_count):
            inset = inset_geometry(
                model,
                radius,
                config.nozzle_width * wall,
                periodic=periodic,
                theta_window=theta_window,
            )
            if inset.is_empty:
                break
            paths["INTERNAL_PERIMETER"].extend(perimeter_paths(inset, "INTERNAL_PERIMETER", radius, config, theta_window))
            last_inset = inset

        solid_source = unary_union([layer.top_surfaces, layer.bottom_surfaces])
        solid_source = solid_source.intersection(last_inset) if not solid_source.is_empty else solid_source
        paths["SOLID_INFILL"].extend(
            line_grid_paths(solid_source, "SOLID_INFILL", radius, config.nozzle_width, config, direction="axial")
        )

        infill_source = inset_geometry(
            last_inset,
            radius,
            config.nozzle_width * 0.5,
            periodic=periodic,
            theta_window=theta_window,
        )
        if not solid_source.is_empty and not infill_source.is_empty:
            infill_source = infill_source.difference(solid_source)
            if not infill_source.is_valid:
                infill_source = shapely.make_valid(infill_source)
            infill_source = infill_source.buffer(0)
        paths["INFILL_RADIAL"].extend(
            line_grid_paths(infill_source, "INFILL_RADIAL", radius, config.infill_spacing, config, direction="radial")
        )
        paths["INFILL_AXIAL"].extend(
            line_grid_paths(infill_source, "INFILL_AXIAL", radius, config.infill_spacing, config, direction="axial")
        )

    support = _clean_polygonal(layer.support)
    interface_source = layer.support_interface if layer.support_interface is not None else GeometryCollection()
    interface = _clean_polygonal(interface_source)
    if not interface.is_empty:
        interface = interface.intersection(support) if not support.is_empty else GeometryCollection()
        interface = _clean_polygonal(interface)
        if not interface.is_empty:
            interface_spacing = config.support_interface_spacing or config.nozzle_width
            interface_direction = layer.support_interface_direction or ("axial" if layer.index % 2 == 0 else "radial")
            if interface_direction == "radial":
                paths["SUPPORT_INTERFACE"].extend(
                    meandered_radial_grid_paths(interface, "SUPPORT_INTERFACE", radius, interface_spacing, config)
                )
            else:
                paths["SUPPORT_INTERFACE"].extend(
                    line_grid_paths(
                        interface,
                        "SUPPORT_INTERFACE",
                        radius,
                        interface_spacing,
                        config,
                        direction=interface_direction,
                    )
                )
            support = _safe_polygon_difference(support, interface)

    if not support.is_empty:
        support_spacing = (config.raft_grid_spacing or config.nozzle_width) if layer.is_raft else config.support_spacing
        if layer.is_raft:
            paths["SUPPORT"].extend(
                line_grid_paths(
                    support,
                    "SUPPORT",
                    radius,
                    support_spacing,
                    config,
                    direction=layer.raft_direction or "radial",
                )
            )
        else:
            radial_support = meandered_radial_grid_paths(support, "SUPPORT", radius, support_spacing, config)
            paths["SUPPORT"].extend(radial_support)
            if radial_support:
                cross_spacing = config.support_cross_spacing or support_spacing * 4.0
                paths["SUPPORT"].extend(
                    line_grid_paths(support, "SUPPORT", radius, cross_spacing, config, direction="axial")
                )

    paths = {
        role: processed
        for role, role_paths in paths.items()
        if (processed := _postprocess_role_paths(role_paths, config))
    }
    return LayerToolpaths(layer_index=layer.index, radius=radius, paths_by_role=paths)


def _init_toolpath_worker(config: SlicerConfig) -> None:
    global _WORKER_CONFIG
    _WORKER_CONFIG = config


def _build_layer_toolpaths_worker(layer: LayerSlice) -> LayerToolpaths:
    if _WORKER_CONFIG is None:
        raise RuntimeError("toolpath worker was not initialized")
    return build_layer_toolpaths(layer, _WORKER_CONFIG)


def gap_support_slices(layers: list[LayerSlice], config: SlicerConfig) -> list[LayerSlice]:
    if not config.generate_gap_support or not layers or any(layer.index < 0 for layer in layers):
        return []

    first_model = next((layer for layer in layers if not layer.occupied.is_empty), None)
    if first_model is None:
        return []

    start = config.bed_radius + config.layer_height * 0.5
    if start >= first_model.radius - 1e-9:
        return []

    pieces = []
    if not first_model.support.is_empty:
        pieces.append(first_model.support)
    if not first_model.occupied.is_empty:
        pieces.append(first_model.occupied)
    footprint = _clean_polygonal(unary_union(pieces)) if pieces else GeometryCollection()
    if footprint.is_empty:
        return []
    raft_footprint = _raft_footprint_with_axial_gap(footprint, config)

    radii = []
    radius = start
    while radius < first_model.radius - 1e-9:
        radii.append(float(radius))
        radius += config.layer_height

    count = len(radii)
    empty = GeometryCollection()
    return [
        LayerSlice(
            index=i - count,
            radius=radius,
            occupied=empty,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=raft_footprint if i < config.raft_layers else footprint,
            theta_window=first_model.theta_window,
            is_full_wrap=first_model.is_full_wrap,
            is_raft=i < config.raft_layers,
            raft_direction=("radial" if i % 2 == 0 else "axial") if i < config.raft_layers else None,
        )
        for i, radius in enumerate(radii)
    ]


def build_toolpaths(layers: list[LayerSlice], config: SlicerConfig) -> list[LayerToolpaths]:
    all_layers = [*gap_support_slices(layers, config), *layers]
    worker_count = min(config.toolpath_workers, len(all_layers))
    if worker_count <= 1:
        return [build_layer_toolpaths(layer, config) for layer in all_layers]

    chunksize = max(1, len(all_layers) // (worker_count * 4))
    with ProcessPoolExecutor(
        max_workers=worker_count,
        initializer=_init_toolpath_worker,
        initargs=(config,),
    ) as executor:
        return list(executor.map(_build_layer_toolpaths_worker, all_layers, chunksize=chunksize))
