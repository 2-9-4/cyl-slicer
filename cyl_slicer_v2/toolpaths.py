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
    "INFILL_RADIAL",
    "INFILL_AXIAL",
]

_WORKER_CONFIG: SlicerConfig | None = None


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
    cleaned = [pts[0]]
    for pt in pts[1:]:
        if cylindrical_distance(cleaned[-1], pt, radius) >= config.min_path_length:
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


def line_grid_paths(geom, role: str, radius: float, spacing: float, config: SlicerConfig, *, direction: str) -> list[Toolpath]:
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
            pts = [(float(t), float(z)) for t, z in line.coords]
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

    if not layer.support.is_empty:
        support_spacing = config.nozzle_width if layer.is_raft else config.support_spacing
        if layer.index >= 0 or layer.is_raft:
            paths["SUPPORT"].extend(
                line_grid_paths(layer.support, "SUPPORT", radius, support_spacing, config, direction="axial")
            )
        paths["SUPPORT"].extend(
            line_grid_paths(layer.support, "SUPPORT", radius, support_spacing, config, direction="radial")
        )

    return LayerToolpaths(layer_index=layer.index, radius=radius, paths_by_role=dict(paths))


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
    theta_min, z_min, theta_max, z_max = footprint.bounds
    raft_footprint = _clean_polygonal(box(theta_min, z_min, theta_max, z_max))

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
