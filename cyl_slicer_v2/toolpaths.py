from __future__ import annotations

import math
from collections import defaultdict

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


def _to_metric(geom, radius: float):
    return shapely.transform(geom, lambda p: p * [radius, 1.0])


def _from_metric(geom, radius: float):
    return shapely.transform(geom, lambda p: p * [1.0 / radius, 1.0])


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


def _ring_to_paths(ring, role: str, radius: float, config: SlicerConfig) -> list[Toolpath]:
    pts = [(float(t), float(z)) for t, z in ring.coords]
    paths = []
    for a, b in zip(pts[:-1], pts[1:]):
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


def perimeter_paths(geom, role: str, radius: float, config: SlicerConfig) -> list[Toolpath]:
    paths = []
    for poly in _polygon_parts(geom):
        paths.extend(_ring_to_paths(poly.exterior, role, radius, config))
        for interior in poly.interiors:
            paths.extend(_ring_to_paths(interior, role, radius, config))
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
        paths["EXTERNAL_PERIMETER"].extend(perimeter_paths(model, "EXTERNAL_PERIMETER", radius, config))
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
            paths["INTERNAL_PERIMETER"].extend(perimeter_paths(inset, "INTERNAL_PERIMETER", radius, config))
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
        paths["SUPPORT"].extend(
            line_grid_paths(layer.support, "SUPPORT", radius, config.support_spacing, config, direction="axial")
        )
        paths["SUPPORT"].extend(
            line_grid_paths(layer.support, "SUPPORT", radius, config.support_spacing, config, direction="radial")
        )

    return LayerToolpaths(layer_index=layer.index, radius=radius, paths_by_role=dict(paths))


def build_toolpaths(layers: list[LayerSlice], config: SlicerConfig) -> list[LayerToolpaths]:
    return [build_layer_toolpaths(layer, config) for layer in layers]
