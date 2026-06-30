from __future__ import annotations

import math

import shapely
from shapely import affinity
from shapely.geometry import GeometryCollection, LineString, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

from .coordinates import TWO_PI


def translate_theta(geom: BaseGeometry, delta: float) -> BaseGeometry:
    if geom.is_empty:
        return geom
    return affinity.translate(geom, xoff=delta)


def _safe_union(geoms) -> BaseGeometry:
    parts = [part for part in geoms if part is not None and not part.is_empty]
    if not parts:
        return GeometryCollection()
    if len(parts) == 1:
        return parts[0]
    try:
        return unary_union(parts)
    except shapely.errors.GEOSException:
        repaired = [shapely.make_valid(part) for part in parts]
        repaired = [part for part in repaired if not part.is_empty]
        if not repaired:
            return GeometryCollection()
        return shapely.union_all(repaired, grid_size=1e-8)


def _safe_intersection(a: BaseGeometry, b: BaseGeometry) -> BaseGeometry:
    if a.is_empty or b.is_empty:
        return GeometryCollection()
    try:
        out = a.intersection(b)
    except shapely.errors.GEOSException:
        out = shapely.intersection(shapely.make_valid(a), shapely.make_valid(b), grid_size=1e-8)
    if not out.is_valid:
        out = shapely.make_valid(out)
    return out


def periodic_copies(geom: BaseGeometry, turns: int = 1) -> BaseGeometry:
    if geom.is_empty:
        return geom
    return _safe_union(translate_theta(geom, k * TWO_PI) for k in range(-turns, turns + 1))


def is_full_wrap(geom: BaseGeometry, samples: int = 96) -> bool:
    if geom.is_empty:
        return False
    z_min, z_max = geom.bounds[1], geom.bounds[3]
    if z_max <= z_min:
        return False
    for i in range(samples):
        theta = (i + 0.5) * TWO_PI / samples
        probe = LineString([(theta, z_min - 1.0), (theta, z_max + 1.0)])
        if geom.intersection(probe).is_empty:
            return False
    return True


def choose_empty_cut(geom: BaseGeometry, samples: int = 144) -> float | None:
    if geom.is_empty:
        return 0.0
    z_min, z_max = geom.bounds[1], geom.bounds[3]
    for i in range(samples):
        theta = (i + 0.5) * TWO_PI / samples
        probe = LineString([(theta, z_min - 1.0), (theta, z_max + 1.0)])
        if geom.intersection(probe).is_empty:
            return theta
    return None


def normalize_periodic(
    geom: BaseGeometry,
    z_min: float,
    z_max: float,
    *,
    cut: float | None = None,
) -> tuple[BaseGeometry, tuple[float, float], bool]:
    if geom.is_empty:
        return GeometryCollection(), (0.0, TWO_PI), False
    base = periodic_copies(geom, turns=1)
    if cut is None:
        cut = choose_empty_cut(base)
    full = cut is None
    lower = 0.0 if full else float(cut)
    upper = lower + TWO_PI
    clipped = _safe_intersection(base, box(lower, z_min, upper, z_max))
    return clipped, (lower, upper), full


def fold_to_canonical(geom: BaseGeometry, z_min: float, z_max: float) -> BaseGeometry:
    if geom.is_empty:
        return geom
    pieces = []
    min_k = math.floor((geom.bounds[0] - TWO_PI) / TWO_PI)
    max_k = math.ceil((geom.bounds[2] + TWO_PI) / TWO_PI)
    clip = box(0.0, z_min, TWO_PI, z_max)
    for k in range(min_k, max_k + 1):
        part = _safe_intersection(translate_theta(geom, -k * TWO_PI), clip)
        if not part.is_empty:
            pieces.append(part)
    return _safe_union(pieces) if pieces else GeometryCollection()
