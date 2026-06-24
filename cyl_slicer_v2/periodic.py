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


def periodic_copies(geom: BaseGeometry, turns: int = 1) -> BaseGeometry:
    if geom.is_empty:
        return geom
    return unary_union([translate_theta(geom, k * TWO_PI) for k in range(-turns, turns + 1)])


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
    clipped = base.intersection(box(lower, z_min, upper, z_max))
    if not clipped.is_valid:
        clipped = shapely.make_valid(clipped)
    return clipped, (lower, upper), full


def fold_to_canonical(geom: BaseGeometry, z_min: float, z_max: float) -> BaseGeometry:
    if geom.is_empty:
        return geom
    pieces = []
    min_k = math.floor((geom.bounds[0] - TWO_PI) / TWO_PI)
    max_k = math.ceil((geom.bounds[2] + TWO_PI) / TWO_PI)
    clip = box(0.0, z_min, TWO_PI, z_max)
    for k in range(min_k, max_k + 1):
        part = translate_theta(geom, -k * TWO_PI).intersection(clip)
        if not part.is_empty:
            pieces.append(part)
    return unary_union(pieces) if pieces else GeometryCollection()
