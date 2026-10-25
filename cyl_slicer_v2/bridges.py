"""Conservative, physical-mm two-ended bridge classification on a cylinder.

Bridge paths are optional because the available mandrel rotates around gravity;
firmware motion and hardware coupons must validate actual bridging performance.
"""
from __future__ import annotations

import math

from shapely.geometry import Point, LineString, GeometryCollection, Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely import affinity

from .periodic import TWO_PI, periodic_copies, fold_to_canonical


def _polygons(g):
    if isinstance(g, Polygon):
        yield g
    elif isinstance(g, MultiPolygon) or hasattr(g, "geoms"):
        for p in g.geoms:
            yield from _polygons(p)


def _span_has_two_anchors(candidate, previous, axis, width, max_span):
    """Validate a narrow, almost rectangular gap with BOTH ends supported.

    Require independent bearing at >half of nine perpendicular samples and
    a straight unsupported corridor. This deliberately rejects rounded,
    branched or one-sided features rather than guessing a bridge direction.
    """
    xmin, ymin, xmax, ymax = candidate.bounds
    span = (xmax - xmin) if axis == "radial" else (ymax-ymin)
    sideways = (ymax-ymin) if axis == "radial" else (xmax-xmin)
    if not width * 1.5 <= span <= max_span or sideways <= width:
        return False
    if candidate.interiors or candidate.area < 0.70 * span * sideways:
        return False
    offset = min(width * 0.25, 0.10)
    supported = 0
    for i in range(9):
        fraction = (i + .5) / 9
        if axis == "radial":
            transverse = ymin + fraction * sideways
            start = (xmin,transverse)
            finish = (xmax,transverse)
            a = (xmin-offset,transverse)
            b = (xmax+offset,transverse)
        else:
            transverse = xmin + fraction * sideways
            start = (transverse,ymin)
            finish = (transverse,ymax)
            a = (transverse,ymin-offset)
            b = (transverse,ymax+offset)
        # The full proposed span must be inside the unsupported candidate
        # and both ends must independently meet previously deposited material.
        if (candidate.buffer(1e-7).covers(LineString([start,finish]))
                and previous.covers(Point(a)) and previous.covers(Point(b))):
            supported += 1
    return supported >= 6


def classify_bridge_regions(
    current,
    previous,
    radius: float,
    bead_width: float,
    max_span_mm: float,
) -> dict[str, object]:
    """Return supported-on-both-ends bridge footprints by print direction.

    The input geometry uses (theta radians, axial mm); all geometric span
    tests use (radius*theta, axial) physical mm, including across 0/2pi.
    The returned regions are folded to canonical theta in [0, 2pi].
    """
    empty = GeometryCollection()
    if (current is None or current.is_empty or previous is None
            or previous.is_empty or radius <= 0 or max_span_mm <= 0):
        return {"axial": empty, "radial": empty}
    source = periodic_copies(current)
    prior = periodic_copies(previous)
    metric_current = affinity.scale(source, xfact=radius, yfact=1, origin=(0,0))
    metric_prior = affinity.scale(prior, xfact=radius, yfact=1, origin=(0,0))
    prior_cover = metric_prior.buffer(bead_width * .5)
    unsup = metric_current.difference(prior_cover)
    directions = {"axial": [], "radial": []}
    for polygon in _polygons(unsup):
        if polygon.is_empty or not polygon.is_valid:
            continue
        possibilities = [
            (axis, (polygon.bounds[2]-polygon.bounds[0]) if axis=="radial"
                   else (polygon.bounds[3]-polygon.bounds[1]))
            for axis in ("radial","axial")
            if _span_has_two_anchors(polygon,prior_cover,axis,bead_width,max_span_mm)
        ]
        if not possibilities:
            continue
        direction = min(possibilities,key=lambda p:p[1])[0]
        directions[direction].append(polygon)

    out = {}
    for direction,parts in directions.items():
        if not parts:
            out[direction] = empty
            continue
        joined = unary_union(parts)
        unwrapped = affinity.scale(joined, xfact=1/radius,yfact=1,origin=(0,0))
        low,high=current.bounds[1],current.bounds[3]
        # Collapse copies only once; never triple-deposit seam-adjacent spans.
        out[direction] = fold_to_canonical(unwrapped,low,high)
    return out
