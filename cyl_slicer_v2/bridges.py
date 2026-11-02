"""Conservative, physical-mm two-ended bridge classification on a cylinder.

Bridge candidates are checked in narrow transverse bands, not assumed to be
printable across an entire irregular unsupported polygon. Unverified regions
retain ordinary support. Actual mandrel gravity/phase still needs coupons.
"""
from __future__ import annotations

import math

from shapely.geometry import Point, LineString, GeometryCollection, Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely import affinity

from .periodic import periodic_copies


def _polygons(g):
    if isinstance(g, Polygon):
        yield g
    elif isinstance(g, MultiPolygon) or hasattr(g, "geoms"):
        for p in g.geoms:
            yield from _polygons(p)


def _lines(geom):
    """All nonempty pieces of a physical-mm scanline intersection."""
    if isinstance(geom, LineString):
        yield geom
    elif hasattr(geom, "geoms"):
        for part in geom.geoms:
            yield from _lines(part)


def _anchored_bridge_lanes(candidate, previous, axis, width, max_span):
    """Conservatively classify *transverse lanes*, not a whole overhang.

    Every printed region must have two independent previous-slab bearings.
    Probe three transverse locations in each narrow band; all must have an
    uninterrupted unsupported chord of printable length and supported ends.
    A band with a missing/one-sided bearing is not converted to bridge.
    Trim the edges of accepted runs by half a bead when adjacent bands fail,
    keeping bridge beads away from the unverified part of an overhang.
    """
    if candidate.interiors or candidate.is_empty or not candidate.is_valid:
        return GeometryCollection()
    xmin, ymin, xmax, ymax = candidate.bounds
    low, high = (ymin, ymax) if axis == "radial" else (xmin, xmax)
    sideways = high - low
    if sideways <= width or width <= 0 or max_span <= 0:
        return GeometryCollection()
    count = max(1, math.ceil(sideways / (0.5 * width)))
    # A fine-textured candidate cannot cause unbounded geometry operations.
    if count > 1024:
        return GeometryCollection()
    pitch = sideways / count
    offset = min(width * 0.25, 0.10)
    approved = []
    for k in range(count):
        for fraction in (0.15, 0.5, 0.85):
            transverse = low + (k + fraction) * pitch
            if axis == "radial":
                probe = LineString([(xmin-width, transverse),
                                    (xmax+width, transverse)])
            else:
                probe = LineString([(transverse, ymin-width),
                                    (transverse, ymax+width)])
            segments = [s for s in _lines(candidate.intersection(probe))
                        if s.length > 1e-8]
            # Reject forks, holes and locally disconnected spans.
            if len(segments) != 1:
                break
            ends = [p[0] if axis == "radial" else p[1]
                    for p in segments[0].coords]
            begin, end = min(ends), max(ends)
            if not width * 1.5 <= end - begin <= max_span:
                break
            if axis == "radial":
                start_bearing = Point(begin-offset, transverse)
                end_bearing = Point(end+offset, transverse)
            else:
                start_bearing = Point(transverse, begin-offset)
                end_bearing = Point(transverse, end+offset)
            if not previous.covers(start_bearing) or not previous.covers(end_bearing):
                break
        else:
            approved.append(k)

    if not approved:
        return GeometryCollection()
    if len(approved) == count:
        # Preserve the previous exact-region behavior for a wholly verified
        # rectangular bridge (no extra edges, quantization or density drift).
        return candidate

    regions = []
    runs = []
    start = last = approved[0]
    for index in approved[1:]:
        if index != last + 1:
            runs.append((start, last))
            start = index
        last = index
    runs.append((start, last))
    for first, last in runs:
        # A .5*w safety margin on either side excludes bead centers near
        # failed bands, while maintaining the original physical endpoints.
        lower = low + first * pitch + width * 0.5
        upper = low + (last + 1) * pitch - width * 0.5
        if upper <= lower:
            continue
        window = (Polygon([(xmin-width, lower), (xmax+width, lower),
                           (xmax+width, upper), (xmin-width, upper)])
                  if axis == "radial"
                  else Polygon([(lower, ymin-width), (upper, ymin-width),
                                (upper, ymax+width), (lower, ymax+width)]))
        clipped = candidate.intersection(window)
        if not clipped.is_empty:
            regions.append(clipped)
    return unary_union(regions) if regions else GeometryCollection()

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
        possibilities = []
        for axis in ("radial", "axial"):
            verified = _anchored_bridge_lanes(
                polygon, prior_cover, axis, bead_width, max_span_mm,
            )
            if not verified.is_empty:
                span = ((polygon.bounds[2]-polygon.bounds[0]) if axis=="radial"
                        else (polygon.bounds[3]-polygon.bounds[1]))
                possibilities.append((axis, verified, span))
        if not possibilities:
            continue
        # Prefer the shorter physical span if both axes can bridge.
        # Favor more verified area only to break equal-span ties. Never
        # replace a short supported bridge with a longer cross-axis span
        # simply because a separate part of the polygon also has bearings.
        direction, footprint, _ = min(
            possibilities, key=lambda p: (p[2], -p[1].area, p[0]),
        )
        directions[direction].append(footprint)

    out = {}
    for direction,parts in directions.items():
        if not parts:
            out[direction] = empty
            continue
        joined = unary_union(parts)
        unwrapped = affinity.scale(joined, xfact=1/radius,yfact=1,origin=(0,0))
        # Intersect the periodic candidates with the EXACT original sampled
        # footprint rather than folding to theta [0,2pi]. The slicer's
        # normalized angular cut can start anywhere in [0,2pi] and extend
        # past 2pi; folding would shift bridges out of that layer's window
        # and create apparent out-of-model regions. This also removes copies
        # without ever triple-printing a seam-adjacent span.
        out[direction] = unwrapped.intersection(current)
    return out
