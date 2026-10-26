"""Measure physical bearing under model walls without changing toolpaths.

Unlike axial/radian fractions, this classifier works in actual unwrapped
cylindrical millimeters and keeps unsupported loop runs continuous at seams.
This is observability only: a rotating mandrel makes gravity/phase important,
so print decisions must not be derived automatically from these numbers.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from shapely import affinity
from shapely.geometry import LineString, Point
from shapely.geometry.base import BaseGeometry

from .coordinates import TWO_PI
from .models import Toolpath
from .periodic import periodic_copies


@dataclass(frozen=True)
class WallSupportStats:
    category: str  # supported | partial | unsupported
    covered_fraction: float
    max_unsupported_run_mm: float
    total_length_mm: float


def _line_parts(geometry):
    if geometry.is_empty:
        return
    if isinstance(geometry, LineString):
        yield geometry
    elif hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            yield from _line_parts(part)


def _longest_uncovered_run(path_line: LineString, uncovered, closed: bool) -> float:
    total = path_line.length
    intervals = []
    for part in _line_parts(uncovered):
        if part.length <= 1e-9:
            continue
        first = path_line.project(Point(part.coords[0]))
        last = path_line.project(Point(part.coords[-1]))
        begin, end = sorted((first, last))
        if end - begin > 1e-9:
            intervals.append((begin, end))
    if not intervals:
        return 0.
    intervals.sort()
    merged = [[*intervals[0]]]
    for begin, end in intervals[1:]:
        if begin <= merged[-1][1] + 1e-7:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([begin, end])
    longest = max(end - begin for begin, end in merged)
    # Shapely splits a closed-loop unsupported run at the arbitrary seam.
    if (closed and len(merged) > 1 and merged[0][0] <= 1e-7
            and merged[-1][1] >= total - 1e-7):
        longest = max(longest, merged[0][1] + total - merged[-1][0])
    return longest


def measure_wall_support(
    path: Toolpath, previous: BaseGeometry, bead_width_mm: float
) -> WallSupportStats:
    """Fraction and longest contiguous unsupported run of a wall centerline.

    Inputs are angle (radians), axial (mm). Previous footprint is repeated
    across the periodic seam and expanded by one half bead width before
    intersecting the actual deposited centerline in R*theta,axial millimeters.
    """
    if (not math.isfinite(path.radius) or path.radius <= 0
            or not math.isfinite(bead_width_mm) or bead_width_mm <= 0):
        raise ValueError("Wall support requires positive finite radius and bead width")
    if len(path.points) < 2:
        return WallSupportStats("supported", 1., 0., 0.)
    phase = TWO_PI * math.floor(path.points[0][0] / TWO_PI)
    physical = [(path.radius * (theta - phase), axial) for theta, axial in path.points]
    wall = LineString(physical)
    total = wall.length
    if total <= 1e-10:
        return WallSupportStats("supported", 1., 0., 0.)
    if previous is None or previous.is_empty:
        return WallSupportStats("unsupported", 0., total, total)
    prior = affinity.scale(periodic_copies(previous), xfact=path.radius,
                           yfact=1., origin=(0., 0.))
    bearing = prior.buffer(bead_width_mm * 0.5)
    uncovered = wall.difference(bearing)
    uncovered_mm = min(total, max(0., uncovered.length))
    fraction = max(0., min(1., 1. - uncovered_mm / total))
    if uncovered_mm <= 1e-8:
        category = "supported"
    elif total - uncovered_mm <= 1e-8:
        category = "unsupported"
    else:
        category = "partial"
    run = _longest_uncovered_run(wall, uncovered, path.closed)
    return WallSupportStats(category, fraction, min(total, max(0., run)), total)
