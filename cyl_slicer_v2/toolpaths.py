from __future__ import annotations

import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

import shapely
from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPolygon, Point, Polygon, box
from shapely.ops import unary_union
from shapely.strtree import STRtree

from .config import SlicerConfig
from .gap_fill import narrow_gap_fill_paths
from .layers import merge_support_body, segment_support_wedges
from .coordinates import TWO_PI, coalesce_collinear_points, cylindrical_distance
from .models import LayerSlice, LayerToolpaths, Toolpath
from .periodic import periodic_copies, translate_theta

ROLE_ORDER = [
    "EXTERNAL_PERIMETER",
    "INTERNAL_PERIMETER",
    "GAP_FILL",
    "BRIDGE",
    "TOP_SKIN",
    "SOLID_INFILL",
    "SUPPORT",
    "SUPPORT_INTERFACE",
    "SUPPORT_BOTTOM_INTERFACE",
    "INFILL_RADIAL",
    "INFILL_AXIAL",
]

_WORKER_CONFIG: SlicerConfig | None = None
PERIMETER_ROLES = frozenset({"EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"})
MIN_GEOMETRY_LENGTH = 1e-8
MIN_OUTPUT_LENGTH_GUARD = 1e-3


def _min_printable_length(config: SlicerConfig) -> float:
    if config.min_path_length <= 0:
        return 0.0
    return config.min_path_length + MIN_OUTPUT_LENGTH_GUARD


def _min_path_length_for_role(role: str, config: SlicerConfig) -> float:
    if role in PERIMETER_ROLES:
        return MIN_GEOMETRY_LENGTH
    return _min_printable_length(config)


def _to_metric(geom, radius: float):
    return shapely.transform(geom, lambda p: p * [radius, 1.0])


def _from_metric(geom, radius: float):
    return shapely.transform(geom, lambda p: p * [1.0 / radius, 1.0])


def _metric_buffer(geom, radius: float, distance: float):
    """Buffer by physical mm in unwrapped theta*R / axial space."""
    return _clean_polygonal(_from_metric(_to_metric(geom, radius).buffer(distance), radius))


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


def _island_path_groups(paths_by_role: dict[str, list[Toolpath]], region) -> list[dict[str, list[Toolpath]]]:
    """Assign each path once to a connected component, joining seam neighbors."""
    parts = _polygon_parts(region)
    if not paths_by_role:
        return []
    if len(parts) <= 1:
        return [paths_by_role]

    tree = STRtree(parts)
    parents = list(range(len(parts)))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    for i, part in enumerate(parts):
        for shift in (0.0, -TWO_PI, TWO_PI):
            shifted = translate_theta(part, shift) if shift else part
            for candidate in tree.query(shifted.envelope.buffer(1e-8)):
                j = int(candidate)
                if root(i) != root(j) and shifted.distance(parts[j]) <= 1e-8:
                    parents[root(j)] = root(i)

    groups: dict[int, dict[str, list[Toolpath]]] = {}
    unassigned: dict[str, list[Toolpath]] = defaultdict(list)
    center = (region.bounds[0] + region.bounds[2]) * 0.5
    for role, paths in paths_by_role.items():
        for path in paths:
            theta, axial = path.points[0]
            theta += round((center - theta) / TWO_PI) * TWO_PI
            owner = None
            for shift in (0.0, -TWO_PI, TWO_PI):
                point = Point(theta + shift, axial)
                for candidate in tree.query(point.buffer(1e-8)):
                    index = int(candidate)
                    if parts[index].distance(point) <= 1e-8:
                        owner = root(index)
                        break
                if owner is not None:
                    break
            if owner is None:
                # Preserve unexpected paths instead of silently dropping material.
                unassigned[role].append(path)
            else:
                groups.setdefault(owner, {}).setdefault(role, []).append(path)
    return [*groups.values(), *([dict(unassigned)] if unassigned else [])]


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


def _raft_footprint_with_axial_gap(footprint, config: SlicerConfig, *, radius: float, theta_window=None):
    _, z_min, _, z_max = footprint.bounds
    # Bridge axial gaps for adhesion. Close only bounded angular gaps: a
    # global bounding box would fill intentional cutouts in a split ring.
    raft = _clean_polygonal(unary_union([
        box(poly.bounds[0], z_min, poly.bounds[2], z_max)
        for poly in _polygon_parts(footprint)
    ]))
    if config.raft_merge_gap_mm > 0:
        angular_gap = config.raft_merge_gap_mm / radius
        intervals = sorted((poly.bounds[0], poly.bounds[2]) for poly in _polygon_parts(raft))
        bridges = [
            box(left_hi, z_min, right_lo, z_max)
            for (_, left_hi), (right_lo, _) in zip(intervals, intervals[1:])
            if 0 < right_lo - left_hi <= angular_gap + 1e-10
        ]
        # The cut is a coordinate boundary, not an island boundary. Fill a
        # small gap straddling it with two pieces in the same theta window.
        if _is_full_theta_window(theta_window):
            lower, upper = theta_window
            first_lo, last_hi = intervals[0][0], intervals[-1][1]
            seam_gap = first_lo + TWO_PI - last_hi
            if (lower <= first_lo <= last_hi <= upper
                    and 0 < seam_gap <= angular_gap + 1e-10):
                if last_hi < upper:
                    bridges.append(box(last_hi, z_min, upper, z_max))
                if first_lo > lower:
                    bridges.append(box(lower, z_min, first_lo, z_max))
        raft = _clean_polygonal(unary_union([raft, *bridges]))
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

    cut = box(raft.bounds[0], cut_lo, raft.bounds[2], cut_hi)
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
    pts = coalesce_collinear_points([(float(t), float(z)) for t, z in ring.coords], radius)
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


def _make_path(role: str, radius: float, pts: list[tuple[float, float]], closed: bool, config: SlicerConfig,
               *, preserve_short: bool = False) -> Toolpath | None:
    if len(pts) < 2:
        return None
    if role in PERIMETER_ROLES:
        cleaned = []
        for pt in pts:
            if not cleaned or not _points_close(cleaned[-1], pt):
                cleaned.append(pt)
        if len(cleaned) < 2:
            return None
        if not any(cylindrical_distance(a, b, radius) > MIN_GEOMETRY_LENGTH for a, b in zip(cleaned, cleaned[1:])):
            return None
        return Toolpath(role=role, radius=radius, points=cleaned, closed=closed)

    # Only simplify collinear geometry within a physical deviation tolerance.
    # Distance-based *vertex* pruning can delete a real corner and silently
    # replace two extruded legs with a diagonal. The minimum printable length
    # applies to the complete feature, not each constituent segment.
    cleaned = coalesce_collinear_points(pts, radius)
    if len(cleaned) < 2:
        return None
    total_length = sum(cylindrical_distance(a, b, radius) for a, b in zip(cleaned, cleaned[1:]))
    if total_length < max(MIN_GEOMETRY_LENGTH, 0.0 if preserve_short else _min_printable_length(config)):
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
    # Whole-turn translation can leave sub-ulp cracks at a nonzero cut.
    # Weld only numerical noise before insetting; otherwise those cracks
    # become real half-bead gaps. This grid is far below physical precision.
    source = shapely.set_precision(periodic_copies(geom, turns=1), 1e-12) if periodic else geom
    metric = _to_metric(source, radius)
    inset = metric.buffer(-distance, join_style="mitre")
    if inset.is_empty:
        return GeometryCollection()
    raw = _from_metric(shapely.make_valid(inset), radius)
    if periodic:
        lower, upper = theta_window if theta_window is not None else (geom.bounds[0], geom.bounds[2])
        raw = raw.intersection(box(lower, geom.bounds[1], upper, geom.bounds[3]))
    return raw


def _registered_axial_grid(radius: float, spacing: float, config: SlicerConfig) -> tuple[float, float, int, int]:
    """Periodic angular pitch, phase, ring count and nested child stride.

    The pitch and phase are fixed across radial slabs. In nested mode a
    power-of-two subset of that same lattice is active on smaller radii;
    pre-existing tracks never move when children become eligible.
    """
    outer = config.axial_grid_outer_radius_mm
    if outer is None:
        raise ValueError("Registered axial grid requires an outer radius; "
                         "use build_toolpaths() or set axial_grid_outer_radius_mm")
    if radius > outer + 1e-7:
        raise ValueError("axial_grid_outer_radius_mm must include all printed radii")
    count = max(1, math.ceil(TWO_PI * outer / spacing - 1e-10))
    step = TWO_PI / count
    stride = 1
    if config.axial_grid_mode == "nested":
        # Adjacent live tracks remain at most spacing mm apart on this
        # layer. On larger radii deterministic midpoints are introduced.
        while stride * 2 < count and step * radius * (stride * 2) <= spacing + 1e-10:
            stride *= 2
    return step, math.radians(config.axial_grid_phase_deg), count, stride


def line_grid_lines(geom, radius: float, spacing: float, *, direction: str,
                    config: SlicerConfig | None = None) -> list[LineString]:
    if geom.is_empty:
        return []
    theta_min, z_min, theta_max, z_max = geom.bounds
    lines = []
    if direction == "axial":
        if config is not None and config.axial_grid_mode != "per-layer":
            step, phase, count, stride = _registered_axial_grid(radius, spacing, config)
            # Integer lattice indices avoid cumulative floating-point
            # phase drift or global realignment when island bounds change.
            first = math.ceil((theta_min - phase) / step - 1e-10)
            last = math.floor((theta_max - phase) / step + 1e-10)
            lower_is_grid = abs((theta_min-phase)/step - round((theta_min-phase)/step)) < 1e-9
            for k in range(first, last + 1):
                # The parent/child subset must repeat every *N* ring
                # tracks. N need not be divisible by stride, so k % stride
                # alone would shift the nested grid by one at the 0/2pi seam.
                if (k % count) % stride != 0:
                    continue
                theta = phase + k * step
                # A complete 2pi closed strip has one physical seam line,
                # not two extrusion passes at its two cut coordinates.
                if (theta_max-theta_min >= TWO_PI - 1e-9
                        and lower_is_grid and abs(theta-theta_max) < 1e-9):
                    continue
                lines.append(LineString([(theta, z_min-spacing), (theta, z_max+spacing)]))
        else:
            # Exact historical behavior for quality and legacy profiles.
            step = spacing / radius
            start = math.floor(theta_min / step) * step
            theta = start
            while theta <= theta_max + step * 0.5:
                lines.append(LineString([(theta, z_min-spacing), (theta, z_max+spacing)]))
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


def line_grid_paths(geom, role: str, radius: float, spacing: float, config: SlicerConfig,
                    *, direction: str, theta_window=None) -> list[Toolpath]:
    result = []
    for line in line_grid_lines(geom, radius, spacing, direction=direction, config=config):
        pts = [(float(t), float(z)) for t, z in line.coords]
        # A tiny clipped end may belong to a long printable pass across the
        # seam. Apply the whole-feature length limit after periodic joining.
        seam_end = _is_full_theta_window(theta_window) and any(
            abs(point[0] - edge) * radius <= 1e-7
            for point in (pts[0], pts[-1]) for edge in theta_window
        )
        path = _make_path(role, radius, pts, False, config, preserve_short=seam_end)
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
    return path.role not in PERIMETER_ROLES


def _join_periodic_seam_paths(paths: list[Toolpath], theta_window) -> list[Toolpath]:
    """Stitch unique collinear continuations across the cut, in unwrapped space.

    No connector is deposited. Only coincident endpoints on opposite window
    edges qualify. Orient fragments into a continuous walk, retaining signed
    theta deltas and reversing bead widths with their segments when needed.
    Existing closed contours and individual full-turn paths stay untouched.
    """
    if not _is_full_theta_window(theta_window):
        return paths
    tolerance = 1e-7  # Physical mm; absorb intersection roundoff only.
    buckets = defaultdict(list)
    endpoints = {}
    for index, path in enumerate(paths):
        if path.closed or len(path.points) < 2:
            continue
        for end in (0, 1):
            point = path.points[0 if end == 0 else -1]
            for side, theta in enumerate(theta_window):
                if abs(point[0] - theta) * path.radius <= tolerance:
                    key = (index, end)
                    endpoints[key] = (side, point[1])
                    buckets[(side, math.floor(point[1] / tolerance))].append(key)
                    break

    contacts = {}
    for key, (side, axial) in endpoints.items():
        path = paths[key[0]]
        bucket = math.floor(axial / tolerance)
        matches = []
        for offset in (-1, 0, 1):
            for other in buckets.get((1 - side, bucket + offset), ()):
                candidate = paths[other[0]]
                if (key[0] != other[0] and path.role == candidate.role
                        and abs(path.radius - candidate.radius) <= 1e-9
                        and abs(axial - endpoints[other][1]) <= tolerance):
                    matches.append(other)
        if len(matches) == 1:
            contacts[key] = matches[0]

    links = {}
    for key, other in contacts.items():
        if contacts.get(other) != key:
            continue  # Ambiguous/branched contacts cannot form a chain.
        path, candidate = paths[key[0]], paths[other[0]]
        if (path.segment_widths is None) != (candidate.segment_widths is None):
            continue
        vectors = []
        for index, end in (key, other):
            points = paths[index].points
            tip, inner = (points[0], points[1]) if end == 0 else (points[-1], points[-2])
            vectors.append(((inner[0] - tip[0]) * path.radius, inner[1] - tip[1]))
        u, v = vectors
        norm = math.hypot(*u) * math.hypot(*v)
        if norm <= 1e-20 or u[0] * v[0] + u[1] * v[1] >= 0:
            continue
        if abs(u[0] * v[1] - u[1] * v[0]) > 1e-5 * norm:
            continue
        links[key] = other

    result = []
    visited = set()
    # Start at free ends before considering any completely periodic cycles.
    starts = [(i, end) for i in range(len(paths)) for end in (0, 1)
              if (i, end) not in links]
    starts.extend((i, 0) for i in range(len(paths)))
    for index, entry in starts:
        if index in visited:
            continue
        original = paths[index]
        points = []
        widths = [] if original.segment_widths is not None else None
        while index not in visited:
            visited.add(index)
            path = paths[index]
            fragment = list(path.points if entry == 0 else reversed(path.points))
            fragment_widths = path.segment_widths
            if fragment_widths is not None and entry == 1:
                fragment_widths = list(reversed(fragment_widths))
            if points:
                shift = round((points[-1][0] - fragment[0][0]) / TWO_PI) * TWO_PI
                fragment = [(theta + shift, axial) for theta, axial in fragment]
                points.extend(fragment[1:])
            else:
                points.extend(fragment)
            if widths is not None:
                widths.extend(fragment_widths)
            next_end = links.get((index, 1 - entry))
            if next_end is None:
                break
            index, entry = next_end
        closed = original.closed or _points_close(points[0], points[-1])
        if closed:
            points[-1] = points[0]
        result.append(replace(original, points=points, segment_widths=widths, closed=closed))
    return result


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
    kept = [path for path in paths if _path_length(path) >= _min_path_length_for_role(path.role, config)]
    if any(path.segment_widths is not None for path in kept):
        # Never join or prune variable-E paths without updating width metadata.
        return kept
    joined = _join_touching_paths(kept)
    if paths and paths[0].role in PERIMETER_ROLES:
        return joined
    return _remove_contained_straight_paths(joined)


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
    if layer.deposited_thickness is not None:
        config = config.for_deposited_thickness(layer.deposited_thickness)
    paths: dict[str, list[Toolpath]] = defaultdict(list)
    model = layer.occupied
    radius = layer.radius

    if not model.is_empty:
        theta_window = layer.theta_window
        # Slanted ribbons can cross the cut without occupying a full angular
        # band at any axial height. Their window still identifies periodic
        # neighbors; inset those together so the cut cannot open a bead gap.
        periodic = layer.is_full_wrap or (
            _is_full_theta_window(theta_window)
            and model.bounds[0] <= theta_window[0] + 1e-8
            and model.bounds[2] >= theta_window[1] - 1e-8
        )
        # Contours are occupied MATERIAL boundaries, not bead CENTERLINES.
        # First wall centers sit half a bead inside outer edges and holes.
        # Remaining walls use the same spacing as dense lines; for a rounded
        # bead this is (area / layer height), not the nominal bead width.
        for wall in range(config.wall_count):
            inset = inset_geometry(
                model,
                radius,
                config.nozzle_width * 0.5 + config.dense_line_spacing * wall,
                periodic=periodic,
                theta_window=theta_window,
            )
            if inset.is_empty:
                break
            role = "EXTERNAL_PERIMETER" if wall == 0 else "INTERNAL_PERIMETER"
            paths[role].extend(perimeter_paths(inset, role, radius, config, theta_window))

        # The last wall covers through w + (wall_count - 1) * dense_pitch.
        # Clip infill to this boundary with independently chosen anchoring
        # overlap; this keeps the rounded cross-section/pitch model coherent.
        # The same rule for solid and sparse avoids retracing wall centerlines
        # with skins, or leaving a sparse-infill anchoring gap.
        overlap = config.nozzle_width * config.infill_wall_overlap_fraction
        infill_source = inset_geometry(
            model,
            radius,
            config.wall_coverage_depth(config.wall_count) - overlap,
            periodic=periodic,
            theta_window=theta_window,
        )
        # Reserve bridge corridors and their bearing overlap before generating
        # ordinary fill. This prevents two intentional extrusions in the same
        # physical area when bridges intersect bottom skins or sparse grids.
        bridge_masks = {}
        if layer.bridge_regions is not None and config.bridge_mode == "print":
            for direction, region in layer.bridge_regions.items():
                if region is None or region.is_empty:
                    continue
                extended = (
                    _metric_buffer(region, radius, config.bridge_anchor_overlap_mm)
                    if config.bridge_anchor_overlap_mm > 0 else region
                )
                bridge_masks[direction] = _clean_polygonal(
                    extended.intersection(infill_source)
                )
        bridge_union = (
            unary_union([part for part in bridge_masks.values() if not part.is_empty])
            if bridge_masks else GeometryCollection()
        )
        for direction, region in bridge_masks.items():
            if not region.is_empty:
                paths["BRIDGE"].extend(line_grid_paths(
                    region, "BRIDGE", radius, config.dense_line_spacing,
                    config, direction=direction, theta_window=theta_window,
                ))
        if not bridge_union.is_empty:
            infill_source = _safe_polygon_difference(infill_source, bridge_union)
        solid_source = unary_union([layer.top_surfaces, layer.bottom_surfaces])
        solid_source = solid_source.intersection(infill_source) if not solid_source.is_empty else solid_source
        solid_direction = (
            "radial" if config.solid_infill_pattern == "alternate" and layer.index % 2 else "axial"
        )
        if config.monotonic_top_skin:
            # Do not fuse exposed skin with the buried bottom/solid regions.
            # Top dominates where both surface masks overlap (thin objects),
            # so the same physical area is not deposited twice.
            top = layer.top_surfaces.intersection(infill_source)
            bottom = layer.bottom_surfaces.intersection(infill_source)
            if not top.is_empty and not bottom.is_empty:
                bottom = bottom.difference(top)
            if not top.is_empty:
                paths["TOP_SKIN"].extend(
                    line_grid_paths(top, "TOP_SKIN", radius, config.dense_line_spacing,
                                    config, direction=solid_direction, theta_window=theta_window)
                )
            if not bottom.is_empty:
                paths["SOLID_INFILL"].extend(
                    line_grid_paths(bottom, "SOLID_INFILL", radius, config.dense_line_spacing,
                                    config, direction=solid_direction, theta_window=theta_window)
                )
        else:
            paths["SOLID_INFILL"].extend(
                line_grid_paths(solid_source, "SOLID_INFILL", radius, config.dense_line_spacing,
                                config, direction=solid_direction, theta_window=theta_window)
            )
        if not solid_source.is_empty and not infill_source.is_empty:
            infill_source = infill_source.difference(solid_source)
            if not infill_source.is_valid:
                infill_source = shapely.make_valid(infill_source)
            infill_source = infill_source.buffer(0)
        # Crosshatching two orthogonal families at spacing S extrudes about
        # twice the line length of one family at S. Halve the one-family
        # spacing when alternating so the deposited volume per *layer*
        # remains comparable, apart from finite polygon clipping effects.
        sparse_spacing = (
            config.infill_spacing * 0.5
            if config.sparse_infill_pattern == "alternate"
            else config.infill_spacing
        )
        if config.sparse_infill_pattern == "both" or layer.index % 2 == 0:
            paths["INFILL_RADIAL"].extend(
                line_grid_paths(infill_source, "INFILL_RADIAL", radius, sparse_spacing, config,
                                direction="radial", theta_window=theta_window)
            )
        if config.sparse_infill_pattern == "both" or layer.index % 2 != 0:
            paths["INFILL_AXIAL"].extend(
                line_grid_paths(infill_source, "INFILL_AXIAL", radius, sparse_spacing, config,
                                direction="axial", theta_window=theta_window)
            )

    if config.gap_fill_enabled and not model.is_empty:
        walls = [path for role in PERIMETER_ROLES for path in paths.get(role, [])]
        paths["GAP_FILL"].extend(narrow_gap_fill_paths(model, walls, radius, config))

    # Also handle externally supplied LayerSlices, and expose the real slots
    # to island grouping and travel containment. Never let cross-ties, dense
    # interfaces or a raft bridge them with the finite deposited bead width.
    support_footprint = segment_support_wedges(layer.support, radius, config)
    support = _clean_polygonal(segment_support_wedges(
        support_footprint, radius, config, bead_clearance=True,
    ))
    interface_source = layer.support_interface if layer.support_interface is not None else GeometryCollection()
    interface = _clean_polygonal(interface_source)
    if not interface.is_empty:
        interface = interface.intersection(support) if not support.is_empty else GeometryCollection()
        interface = _clean_polygonal(interface)
        if not interface.is_empty:
            interface_spacing = config.support_interface_spacing or config.dense_line_spacing
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
                        theta_window=layer.theta_window,
                    )
                )
            support = _safe_polygon_difference(support, interface)

    # Top and bottom support contact masks can coexist on the same layer.
    # Give the upper/overhang-facing interface precedence where they overlap;
    # the two dense roles must never deposit twice in the same footprint.
    bottom_source = (layer.support_bottom_interface if
                     layer.support_bottom_interface is not None else GeometryCollection())
    bottom = _clean_polygonal(bottom_source)
    if not bottom.is_empty:
        bottom = _clean_polygonal(bottom.intersection(support))
        if not bottom.is_empty:
            bottom_spacing = (config.support_bottom_interface_spacing
                              or config.support_interface_spacing
                              or config.dense_line_spacing)
            bottom_direction = "axial" if layer.index % 2 == 0 else "radial"
            if bottom_direction == "radial":
                paths["SUPPORT_BOTTOM_INTERFACE"].extend(
                    meandered_radial_grid_paths(
                        bottom, "SUPPORT_BOTTOM_INTERFACE", radius,
                        bottom_spacing, config,
                    )
                )
            else:
                paths["SUPPORT_BOTTOM_INTERFACE"].extend(
                    line_grid_paths(bottom, "SUPPORT_BOTTOM_INTERFACE", radius,
                                    bottom_spacing, config, direction="axial",
                                    theta_window=layer.theta_window)
                )
            support = _safe_polygon_difference(support, bottom)

    if not support.is_empty:
        support_spacing = (config.raft_grid_spacing or config.dense_line_spacing) if layer.is_raft else config.support_spacing
        if layer.is_raft:
            paths["SUPPORT"].extend(
                line_grid_paths(
                    support,
                    "SUPPORT",
                    radius,
                    support_spacing,
                    config,
                    direction=layer.raft_direction or "radial",
                    theta_window=layer.theta_window,
                )
            )
        else:
            radial_support = meandered_radial_grid_paths(support, "SUPPORT", radius, support_spacing, config)
            paths["SUPPORT"].extend(radial_support)
            if radial_support:
                cross_spacing = config.support_cross_spacing or support_spacing * 4.0
                paths["SUPPORT"].extend(
                    line_grid_paths(support, "SUPPORT", radius, cross_spacing, config,
                                    direction="axial", theta_window=layer.theta_window)
                )

    paths = {
        role: processed
        for role, role_paths in paths.items()
        if (processed := _postprocess_role_paths(
            _join_periodic_seam_paths(role_paths, layer.theta_window), config,
        ))
    }
    travel_region = (
        model if support_footprint.is_empty else support_footprint if model.is_empty
        else unary_union([model, support_footprint])
    )
    support_roles = {"SUPPORT", "SUPPORT_INTERFACE", "SUPPORT_BOTTOM_INTERFACE"}
    model_paths = {role: items for role, items in paths.items() if role not in support_roles}
    support_paths = {role: items for role, items in paths.items() if role in support_roles}
    islands = _island_path_groups(model_paths, model) + _island_path_groups(support_paths, support_footprint)
    # Slicing/intersection geometry remains at the slab midpoint. Nozzle Z
    # must be at its OUTER deposited surface, not at the sampled midpoint.
    return LayerToolpaths(layer_index=layer.index, radius=radius, paths_by_role=paths,
                          travel_region=travel_region, island_paths=islands,
                          nozzle_radius=radius + config.layer_height * 0.5,
                          deposited_thickness=config.layer_height,
                          previous_occupied=layer.previous_occupied)


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

    # The off-bed foundation is support just as much as a column below an
    # overhang.  In particular, it must not silently bypass the configured
    # contact clearance merely because the first model layer has no preceding
    # model slice.  Compare deposited faces here, as propagate_support() does
    # for ordinary supports.  The legacy layer count is expressed in nominal
    # radial slabs; an explicit physical clearance takes precedence.
    model_height = (first_model.deposited_thickness
                    if first_model.deposited_thickness is not None
                    else config.layer_height)
    model_inner = first_model.radius - model_height * 0.5
    requested_gap = (config.support_contact_gap_mm
                     if config.support_contact_gap_mm is not None
                     else config.effective_support_gap_layers * config.layer_height)
    support_outer_limit = model_inner - requested_gap

    if config.radial_layer_mode == "complete-terminal":
        # The model is sampled at its slab midpoint, so foundation slabs must
        # stop at the physical contact limit derived from its INNER DEPOSITED
        # FACE. An ordinary fixed-height loop based on sample centers can
        # overlap that limit on partially filled terminal slabs.
        raw_span = model_inner - config.bed_radius
        if raw_span <= 1e-9:
            return []
        if raw_span < (config.first_layer_height or config.layer_height * 0.5) - 1e-9:
            raise ValueError(
                "Gap-support foundation is thinner than half the requested "
                "radial layer height; decrease --layer-height or adjust model radius"
            )
        span = support_outer_limit - config.bed_radius
        # A deliberately requested contact gap can consume all, or all but a
        # sub-half-height sliver, of the available foundation.  Do not turn
        # that request into an unsafe thin support bead.
        if span < (config.first_layer_height or config.layer_height * 0.5) - 1e-9:
            return []
        h = config.layer_height
        full_count = int(math.floor(span / h + 1e-10))
        remainder = span - full_count * h
        if abs(remainder) < 1e-9:
            remainder = 0.0
        heights = [h] * full_count
        if remainder > 0:
            if not heights:
                heights = [remainder]
            elif remainder < h * 0.5:
                # Rebalance the upper support layers rather than extruding
                # a sliver or overlapping the first model slab.
                shared = (heights.pop() + remainder) * 0.5
                heights.extend((shared, shared))
            else:
                heights.append(remainder)
        support_slabs = []
        edge = config.bed_radius
        for height in heights:
            support_slabs.append((float(edge + 0.5 * height), float(height)))
            edge += height
    else:
        # Fixed-height foundations keep their bed-anchored schedule, but only
        # emit slabs whose outer deposited face is below the required contact
        # limit.  This can leave a little more than the requested clearance
        # when the model and bed schedules are not phase-aligned, never less.
        start = config.bed_radius + config.layer_height * 0.5
        if start + config.layer_height * 0.5 > support_outer_limit + 1e-9 and config.first_layer_height is None:
            return []
        support_slabs = []
        radius = start
        while radius + config.layer_height * 0.5 <= support_outer_limit + 1e-9:
            support_slabs.append((float(radius), None))
            radius += config.layer_height

    if config.first_layer_height is not None:
        # Foundation begins at the bed, independently of the model schedule.
        first = config.first_layer_height
        span = support_outer_limit - config.bed_radius
        if span < first - 1e-9:
            return []
        support_slabs = [(float(config.bed_radius + first * 0.5), first)]
        remaining = max(0.0, span - first)
        h = config.layer_height
        count = int(math.floor((remaining + 1e-10) / h))
        heights = [h] * count
        remainder = remaining - count * h
        if config.radial_layer_mode == "complete-terminal" and remainder > 1e-9:
            if remainder >= h * 0.5 - 1e-9:
                heights.append(remainder)
            elif heights:
                shared = (heights.pop() + remainder) * 0.5
                heights.extend((shared, shared))
        edge = config.bed_radius + first
        for height in heights:
            support_slabs.append((float(edge + height * 0.5), height))
            edge += height

    pieces = []
    if not first_model.support.is_empty:
        pieces.append(first_model.support)
    if not first_model.occupied.is_empty:
        pieces.append(first_model.occupied)
    footprint = _clean_polygonal(unary_union(pieces)) if pieces else GeometryCollection()
    if footprint.is_empty:
        return []
    raft_footprint = _raft_footprint_with_axial_gap(
        footprint, config, radius=support_slabs[0][0], theta_window=first_model.theta_window,
    )

    count = len(support_slabs)
    empty = GeometryCollection()
    return [
        LayerSlice(
            index=i - count,
            radius=radius,
            occupied=empty,
            top_surfaces=empty,
            bottom_surfaces=empty,
            needs_support=empty,
            support=segment_support_wedges(
                raft_footprint if i < config.raft_layers
                else merge_support_body(footprint, radius, config), radius, config,
            ),
            theta_window=first_model.theta_window,
            is_full_wrap=first_model.is_full_wrap,
            is_raft=i < config.raft_layers,
            raft_direction=("radial" if i % 2 == 0 else "axial") if i < config.raft_layers else None,
            deposited_thickness=thickness,
        )
        for i, (radius, thickness) in enumerate(support_slabs)
    ]


def build_toolpaths(layers: list[LayerSlice], config: SlicerConfig) -> list[LayerToolpaths]:
    all_layers = [*gap_support_slices(layers, config), *layers]
    if config.axial_grid_mode != "per-layer" and config.axial_grid_outer_radius_mm is None and all_layers:
        # One radius shared by all workers, model islands, support, and raft.
        # Use the physical deposited outer face, not a local model/island
        # bounding box or the sample radius of a truncated worker batch.
        outer = max(layer.radius + 0.5 * (
            layer.deposited_thickness if layer.deposited_thickness is not None
            else config.layer_height) for layer in all_layers)
        config = replace(config, axial_grid_outer_radius_mm=outer)
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
