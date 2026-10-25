from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import itertools
import math
import warnings

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import GeometryCollection, LineString, MultiLineString, Point, Polygon, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import polygonize, unary_union
import trimesh

from .config import SlicerConfig
from .coordinates import TWO_PI
from .mesh_io import axial_bounds
from .periodic import normalize_periodic

_WORKER_MESH: trimesh.Trimesh | None = None
_WORKER_CONFIG: SlicerConfig | None = None
_WORKER_Z_BOUNDS: tuple[float, float] | None = None
_MIXED_FACE_SOLID_THRESHOLD = 0.35
_INTERSECTION_POINT_WELD_TOL = 1e-5
_POINT_WELD_NEIGHBOR_OFFSETS = tuple(itertools.product((-1, 0, 1), repeat=3))


def _cart_to_theta_z(points: np.ndarray) -> np.ndarray:
    theta = np.mod(np.arctan2(points[:, 1], points[:, 0]), 2.0 * np.pi)
    return np.column_stack([theta, points[:, 2]])


def _restore_loop_walk(coords: np.ndarray) -> np.ndarray:
    if len(coords) <= 1:
        return coords
    out = coords.copy()
    out[:, 0] = np.unwrap(out[:, 0])
    if abs(float(out[-1, 0] - out[0, 0])) < 1e-8 and abs(float(out[-1, 1] - out[0, 1])) < 1e-8:
        out[-1] = out[0]
    return out


def _weld_intersection_graph(
    points: np.ndarray,
    edges: list[tuple[int, int]],
    tol: float = _INTERSECTION_POINT_WELD_TOL,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    if len(points) == 0 or not edges:
        return points, edges

    inv_tol = 1.0 / tol
    tol_sq = tol * tol
    buckets: dict[tuple[int, int, int], list[int]] = {}
    old_to_new = np.empty(len(points), dtype=np.int64)
    welded: list[np.ndarray] = []

    for old_index, point in enumerate(points):
        key = (
            math.floor(float(point[0]) * inv_tol),
            math.floor(float(point[1]) * inv_tol),
            math.floor(float(point[2]) * inv_tol),
        )
        match = None
        for offset in _POINT_WELD_NEIGHBOR_OFFSETS:
            neighbor_key = (key[0] + offset[0], key[1] + offset[1], key[2] + offset[2])
            for candidate in buckets.get(neighbor_key, ()):
                delta = welded[candidate] - point
                if float(np.dot(delta, delta)) <= tol_sq:
                    match = candidate
                    break
            if match is not None:
                break

        if match is None:
            match = len(welded)
            welded.append(point)
            buckets.setdefault(key, []).append(match)
        old_to_new[old_index] = match

    welded_edges = []
    seen_edges = set()
    for start, end in edges:
        welded_start = int(old_to_new[start])
        welded_end = int(old_to_new[end])
        if welded_start == welded_end:
            continue
        edge = (welded_start, welded_end)
        if edge in seen_edges:
            continue
        seen_edges.add(edge)
        welded_edges.append(edge)

    return np.asarray(welded, dtype=float), welded_edges


def _polygon_parts(geom: BaseGeometry) -> list[Polygon]:
    if geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if hasattr(geom, "geoms"):
        result = []
        for part in geom.geoms:
            result.extend(_polygon_parts(part))
        return result
    return []


def _line_parts(geom: BaseGeometry) -> list[LineString]:
    if geom.is_empty:
        return []
    if isinstance(geom, LineString):
        return [geom]
    if isinstance(geom, MultiLineString):
        return list(geom.geoms)
    if hasattr(geom, "geoms"):
        result = []
        for part in geom.geoms:
            result.extend(_line_parts(part))
        return result
    return []


def _is_full_span_loop(coords: np.ndarray, period: float = 2.0 * np.pi, tol: float = 0.05) -> bool:
    if len(coords) < 4:
        return False
    theta = coords[:, 0]
    span = float(np.max(theta) - np.min(theta))
    endpoint_span = abs(float(theta[-1] - theta[0]))
    return abs(span - period) <= tol or abs(endpoint_span - period) <= tol


def _full_span_trace(coords: np.ndarray, z_min: float, z_max: float, period: float = 2.0 * np.pi) -> LineString | None:
    if len(coords) < 4:
        return None
    trace_coords = np.asarray(coords, dtype=float)
    if trace_coords[-1, 0] < trace_coords[0, 0]:
        trace_coords = trace_coords[::-1].copy()

    z_pad = max(1.0, abs(z_max - z_min))
    clip = box(0.0, z_min - z_pad, period, z_max + z_pad)
    fragments = []
    theta_min = float(np.min(trace_coords[:, 0]))
    theta_max = float(np.max(trace_coords[:, 0]))
    k_min = int(np.floor((0.0 - theta_max) / period)) - 1
    k_max = int(np.ceil((period - theta_min) / period)) + 1
    for shift in range(k_min, k_max + 1):
        shifted = trace_coords.copy()
        shifted[:, 0] += shift * period
        inter = LineString(shifted).intersection(clip)
        fragments.extend(part for part in _line_parts(inter) if part.length > 1e-9)
    if not fragments:
        return None
    pieces = []
    for fragment in fragments:
        arr = np.asarray(fragment.coords, dtype=float)
        if len(arr) < 2:
            continue
        if arr[-1, 0] < arr[0, 0]:
            arr = arr[::-1].copy()
        pieces.append(arr)
    if not pieces:
        return None
    pieces.sort(key=lambda arr: (arr[0, 0], arr[-1, 0]))
    stitched = [pieces[0]]
    for arr in pieces[1:]:
        if np.linalg.norm(stitched[-1][-1] - arr[0]) <= 1e-5:
            stitched.append(arr[1:])
        else:
            stitched.append(arr)
    trace_arr = np.vstack(stitched)
    if len(trace_arr) < 2:
        return None
    if trace_arr[0, 0] > 0.05 or trace_arr[-1, 0] < period - 0.05:
        return None
    trace_arr[0, 0] = 0.0
    trace_arr[-1, 0] = period
    return LineString(trace_arr)


def _full_span_strips(walks: list[np.ndarray], z_min: float, z_max: float, period: float = 2.0 * np.pi) -> list[Polygon]:
    traces = []
    for walk in walks:
        if not _is_full_span_loop(walk, period=period):
            continue
        trace = _full_span_trace(walk, z_min, z_max, period=period)
        if trace is None:
            continue
        coords = np.asarray(trace.coords, dtype=float)
        finite_z = coords[:, 1]
        traces.append((float(np.median(finite_z)), coords))
    if len(traces) < 2:
        return []

    traces.sort(key=lambda item: item[0])
    strips = []
    for (_, lower), (_, upper) in zip(traces[0::2], traces[1::2]):
        ring = np.vstack([lower, upper[::-1], lower[:1]])
        try:
            geom = shapely.make_valid(Polygon(ring))
        except (ValueError, shapely.errors.GEOSException):
            continue
        for part in _polygon_parts(geom):
            if part.area > 1e-9:
                strips.append(part)
    return strips


def _point_is_solid(mesh: trimesh.Trimesh, radius: float, theta: float, axial: float) -> bool:
    point = np.array([[radius * np.cos(theta), radius * np.sin(theta), axial]], dtype=float)
    return bool(mesh.contains(point)[0])


def _sample_face_solid_fraction(
    mesh: trimesh.Trimesh,
    radius: float,
    face: Polygon,
    samples_per_axis: int,
) -> tuple[float, int]:
    theta_min, z_min, theta_max, z_max = face.bounds
    solid = 0
    total = 0
    for theta_frac in np.linspace(0.1, 0.9, samples_per_axis):
        for z_frac in np.linspace(0.1, 0.9, samples_per_axis):
            theta = theta_min + (theta_max - theta_min) * float(theta_frac)
            axial = z_min + (z_max - z_min) * float(z_frac)
            p = Point(theta, axial)
            if not face.covers(p):
                continue
            total += 1
            if _point_is_solid(mesh, radius, theta, axial):
                solid += 1
    return (solid / total, total) if total else (0.0, 0)


def _vote_mixed_face(mesh: trimesh.Trimesh, radius: float, face: Polygon, config: SlicerConfig) -> BaseGeometry:
    refine = config.intersection_refine_factor
    theta_step = config.theta_step / refine
    z_step = config.axial_step / refine
    theta_min, z_min, theta_max, z_max = face.bounds
    theta_start = np.floor(theta_min / theta_step) * theta_step + theta_step * 0.5
    z_start = np.floor(z_min / z_step) * z_step + z_step * 0.5
    theta_values = np.arange(theta_start, theta_max + theta_step * 0.5, theta_step)
    z_values = np.arange(z_start, z_max + z_step * 0.5, z_step)
    solid = 0
    total = 0
    for theta in theta_values:
        for axial in z_values:
            p = Point(float(theta), float(axial))
            if not face.contains(p):
                continue
            total += 1
            if _point_is_solid(mesh, radius, float(theta), float(axial)):
                solid += 1
    if total == 0:
        rep = face.representative_point()
        return face if _point_is_solid(mesh, radius, float(rep.x), float(rep.y)) else GeometryCollection()
    return face if solid / total >= _MIXED_FACE_SOLID_THRESHOLD else GeometryCollection()


def _vote_face(
    mesh: trimesh.Trimesh,
    radius: float,
    face: Polygon,
    samples_per_axis: int = 7,
    threshold: float = 0.5,
) -> BaseGeometry:
    fraction, total = _sample_face_solid_fraction(mesh, radius, face, samples_per_axis)
    if total:
        return face if fraction >= threshold else GeometryCollection()
    rep = face.representative_point()
    return face if _point_is_solid(mesh, radius, float(rep.x), float(rep.y)) else GeometryCollection()


def _classify_face(mesh: trimesh.Trimesh, radius: float, face: Polygon, config: SlicerConfig) -> BaseGeometry:
    theta_min, z_min, theta_max, z_max = face.bounds
    samples = [
        face.representative_point(),
        Point((theta_min + theta_max) * 0.5, (z_min + z_max) * 0.5),
    ]
    for theta_frac in (0.25, 0.75):
        for z_frac in (0.25, 0.75):
            samples.append(Point(theta_min + (theta_max - theta_min) * theta_frac, z_min + (z_max - z_min) * z_frac))
    inside = []
    for sample in samples:
        if face.covers(sample):
            inside.append(_point_is_solid(mesh, radius, float(sample.x), float(sample.y)))
    if inside and all(inside):
        if theta_max - theta_min >= TWO_PI * 1.5:
            return _vote_face(mesh, radius, face)
        return face
    if inside and not any(inside):
        return GeometryCollection()
    if config.intersection_refine_factor > 0:
        return _vote_mixed_face(mesh, radius, face, config)
    if inside:
        return _vote_face(mesh, radius, face, samples_per_axis=5)
    rep = face.representative_point()
    return face if _point_is_solid(mesh, radius, float(rep.x), float(rep.y)) else GeometryCollection()


def _has_full_width_horizontal_span(geom: BaseGeometry, samples: int = 32) -> bool:
    if geom.is_empty:
        return False
    z_min = geom.bounds[1]
    z_max = geom.bounds[3]
    if z_max <= z_min:
        return False
    for z in np.linspace(z_min, z_max, samples):
        line = LineString([(0.0, float(z)), (2.0 * np.pi, float(z))])
        inter = geom.intersection(line)
        if inter.length >= 2.0 * np.pi * 0.99:
            return True
    return False


def _loops_from_directed_edges(points: np.ndarray, edges: list[tuple[int, int]]) -> list[np.ndarray]:
    next_of: dict[int, int | None] = {}
    for start, end in edges:
        if start == end:
            continue
        if start in next_of:
            next_of[start] = None
        else:
            next_of[start] = end

    loops = []
    state: dict[int, int] = {}
    for start in next_of:
        if state.get(start, 0) != 0:
            continue
        chain = []
        position = {}
        current = start
        while next_of.get(current) is not None and state.get(current, 0) == 0 and current not in position:
            state[current] = 1
            position[current] = len(chain)
            chain.append(current)
            current = next_of[current]  # type: ignore[assignment]
        if current in position:
            cycle = chain[position[current] :]
            if len(cycle) >= 3:
                loops.append(points[cycle + [cycle[0]]])
        for node in chain:
            state[node] = 2
    return loops


def cylindrical_intersection_segments(
    mesh: trimesh.Trimesh,
    radius: float,
    max_theta_step: float | None = None,
    max_axial_step: float | None = None,
    max_curve_error_mm: float | None = None,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    # Backward compatible explicit angular/axial caps. When supplied, the
    # physical deviation limit is independent and additionally constrains the
    # curvature between exact edge/cylinder intersections.
    if max_curve_error_mm is not None and (
        not math.isfinite(max_curve_error_mm) or max_curve_error_mm <= 0
    ):
        raise ValueError("max_curve_error_mm must be finite and positive")
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    vertex_radii = np.hypot(vertices[:, 0], vertices[:, 1])
    points: list[np.ndarray] = []
    edge_cache: dict[tuple[int, int], int] = {}
    edges: list[tuple[int, int]] = []

    def crossing_index(a: int, b: int) -> int:
        key = (min(a, b), max(a, b))
        if key not in edge_cache:
            start = vertices[a]
            delta = vertices[b] - start
            qa = float(delta[0] * delta[0] + delta[1] * delta[1])
            qb = float(2.0 * (start[0] * delta[0] + start[1] * delta[1]))
            qc = float(start[0] * start[0] + start[1] * start[1] - radius * radius)

            if qa <= 1e-18:
                t = 0.0
            else:
                disc = max(0.0, qb * qb - 4.0 * qa * qc)
                sqrt_disc = float(np.sqrt(disc))
                candidates = [
                    root
                    for root in (
                        (-qb - sqrt_disc) / (2.0 * qa),
                        (-qb + sqrt_disc) / (2.0 * qa),
                    )
                    if -1e-9 <= root <= 1.0 + 1e-9
                ]
                if candidates:
                    t = min(max(candidates[0], 0.0), 1.0)
                else:
                    ra = vertex_radii[a]
                    rb = vertex_radii[b]
                    denom = rb - ra
                    t = 0.0 if abs(denom) < 1e-12 else (radius - ra) / denom
                    t = min(max(float(t), 0.0), 1.0)

            edge_cache[key] = len(points)
            points.append(start + t * delta)
        return edge_cache[key]

    def append_curve_edges(from_idx: int, to_idx: int, face_vertices: np.ndarray) -> None:
        p0 = points[from_idx]
        p1 = points[to_idx]
        theta0 = math.atan2(float(p0[1]), float(p0[0]))
        theta1 = math.atan2(float(p1[1]), float(p1[0]))
        delta_theta = (theta1 - theta0 + math.pi) % (2.0 * math.pi) - math.pi

        normal = np.cross(face_vertices[1] - face_vertices[0], face_vertices[2] - face_vertices[0])
        offset = -float(np.dot(normal, face_vertices[0]))
        # For a nonvertical plane n·(R cos(t),R sin(t),z) + d = 0,
        # z(t) is sinusoidal. Its second derivative magnitude is at most
        # R*hypot(nx,ny)/abs(nz) (mm/radian^2). Linear interpolation error
        # of a twice differentiable function is <= M*dt^2/8.
        # Because u = R*t is linear and z is the only nonlinear unwrapped
        # component, the same bound applies to physical surface distance.
        curvature_bound = (
            radius * math.hypot(float(normal[0]), float(normal[1])) / abs(float(normal[2]))
            if abs(float(normal[2])) > 1e-12 else 0.0
        )
        steps = 1
        if max_curve_error_mm is not None and curvature_bound > 0.0:
            needed = abs(delta_theta) * math.sqrt(curvature_bound / (8.0 * max_curve_error_mm))
            steps = max(steps, int(math.ceil(needed)))
        if max_theta_step is not None and max_theta_step > 0:
            steps = max(steps, int(math.ceil(abs(delta_theta) / max_theta_step)))
        # No artificial subdivision of genuinely straight axial lines in
        # automatic error-bounded mode; otherwise keep the older axial cap.
        if max_axial_step is not None and max_axial_step > 0 and (
            max_curve_error_mm is None or curvature_bound > 1e-12
        ):
            steps = max(steps, int(math.ceil(abs(float(p1[2] - p0[2])) / max_axial_step)))
        if steps > 65536:
            raise ValueError("Intersection curve requires >65536 segments: check face degeneracy/tolerance")
        if steps <= 1:
            edges.append((from_idx, to_idx))
            return
        chain = [from_idx]
        for step in range(1, steps):
            frac = step / steps
            if abs(float(normal[2])) > 1e-12:
                theta = theta0 + delta_theta * frac
                point = np.array(
                    [
                        radius * math.cos(theta),
                        radius * math.sin(theta),
                        -(
                            float(normal[0]) * radius * math.cos(theta)
                            + float(normal[1]) * radius * math.sin(theta)
                            + offset
                        )
                        / float(normal[2]),
                    ],
                    dtype=float,
                )
            else:
                point = p0 + (p1 - p0) * frac
                radial = math.hypot(float(point[0]), float(point[1]))
                if radial > 1e-12:
                    point = point.copy()
                    point[:2] *= radius / radial

            chain.append(len(points))
            points.append(point)
        chain.append(to_idx)
        edges.extend(zip(chain[:-1], chain[1:]))

    for v0, v1, v2 in faces:
        inside = vertex_radii[[v0, v1, v2]] < radius
        if bool(inside[0]) == bool(inside[1]) == bool(inside[2]):
            continue

        from_idx = None
        to_idx = None
        for a, b in ((v0, v1), (v1, v2), (v2, v0)):
            a_inside = vertex_radii[a] < radius
            b_inside = vertex_radii[b] < radius
            if a_inside and not b_inside:
                from_idx = crossing_index(a, b)
            elif not a_inside and b_inside:
                to_idx = crossing_index(a, b)

        if from_idx is not None and to_idx is not None and from_idx != to_idx:
            append_curve_edges(from_idx, to_idx, vertices[[v0, v1, v2]])

    return _weld_intersection_graph(np.asarray(points, dtype=float), edges)


def _sample_intersection_layer_once(
    mesh: trimesh.Trimesh,
    radius: float,
    config: SlicerConfig,
    z_bounds: tuple[float, float] | None = None,
):
    z_min, z_max = z_bounds if z_bounds is not None else axial_bounds(mesh, config.axial_margin)
    clip_z_min, clip_z_max = axial_bounds(mesh, 0.0)
    points_3d, edges = cylindrical_intersection_segments(
        mesh,
        radius,
        # Sampling density follows the physical error of the actual
        # face/cylinder curve, not a global 2-degree angular grid.
        max_axial_step=config.max_axial_segment,
        max_curve_error_mm=config.intersection_curve_error_mm,
    )
    if len(points_3d) == 0 or not edges:
        return GeometryCollection(), (0.0, 2.0 * np.pi), False

    theta_z = _cart_to_theta_z(points_3d)
    loops = _loops_from_directed_edges(theta_z, edges)
    walks = []
    lines = []
    for loop in loops:
        walk = _restore_loop_walk(loop)
        if len(walk) < 4:
            continue
        walks.append(walk)
        lines.append(LineString(walk))

    if not lines:
        return GeometryCollection(), (0.0, 2.0 * np.pi), False

    period = TWO_PI
    loop_span = max(float(np.max(walk[:, 0]) - np.min(walk[:, 0])) for walk in walks)
    if loop_span > period * 1.5:
        theta_min = min(float(line.bounds[0]) for line in lines)
        theta_max = max(float(line.bounds[2]) for line in lines)
        unwrapped_frame = box(theta_min - period, z_min, theta_max + period, z_max)
        arrangement_lines = []
        for line in lines:
            inter = line.intersection(unwrapped_frame)
            arrangement_lines.extend(part for part in _line_parts(inter) if part.length > 1e-12)

        full_span_strips = _full_span_strips(walks, z_min, z_max, period=period)
        strip_boundaries = []
        for strip in full_span_strips:
            shift_min = int(np.floor((unwrapped_frame.bounds[0] - strip.bounds[2]) / period)) - 1
            shift_max = int(np.ceil((unwrapped_frame.bounds[2] - strip.bounds[0]) / period)) + 1
            for shift in range(shift_min, shift_max + 1):
                shifted = affinity.translate(strip.boundary, xoff=shift * period)
                inter = shifted.intersection(unwrapped_frame)
                strip_boundaries.extend(part for part in _line_parts(inter) if part.length > 1e-12)
        arrangement = unary_union([unwrapped_frame.boundary, *arrangement_lines, *strip_boundaries])
        polygons = []
        for face in polygonize(arrangement):
            if not isinstance(face, Polygon) or face.area <= 1e-12:
                continue
            classified = _classify_face(mesh, radius, face, config)
            if not classified.is_empty:
                clipped = classified.intersection(unwrapped_frame)
                if not clipped.is_empty:
                    polygons.append(clipped)

        if not polygons:
            return GeometryCollection(), (0.0, 2.0 * np.pi), False

        geom = unary_union(polygons)
    else:
        central_frame = box(0.0, z_min, period, z_max)
        expanded_frame = box(-period, z_min, 2.0 * period, z_max)
        expanded_lines = []
        for line in lines:
            for shift in (-2, -1, 0, 1, 2):
                shifted = affinity.translate(line, xoff=shift * period)
                inter = shifted.intersection(expanded_frame)
                expanded_lines.extend(part for part in _line_parts(inter) if part.length > 1e-12)
        full_span_strips = _full_span_strips(walks, z_min, z_max, period=period)
        strip_boundaries = []
        for strip in full_span_strips:
            for shift in (-1, 0, 1):
                shifted = affinity.translate(strip.boundary, xoff=shift * period)
                inter = shifted.intersection(expanded_frame)
                strip_boundaries.extend(part for part in _line_parts(inter) if part.length > 1e-12)
        arrangement = unary_union([expanded_frame.boundary, *expanded_lines, *strip_boundaries])
        polygons = []
        for face in polygonize(arrangement):
            if not isinstance(face, Polygon) or face.area <= 1e-12:
                continue
            if face.disjoint(central_frame):
                continue
            classified = _classify_face(mesh, radius, face, config)
            if not classified.is_empty:
                clipped = classified.intersection(central_frame)
                if not clipped.is_empty:
                    polygons.append(clipped)

        if not polygons:
            return GeometryCollection(), (0.0, 2.0 * np.pi), False

        geom = unary_union(polygons)
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    geom = geom.buffer(0)
    geom, window, _ = normalize_periodic(geom, clip_z_min, clip_z_max)
    return geom, window, _has_full_width_horizontal_span(geom)


def _full_span_disagrees_with_mesh(
    mesh: trimesh.Trimesh,
    radius: float,
    geom: BaseGeometry,
    *,
    z_samples: int = 9,
    theta_samples: int = 48,
    min_solid_fraction: float = 0.90,
) -> bool:
    if geom.is_empty:
        return False
    z_min, z_max = geom.bounds[1], geom.bounds[3]
    if z_max <= z_min:
        return False

    theta_values = np.linspace(0.0, TWO_PI, theta_samples, endpoint=False) + math.pi / theta_samples
    for z in np.linspace(z_min, z_max, z_samples):
        probe = LineString([(0.0, float(z)), (TWO_PI, float(z))])
        if geom.intersection(probe).length < TWO_PI * 0.99:
            continue
        points = np.column_stack(
            [
                radius * np.cos(theta_values),
                radius * np.sin(theta_values),
                np.full_like(theta_values, float(z)),
            ]
        )
        if float(np.mean(mesh.contains(points))) < min_solid_fraction:
            return True
    return False


def _looks_like_degenerate_full_span(
    mesh: trimesh.Trimesh,
    radius: float,
    geom: BaseGeometry,
    full: bool,
    config: SlicerConfig,
) -> bool:
    if not full or geom.is_empty:
        return False

    mesh_z_min, mesh_z_max = axial_bounds(mesh, 0.0)
    z_tol = max(config.layer_height, config.nozzle_width, 1e-3)
    if geom.bounds[1] < mesh_z_min - z_tol or geom.bounds[3] > mesh_z_max + z_tol:
        return True

    return _full_span_disagrees_with_mesh(mesh, radius, geom)


def _degenerate_radius_offsets(config: SlicerConfig) -> tuple[float, ...]:
    unit = max(min(config.layer_height, config.nozzle_width), 1e-6)
    offsets: list[float] = []
    for scale in (1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 2.5e-2):
        delta = unit * scale
        offsets.extend((delta, -delta))
    for scale in (5e-2, 1e-1):
        delta = config.layer_height * scale
        offsets.extend((delta, -delta))

    deduped = []
    seen = set()
    for offset in offsets:
        key = round(offset, 12)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(offset)
    return tuple(deduped)


def _recovery_compatible(a, b, radius: float, error_mm: float) -> tuple[bool, float]:
    """Require the two sides of the requested radius to agree geometrically.

    The center sample is already known to be degenerate and cannot act as a
    valid comparison geometry. Sampling on both sides bounds the uncertainty
    of substituting a nearby radius. Topology (islands + holes) must match and
    physical, unwrapped boundary deviation must fit the slice error budget.
    """
    ga, gb = a[0], b[0]
    if ga.is_empty or gb.is_empty or bool(a[2]) != bool(b[2]):
        return False, math.inf
    pa, pb = _polygon_parts(ga), _polygon_parts(gb)
    signature = lambda parts: (len(parts), sum(len(poly.interiors) for poly in parts))
    if signature(pa) != signature(pb):
        return False, math.inf
    ma = shapely.transform(ga.boundary, lambda xy: xy * [radius, 1.0])
    mb = shapely.transform(gb.boundary, lambda xy: xy * [radius, 1.0])
    deviation = float(ma.hausdorff_distance(mb))
    return math.isfinite(deviation) and deviation <= 2.0 * error_mm, deviation


def sample_intersection_layer(
    mesh: trimesh.Trimesh,
    radius: float,
    config: SlicerConfig,
    z_bounds: tuple[float, float] | None = None,
):
    first_error = None
    try:
        geom, window, full = _sample_intersection_layer_once(mesh, radius, config, z_bounds)
    except shapely.errors.GEOSException as exc:
        first_error = exc
        geom, window, full = GeometryCollection(), (0.0, TWO_PI), False
    else:
        if not _looks_like_degenerate_full_span(mesh, radius, geom, full, config):
            return geom, window, full

    # Restrict radius substitution by both sampling precision and deposited
    # slab thickness. A 10%-layer offset is not automatically acceptable if
    # the geometry can vary sharply over that radial distance.
    max_shift = 0.1 * config.layer_height
    if config.intersection_strict_recovery:
        max_shift = min(max_shift, config.intersection_curve_error_mm)
    for delta in _degenerate_radius_offsets(config):
        if abs(delta) > max_shift + 1e-12:
            continue
        try:
            candidate = _sample_intersection_layer_once(mesh, radius + delta, config, z_bounds)
        except shapely.errors.GEOSException:
            continue
        if _looks_like_degenerate_full_span(mesh, radius + delta, candidate[0], candidate[2], config):
            continue

        # A true radial tangent can change topology at the requested radius.
        # Sampling the opposite side may therefore disagree even when the
        # chosen one-sided geometry is physically meaningful. In ordinary
        # mode preserve that legacy behavior, but explicitly WARN of the
        # unverified substitution. In strict mode require two-sided consensus.
        compatible, deviation = False, math.inf
        try:
            opposite = _sample_intersection_layer_once(mesh, radius - delta, config, z_bounds)
        except shapely.errors.GEOSException:
            pass
        else:
            if not _looks_like_degenerate_full_span(mesh, radius - delta, opposite[0], opposite[2], config):
                compatible, deviation = _recovery_compatible(
                    candidate, opposite, radius, config.intersection_curve_error_mm
                )
        if not compatible and config.intersection_strict_recovery:
            continue
        status = ("two-sided verified" if compatible
                  else "ONE-SIDED, topology or dimensional disagreement unverified")
        if abs(delta) > config.intersection_curve_error_mm:
            status += "; SHIFT EXCEEDS CURVE-ERROR BUDGET"
        warnings.warn(
            f"Intersection layer R={radius:.6f} mm recovered from R={radius + delta:.6f} mm "
            f"(radial shift={delta:+.6g} mm, {status}, "
            f"boundary deviation={deviation:.6g} mm)",
            RuntimeWarning,
            stacklevel=2,
        )
        return candidate

    # Never silently return geometry already flagged as degenerate.
    # Both modes require a valid, bounded shifted sample; strict mode also
    # requires two-sided topological/dimensional consensus.
    if first_error is not None:
        raise RuntimeError(
            f"Intersection layer R={radius:.6f} mm failed bounded two-sided recovery"
        ) from first_error
    raise RuntimeError(
        f"Intersection layer R={radius:.6f} mm appears degenerate and no nearby "
        f"geometry agrees within {config.intersection_curve_error_mm:g} mm"
    )


def _init_worker(mesh: trimesh.Trimesh, config: SlicerConfig, z_bounds: tuple[float, float]) -> None:
    global _WORKER_MESH, _WORKER_CONFIG, _WORKER_Z_BOUNDS
    _WORKER_MESH = mesh
    _WORKER_CONFIG = config
    _WORKER_Z_BOUNDS = z_bounds


def _sample_worker(item: tuple[int, float]):
    if _WORKER_MESH is None or _WORKER_CONFIG is None or _WORKER_Z_BOUNDS is None:
        raise RuntimeError("intersection worker was not initialized")
    index, radius = item
    geom, window, full = sample_intersection_layer(_WORKER_MESH, radius, _WORKER_CONFIG, _WORKER_Z_BOUNDS)
    return index, geom, window, full


def sample_intersection_layers(mesh: trimesh.Trimesh, config: SlicerConfig, radii: list[float]):
    z_bounds = axial_bounds(mesh, config.axial_margin)
    worker_count = min(config.slice_workers, len(radii))
    if worker_count <= 1:
        return [
            (i, *sample_intersection_layer(mesh, radius, config, z_bounds))
            for i, radius in enumerate(radii)
        ]

    chunksize = max(1, len(radii) // (worker_count * 4))
    with ProcessPoolExecutor(
        max_workers=worker_count,
        initializer=_init_worker,
        initargs=(mesh, config, z_bounds),
    ) as executor:
        return list(executor.map(_sample_worker, enumerate(radii), chunksize=chunksize))
