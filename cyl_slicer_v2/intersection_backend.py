from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

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
from .periodic import fold_to_canonical

_WORKER_MESH: trimesh.Trimesh | None = None
_WORKER_CONFIG: SlicerConfig | None = None
_WORKER_Z_BOUNDS: tuple[float, float] | None = None
_MIXED_FACE_SOLID_THRESHOLD = 0.35


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


def _solid_mask(mesh: trimesh.Trimesh, radius: float, theta_z: list[tuple[float, float]]) -> np.ndarray:
    if not theta_z:
        return np.zeros(0, dtype=bool)
    arr = np.asarray(theta_z, dtype=float)
    points = np.column_stack(
        [
            radius * np.cos(arr[:, 0]),
            radius * np.sin(arr[:, 0]),
            arr[:, 1],
        ]
    )
    return np.asarray(mesh.contains(points), dtype=bool)


def _classify_face_cells(mesh: trimesh.Trimesh, radius: float, face: Polygon, config: SlicerConfig) -> BaseGeometry:
    theta_step = config.theta_step
    z_step = config.axial_step
    theta_min, z_min, theta_max, z_max = face.bounds
    theta = np.floor(theta_min / theta_step) * theta_step
    cell_bounds = []
    probes = []
    while theta < theta_max - 1e-12:
        next_theta = min(theta + theta_step, theta_max)
        z = np.floor(z_min / z_step) * z_step
        while z < z_max - 1e-12:
            next_z = min(z + z_step, z_max)
            probe = Point((theta + next_theta) * 0.5, (z + next_z) * 0.5)
            if face.covers(probe):
                cell_bounds.append((theta, z, next_theta, next_z))
                probes.append((float(probe.x), float(probe.y)))
            z += z_step
        theta += theta_step

    if not cell_bounds:
        probe = face.representative_point()
        return face if _point_is_solid(mesh, radius, float(probe.x), float(probe.y)) else GeometryCollection()
    solid = _solid_mask(mesh, radius, probes)
    kept = []
    for (theta0, z0, theta1, z1), is_solid in zip(cell_bounds, solid):
        if not is_solid:
            continue
        cell = face.intersection(box(theta0, z0, theta1, z1))
        if not cell.is_empty and cell.area > 1e-12:
            probe = cell.representative_point()
            if _point_is_solid(mesh, radius, float(probe.x), float(probe.y)):
                kept.append(cell)
    if not kept:
        return GeometryCollection()
    geom = unary_union(kept)
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    return geom.buffer(0)


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
    if inside and all(inside) and len(inside) >= 3:
        return face
    if inside and not any(inside) and len(inside) >= 3:
        return GeometryCollection()
    if config.intersection_refine_factor > 0:
        return _vote_mixed_face(mesh, radius, face, config)
    if inside and len(inside) < 3:
        return _classify_face_cells(mesh, radius, face, config)
    if inside:
        fraction, total = _sample_face_solid_fraction(mesh, radius, face, 5)
        if total:
            return face if fraction >= 0.5 else GeometryCollection()
        return face if sum(inside) / len(inside) >= 0.5 else GeometryCollection()
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


def cylindrical_intersection_segments(mesh: trimesh.Trimesh, radius: float) -> tuple[np.ndarray, list[tuple[int, int]]]:
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    vertex_radii = np.hypot(vertices[:, 0], vertices[:, 1])
    points: list[np.ndarray] = []
    edge_cache: dict[tuple[int, int], int] = {}
    edges: list[tuple[int, int]] = []

    def crossing_index(a: int, b: int) -> int:
        key = (min(a, b), max(a, b))
        if key not in edge_cache:
            ra = vertex_radii[a]
            rb = vertex_radii[b]
            denom = rb - ra
            if abs(denom) < 1e-12:
                edge_cache[key] = len(points)
                points.append(vertices[a].copy())
            else:
                t = (radius - ra) / denom
                edge_cache[key] = len(points)
                points.append(vertices[a] + t * (vertices[b] - vertices[a]))
        return edge_cache[key]

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
            edges.append((from_idx, to_idx))

    return np.asarray(points, dtype=float), edges


def sample_intersection_layer(
    mesh: trimesh.Trimesh,
    radius: float,
    config: SlicerConfig,
    z_bounds: tuple[float, float] | None = None,
):
    z_min, z_max = z_bounds if z_bounds is not None else axial_bounds(mesh, config.axial_margin)
    points_3d, edges = cylindrical_intersection_segments(mesh, radius)
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

        raw = unary_union(polygons)
        if not raw.is_valid:
            raw = shapely.make_valid(raw)
        raw = raw.buffer(0)
        geom = fold_to_canonical(raw, z_min, z_max)
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
    return geom, (0.0, 2.0 * np.pi), _has_full_width_horizontal_span(geom)


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
