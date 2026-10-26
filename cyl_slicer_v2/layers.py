from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import math

import numpy as np
import shapely
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, box
from shapely.ops import unary_union
import trimesh

from .bridges import classify_bridge_regions
from .config import SlicerConfig
from .coordinates import TWO_PI
from .intersection_backend import sample_intersection_layers
from .mesh_io import radial_bounds
from .models import LayerSlice
from .occupancy_backends import MeshContainsBackend, OccupancyBackend
from .occupancy import axial_edges, sample_layer_with_edges, theta_edges
from .periodic import translate_theta

_WORKER_BACKEND: OccupancyBackend | None = None
_WORKER_CONFIG: SlicerConfig | None = None
_WORKER_T_EDGES: np.ndarray | None = None
_WORKER_Z_EDGES: np.ndarray | None = None


def _clean(geom):
    if geom.is_empty:
        return GeometryCollection()
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    # A zero-width buffer can discard valid components of a MultiPolygon when
    # very thin components are present. Preserve valid surfaces; only strip non-area
    # pieces that make_valid or a polygon intersection may have produced.
    if not isinstance(geom, (Polygon, MultiPolygon)):
        geom = _safe_union(_polygon_parts(geom))
    if geom.is_empty or geom.area <= 1e-12:
        return GeometryCollection()
    return geom


def _safe_union(geoms):
    parts = [part for part in geoms if part is not None and not part.is_empty]
    if not parts:
        return GeometryCollection()
    if len(parts) == 1:
        return _clean(parts[0])
    try:
        return _clean(unary_union(parts))
    except shapely.errors.GEOSException:
        repaired = []
        for part in parts:
            part = _clean(part)
            if not part.is_empty:
                repaired.append(part)
        if not repaired:
            return GeometryCollection()
        return _clean(shapely.union_all(repaired, grid_size=1e-8))


def _periodic_union(geom, turns: int = 1):
    if geom is None or geom.is_empty:
        return GeometryCollection()
    shifts = (k * TWO_PI for k in range(-turns, turns + 1))
    return _safe_union(translate_theta(geom, shift) for shift in shifts)


def _safe_intersection(a, b):
    if a is None or b is None or a.is_empty or b.is_empty:
        return GeometryCollection()
    try:
        return _clean(a.intersection(b))
    except shapely.errors.GEOSException:
        left = _clean(a)
        right = _clean(b)
        if left.is_empty or right.is_empty:
            return GeometryCollection()
        return _clean(shapely.intersection(left, right, grid_size=1e-8))


def _safe_difference(a, b):
    if a is None or a.is_empty:
        return GeometryCollection()
    if b is None or b.is_empty:
        return _clean(a)
    try:
        return _clean(a.difference(b))
    except shapely.errors.GEOSException:
        left = _clean(a)
        right = _clean(b)
        if left.is_empty:
            return GeometryCollection()
        if right.is_empty:
            return left
        return _clean(shapely.difference(left, right, grid_size=1e-8))


def _polygon_parts(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return list(geom.geoms)
    if hasattr(geom, "geoms"):
        result = []
        for part in geom.geoms:
            result.extend(_polygon_parts(part))
        return result
    return []


def _to_metric(geom, radius: float):
    scale = max(radius, 1e-9)
    return shapely.transform(geom, lambda p: p * [scale, 1.0])


def _from_metric(geom, radius: float):
    scale = max(radius, 1e-9)
    return shapely.transform(geom, lambda p: p * [1.0 / scale, 1.0])


def _metric_buffer(geom, radius: float, distance: float):
    metric = _to_metric(geom, radius)
    buffered = metric.buffer(distance)
    if buffered.is_empty:
        return GeometryCollection()
    raw = _from_metric(buffered, radius)
    return _clean(raw)


def _metric_area(geom, radius: float) -> float:
    if geom is None or geom.is_empty:
        return 0.0
    return float(_to_metric(geom, radius).area)


def _fold_metric_result(raw):
    raw = _clean(raw)
    if raw.is_empty:
        return GeometryCollection()
    z_min, z_max = raw.bounds[1], raw.bounds[3]
    clip = box(0.0, z_min, TWO_PI, z_max)
    min_k = math.floor((raw.bounds[0] - TWO_PI) / TWO_PI)
    max_k = math.ceil((raw.bounds[2] + TWO_PI) / TWO_PI)
    pieces = []
    for k in range(min_k, max_k + 1):
        part = _safe_intersection(translate_theta(raw, -k * TWO_PI), clip)
        if not part.is_empty:
            pieces.append(part)
    return _safe_union(pieces)


def _metric_region_filter(
    geom,
    radius: float,
    *,
    min_area_mm2: float,
    promote_area_mm2: float = 0.0,
    merge_gap_mm: float = 0.0,
):
    if geom is None or geom.is_empty:
        return GeometryCollection()

    geom = _clean(geom)
    if geom.is_empty:
        return GeometryCollection()

    distance_margin = max(
        float(np.sqrt(promote_area_mm2 / np.pi)) if promote_area_mm2 > 0 else 0.0,
        merge_gap_mm * 0.5,
    )
    theta_margin = distance_margin / max(radius, 1e-9) + 1e-9
    theta_min, _, theta_max, _ = geom.bounds
    can_cross_boundary = distance_margin > 0
    reaches_lower_boundary = theta_min <= theta_margin
    reaches_upper_boundary = theta_max >= TWO_PI - theta_margin
    out_of_range = theta_min < -theta_margin or theta_max > TWO_PI + theta_margin
    needs_periodic_neighbors = (
        out_of_range
        or (can_cross_boundary and reaches_lower_boundary and reaches_upper_boundary)
    )
    needs_fold = out_of_range or (can_cross_boundary and (reaches_lower_boundary or reaches_upper_boundary))
    shifts = (-TWO_PI, 0.0, TWO_PI) if needs_periodic_neighbors else (0.0,)

    kept = []
    for shift in shifts:
        source = translate_theta(geom, shift) if shift else geom
        if source.is_empty:
            continue
        metric = _to_metric(source, radius)
        for part in _polygon_parts(metric):
            if part.is_empty:
                continue
            if not part.is_valid:
                part = shapely.make_valid(part)
            part = part.buffer(0)
            if part.is_empty:
                continue
            area = float(part.area)
            if area < min_area_mm2:
                continue
            if promote_area_mm2 > area:
                target_radius = float(np.sqrt(promote_area_mm2 / np.pi))
                current_radius = float(np.sqrt(max(area, 1e-12) / np.pi))
                part = part.buffer(target_radius - current_radius)
                if part.is_empty:
                    continue
            kept.append(part)

    if not kept:
        return GeometryCollection()

    metric = _safe_union(kept)
    if merge_gap_mm > 0:
        metric = _clean(metric.buffer(merge_gap_mm * 0.5).buffer(-merge_gap_mm * 0.5))
    if metric.is_empty:
        return GeometryCollection()
    raw = _from_metric(metric, radius)
    return _fold_metric_result(raw) if needs_fold else _clean(raw)


def _bearing_regions(geom, radius: float, config: SlicerConfig):
    return _metric_region_filter(
        geom,
        radius,
        min_area_mm2=config.support_bearing_min_area,
    )


def _support_regions(geom, radius: float, config: SlicerConfig, *, promote: bool = True):
    return _metric_region_filter(
        geom,
        radius,
        min_area_mm2=config.support_min_island_area,
        promote_area_mm2=config.support_min_region_area if promote else 0.0,
        merge_gap_mm=config.support_merge_gap,
    )


def _support_interface_regions(geom, radius: float, config: SlicerConfig):
    return _metric_region_filter(
        geom,
        radius,
        min_area_mm2=config.support_min_island_area,
    )


def _unsupported_from_previous(geom, previous, radius: float, config: SlicerConfig):
    if geom.is_empty:
        return GeometryCollection()
    if previous.is_empty:
        return _clean(geom)
    previous_cover = _metric_buffer(
        _periodic_union(previous, turns=1),
        radius,
        config.nozzle_width * 0.5,
    )
    uncovered = _safe_difference(geom, previous_cover)
    if _metric_area(uncovered, radius) < max(config.support_min_island_area, 1e-9):
        return GeometryCollection()

    covered = _safe_intersection(geom, previous_cover)
    bearing_cover = _bearing_regions(covered, radius, config)
    return _safe_difference(geom, bearing_cover)


def _needs_support_from_previous(geom, previous, radius: float, config: SlicerConfig):
    return _support_regions(_unsupported_from_previous(geom, previous, radius, config), radius, config)


def _merge_support_carry(items, radius: float, config: SlicerConfig):
    grouped = {}
    for geom, interface_target, gap_layers, interface_layers in items:
        if geom is None or geom.is_empty:
            continue
        grouped.setdefault((gap_layers, interface_layers), []).append((geom, interface_target))

    merged = []
    for (gap_layers, interface_layers), items in grouped.items():
        geoms = [geom for geom, _ in items]
        geom = _support_regions(_safe_union(geoms), radius, config, promote=False)
        if not geom.is_empty:
            interface_target = _safe_union(
                target for _, target in items if target is not None and not target.is_empty
            )
            merged.append((geom, interface_target, gap_layers, interface_layers))
    return merged


def propagate_support(layers: list[LayerSlice], config: SlicerConfig) -> None:
    carry = []
    for layer in reversed(layers):
        support_pieces = []
        interface_pieces = []
        interface_direction = None
        next_carry = []
        for carried_geom, interface_target, gap_layers, interface_layers in carry:
            occupied = _periodic_union(layer.occupied, turns=1) if not layer.occupied.is_empty else GeometryCollection()
            if occupied.is_empty:
                ungrounded = carried_geom
            else:
                contact = _safe_intersection(carried_geom, occupied)
                grounded = _bearing_regions(contact, layer.radius, config)
                ungrounded = _safe_difference(carried_geom, grounded)
            carried_support = _support_regions(ungrounded, layer.radius, config, promote=False)

            printable_support = carried_support
            if not occupied.is_empty and not carried_support.is_empty:
                # Same-layer support is removed inside the model, optionally
                # with a physical-mm buffer for easy-release lateral clearance.
                # Carry remains UNbuffered so lower support foundations are
                # not eroded or silently disconnected.
                clearance = (
                    _metric_buffer(occupied, layer.radius, config.support_lateral_gap_mm)
                    if config.support_lateral_gap_mm > 0 else occupied
                )
                printable_support = _fold_metric_result(
                    _safe_difference(carried_support, clearance)
                )
                printable_support = _support_regions(
                    printable_support, layer.radius, config, promote=False
                )

            if gap_layers <= 0 and not printable_support.is_empty:
                printable_contact = _safe_intersection(printable_support, interface_target)
                printable_contact = _support_interface_regions(printable_contact, layer.radius, config)
                if printable_contact.is_empty:
                    if not carried_support.is_empty:
                        next_gap = max(0, gap_layers - 1)
                        next_interface = interface_layers if gap_layers > 0 else max(0, interface_layers - 1)
                        next_carry.append((carried_support, interface_target, next_gap, next_interface))
                    continue

                support_pieces.append(printable_support)
                if interface_layers > 0:
                    interface_pieces.append(printable_contact)
                    interface_index = config.support_interface_layers - interface_layers
                    direction = "axial" if interface_index % 2 == 0 else "radial"
                    if interface_direction is None:
                        interface_direction = direction

            if not carried_support.is_empty:
                next_gap = max(0, gap_layers - 1)
                next_interface = interface_layers if gap_layers > 0 else max(0, interface_layers - 1)
                next_carry.append((carried_support, interface_target, next_gap, next_interface))

        layer.support = _safe_union(support_pieces)
        layer.support_interface = _safe_union(interface_pieces)
        layer.support_interface_direction = interface_direction

        if not layer.needs_support.is_empty:
            interface_target = (
                layer.support_interface_target
                if layer.support_interface_target is not None and not layer.support_interface_target.is_empty
                else layer.needs_support
            )
            next_carry.append(
                (
                    layer.needs_support,
                    interface_target,
                    config.effective_support_gap_layers,
                    config.support_interface_layers,
                )
            )
        carry = _merge_support_carry(next_carry, layer.radius, config)


def _assign_printable_bridges(layers: list[LayerSlice], config: SlicerConfig) -> None:
    """Replace only confirmed two-ended short bridge support with bridge fills."""
    if config.bridge_mode != "print":
        return
    for previous, layer in zip(layers, layers[1:]):
        regions = classify_bridge_regions(
            layer.occupied, previous.occupied, layer.radius,
            config.nozzle_width, config.bridge_max_span_mm,
        )
        parts = [part for part in regions.values() if not part.is_empty]
        if not parts:
            continue
        mask = _safe_union(parts)
        layer.bridge_regions = regions
        if not layer.needs_support.is_empty:
            layer.needs_support = _safe_difference(layer.needs_support, mask)
        if layer.support_interface_target is not None and not layer.support_interface_target.is_empty:
            layer.support_interface_target = _safe_difference(
                layer.support_interface_target, mask
            )


def propagate_solid_surfaces(layers: list[LayerSlice], config: SlicerConfig) -> None:
    top_sources = [
        (index, layer.top_surfaces)
        for index, layer in enumerate(layers)
        if not layer.top_surfaces.is_empty
    ]
    bottom_sources = [
        (index, layer.bottom_surfaces)
        for index, layer in enumerate(layers)
        if not layer.bottom_surfaces.is_empty
    ]

    for index, surface in top_sources:
        carry = surface
        for depth in range(1, config.top_layers):
            target_index = index - depth
            if target_index < 0:
                break
            target = layers[target_index]
            carry = _safe_intersection(carry, target.occupied)
            if carry.is_empty:
                break
            target.top_surfaces = _safe_union([target.top_surfaces, carry])

    for index, surface in bottom_sources:
        carry = surface
        for depth in range(1, config.bottom_layers):
            target_index = index + depth
            if target_index >= len(layers):
                break
            target = layers[target_index]
            carry = _safe_intersection(carry, target.occupied)
            if carry.is_empty:
                break
            target.bottom_surfaces = _safe_union([target.bottom_surfaces, carry])


def layer_slabs(mesh: trimesh.Trimesh, config: SlicerConfig) -> list[tuple[float, float]]:
    """Sample midpoint radius and physical thickness for each radial slab.

    The original fixed schedule remains byte-for-byte equivalent by default.
    The opt-in complete-terminal schedule covers the exact requested radial
    interval; its final two slabs share a tiny remainder rather than creating
    an unreliable, arbitrarily thin extrusion.
    """
    r_min_mesh, r_max_mesh = radial_bounds(mesh)
    r_min = config.min_radius if config.min_radius is not None else r_min_mesh
    r_max = config.max_radius if config.max_radius is not None else r_max_mesh
    r_min = max(r_min, config.bed_radius)
    h = config.layer_height
    print(f"layer radii goes up {r_min:.2f}->{r_max:.2f}")
    if config.radial_layer_mode == "fixed":
        start = r_min + h * 0.5
        stop = r_max - h * 0.5
        if stop < start:
            return []
        count = int(np.floor((stop - start) / h)) + 1
        return [(float(start + i * h), h) for i in range(count)]

    span = r_max - r_min
    if span <= 0:
        return []
    if span < h * 0.5 - 1e-9:
        raise ValueError("Complete-terminal radial range is less than half a layer; "
                         "reduce --layer-height")
    full_count = int(math.floor(span / h + 1e-10))
    remainder = span - full_count * h
    if abs(remainder) < 1e-9:
        remainder = 0.0
    heights = [h] * full_count
    if remainder > 0:
        if not heights:
            heights.append(remainder)
        elif remainder < h * 0.5:
            shared = (heights.pop() + remainder) * 0.5
            heights.extend((shared, shared))
        else:
            heights.append(remainder)
    edge = r_min
    slabs = []
    for height in heights:
        slabs.append((float(edge + 0.5 * height), float(height)))
        edge += height
    return slabs


def layer_radii(mesh: trimesh.Trimesh, config: SlicerConfig) -> list[float]:
    return [radius for radius, _ in layer_slabs(mesh, config)]


def _init_slice_worker(backend: OccupancyBackend, config: SlicerConfig, t_edges: np.ndarray, z_edges: np.ndarray) -> None:
    global _WORKER_BACKEND, _WORKER_CONFIG, _WORKER_T_EDGES, _WORKER_Z_EDGES
    _WORKER_BACKEND = backend
    _WORKER_CONFIG = config
    _WORKER_T_EDGES = t_edges
    _WORKER_Z_EDGES = z_edges


def _sample_layer_worker(item: tuple[int, float]):
    if _WORKER_BACKEND is None or _WORKER_CONFIG is None or _WORKER_T_EDGES is None or _WORKER_Z_EDGES is None:
        raise RuntimeError("slice worker was not initialized")
    index, radius = item
    geom, window, full = sample_layer_with_edges(
        _WORKER_BACKEND,
        radius,
        _WORKER_CONFIG,
        _WORKER_T_EDGES,
        _WORKER_Z_EDGES,
    )
    return index, geom, window, full


def _sample_layers(mesh: trimesh.Trimesh, config: SlicerConfig, radii: list[float], backend: OccupancyBackend):
    t_edges = theta_edges(config)
    z_edges = axial_edges(mesh, config)
    worker_count = min(config.slice_workers, len(radii))
    if worker_count <= 1:
        return [
            (i, *sample_layer_with_edges(backend, radius, config, t_edges, z_edges))
            for i, radius in enumerate(radii)
        ]

    chunksize = max(1, len(radii) // (worker_count * 4))
    with ProcessPoolExecutor(
        max_workers=worker_count,
        initializer=_init_slice_worker,
        initargs=(backend, config, t_edges, z_edges),
    ) as executor:
        return list(executor.map(_sample_layer_worker, enumerate(radii), chunksize=chunksize))


def slice_mesh(
    mesh: trimesh.Trimesh,
    config: SlicerConfig,
    max_layers: int | None = None,
    occupancy_backend: OccupancyBackend | None = None,
) -> list[LayerSlice]:
    config.validate()
    backend = occupancy_backend or MeshContainsBackend(mesh)
    slabs = layer_slabs(mesh, config)
    if max_layers is not None:
        slabs = slabs[:max_layers]
    radii = [radius for radius, _ in slabs]

    occupied = [GeometryCollection() for _ in radii]
    windows = [(0.0, 0.0) for _ in radii]
    full_flags = [False for _ in radii]
    if config.layer_backend == "intersection":
        sampled = sample_intersection_layers(mesh, config, radii)
    else:
        sampled = _sample_layers(mesh, config, radii, backend)
    for i, geom, window, full in sampled:
        occupied[i] = _clean(geom)
        windows[i] = window
        full_flags[i] = full

    layers: list[LayerSlice] = []
    for i, geom in enumerate(occupied):
        prev = occupied[i - 1] if i > 0 else GeometryCollection()
        nxt = occupied[i + 1] if i + 1 < len(occupied) else GeometryCollection()
        top = _clean(geom.difference(nxt)) if not geom.is_empty else GeometryCollection()
        bottom = _clean(geom.difference(prev)) if not geom.is_empty else GeometryCollection()
        unsupported = _unsupported_from_previous(geom, prev, radii[i], config) if i > 0 else GeometryCollection()
        needs = _support_regions(unsupported, radii[i], config) if i > 0 else GeometryCollection()
        interface_target = _support_interface_regions(unsupported, radii[i], config) if i > 0 else GeometryCollection()
        layers.append(
            LayerSlice(
                index=i,
                radius=radii[i],
                occupied=geom,
                top_surfaces=top,
                bottom_surfaces=bottom,
                needs_support=needs,
                support=GeometryCollection(),
                theta_window=windows[i],
                is_full_wrap=full_flags[i],
                support_interface_target=interface_target,
                previous_occupied=(prev if i > 0 else None),
                deposited_thickness=(slabs[i][1] if config.radial_layer_mode == "complete-terminal" else None),
            )
        )

    _assign_printable_bridges(layers, config)
    propagate_solid_surfaces(layers, config)
    propagate_support(layers, config)

    return layers
