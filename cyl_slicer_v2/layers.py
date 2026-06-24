from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
import shapely
from shapely.geometry import GeometryCollection
from shapely.ops import unary_union
import trimesh

from .config import SlicerConfig
from .mesh_io import radial_bounds
from .models import LayerSlice
from .occupancy_backends import MeshContainsBackend, OccupancyBackend
from .occupancy import axial_edges, sample_layer_with_edges, theta_edges

_WORKER_BACKEND: OccupancyBackend | None = None
_WORKER_CONFIG: SlicerConfig | None = None
_WORKER_T_EDGES: np.ndarray | None = None
_WORKER_Z_EDGES: np.ndarray | None = None


def _clean(geom):
    if geom.is_empty:
        return GeometryCollection()
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    return geom


def layer_radii(mesh: trimesh.Trimesh, config: SlicerConfig) -> list[float]:
    r_min_mesh, r_max_mesh = radial_bounds(mesh)
    r_min = config.min_radius if config.min_radius is not None else r_min_mesh
    r_max = config.max_radius if config.max_radius is not None else r_max_mesh
    r_min = max(r_min, config.bed_radius)
    start = r_min + config.layer_height * 0.5
    stop = r_max - config.layer_height * 0.5
    if stop < start:
        return []
    count = int(np.floor((stop - start) / config.layer_height)) + 1
    return [float(start + i * config.layer_height) for i in range(count)]


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
    radii = layer_radii(mesh, config)
    if max_layers is not None:
        radii = radii[:max_layers]

    occupied = [GeometryCollection() for _ in radii]
    windows = [(0.0, 0.0) for _ in radii]
    full_flags = [False for _ in radii]
    for i, geom, window, full in _sample_layers(mesh, config, radii, backend):
        occupied[i] = _clean(geom)
        windows[i] = window
        full_flags[i] = full

    first_non_empty = next((i for i, geom in enumerate(occupied) if not geom.is_empty), None)
    layers: list[LayerSlice] = []
    for i, geom in enumerate(occupied):
        prev = occupied[i - 1] if i > 0 else GeometryCollection()
        nxt = occupied[i + 1] if i + 1 < len(occupied) else GeometryCollection()
        top = _clean(geom.difference(nxt)) if not geom.is_empty else GeometryCollection()
        bottom = _clean(geom.difference(prev)) if not geom.is_empty else GeometryCollection()
        if first_non_empty is None or i <= first_non_empty or geom.is_empty:
            needs = GeometryCollection()
        else:
            needs = _clean(geom.difference(prev))
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
            )
        )

    carry = GeometryCollection()
    for layer in reversed(layers):
        support_here = _clean(carry.difference(layer.occupied)) if not carry.is_empty else GeometryCollection()
        layer.support = support_here
        pieces = []
        if not support_here.is_empty:
            pieces.append(support_here)
        if not layer.needs_support.is_empty:
            pieces.append(layer.needs_support)
        carry = _clean(unary_union(pieces)) if pieces else GeometryCollection()

    return layers
