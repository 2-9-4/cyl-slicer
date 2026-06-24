from __future__ import annotations

import numpy as np
import shapely
from shapely.geometry import GeometryCollection
from shapely.ops import unary_union
import trimesh

from .config import SlicerConfig
from .mesh_io import radial_bounds
from .models import LayerSlice
from .occupancy import sample_layer


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


def slice_mesh(mesh: trimesh.Trimesh, config: SlicerConfig, max_layers: int | None = None) -> list[LayerSlice]:
    config.validate()
    radii = layer_radii(mesh, config)
    if max_layers is not None:
        radii = radii[:max_layers]

    occupied = []
    windows = []
    full_flags = []
    for radius in radii:
        geom, window, full = sample_layer(mesh, radius, config)
        occupied.append(_clean(geom))
        windows.append(window)
        full_flags.append(full)

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
