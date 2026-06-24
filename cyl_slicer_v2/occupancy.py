from __future__ import annotations

import numpy as np
import shapely
from shapely.geometry import GeometryCollection, box
from shapely.ops import unary_union
import trimesh

from .config import SlicerConfig
from .coordinates import cylindrical_to_cartesian
from .mesh_io import axial_bounds
from .periodic import normalize_periodic


def axial_edges(mesh: trimesh.Trimesh, config: SlicerConfig) -> np.ndarray:
    z_min, z_max = axial_bounds(mesh, config.axial_margin)
    count = max(1, int(np.ceil((z_max - z_min) / config.axial_step)))
    return np.linspace(z_min, z_min + count * config.axial_step, count + 1)


def theta_edges(config: SlicerConfig) -> np.ndarray:
    step = config.theta_step
    n = int(np.ceil((2.0 * np.pi) / step))
    step = (2.0 * np.pi) / n
    start = -config.guard_columns * step
    stop = 2.0 * np.pi + config.guard_columns * step
    return np.arange(start, stop + step * 0.5, step)


def occupancy_mask(mesh: trimesh.Trimesh, radius: float, t_edges: np.ndarray, z_edges: np.ndarray) -> np.ndarray:
    t_centers = (t_edges[:-1] + t_edges[1:]) * 0.5
    z_centers = (z_edges[:-1] + z_edges[1:]) * 0.5
    tt, zz = np.meshgrid(t_centers, z_centers, indexing="ij")
    points = cylindrical_to_cartesian(tt.ravel(), zz.ravel(), radius)
    inside = mesh.contains(points)
    return inside.reshape((len(t_centers), len(z_centers)))


def mask_to_geometry(mask: np.ndarray, t_edges: np.ndarray, z_edges: np.ndarray) -> shapely.Geometry:
    cells = []
    for ti, zi in np.argwhere(mask):
        cells.append(box(t_edges[ti], z_edges[zi], t_edges[ti + 1], z_edges[zi + 1]))
    if not cells:
        return GeometryCollection()
    geom = unary_union(cells)
    if not geom.is_valid:
        geom = shapely.make_valid(geom)
    return geom


def sample_layer(mesh: trimesh.Trimesh, radius: float, config: SlicerConfig) -> tuple[shapely.Geometry, tuple[float, float], bool]:
    z_edges = axial_edges(mesh, config)
    t_edges = theta_edges(config)
    mask = occupancy_mask(mesh, radius, t_edges, z_edges)
    raw = mask_to_geometry(mask, t_edges, z_edges)
    return normalize_periodic(raw, float(z_edges[0]), float(z_edges[-1]))
