from __future__ import annotations

from pathlib import Path

import numpy as np
import trimesh


def load_mesh(path: str | Path) -> trimesh.Trimesh:
    mesh = trimesh.load_mesh(str(path), force="mesh")
    if not isinstance(mesh, trimesh.Trimesh):
        mesh = trimesh.util.concatenate(tuple(mesh.geometry.values()))
    mesh = mesh.copy()
    mesh.update_faces(mesh.unique_faces())
    mesh.remove_unreferenced_vertices()
    if not mesh.is_watertight:
        trimesh.repair.fix_normals(mesh)
    return mesh


def radial_bounds(mesh: trimesh.Trimesh) -> tuple[float, float]:
    radii = np.hypot(mesh.vertices[:, 0], mesh.vertices[:, 1])
    return float(radii.min()), float(radii.max())


def axial_bounds(mesh: trimesh.Trimesh, margin: float = 0.0) -> tuple[float, float]:
    z_min = float(mesh.vertices[:, 2].min()) - margin
    z_max = float(mesh.vertices[:, 2].max()) + margin
    return z_min, z_max
