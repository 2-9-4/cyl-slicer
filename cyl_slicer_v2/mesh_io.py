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
    """True radial extrema of all mesh triangles about the global Z axis.

    Max distance is attained at a vertex (convex norm). Min distance may be
    attained on a projected edge, or be zero inside a triangle's XY projection.
    Looking only at vertices can skip the first occupied radial slabs.
    """
    radii = np.hypot(mesh.vertices[:, 0], mesh.vertices[:, 1])
    vertex_min, maximum = float(radii.min()), float(radii.max())
    if len(mesh.faces) == 0:
        return vertex_min, maximum

    triangles = mesh.vertices[mesh.faces, :2]
    edges = np.roll(triangles, -1, axis=1) - triangles
    lengths_sq = np.sum(edges * edges, axis=2)
    # Nearest points to origin on all three closed triangle edges. This also
    # handles zero-length edges and triangles degenerate in XY projection.
    fractions = np.divide(
        -np.sum(triangles * edges, axis=2),
        lengths_sq,
        out=np.zeros_like(lengths_sq),
        where=lengths_sq > 0.0,
    )
    closest = triangles + np.clip(fractions, 0.0, 1.0)[..., None] * edges
    face_min_sq = np.min(np.sum(closest * closest, axis=2), axis=1)

    # A nondegenerate projected triangle containing the axis reaches radius
    # zero in its interior, even if none of its edges comes close to it.
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    cross_ab = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
    cross_bc = b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0]
    cross_ca = c[:, 0] * a[:, 1] - c[:, 1] * a[:, 0]
    same_sign = (
        ((cross_ab >= 0) & (cross_bc >= 0) & (cross_ca >= 0))
        | ((cross_ab <= 0) & (cross_bc <= 0) & (cross_ca <= 0))
    )
    twice_area = (
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
        - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    )
    side_scale_sq = np.max(lengths_sq, axis=1)
    nondegenerate = np.abs(twice_area) > 32.0 * np.finfo(float).eps * side_scale_sq
    face_min_sq[same_sign & nondegenerate] = 0.0
    return float(np.sqrt(face_min_sq.min())), maximum


def axial_bounds(mesh: trimesh.Trimesh, margin: float = 0.0) -> tuple[float, float]:
    z_min = float(mesh.vertices[:, 2].min()) - margin
    z_max = float(mesh.vertices[:, 2].max()) + margin
    return z_min, z_max
