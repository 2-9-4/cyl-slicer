from __future__ import annotations

import math

import numpy as np
import trimesh
from shapely.geometry import Point


def _positive_z(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    mesh = mesh.copy()
    mesh.apply_translation([0.0, 0.0, -float(mesh.bounds[0, 2])])
    return mesh


def solid_box(radius: float = 20.0, size: tuple[float, float, float] = (5.0, 5.0, 10.0)) -> trimesh.Trimesh:
    mesh = trimesh.creation.box(extents=size)
    mesh.apply_translation([radius, 0.0, size[2] / 2.0])
    return mesh


def offset_box_pair(radius: float = 20.0) -> trimesh.Trimesh:
    a = solid_box(radius=radius, size=(4.0, 4.0, 8.0))
    b = trimesh.creation.box(extents=(4.0, 4.0, 8.0))
    angle = math.radians(130.0)
    b.apply_translation([radius * math.cos(angle), radius * math.sin(angle), 4.0])
    return trimesh.util.concatenate([a, b])


def tube(inner_radius: float = 10.0, outer_radius: float = 12.0, height: float = 12.0) -> trimesh.Trimesh:
    return _positive_z(trimesh.creation.annulus(inner_radius, outer_radius, height=height, sections=72))


def thin_ring(inner_radius: float = 10.0, outer_radius: float = 10.8, height: float = 8.0) -> trimesh.Trimesh:
    return _positive_z(trimesh.creation.annulus(inner_radius, outer_radius, height=height, sections=72))


def helix_like_ribbon(radius: float = 18.0, height: float = 16.0, turns: float = 1.4) -> trimesh.Trimesh:
    ts = np.linspace(-0.35, turns * 2.0 * math.pi + 0.35, 80)
    path = np.column_stack(
        [
            radius * np.cos(ts),
            radius * np.sin(ts),
            (ts - ts.min()) / (ts.max() - ts.min()) * height,
        ]
    )
    profile = Point(0, 0).buffer(0.9, quad_segs=8)
    return trimesh.creation.sweep_polygon(profile, path)


def overhang_step() -> trimesh.Trimesh:
    post = trimesh.creation.cylinder(radius=6.0, height=12.0, sections=48)
    post = _positive_z(post)
    shelf = trimesh.creation.annulus(10.0, 14.0, height=4.0, sections=72)
    shelf.apply_translation([0.0, 0.0, 10.0])
    return trimesh.util.concatenate([post, shelf])


FIXTURES = {
    "solid_box": solid_box,
    "offset_box_pair": offset_box_pair,
    "tube": tube,
    "thin_ring": thin_ring,
    "helix_like_ribbon": helix_like_ribbon,
    "overhang_step": overhang_step,
}


def make_fixture(name: str) -> trimesh.Trimesh:
    try:
        return FIXTURES[name]()
    except KeyError as exc:
        known = ", ".join(sorted(FIXTURES))
        raise ValueError(f"unknown fixture {name!r}; expected one of: {known}") from exc
