"""Conservative radial topography queries for future conformal top-skin printing.

This module is geometry/analysis only: it does not authorize variable-Z
extrusion. The printhead's swept envelope, layer ownership and bead packing
must be verified separately before nonplanar paths are emitted.

Mesh Cartesian axes: XY = radial plane, Z = machine axial X.
Machine Y = theta in *radians* inside the slicer (degrees in Marlin G-code).
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import trimesh


TAU = 2.0 * math.pi


@dataclass(frozen=True)
class RadialSurfaceSample:
    theta: float
    axial: float
    radius: float
    face_index: int
    outward_normal_radial: float


class RadialSurfaceSampler:
    """Intersect an outward radial ray with exact mesh triangles.

    Only a unique radial outward-facing *outermost* triangle sheet is
    eligible. A hollow shell may have an inner inward-facing wall; two
    separate outward-facing sheets at different radii are ambiguous and
    must not be connected by a conformal finish pass.
    """

    def __init__(self, mesh: trimesh.Trimesh, *,
                 min_radial_normal: float = 0.03,
                 radius_merge_tolerance_mm: float = 1e-5):
        if not (0 < min_radial_normal <= 1 and math.isfinite(min_radial_normal)):
            raise ValueError("min_radial_normal must be finite within (0, 1]")
        if not (radius_merge_tolerance_mm > 0
                and math.isfinite(radius_merge_tolerance_mm)):
            raise ValueError("radius_merge_tolerance_mm must be finite and positive")
        tri = np.asarray(mesh.triangles, dtype=np.float64)
        if tri.ndim != 3 or tri.shape[1:] != (3, 3) or len(tri) == 0:
            raise ValueError("mesh must contain triangles")
        if not np.isfinite(tri).all():
            raise ValueError("mesh triangles must be finite")
        e0, e1 = tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]
        normal = np.cross(e0, e1)
        norm = np.linalg.norm(normal, axis=1)
        d00 = np.einsum("ij,ij->i", e0, e0)
        d01 = np.einsum("ij,ij->i", e0, e1)
        d11 = np.einsum("ij,ij->i", e1, e1)
        det = d00 * d11 - d01 * d01
        keep = (norm > 1e-12) & (det > 1e-22)
        self.vertices = tri[keep, 0]
        self.e0, self.e1 = e0[keep], e1[keep]
        self.normal = normal[keep] / norm[keep, None]
        self.d00, self.d01, self.d11 = d00[keep], d01[keep], d11[keep]
        self.det = det[keep]
        self.face_indices = np.flatnonzero(keep)
        self.axial_lo = tri[keep, :, 2].min(axis=1)
        self.axial_hi = tri[keep, :, 2].max(axis=1)
        self.min_radial_normal = min_radial_normal
        self.radius_merge_tolerance_mm = radius_merge_tolerance_mm

    def sample(self, theta: float, axial: float) -> RadialSurfaceSample | None:
        """Return an unambiguous, outward-facing outer radial surface hit.

        Reject grazing/degenerate facets, radial self-occlusions and holes.
        Equal-radius hits on shared mesh edges are merged deterministically.
        This only certifies an individual ray, not collision clearance.
        """
        if not (math.isfinite(theta) and math.isfinite(axial)):
            raise ValueError("sample coordinates must be finite")
        c, s = math.cos(theta), math.sin(theta)
        mask = ((axial >= self.axial_lo - 1e-9)
                & (axial <= self.axial_hi + 1e-9))
        if not np.any(mask):
            return None
        indices = np.flatnonzero(mask)
        n = self.normal[indices]
        dot = n[:, 0] * c + n[:, 1] * s
        valid = np.abs(dot) >= self.min_radial_normal
        indices, n, dot = indices[valid], n[valid], dot[valid]
        if not len(indices):
            return None
        v0 = self.vertices[indices]
        rho = (np.einsum("ij,ij->i", n, v0) - n[:, 2] * axial) / dot
        valid = np.isfinite(rho) & (rho > 1e-8)
        indices, dot, rho = indices[valid], dot[valid], rho[valid]
        if not len(indices):
            return None
        v0 = self.vertices[indices]
        p = np.column_stack((rho * c, rho * s, np.full_like(rho, axial))) - v0
        d20 = np.einsum("ij,ij->i", p, self.e0[indices])
        d21 = np.einsum("ij,ij->i", p, self.e1[indices])
        v = (self.d11[indices] * d20 - self.d01[indices] * d21) / self.det[indices]
        w = (self.d00[indices] * d21 - self.d01[indices] * d20) / self.det[indices]
        u = 1.0 - v - w
        valid = ((u >= -1e-8) & (v >= -1e-8) & (w >= -1e-8)
                 & (u <= 1 + 1e-8) & (v <= 1 + 1e-8) & (w <= 1 + 1e-8))
        indices, dot, rho = indices[valid], dot[valid], rho[valid]
        if not len(indices):
            return None
        order = np.argsort(rho)[::-1]
        indices, dot, rho = indices[order], dot[order], rho[order]
        # If the visible outermost sheet faces inward, its outward winding
        # is invalid or the ray has encountered a forbidden geometry.
        if dot[0] <= self.min_radial_normal:
            return None
        # Multiple distinct radially outward sheets along one ray cannot
        # be represented by a single-valued r(theta, X) finish surface.
        for distance, orientation in zip(rho[1:], dot[1:]):
            if distance < rho[0] - self.radius_merge_tolerance_mm:
                if orientation > self.min_radial_normal:
                    return None
            elif orientation < -self.min_radial_normal:
                return None  # Zero-thickness/double-sided surface.
        return RadialSurfaceSample(float(theta), float(axial), float(rho[0]),
                                   int(self.face_indices[indices[0]]),
                                   float(dot[0]))


def segment_surface_length(a: RadialSurfaceSample, b: RadialSurfaceSample) -> float:
    """Short-move 3D cylindrical arc metric (X, R*theta, R), in millimeters."""
    radius = (a.radius + b.radius) * 0.5
    return math.sqrt((b.axial - a.axial)**2
                     + (radius * (b.theta - a.theta))**2
                     + (b.radius - a.radius)**2)


def sample_surface_ribbon(
    sampler: RadialSurfaceSampler,
    centers: Sequence[tuple[float, float]],
    *,
    max_segment_mm: float = 0.5,
    max_radial_chord_error_mm: float = 0.015,
    max_radial_slope: float = 0.3,
    max_total_radial_excursion_mm: float = 1.0,
    max_points: int = 4096,
) -> list[RadialSurfaceSample]:
    """Adaptive candidate only; never emit unapproved extrusion from this.

    Input centers are in *continuous* radians theta and axial mm. Split
    edges to a maximum physical travel length and bounded midpoint radial
    interpolation error; reject holes, multiple surfaces, excessive slope,
    total radial excursion or subdivision budget. These samples alone do
    not prove freedom from sub-resolution holes or toolhead collisions.
    """
    for name, value in (
        ("max_segment_mm", max_segment_mm),
        ("max_radial_chord_error_mm", max_radial_chord_error_mm),
        ("max_radial_slope", max_radial_slope),
        ("max_total_radial_excursion_mm", max_total_radial_excursion_mm),
    ):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not isinstance(max_points, int) or max_points < 3:
        raise ValueError("max_points must be at least 3")
    if len(centers) < 2:
        raise ValueError("a ribbon requires at least two XY positions")
    initial = []
    for theta, axial in centers:
        hit = sampler.sample(float(theta), float(axial))
        if hit is None:
            raise ValueError("no unique outward radial top surface at ribbon point")
        initial.append(hit)
    out = [initial[0]]

    def refine(a: RadialSurfaceSample, b: RadialSurfaceSample, depth: int) -> None:
        if len(out) >= max_points or depth >= 32:
            raise ValueError("nonplanar ribbon exceeded safe subdivision budget")
        length = segment_surface_length(a, b)
        if length <= 1e-10:
            raise ValueError("nonplanar ribbon has duplicate/zero-length centers")
        middle = sampler.sample((a.theta + b.theta) * .5,
                                (a.axial + b.axial) * .5)
        if middle is None:
            raise ValueError("nonplanar ribbon crosses a hole or ambiguous surface")
        radial_error = abs(middle.radius - (a.radius + b.radius) * .5)
        if length > max_segment_mm + 1e-9 or radial_error > max_radial_chord_error_mm:
            refine(a, middle, depth + 1)
            refine(middle, b, depth + 1)
            return
        planar_length = math.hypot(b.axial - a.axial,
                                   (a.radius + b.radius) * .5 * (b.theta - a.theta))
        if planar_length < 1e-10 or abs(b.radius - a.radius) > max_radial_slope * planar_length + 1e-9:
            raise ValueError("nonplanar ribbon slope exceeds limit")
        # The actual midpoint slope can differ substantially from the
        # endpoint slope even when its radial chord error is small.
        for left, right in ((a, middle), (middle, b)):
            transverse = math.hypot(right.axial - left.axial,
                                    (left.radius + right.radius) * .5 *
                                    (right.theta - left.theta))
            if transverse <= 1e-10 or abs(right.radius - left.radius) > max_radial_slope * transverse + 1e-9:
                raise ValueError("nonplanar ribbon midpoint slope exceeds limit")
        out.append(b)

    for a, b in zip(initial, initial[1:]):
        refine(a, b, 0)
    radii = [p.radius for p in out]
    if max(radii) - min(radii) > max_total_radial_excursion_mm + 1e-9:
        raise ValueError("nonplanar ribbon radial excursion exceeds limit")
    return out
