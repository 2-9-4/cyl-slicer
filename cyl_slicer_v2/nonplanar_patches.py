"""Conservative candidate-region discovery for cylindrical conformal finishing.

No G-code emission or printability approval.  A sampled cell can be a
candidate only if its nine surface probes are unambiguous and its *exact
mesh feature edges* cannot cross its theta/axial footprint. This prevents
sub-grid topological holes in an open mesh from being silently bridged.
Closed-mesh interior occlusion still needs swept-volume and layer checks.
"""
from __future__ import annotations

from collections import Counter, deque
import math

import numpy as np
import trimesh
from shapely.geometry import box
from shapely.strtree import STRtree

from .nonplanar_surface import RadialSurfaceSampler, RadialSurfaceSample, TAU


def _radial_feature_boxes(mesh: trimesh.Trimesh, theta_min: float,
                          theta_max: float, *, min_radial_normal: float,
                          sharp_dihedral_deg: float) -> list:
    """Conservative boxes for open, nonmanifold and radial-fold mesh edges.

    Match coincident Cartesian vertices at 1e-7 mm before edge incidence
    counting. The theta image of a *straight Cartesian edge not passing the
    axis* is a monotone angular interval on its shorter arc. Enclosing its
    whole (theta, X) image in a rectangle may reject extra good cells,
    but cannot skip a tiny hole between sample points.
    """
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.intp)
    if not len(faces):
        return []
    if not np.isfinite(vertices).all():
        raise ValueError("mesh vertices must be finite")
    # A tessellated full turn often duplicates its seam vertices at theta
    # 0 and 2pi. Weld by position for topology inspection only; do not
    # alter the source mesh, its face IDs or emitted geometry.
    _, canonical = np.unique(np.round(vertices, decimals=7), axis=0,
                             return_inverse=True)
    welded_faces = canonical[faces]
    edges = np.concatenate((welded_faces[:, [0, 1]],
                            welded_faces[:, [1, 2]],
                            welded_faces[:, [2, 0]]), axis=0)
    owners = np.tile(np.arange(len(faces)), 3)
    edges.sort(axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    unique, starts, counts = np.unique(edges[order], axis=0,
                                       return_index=True, return_counts=True)
    owners = owners[order]
    # Indexed by canonical ID; no risk of choosing a differently positioned
    # duplicate after weld at the specified numerical tolerance.
    coords = np.zeros((int(canonical.max()) + 1, 3), dtype=float)
    coords[canonical] = vertices
    normal = np.asarray(mesh.face_normals, dtype=float)
    cos_sharp = math.cos(math.radians(sharp_dihedral_deg))
    features = []
    for pair, start, count in zip(unique, starts, counts):
        a, b = coords[pair]
        mid = (a + b) * .5
        radius = math.hypot(mid[0], mid[1])
        selected = count != 2
        if not selected:
            adjacent = normal[owners[start:start+2]]
            if radius <= 1e-8:
                selected = True
            else:
                radial = mid[:2] / radius
                rdot = adjacent[:, :2] @ radial
                selected = (np.any(rdot <= min_radial_normal)
                            or float(np.dot(adjacent[0], adjacent[1])) < cos_sharp)
        if not selected:
            continue

        z0, z1 = sorted((float(a[2]), float(b[2])))
        # Passing through the build axis means angular projection cannot
        # be bounded by a short interval: reject every angular column
        # intersecting its axial envelope.
        v = b[:2] - a[:2]
        vv = float(np.dot(v, v))
        if vv < 1e-16:
            near_axis = radius < 1e-8
        else:
            t = max(0.0, min(1.0, -float(np.dot(a[:2], v)) / vv))
            near_axis = float(np.linalg.norm(a[:2] + t * v)) < 1e-8
        if near_axis:
            features.append(box(theta_min, z0, theta_max, z1))
            continue

        ta, tb = math.atan2(a[1], a[0]), math.atan2(b[1], b[0])
        delta = (tb - ta + math.pi) % TAU - math.pi
        # Antipodal radial directions need more than an endpoint test.
        if abs(abs(delta) - math.pi) < 1e-8:
            features.append(box(theta_min, z0, theta_max, z1))
            continue
        lo, hi = sorted((ta, ta + delta))
        first_shift = math.floor((theta_min - hi) / TAU)
        last_shift = math.ceil((theta_max - lo) / TAU)
        for shift in range(first_shift, last_shift + 1):
            left, right = lo + shift * TAU, hi + shift * TAU
            if left <= theta_max + 1e-12 and right >= theta_min - 1e-12:
                features.append(box(max(left, theta_min), z0,
                                    min(right, theta_max), z1))
    return features


def _grid_samples(sampler: RadialSurfaceSampler, theta_min: float,
                  theta_max: float, axial_min: float, axial_max: float,
                  theta_cells: int, axial_cells: int
                  ) -> list[list[RadialSurfaceSample | None]]:
    """Cache each corner/edge/center probe once for all neighboring cells."""
    step_t = (theta_max - theta_min) / (2 * theta_cells)
    step_x = (axial_max - axial_min) / (2 * axial_cells)
    return [
        [sampler.sample(theta_min + i * step_t, axial_min + j * step_x)
         for i in range(2 * theta_cells + 1)]
        for j in range(2 * axial_cells + 1)
    ]


def scan_surface_patches(
    mesh: trimesh.Trimesh,
    *,
    theta_min: float = 0.0,
    theta_max: float = TAU,
    axial_min: float,
    axial_max: float,
    theta_cells: int = 72,
    axial_cells: int = 32,
    max_slope: float = 0.3,
    max_chord_error_mm: float = 0.03,
    max_cell_excursion_mm: float = 0.5,
    sharp_dihedral_deg: float = 55.0,
    max_cells: int = 12000,
    min_patch_cells: int = 1,
) -> dict:
    """Discover connected *candidate* surface patches in a bounded grid.

    Cell vertices, four edge midpoints and center must each have one unique
    outward radial surface. Reject mesh boundary/crease edge crossings
    independently of the point grid. Reject local excessive slope, bilinear
    height discrepancy or large radial excursion. Wrap patch connectivity
    across one full revolution, never duplicate seam deposition.

    A candidate is not a continuous swept-nozzle clearance certificate.
    """
    values = (theta_min, theta_max, axial_min, axial_max, max_slope,
              max_chord_error_mm, max_cell_excursion_mm, sharp_dihedral_deg)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("scan limits and thresholds must be finite")
    if (not 0 < theta_max - theta_min <= TAU + 1e-9
            or axial_min >= axial_max):
        raise ValueError("scan needs a positive axial interval and theta span <= 2pi")
    if (max_slope <= 0 or max_chord_error_mm <= 0 or
            max_cell_excursion_mm <= 0 or
            not 0 < sharp_dihedral_deg <= 180):
        raise ValueError("scan safety thresholds must be positive and valid")
    if (type(theta_cells) is not int or type(axial_cells) is not int
            or type(max_cells) is not int or type(min_patch_cells) is not int
            or theta_cells < 1 or axial_cells < 1
            or min_patch_cells < 1 or max_cells < 1
            or theta_cells * axial_cells > max_cells):
        raise ValueError("scan grid dimensions exceed validated cell budget")

    sampler = RadialSurfaceSampler(mesh)
    hits = _grid_samples(sampler, theta_min, theta_max, axial_min, axial_max,
                         theta_cells, axial_cells)
    feature_boxes = _radial_feature_boxes(
        mesh, theta_min, theta_max,
        min_radial_normal=sampler.min_radial_normal,
        sharp_dihedral_deg=sharp_dihedral_deg,
    )
    tree = STRtree(feature_boxes) if feature_boxes else None
    dt = (theta_max - theta_min) / theta_cells
    dx = (axial_max - axial_min) / axial_cells
    mask = np.zeros((axial_cells, theta_cells), dtype=bool)
    radii = np.full_like(mask, np.nan, dtype=float)
    min_normals = np.full_like(mask, np.nan, dtype=float)
    rejections = Counter()
    full_turn = abs(theta_max - theta_min - TAU) <= 1e-8

    for j in range(axial_cells):
        for i in range(theta_cells):
            cell = [hits[2*j + dy][2*i + di] for dy in range(3)
                    for di in range(3)]
            if any(v is None for v in cell):
                rejections["missing_or_ambiguous_surface"] += 1
                continue
            if tree is not None:
                footprint = box(theta_min + i*dt, axial_min + j*dx,
                                theta_min + (i+1)*dt, axial_min + (j+1)*dx)
                if len(tree.query(footprint)):
                    rejections["mesh_boundary_or_radial_crease"] += 1
                    continue
            r = np.asarray([point.radius for point in cell])
            if float(r.max() - r.min()) > max_cell_excursion_mm + 1e-9:
                rejections["radial_excursion"] += 1
                continue
            rgrid = r.reshape(3, 3)
            # Each center or edge midpoint must agree with adjacent exact
            # triangle-intersection endpoints within requested chord bound.
            discrepancies = [
                abs(rgrid[0,1] - (rgrid[0,0]+rgrid[0,2])*.5),
                abs(rgrid[2,1] - (rgrid[2,0]+rgrid[2,2])*.5),
                abs(rgrid[1,0] - (rgrid[0,0]+rgrid[2,0])*.5),
                abs(rgrid[1,2] - (rgrid[0,2]+rgrid[2,2])*.5),
                abs(rgrid[1,1] - .25*(rgrid[0,0]+rgrid[0,2]
                                      +rgrid[2,0]+rgrid[2,2])),
            ]
            if max(discrepancies) > max_chord_error_mm + 1e-9:
                rejections["radial_chord_error"] += 1
                continue

            failed_slope = False
            for row in range(3):
                for col in range(2):
                    a, b = rgrid[row,col:col+2]
                    if abs(b-a) > max_slope * (a+b)*.5*dt*.5 + 1e-9:
                        failed_slope = True
            for row in range(2):
                for col in range(3):
                    a, b = rgrid[row:row+2,col]
                    if abs(b-a) > max_slope*dx*.5 + 1e-9:
                        failed_slope = True
            if failed_slope:
                rejections["radial_slope"] += 1
                continue
            mask[j,i] = True
            radii[j,i] = float(np.mean(r))
            min_normals[j,i] = min(point.outward_normal_radial for point in cell)

    visited = np.zeros_like(mask, dtype=bool)
    patches = []
    for j in range(axial_cells):
        for i in range(theta_cells):
            if not mask[j,i] or visited[j,i]:
                continue
            queue = deque([(i,j)])
            visited[j,i] = True
            members = []
            while queue:
                col,row = queue.popleft()
                members.append((col,row))
                for ni,nj in ((col-1,row),(col+1,row),
                              (col,row-1),(col,row+1)):
                    if full_turn:
                        ni %= theta_cells
                    if (0 <= ni < theta_cells and 0 <= nj < axial_cells
                            and mask[nj,ni] and not visited[nj,ni]):
                        visited[nj,ni] = True
                        queue.append((ni,nj))
            if len(members) < min_patch_cells:
                rejections["below_min_patch_size"] += len(members)
                for col,row in members:
                    mask[row,col] = False
                continue
            angles = [p[0] for p in members]
            axials = [p[1] for p in members]
            source_radius = [float(radii[row,col]) for col,row in members]
            patches.append({
                "id": len(patches)+1,
                "cell_count": len(members),
                "theta_range_deg": [
                    math.degrees(theta_min + min(angles)*dt),
                    math.degrees(theta_min + (max(angles)+1)*dt)],
                "axial_range_mm": [
                    axial_min + min(axials)*dx,
                    axial_min + (max(axials)+1)*dx],
                "crosses_periodic_seam": (full_turn and 0 in angles
                                           and theta_cells-1 in angles
                                           and theta_cells > 1),
                "min_radius_mm": min(source_radius),
                "max_radius_mm": max(source_radius),
                "min_radial_normal": float(min(min_normals[row,col]
                                                for col,row in members)),
                "approx_unwrapped_area_mm2": float(sum(
                    dt*dx*radii[row,col] for col,row in members)),
            })
    patches.sort(key=lambda p: (-p["cell_count"], p["id"]))
    return {
        "status": "geometry_candidates_only",
        "nonplanar_extrusion_enabled": False,
        "collision_certified": False,
        "warnings": [
            "Cell grid and feature-edge tests are conservative but not a complete continuous-surface proof.",
            "No nozzle swept-volume, shell ownership, bead support, or machine-planner validation.",
        ],
        "theta_range_deg": [math.degrees(theta_min), math.degrees(theta_max)],
        "axial_range_mm": [axial_min, axial_max],
        "theta_cells": theta_cells,
        "axial_cells": axial_cells,
        "accepted_cells": int(mask.sum()),
        "rejected_cells": int(mask.size-mask.sum()),
        "rejection_counts": dict(sorted(rejections.items())),
        "candidate_mask": ["".join("1" if flag else "0" for flag in row)
                           for row in mask],
        "patches": patches,
    }
