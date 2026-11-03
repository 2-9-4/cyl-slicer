"""Conservative continuous head-section vs mesh clearance bounds.

Analysis ONLY: this proves absence of contact between *specified cylindrical
sections* and the provided triangle surfaces, over continuous native-linear
(theta, axial, radius) movements. It does NOT include the nozzle tip's
intentional contact, heater block shapes that were not declared, mandrel,
previous/future deposited layers, machine travel limits, or extrusion. A
successful section audit is NEVER full-printer collision certification.

At a machine pose (theta, axial, radius), transform mesh CAD Cartesian
coordinates (u,v,X) into the radial head frame:
  tangential = -u*sin(theta) + v*cos(theta)
  axial = X - machine_X
  radial = u*cos(theta) + v*sin(theta) - machine_radius

The head sections are SOLID circular cylinders with their axis radial.
To bound the *entire* motion of every triangle, evaluate its spatial box at
the midpoint of each time interval and enlarge it using exact Lipschitz
bounds from dtheta, dX, dr and the triangle's maximum radial vertex norm.
A separated box certifies that triangle cannot contact the section at any
point in the continuous interval. Ambiguous regions are recursively split
in time and across the triangle; unresolved regions fail closed. Witnessed
mesh points within a head section can provide definite collision evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
import trimesh


@dataclass(frozen=True)
class RadialHeadSection:
    """Cylindrical disk of head material in local radial-above-tip coordinates."""

    min_height_mm: float
    max_height_mm: float
    radius_mm: float

    def validate(self) -> None:
        values = (self.min_height_mm, self.max_height_mm, self.radius_mm)
        if (not all(math.isfinite(v) for v in values)
                or self.min_height_mm < 0
                or self.max_height_mm <= self.min_height_mm
                or self.radius_mm <= 0):
            raise ValueError("head section heights and radius must be finite and positive")


def cone_head_sections(
    *,
    flat_tip_radius_mm: float,
    cone_half_angle_deg: float,
    cone_height_mm: float,
    section_count: int = 12,
    min_height_mm: float = 0.1,
) -> tuple[RadialHeadSection, ...]:
    """Conservatively envelope a conical head *above* the excluded contact tip.

    Each radial band uses the cone's radius at its UPPER height, so the
    union contains the physical cone in that band. This is safe for
    clear/no-clear classification, with possible false positives.
    """
    values = (flat_tip_radius_mm, cone_half_angle_deg, cone_height_mm,
              min_height_mm)
    if (not all(math.isfinite(v) for v in values)
            or flat_tip_radius_mm <= 0
            or not 0 < cone_half_angle_deg < 85
            or cone_height_mm <= min_height_mm
            or min_height_mm <= 0
            or type(section_count) is not int
            or not 1 <= section_count <= 128):
        raise ValueError("invalid cone envelope dimensions or section count")
    tan = math.tan(math.radians(cone_half_angle_deg))
    step = (cone_height_mm - min_height_mm) / section_count
    return tuple(
        RadialHeadSection(min_height_mm + i*step,
                          min_height_mm + (i+1)*step,
                          flat_tip_radius_mm + (min_height_mm + (i+1)*step)*tan)
        for i in range(section_count)
    )


def audit_continuous_head_sweep(
    mesh: trimesh.Trimesh,
    positions: Sequence[tuple[float, float, float]],
    sections: Sequence[RadialHeadSection],
    *,
    radial_margin_mm: float = 0.01,
    lateral_margin_mm: float = 0.01,
    max_time_depth: int = 8,
    max_face_depth: int = 5,
    max_nodes: int = 200_000,
) -> dict:
    """Conservatively check continuous modeled sections against ALL triangles.

    Positions are (theta radians, axial mm, absolute geometric radius mm)
    and consecutive positions interpolate linearly in native machine axes.
    The source mesh must use Cartesian XY radial plane and Z axial. Returns
    a clear boolean ONLY for supplied body sections vs this specific mesh.
    Never convert that into a claim about the full head, print sequence,
    mandrel or actual slicer G-code.
    """
    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError("mesh must be a trimesh.Trimesh")
    if not positions or not sections:
        raise ValueError("at least one pose and one head section required")
    if (type(max_time_depth) is not int or not 0 <= max_time_depth <= 16
            or type(max_face_depth) is not int or not 0 <= max_face_depth <= 10
            or type(max_nodes) is not int or max_nodes < 1):
        raise ValueError("invalid sweep recursion limits or node budget")
    if (not math.isfinite(radial_margin_mm)
            or not math.isfinite(lateral_margin_mm)
            or radial_margin_mm < 0 or lateral_margin_mm < 0):
        raise ValueError("clearance margins must be finite and nonnegative")
    for section in sections:
        if not isinstance(section, RadialHeadSection):
            raise ValueError("head sections must be RadialHeadSection")
        section.validate()
    poses = np.asarray(positions, dtype=np.float64)
    if (poses.ndim != 2 or poses.shape[1] != 3 or not np.isfinite(poses).all()
            or (poses[:,2] <= 0).any()):
        raise ValueError("poses must be finite (theta, axial, positive radius)")
    triangles = np.asarray(mesh.triangles, dtype=np.float64)
    if triangles.ndim != 3 or triangles.shape[1:] != (3,3) or not len(triangles):
        raise ValueError("mesh must contain triangles")
    if not np.isfinite(triangles).all():
        raise ValueError("mesh triangles must be finite")

    # Source triangles may be arbitrary size/triangulation. Do not discard
    # a large facet simply because none of its VERTICES crosses the tool.
    sections = tuple(sections)
    nodes = 0
    first_unknown = None
    unresolved = 0
    tested_faces = 0

    def overlap_bounds(tris: np.ndarray, pose0: np.ndarray, pose1: np.ndarray,
                       t0: float, t1: float) -> np.ndarray:
        """Vectorized conservative triangle/section overlap mask."""
        midpoint = (t0+t1)*0.5
        p = pose0 + (pose1-pose0)*midpoint
        c, s = math.cos(float(p[0])), math.sin(float(p[0]))
        u, v, axial = tris[:,:,0], tris[:,:,1], tris[:,:,2]
        tangential = -u*s + v*c
        local_axial = axial - p[1]
        local_radial = u*c + v*s - p[2]

        angular_speed_bound = np.hypot(u,v).max(axis=1) * abs(pose1[0]-pose0[0])
        halfwidth = (t1-t0)*.5
        tang_exp = angular_speed_bound * halfwidth + lateral_margin_mm
        ax_exp = abs(pose1[1]-pose0[1])*halfwidth + lateral_margin_mm
        rad_exp = (angular_speed_bound + abs(pose1[2]-pose0[2]))*halfwidth + radial_margin_mm
        tlo, thi = tangential.min(axis=1)-tang_exp, tangential.max(axis=1)+tang_exp
        xlo, xhi = local_axial.min(axis=1)-ax_exp, local_axial.max(axis=1)+ax_exp
        rlo, rhi = local_radial.min(axis=1)-rad_exp, local_radial.max(axis=1)+rad_exp

        # Minimum lateral distance between any points in each XY rectangle
        # and the nozzle axis. The triangle lies INSIDE this rectangle.
        dt = np.maximum(np.maximum(tlo,-thi),0)
        dx = np.maximum(np.maximum(xlo,-xhi),0)
        min_lateral = np.hypot(dt,dx)
        return np.logical_or.reduce([
            (min_lateral <= section.radius_mm)
            & (rhi >= section.min_height_mm)
            & (rlo <= section.max_height_mm)
            for section in sections
        ])

    def witnessed(triangle: np.ndarray, p0: np.ndarray, p1: np.ndarray,
                  t0: float, t1: float):
        """A real point strictly in the modeled solid is definitive evidence."""
        mid = (t0+t1)*.5
        p = p0+(p1-p0)*mid
        c,s = math.cos(float(p[0])),math.sin(float(p[0]))
        a,b,d = triangle
        samples = np.stack((a,b,d,(a+b)*.5,(b+d)*.5,(d+a)*.5,
                            (a+b+d)/3),axis=0)
        t = -samples[:,0]*s + samples[:,1]*c
        x = samples[:,2]-p[1]
        r = samples[:,0]*c + samples[:,1]*s - p[2]
        laterals = np.hypot(t,x)
        for i,section in enumerate(sections):
            # Keep witness separate from conservative bounding margins.
            inside = np.where(
                (laterals < section.radius_mm - 1e-8)
                & (r > section.min_height_mm + 1e-8)
                & (r < section.max_height_mm - 1e-8)
            )[0]
            if len(inside):
                hit = int(inside[0])
                return {
                    "section_index": i, "path_fraction": float(mid),
                    "head_local_tangential_mm": float(t[hit]),
                    "head_local_axial_mm": float(x[hit]),
                    "head_local_radial_mm": float(r[hit]),
                }
        return None

    def subdivide(triangle: np.ndarray) -> tuple[np.ndarray, ...]:
        a,b,c=triangle
        ab=(a+b)*.5
        bc=(b+c)*.5
        ca=(c+a)*.5
        return (
            np.stack((a,ab,ca)),
            np.stack((ab,b,bc)),
            np.stack((ca,bc,c)),
            np.stack((ab,bc,ca)),
        )

    class Budget(Exception):
        pass

    class Hit(Exception):
        def __init__(self, detail: dict):
            self.detail = detail

    def verify(triangle: np.ndarray, p0: np.ndarray, p1: np.ndarray,
               t0: float, t1: float, tdepth: int, fdepth: int, face: int,
               segment: int) -> None:
        nonlocal nodes, unresolved, first_unknown
        nodes += 1
        if nodes > max_nodes:
            raise Budget()
        if not overlap_bounds(triangle[None,:,:],p0,p1,t0,t1)[0]:
            return
        hit = witnessed(triangle,p0,p1,t0,t1)
        if hit is not None:
            hit.update({"face_index":face,"segment_index":segment})
            raise Hit(hit)

        # Subdivide the source triangle if its AABB is much larger than
        # the printhead; then subdivide time if native motion is significant.
        # Both dimensions have independent bounded refinement budgets.
        extent = float(np.linalg.norm(np.ptp(triangle,axis=0)))
        long_move = (t1-t0)*float(np.linalg.norm(p1-p0))
        if fdepth < max_face_depth and extent > .5*min(s.radius_mm for s in sections):
            for part in subdivide(triangle):
                verify(part,p0,p1,t0,t1,tdepth,fdepth+1,face,segment)
            return
        if tdepth < max_time_depth and long_move > 1e-8:
            center=(t0+t1)*.5
            verify(triangle,p0,p1,t0,center,tdepth+1,fdepth,face,segment)
            verify(triangle,p0,p1,center,t1,tdepth+1,fdepth,face,segment)
            return
        if fdepth < max_face_depth:
            for part in subdivide(triangle):
                verify(part,p0,p1,t0,t1,tdepth,fdepth+1,face,segment)
            return
        unresolved += 1
        if first_unknown is None:
            first_unknown={"segment_index":segment, "face_index":face,
                           "time_fraction_range":[t0,t1]}

    try:
        for segment in range(max(1,len(poses)-1)):
            p0=poses[segment]
            p1=poses[min(segment+1,len(poses)-1)]
            candidate=overlap_bounds(triangles,p0,p1,0.,1.)
            indices=np.flatnonzero(candidate)
            tested_faces += len(triangles)
            for face in indices:
                verify(triangles[face],p0,p1,0.,1.,0,0,int(face),segment)
    except Hit as exc:
        return {
            "status":"witnessed_interference",
            "modeled_sections_clear":False,
            "full_head_collision_certified":False,
            "mesh_watertight":bool(mesh.is_watertight),
            "witness":exc.detail,
            "checked_nodes":nodes,
            "checked_triangle_roots":tested_faces,
        }
    except Budget:
        return {
            "status":"budget_exceeded",
            "modeled_sections_clear":False,
            "full_head_collision_certified":False,
            "mesh_watertight":bool(mesh.is_watertight),
            "checked_nodes":nodes,
            "checked_triangle_roots":tested_faces,
            "first_uncertain":first_unknown,
        }
    return {
        "status":"unresolved_geometry" if unresolved else "modeled_sections_clear",
        "modeled_sections_clear":unresolved==0,
        "full_head_collision_certified":False,
        "warning":"Only declared cylindrical head sections vs supplied triangle mesh are bounded; nozzle tip/contact, mandrel and printed layers are NOT certified.",
        "mesh_watertight":bool(mesh.is_watertight),
        "checked_nodes":nodes,
        "checked_triangle_roots":tested_faces,
        "unresolved_regions":unresolved,
        "first_uncertain":first_unknown,
    }
