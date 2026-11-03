"""Topology-clipped axial conformal finish RIBBON candidates (no G-code).

Construct optional, connected surface-following axial rasters from the
conservative grid of scan_surface_patches(). Each ribbon must stay inside a
contiguous run of accepted scan cells in the same angular column; scan
rejections, including sub-grid mesh boundaries, split the raster.

These are centerline *proposals*, not a printable top skin: actual
surface-metric bead footprint/spacing, perimeter ownership, conformal
shell replacement, support, collision and extrusion remain unresolved.
"""
from __future__ import annotations

from collections import Counter
import math

import trimesh

from .mesh_io import radial_bounds
from .nonplanar_surface import (
    RadialSurfaceSampler, sample_surface_ribbon, segment_surface_length,
)


def propose_axial_conformal_rasters(
    mesh: trimesh.Trimesh,
    scan: dict,
    *,
    nominal_surface_pitch_mm: float = 0.42,
    allowed_surface_slope: float = 0.3,
    edge_inset_mm: float = 0.21,
    min_path_length_mm: float = 0.5,
    max_segment_mm: float = 0.5,
    max_radial_chord_error_mm: float = 0.015,
    max_radial_excursion_mm: float = 1.0,
    max_paths: int = 2048,
    max_total_points: int = 100_000,
    include_points: bool = False,
) -> dict:
    """Clip periodic theta-registered axial rasters to candidate scan cells.

    The theta track lattice is globally uniform, not regenerated per
    angular cell/island. A global outer mesh radius is used to pick a
    conservative projected angular pitch. The sqrt(1+slope²) factor
    offsets the *maximum configured slope*, but the actual geodesic
    surface spacing/normal bead coverage is not yet certified. No
    G-code/retraction/travel/skin material is created or removed.
    """
    values=(nominal_surface_pitch_mm,allowed_surface_slope,edge_inset_mm,
            min_path_length_mm,max_segment_mm,max_radial_chord_error_mm,
            max_radial_excursion_mm)
    if (not all(math.isfinite(v) for v in values)
            or nominal_surface_pitch_mm <= 0
            or allowed_surface_slope <= 0
            or edge_inset_mm < 0
            or min_path_length_mm <= 0
            or max_segment_mm <= 0
            or max_radial_chord_error_mm <= 0
            or max_radial_excursion_mm <= 0):
        raise ValueError("raster spacing, safety and path thresholds must be finite and valid")
    if (type(max_paths) is not int or max_paths < 1
            or type(max_total_points) is not int or max_total_points < 2):
        raise ValueError("raster path and sample budgets must be finite and positive")
    if scan.get("status") != "geometry_candidates_only":
        raise ValueError("raster preview requires automatic candidate patch scan")
    if scan.get("nonplanar_extrusion_enabled") is not False:
        raise ValueError("raster preview refuses non-analysis scan inputs")
    mask=scan.get("candidate_mask")
    theta_range=scan.get("theta_range_deg")
    axial_range=scan.get("axial_range_mm")
    columns=scan.get("theta_cells")
    rows=scan.get("axial_cells")
    if (type(columns) is not int or type(rows) is not int
            or columns < 1 or rows < 1
            or not isinstance(mask,list) or len(mask)!=rows
            or any(type(row) is not str or len(row)!=columns
                   or any(c not in "01" for c in row) for row in mask)
            or not isinstance(theta_range,list) or len(theta_range)!=2
            or not isinstance(axial_range,list) or len(axial_range)!=2
            or not all(math.isfinite(v) for v in theta_range+axial_range)
            or theta_range[1] <= theta_range[0]
            or theta_range[1]-theta_range[0] > 360.0000001
            or axial_range[1] <= axial_range[0]):
        raise ValueError("candidate scan mask/metadata is invalid")

    theta_min,theta_max=map(math.radians,theta_range)
    axial_min,axial_max=axial_range
    delta_theta=(theta_max-theta_min)/columns
    delta_axial=(axial_max-axial_min)/rows
    outer_radius=radial_bounds(mesh)[1]
    if not (math.isfinite(outer_radius) and outer_radius>0):
        raise ValueError("mesh has no positive finite radial extent")
    projected_pitch=nominal_surface_pitch_mm / math.sqrt(1+allowed_surface_slope**2)
    track_count=max(1,math.ceil(
        (theta_max-theta_min)*outer_radius / projected_pitch - 1e-10))
    if track_count*rows>max_paths*max(rows,1)*64:
        raise ValueError("angular raster lattice exceeds finite path budget")
    angular_step=(theta_max-theta_min)/track_count

    sampler=RadialSurfaceSampler(mesh)
    candidate=[]
    skip=Counter()
    num_points=0
    for k in range(track_count):
        theta=theta_min+(k+0.5)*angular_step
        column=min(columns-1,max(0,int((theta-theta_min)/delta_theta)))
        active=[mask[row][column]=="1" for row in range(rows)]
        j=0
        while j<rows:
            if not active[j]:
                j+=1
                continue
            start=j
            while j<rows and active[j]:
                j+=1
            stop=j # exclusive
            a=axial_min+start*delta_axial+edge_inset_mm
            b=axial_min+stop*delta_axial-edge_inset_mm
            if b-a < min_path_length_mm:
                skip["short_candidate"]+=1
                continue
            if len(candidate)>=max_paths:
                raise ValueError("raster preview exceeds maximum path budget")
            try:
                line=sample_surface_ribbon(
                    sampler,[(theta,a),(theta,b)],
                    max_segment_mm=max_segment_mm,
                    max_radial_chord_error_mm=max_radial_chord_error_mm,
                    max_radial_slope=allowed_surface_slope,
                    max_total_radial_excursion_mm=max_radial_excursion_mm,
                    max_points=min(4096,max_total_points),
                )
            except ValueError:
                skip["failed_surface_ribbon"]+=1
                continue
            num_points+=len(line)
            if num_points>max_total_points:
                raise ValueError("raster preview exceeds maximum point budget")
            info={
                "id":len(candidate)+1,
                "angular_track_index":k,
                "theta_deg":math.degrees(theta),
                "scan_column":column,
                "axial_row_span":[start,stop],
                "axial_range_mm":[a,b],
                "point_count":len(line),
                "physical_length_mm":sum(
                    segment_surface_length(p,q) for p,q in zip(line,line[1:])),
                "min_radius_mm":min(point.radius for point in line),
                "max_radius_mm":max(point.radius for point in line),
            }
            if include_points:
                info["points"]=[
                    {"theta_deg":math.degrees(p.theta),
                     "axial_mm":p.axial,
                     "radius_mm":p.radius}
                    for p in line
                ]
            candidate.append(info)
    return {
        "status":"candidate_axial_rasters_only",
        "nonplanar_extrusion_enabled":False,
        "actual_surface_coverage_certified":False,
        "full_head_collision_certified":False,
        "warning":"No bead/skin ownership, swept collision, support, surface-metric coverage or G-code validation.",
        "global_angular_track_count":track_count,
        "nominal_surface_pitch_mm":nominal_surface_pitch_mm,
        "projected_angular_pitch_at_outer_radius_mm":outer_radius*angular_step,
        "candidate_count":len(candidate),
        "total_centerline_length_mm":sum(p["physical_length_mm"] for p in candidate),
        "point_count":num_points,
        "skip_reasons":dict(sorted(skip.items())),
        "paths":candidate,
    }
