"""Read-only exclusive top-skin replacement ownership audit.

This audit does NOT remove conventional top-skin paths or authorize
nonplanar extrusion. Its accepted output identifies precisely one
previously generated MONOTONIC TOP_SKIN toolpath whose complete 2D
centerline is matched by a proposed conformal replacement. It tests
exterior top-layer exposure, solid backing, radial slab eligibility,
reserved vs protected role centerlines, and approximate volume
conservation before suggesting future ownership transfer.

The deliberately strict "whole toolpath" contract avoids partial
clipping, which would accidentally drop or double-print portions of
an existing SOLID_INFILL, perimeter, bridge or support path.

In future, a live renderer may consume a separate validated ownership
transfer and must remove that toolpath from BOTH paths_by_role and
island_paths; never mutate one alias while leaving the other intact.
"""
from __future__ import annotations

import math
from typing import Sequence

from shapely import affinity
from shapely.geometry import GeometryCollection, LineString
from shapely.geometry.base import BaseGeometry

from .config import SlicerConfig
from .coordinates import cylindrical_distance
from .models import LayerSlice, LayerToolpaths, Toolpath
from .nonplanar_motion import evaluate_conformal_move
from .periodic import periodic_copies


TAU = 2 * math.pi


def _metric_geometry(geom: BaseGeometry, radius: float) -> BaseGeometry:
    return affinity.scale(geom, xfact=radius, yfact=1.0, origin=(0, 0))


def _metric_line(points: Sequence[tuple[float, float]], radius: float,
                 *, near_theta: float | None = None) -> LineString:
    """Map an unwrapped constant-period path into physical theta*R, axial."""
    if len(points) < 2:
        raise ValueError("ownership paths need at least two points")
    if not all(math.isfinite(t) and math.isfinite(x) for t,x in points):
        raise ValueError("ownership centerline coordinates must be finite")
    shift=0.0
    if near_theta is not None:
        shift=round((near_theta - points[0][0]) / TAU) * TAU
    unwrapped=[(points[0][0]+shift,points[0][1])]
    prev=points[0][0]
    for theta,axial in points[1:]:
        delta=theta-prev
        if abs(delta)>math.pi+1e-9:
            raise ValueError("ownership cannot resolve a single segment exceeding half a turn")
        unwrapped.append((unwrapped[-1][0]+delta,axial))
        prev=theta
    return LineString([(theta*radius,axial) for theta,axial in unwrapped])


def _physical_volume(path: Toolpath, config: SlicerConfig) -> float:
    widths=path.segment_widths
    if widths is not None and len(widths)!=len(path.points)-1:
        raise ValueError("legacy variable-width path is internally inconsistent")
    volume=0.0
    for i,(a,b) in enumerate(zip(path.points,path.points[1:])):
        width=widths[i] if widths is not None else config.nozzle_width
        height=config.layer_height
        if width<height-1e-9 or width>config.nozzle_width+1e-9:
            raise ValueError("unprintable legacy skin bead width")
        area=(height*(width-height)+math.pi*height**2/4
              if config.bead_model=="rounded" else width*height)
        volume+=area*cylindrical_distance(a,b,path.radius)
    return volume


def audit_conformal_top_skin_ownership(
    layers: Sequence[LayerSlice],
    generated: Sequence[LayerToolpaths],
    candidate: Toolpath,
    config: SlicerConfig,
    *,
    target_layer_index: int,
    match_tolerance_mm: float = 0.025,
    protected_role_clearance_mm: float = 0.025,
    radial_face_tolerance_mm: float = 0.08,
    max_volume_difference_fraction: float = 0.05,
    minimum_solid_backing_layers: int = 1,
    slices_cover_full_radial_range: bool = False,
) -> dict:
    """Audit whether a candidate may *eventually* replace one whole top-skin path.

    A positive result is an *exclusive ownership reservation candidate*,
    not an approval to mutate slices or emit G-code. Only monotonic top
    skins are distinguishable from underlying bottom solid infill at this
    stage; avoid removing a shared SOLID_INFILL path accidentally.
    """
    result={
        "status":"rejected",
        "reason":None,
        "candidate_only":True,
        "nonplanar_extrusion_enabled":False,
        "conventional_toolpaths_mutated":False,
        "full_head_collision_certified":False,
        "target_layer_index":target_layer_index,
        "reserved_legacy_paths":[],
    }
    def reject(reason: str) -> dict:
        result["reason"]=reason
        return result

    config.validate()
    for v in (match_tolerance_mm,protected_role_clearance_mm,
              radial_face_tolerance_mm,max_volume_difference_fraction):
        if not math.isfinite(v) or v<0:
            raise ValueError("ownership tolerances must be finite and nonnegative")
    if (max_volume_difference_fraction>=1
            or type(minimum_solid_backing_layers) is not int
            or minimum_solid_backing_layers<1):
        raise ValueError("ownership volume and backing limits are invalid")
    if type(target_layer_index) is not int:
        raise ValueError("target layer index must be an integer")
    if not config.monotonic_top_skin:
        return reject("requires_separate_monotonic_TOP_SKIN_role")
    if candidate.role!="TOP_SKIN" or candidate.closed:
        return reject("conformal_candidate_must_be_open_TOP_SKIN")
    if candidate.segment_widths is not None:
        return reject("variable_planar_widths_cannot_mix_with_conformal_bead_areas")
    if (len(candidate.points)<2
            or candidate.point_radii is None or candidate.segment_areas_mm2 is None
            or len(candidate.point_radii)!=len(candidate.points)
            or len(candidate.segment_areas_mm2)!=len(candidate.points)-1):
        return reject("missing_pointwise_radii_or_segment_bead_areas")
    if (any(not math.isfinite(r) or r<config.bed_radius
            for r in candidate.point_radii)
            or any(not math.isfinite(a) or a<=0 for a in candidate.segment_areas_mm2)):
        return reject("invalid_conformal_radii_or_bead_areas")
    if not layers or not generated:
        return reject("complete_layer_and_toolpath_sets_are_required")
    lookup={layer.index:layer for layer in layers}
    generated_by_index={layer.layer_index:layer for layer in generated}
    if (len(lookup)!=len(layers) or len(generated_by_index)!=len(generated)
            or target_layer_index not in lookup
            or target_layer_index not in generated_by_index):
        return reject("missing_or_duplicate_target_layer")
    target=lookup[target_layer_index]
    target_paths=generated_by_index[target_layer_index]
    if target.is_raft or target.occupied.is_empty or target.top_surfaces.is_empty:
        return reject("target_has_no_model_top_surface")

    next_layer=lookup.get(target_layer_index+1)
    if next_layer is None and not slices_cover_full_radial_range:
        return reject("outer_layer_occupancy_is_unknown_for_truncated_job")
    previous_layers=[]
    for layer_index in range(target_layer_index-minimum_solid_backing_layers,
                             target_layer_index):
        below=lookup.get(layer_index)
        if below is None or below.occupied.is_empty:
            return reject("required_solid_backing_layer_is_missing")
        previous_layers.append(below)

    thickness=(target.deposited_thickness if target.deposited_thickness is not None
               else config.layer_height)
    slab_config=config.for_deposited_thickness(thickness)
    top_face=target.radius+thickness*.5
    if any(abs(r-top_face)>radial_face_tolerance_mm+1e-9
           for r in candidate.point_radii):
        return reject("candidate_radius_outside_single_deposited_top_slab")
    if not all(math.isfinite(t) and math.isfinite(x)
               for t,x in candidate.points):
        return reject("candidate_positions_are_nonfinite")

    reference_radius=target.radius
    line=_metric_line(candidate.points,reference_radius)
    if line.length<1e-7 or not line.is_simple:
        return reject("candidate_centerline_is_zero_length_or_self_intersecting")
    # A purely projected XY bead corridor; later full 3D swept-volume
    # ownership must confirm every radius and bead section along this line.
    corridor=line.buffer(config.nozzle_width*.5,cap_style=2,join_style=2)
    top=periodic_copies(target.top_surfaces)
    exposed=periodic_copies(target.occupied)
    if next_layer is not None and not next_layer.occupied.is_empty:
        exposed=exposed.difference(periodic_copies(next_layer.occupied))
    permitted=_metric_geometry(top.intersection(exposed),reference_radius)
    if corridor.difference(permitted).area>1e-7:
        return reject("candidate_bead_corridor_is_not_fully_exposed_top_skin")
    if not target.bottom_surfaces.is_empty:
        bottom=_metric_geometry(periodic_copies(target.bottom_surfaces),reference_radius)
        if corridor.intersects(bottom):
            return reject("candidate_overlaps_undifferentiated_bottom_solid")
    for previous in previous_layers:
        support=_metric_geometry(periodic_copies(previous.occupied),
                                 reference_radius)
        if corridor.difference(support).area>1e-7:
            return reject("candidate_lacks_full_solid_backing")

    matches=[]
    legacy_paths=target_paths.paths_by_role.get("TOP_SKIN",[])
    near_theta=candidate.points[0][0]
    for index,old in enumerate(legacy_paths):
        old_line=_metric_line(old.points,reference_radius,near_theta=near_theta)
        distance=old_line.distance(line)
        if distance<=match_tolerance_mm+1e-9:
            if old_line.hausdorff_distance(line)>match_tolerance_mm+1e-9:
                return reject("partial_legacy_TOP_SKIN_overlap")
            matches.append((index,old))
    if len(matches)!=1:
        return reject("must_match_exactly_one_whole_legacy_TOP_SKIN_path")

    for role,paths in target_paths.paths_by_role.items():
        if role=="TOP_SKIN":
            continue
        for protected in paths:
            protected_line=_metric_line(protected.points,reference_radius,
                                        near_theta=near_theta)
            if protected_line.distance(line)<=protected_role_clearance_mm+1e-9:
                return reject("conformal_candidate_conflicts_with_protected_"+role)

    legacy_index,old=matches[0]
    old_volume=_physical_volume(old,slab_config)
    new_volume=0.0
    for i,((theta,axial),(next_theta,next_axial)) in enumerate(
        zip(candidate.points,candidate.points[1:])
    ):
        move=evaluate_conformal_move(
            (theta,axial,candidate.point_radii[i]),
            (next_theta,next_axial,candidate.point_radii[i+1]),
            physical_print_feed_mm_per_min=config.feed_mm_per_min,
            bead_area_mm2=candidate.segment_areas_mm2[i],
            filament_diameter_mm=config.filament_diameter,
        )
        new_volume+=move.length_mm*candidate.segment_areas_mm2[i]
    if old_volume<=0 or (abs(new_volume-old_volume)
                         >old_volume*max_volume_difference_fraction+1e-9):
        result.update({
            "original_volume_mm3":old_volume,
            "conformal_volume_mm3":new_volume,
        })
        return reject("conformal_volume_does_not_replace_legacy_volume_within_budget")

    result.update({
        "status":"exclusive_legacy_TOP_SKIN_replacement_candidate",
        "reason":None,
        "reserved_legacy_paths":[{
            "layer_index":target_layer_index,
            "role":"TOP_SKIN",
            "role_path_index":legacy_index,
        }],
        "original_volume_mm3":old_volume,
        "conformal_volume_mm3":new_volume,
        "relative_volume_difference":abs(new_volume-old_volume)/old_volume,
        "backing_layers_checked":len(previous_layers),
        "unresolved_requirements":[
            "3D bead-volume overlap and top-shell contact geometry",
            "complete measured printhead and mandrel collision clearance",
            "machine-safe travels, acceleration, retracts and rotary rebasing",
            "atomic replacement in paths_by_role AND island_paths",
            "actual material and machine print coupon validation",
        ],
    })
    return result
