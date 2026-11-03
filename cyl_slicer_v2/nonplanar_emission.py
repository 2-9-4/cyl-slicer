"""NON-PRINTABLE conformal XYZ G-code *motion* previews.

This module exposes the missing native X/Y(degrees)/Z(mm)/E/F conversion
for radius-varying candidate Toolpaths, but DOES NOT connect the result
to the production emitter. It does not generate homing, heatup, retract,
fan, start/end travel, printer coordinate setup or a complete job.

G-code-like commands are for inspection and regression tests only. A
true print exporter must also OWN and replace the conventional top skin
instead of overprinting it; clear the entire physical printhead,
mandrel and deposited layers; and validate backing support, collision,
travel/retraction/acceleration and firmware precision.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .config import SlicerConfig
from .models import Toolpath
from .nonplanar_plan import plan_conformal_motion
from .nonplanar_surface import RadialSurfaceSample


_NATIVE_Y_BUDGET_DEGREES = 3599.0


@dataclass(frozen=True)
class ConformalMotionPreview:
    status: str
    nonplanar_extrusion_enabled: bool
    commands: tuple[str, ...]
    segment_count: int
    source_segment_count: int
    exact_physical_length_mm: float
    ideal_filament_length_mm: float
    emitted_filament_length_mm: float
    maximum_local_e_error_fraction: float
    min_radius_mm: float
    max_radius_mm: float


def preview_conformal_toolpath_xyz(
    path: Toolpath,
    config: SlicerConfig,
    *,
    output_x_offset_mm: float = 0.0,
    max_local_e_error_fraction: float = 0.01,
    max_segments: int = 8192,
) -> ConformalMotionPreview:
    """Render one UNOWNED toolpath into isolated, *unexecutable* motion lines.

    All points are (theta-radians, axial-mm) with matching per-point
    *absolute* radius in mm. Each source segment has an explicit physical
    bead cross-sectional area; cannot silently reuse a planar layer's
    E/mm. Native-linear XYZ geometry is subdivided by the independent
    volume-variation estimator. X/Y/Z coordinates are absolute *relative
    to an unspecified preview origin*, with rotary Y=0 at the first
    theta position. No start travel or origin setup is emitted.

    Rounded XYZ/E/F are checked: never emit a positive-E zero-distance
    block, and cap volumetric flow against actual rounded native values.
    """
    config.validate()
    if path.role != "TOP_SKIN" or path.closed:
        raise ValueError("preview only supports open, unowned TOP_SKIN centerlines")
    if len(path.points) < 2:
        raise ValueError("conformal path needs at least two centerline positions")
    if path.segment_widths is not None:
        raise ValueError("conformal segment areas cannot be combined with planar segment_widths")
    radii = path.point_radii
    areas = path.segment_areas_mm2
    if (radii is None or areas is None
            or len(radii) != len(path.points)
            or len(areas) != len(path.points)-1):
        raise ValueError("conformal path requires per-point radii and per-segment bead areas")
    if not math.isfinite(output_x_offset_mm):
        raise ValueError("output X offset must be finite")
    if (not all(math.isfinite(r) and r >= config.bed_radius
                for r in radii)
            or not all(math.isfinite(a) and a > 0 for a in areas)
            or not all(math.isfinite(t) and math.isfinite(x)
                       for t,x in path.points)):
        raise ValueError("conformal points, radii and bead areas must be finite, positive and above the mandrel")
    if type(max_segments) is not int or max_segments <= 0:
        raise ValueError("conformal preview requires positive segment budget")

    phase = path.points[0][0]
    inputs = [
        (float(theta),float(axial),float(radius))
        for (theta,axial),radius in zip(path.points,radii)
    ]
    for a,b in zip(inputs,inputs[1:]):
        if a == b:
            raise ValueError("duplicate nonplanar endpoints cannot extrude")
    if any(abs(math.degrees(point[0]-phase)) > _NATIVE_Y_BUDGET_DEGREES
           for point in inputs):
        raise ValueError("preview rotary motion exceeds native-Y precision budget; output requires coordinated Marlin rebasing")

    total_physical=0.0
    total_e=0.0
    split_count=0
    maximum_error=0.0
    lines=[
        ";CONFORMAL_XYZ_MOTION_PREVIEW_ONLY -- NOT A PRINT JOB",
        ";NO SHELL OWNERSHIP, SAFE TRAVEL, NOZZLE CLEARANCE OR PRINTER ORIGIN HAS BEEN VERIFIED",
        ";Y VALUES ARE IN DEGREES RELATIVE TO THE FIRST CANDIDATE THETA",
    ]
    emitted_e=0.0
    previous_native = (
        round(inputs[0][1]+output_x_offset_mm,3),
        round(0.0,3),
        round(inputs[0][2]-config.bed_radius,4),
    )
    filament_area=math.pi*(config.filament_diameter/2)**2
    max_volume_flow=(config.max_volumetric_flow_mm3_s
                     if config.max_volumetric_flow_mm3_s > 0 else None)
    for i,((start,end),area) in enumerate(zip(zip(inputs,inputs[1:]),areas)):
        remaining=max_segments-split_count
        if remaining <= 0:
            raise ValueError("conformal XYZ preview exceeded total segment budget")
        samples=[
            RadialSurfaceSample(start[0],start[1],start[2],-1,1.0),
            RadialSurfaceSample(end[0],end[1],end[2],-1,1.0),
        ]
        part=plan_conformal_motion(
            samples,
            physical_feed_mm_per_min=config.feed_mm_per_min,
            bead_area_mm2=area,
            filament_diameter_mm=config.filament_diameter,
            max_local_volume_error_fraction=max_local_e_error_fraction,
            max_volumetric_flow_mm3_s=max_volume_flow,
            max_segments=remaining,
        )
        split_count+=len(part.segments)
        total_physical+=part.physical_length_mm
        maximum_error=max(maximum_error,part.largest_volume_error_fraction)
        for seg in part.segments:
            tip=seg.end
            native=(
                round(tip[1]+output_x_offset_mm,3),
                round(math.degrees(tip[0]-phase),3),
                round(tip[2]-config.bed_radius,4),
            )
            native_distance=math.dist(previous_native,native)
            if native_distance <= 1e-10:
                raise ValueError("conformal motion rounds to a stationary positive-extrusion step")
            # E is computed from the exact integrated physical 3D arc, then
            # multiplied by the existing machine flow correction. Preserve
            # the config's calibrated extrusion-per-mm ratio vs its ideal
            # model, including optional nonstandard filament loading.
            ideal_per_mm=config.bead_area/filament_area
            multiplier=(config.extrusion_per_mm / ideal_per_mm
                        *config.flow_percent/100.0)
            total_e += seg.metrics.filament_mm * multiplier
            rounded_e=round(total_e,5)
            if rounded_e <= emitted_e:
                raise ValueError("conformal extrusion is below 0.00001 mm E precision")
            # Native F is already bounded by peak true physical XYZ speed.
            # Conservative floor after output rounding avoids pushing F up.
            feed=math.floor(seg.metrics.marlin_feed_mm_per_min*10)/10
            if max_volume_flow is not None:
                actual_e=rounded_e-emitted_e
                native_flow_cap=(max_volume_flow*60*native_distance
                                 /(actual_e*filament_area))
                feed=min(feed,math.floor(native_flow_cap*10)/10)
            if feed <= 0 or not math.isfinite(feed):
                raise ValueError("conformal feed below Marlin's 0.1 native units/min resolution")
            lines.append(
                f"G1 X{native[0]:.3f} Y{native[1]:.3f} Z{native[2]:.4f}"
                f" E{rounded_e:.5f} F{feed:.1f}"
            )
            previous_native=native
            emitted_e=rounded_e

    return ConformalMotionPreview(
        status="uncertified_motion_preview_only",
        nonplanar_extrusion_enabled=False,
        commands=tuple(lines),
        segment_count=split_count,
        source_segment_count=len(inputs)-1,
        exact_physical_length_mm=total_physical,
        ideal_filament_length_mm=total_e,
        emitted_filament_length_mm=emitted_e,
        maximum_local_e_error_fraction=maximum_error,
        min_radius_mm=min(radii),
        max_radius_mm=max(radii),
    )
