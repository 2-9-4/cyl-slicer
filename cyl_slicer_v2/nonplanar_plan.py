"""Volume-consistent, bounded segmentation of candidate nonplanar XYZ moves.

This is a *planning and diagnostic* module, NEVER print-ready G-code.
It preserves each input machine-linear (theta, axial, radius) trajectory
exactly while subdividing intervals where Marlin's native-linear E motion
would vary materially in deposited plastic per physical millimeter.

Upstream surface samples must already obey mesh chord/slope/patch checks;
splitting alone cannot certify nozzle collision, shell ownership or support.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

from .nonplanar_motion import ConformalMove, evaluate_conformal_move
from .nonplanar_surface import RadialSurfaceSample


@dataclass(frozen=True)
class ConformalPlannedSegment:
    start: tuple[float, float, float]
    end: tuple[float, float, float]
    metrics: ConformalMove
    cumulative_e_start_mm: float
    cumulative_e_end_mm: float


@dataclass(frozen=True)
class ConformalMotionPlan:
    segments: tuple[ConformalPlannedSegment, ...]
    physical_length_mm: float
    filament_length_mm: float
    largest_volume_error_fraction: float
    split_count: int
    input_segment_count: int


def plan_conformal_motion(
    samples: Sequence[RadialSurfaceSample],
    *,
    physical_feed_mm_per_min: float,
    bead_area_mm2: float,
    filament_diameter_mm: float,
    max_local_volume_error_fraction: float = 0.01,
    max_volumetric_flow_mm3_s: float | None = None,
    max_segments: int = 8192,
    max_split_depth: int = 20,
) -> ConformalMotionPlan:
    """Subdivide native-linear polar moves, preserving their physical sweep.

    An accepted segment's native-linear E can at worst deviate from its
    desired bead-area/physical-length by the caller's relative tolerance.
    E is computed by integrating actual 3D cylindrical length. A split
    uses midpoint of ALL THREE machine coordinates, including radial R,
    so it does *not* change the original straight-native motion geometry.
    This is different from resampling an actual mesh: do that in the
    upstream topography/chord-error stage, before planning physical E.

    No Marlin commands, retracts, z-hops or machine actions are emitted.
    """
    if (type(max_segments) is not int or max_segments < 1
            or type(max_split_depth) is not int or not 0 <= max_split_depth <= 28):
        raise ValueError("motion plan needs a finite positive segment/depth budget")
    if len(samples) < 2:
        raise ValueError("motion plan needs at least two surface samples")
    for point in samples:
        if not isinstance(point, RadialSurfaceSample):
            raise ValueError("motion plan requires RadialSurfaceSample positions")
        if not all(math.isfinite(x) for x in
                   (point.theta, point.axial, point.radius)) or point.radius <= 0:
            raise ValueError("surface positions must be finite with positive radius")
    for a,b in zip(samples,samples[1:]):
        if (a.theta == b.theta and a.axial == b.axial
                and a.radius == b.radius):
            raise ValueError("duplicated surface sample has no machine travel")
    initial_count = len(samples)-1
    if initial_count > max_segments:
        raise ValueError("motion plan exceeds maximum segment count")

    raw: list[tuple[tuple[float, float, float],
                    tuple[float, float, float], ConformalMove]] = []
    split_count = 0

    def evaluate(a, b):
        return evaluate_conformal_move(
            a, b, physical_print_feed_mm_per_min=physical_feed_mm_per_min,
            bead_area_mm2=bead_area_mm2,
            filament_diameter_mm=filament_diameter_mm,
            max_local_volume_error_fraction=max_local_volume_error_fraction,
            max_volumetric_flow_mm3_s=max_volumetric_flow_mm3_s,
        )

    def subdivide(a, b, depth):
        nonlocal split_count
        if len(raw)+1 > max_segments:
            raise ValueError("conformal E segmentation exceeds maximum segment budget")
        metric=evaluate(a,b)
        if metric.needs_shorter_segments:
            if depth >= max_split_depth:
                raise ValueError("conformal E segmentation could not meet the requested local volume accuracy")
            mid=tuple((v+w)*.5 for v,w in zip(a,b))
            if mid == a or mid == b:
                raise ValueError("conformal E segmentation reached numerical precision limit")
            split_count += 1
            subdivide(a,mid,depth+1)
            subdivide(mid,b,depth+1)
            return
        raw.append((a,b,metric))

    for a,b in zip(samples,samples[1:]):
        start=(a.theta,a.axial,a.radius)
        end=(b.theta,b.axial,b.radius)
        subdivide(start,end,0)
    # Recursive subtrees may each temporarily fit the budget, but the
    # combined plan must also respect the global cap.
    if len(raw) > max_segments:
        raise ValueError("conformal E segmentation exceeds maximum segment budget")

    cumulative=0.0
    segments=[]
    for a,b,m in raw:
        begin=cumulative
        cumulative+=m.filament_mm
        segments.append(ConformalPlannedSegment(
            start=a,end=b,metrics=m,
            cumulative_e_start_mm=begin,
            cumulative_e_end_mm=cumulative,
        ))
    return ConformalMotionPlan(
        segments=tuple(segments),
        physical_length_mm=sum(segment.metrics.length_mm
                               for segment in segments),
        filament_length_mm=cumulative,
        largest_volume_error_fraction=max(
            segment.metrics.max_local_volume_error_fraction
            for segment in segments
        ),
        split_count=split_count,
        input_segment_count=initial_count,
    )
