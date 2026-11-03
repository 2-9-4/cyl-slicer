"""Exact cylindrical 3D motion metrics for future conformal extrusion.

This is an analysis kernel ONLY. The existing constant-radius toolpath/G-code
pipeline is unchanged. A future live XYZ emitter must also implement complete
nozzle swept clearance, support, layer ownership and volumetric calibration.

Machine axes: X=axial mm, Y=rotary degrees, Z=radial mm above bed cylinder.
The slicer describes theta in radians; it may cross +/-2*pi continuously.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


# Positive roots and weights of eighth-order Gauss-Legendre quadrature.
# The exact physical speed for a *linear move in machine coordinates* is
# sqrt(dx^2 + dr^2 + (r(t)*dtheta)^2), with r(t) = r0 + t*(r1-r0).
_NODES = (0.1834346424956498, 0.5255324099163290,
          0.7966664774136267, 0.9602898564975363)
_WEIGHTS = (0.3626837833783620, 0.3137066458778873,
            0.2223810344533745, 0.1012285362903763)


@dataclass(frozen=True)
class ConformalMove:
    """Numerical characteristics; never interpret this as an approved move."""

    length_mm: float
    marlin_native_length: float
    marlin_feed_mm_per_min: float
    filament_mm: float
    min_physical_speed_mm_per_min: float
    max_physical_speed_mm_per_min: float
    max_local_volume_error_fraction: float
    needs_shorter_segments: bool
    volumetric_rate_mm3_s: float


def _physical_length(dx: float, dtheta: float, r0: float, r1: float) -> float:
    """Integrate physical arc length for linear simultaneous X,Y,Z travel.

    In the cylindrical radial geometry, physical speed is *not constant*
    during native linear XYZ interpolation if both theta and radius change.
    This is a Gauss-Legendre integral of the exact Cartesian velocity, not
    a chord shortcut or use of one arbitrarily chosen sample radius.
    """
    dr = r1 - r0
    squared = dx * dx + dr * dr
    if dtheta == 0:
        return math.sqrt(squared)
    if dr == 0:
        return math.hypot(dx, r0 * dtheta)
    total = 0.0
    for node, weight in zip(_NODES, _WEIGHTS):
        lo = r0 + (.5 - node*.5)*dr
        hi = r0 + (.5 + node*.5)*dr
        total += weight*(math.sqrt(squared + (lo*dtheta)**2)
                         + math.sqrt(squared + (hi*dtheta)**2))
    return .5*total


def evaluate_conformal_move(
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    *,
    physical_print_feed_mm_per_min: float,
    bead_area_mm2: float,
    filament_diameter_mm: float,
    max_local_volume_error_fraction: float = 0.01,
    max_volumetric_flow_mm3_s: float | None = None,
) -> ConformalMove:
    """Quantify one continuous (theta-radians, axial-mm, radius-mm) move.

    Native Marlin planner distance includes Y in DEGREES, not r*dtheta.
    Bound native planner F against the *largest instantaneous* physical
    speed anywhere along the polar move, and optionally average flow.
    E is based on integrated true 3D cylindrical length and a fixed
    per-segment bead area. Since E advances linearly with native movement,
    its local volume/physical-mm changes as r(t) varies. Flag moves with
    excessive variation for splitting into shorter segments: NEVER silently
    assume a constant bead width on long simultaneous Y/Z moves.
    """
    if len(start) != 3 or len(end) != 3:
        raise ValueError("conformal move endpoints must be theta, axial, radius")
    values = (*start, *end, physical_print_feed_mm_per_min,
              bead_area_mm2, filament_diameter_mm,
              max_local_volume_error_fraction)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("conformal move parameters must be finite")
    if (start[2] <= 0 or end[2] <= 0 or
            physical_print_feed_mm_per_min <= 0 or bead_area_mm2 <= 0
            or filament_diameter_mm <= 0 or
            not 0 <= max_local_volume_error_fraction < 1):
        raise ValueError("conformal radius, feed, bead and filament must be positive")
    if (max_volumetric_flow_mm3_s is not None
            and (not math.isfinite(max_volumetric_flow_mm3_s)
                 or max_volumetric_flow_mm3_s <= 0)):
        raise ValueError("volumetric ceiling must be finite and positive")

    dtheta = end[0] - start[0]
    dx = end[1] - start[1]
    dr = end[2] - start[2]
    native_length = math.sqrt(dx*dx + math.degrees(dtheta)**2 + dr*dr)
    if native_length <= 1e-12:
        raise ValueError("conformal move must have finite nonzero machine travel")
    length = _physical_length(dx,dtheta,start[2],end[2])
    if not (math.isfinite(length) and length > 0):
        raise ValueError("conformal physical length is invalid")
    squared = dx*dx + dr*dr
    minimum = math.sqrt(squared + (min(start[2],end[2])*dtheta)**2)
    maximum = math.sqrt(squared + (max(start[2],end[2])*dtheta)**2)
    # At requested F, this is the native feed that cannot exceed F at
    # *either* radial endpoint (or anywhere between them).
    feed = physical_print_feed_mm_per_min * native_length / maximum
    if max_volumetric_flow_mm3_s is not None:
        # E (and volume) is linear over native planner travel, so volume/s
        # is area*length*feed/(60*native_length).
        feed = min(feed, 60*max_volumetric_flow_mm3_s
                   * native_length / (bead_area_mm2*length))
    mean = feed*length/native_length
    local_max = feed*maximum/native_length
    local_min = feed*minimum/native_length
    filament_area = math.pi*(filament_diameter_mm*.5)**2
    e = bead_area_mm2*length/filament_area
    # Difference between commanded E/mm of actual plastic and the ideal
    # constant bead E/mm. The physical-speed distribution is monotone in
    # positive radius; its maximum proportional error occurs at endpoints.
    variation = max(abs(length/minimum-1), abs(length/maximum-1))
    return ConformalMove(
        length_mm=length,
        marlin_native_length=native_length,
        marlin_feed_mm_per_min=feed,
        filament_mm=e,
        min_physical_speed_mm_per_min=local_min,
        max_physical_speed_mm_per_min=local_max,
        max_local_volume_error_fraction=variation,
        needs_shorter_segments=variation > max_local_volume_error_fraction,
        volumetric_rate_mm3_s=bead_area_mm2*mean/60,
    )
