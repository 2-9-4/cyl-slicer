"""Sample possible radial printhead interference around conformal ribbon points.

This is a NECESSARY-CHECK-ONLY feasibility diagnostic, not a certificate of
swept-volume safety and NEVER permission to emit live nonplanar G-code.
The modeled tip is a flat disk and conical exterior only; the actual
heater block, sock, gantry, bed, printed neighboring layers and continuous
motion between samples require separate swept-shape validation.

In the rotating world frame at the tool point:
  head-local tangential = neighbor_radius * sin(neighbor_theta - tool_theta)
  head-local axial      = neighbor_axial - tool_axial
  head-local radial     = neighbor_radius * cos(delta_theta) - tool_radius
The cosine is crucial: subtracting neighbor radius directly would give
wrong interference near the rotating cylinder's tangent plane.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from .nonplanar_surface import RadialSurfaceSampler, segment_surface_length


@dataclass(frozen=True)
class NozzleTipProfile:
    flat_tip_radius_mm: float = 0.2
    cone_half_angle_deg: float = 35.0
    scan_radius_mm: float = 1.5
    neighbor_step_mm: float = 0.2
    path_step_mm: float = 0.4
    contact_allowance_mm: float = 0.025
    safety_margin_mm: float = 0.01
    max_total_probes: int = 250_000

    def validate(self) -> None:
        values = (self.flat_tip_radius_mm, self.cone_half_angle_deg,
                  self.scan_radius_mm, self.neighbor_step_mm,
                  self.path_step_mm, self.contact_allowance_mm,
                  self.safety_margin_mm)
        if not all(math.isfinite(v) for v in values):
            raise ValueError("nozzle profile must be finite")
        if (not 0 < self.flat_tip_radius_mm < self.scan_radius_mm
                or not 0 < self.cone_half_angle_deg < 85
                or self.neighbor_step_mm <= 0
                or self.path_step_mm <= 0
                or self.contact_allowance_mm < 0
                or self.safety_margin_mm < 0
                or type(self.max_total_probes) is not int
                or self.max_total_probes < 1):
            raise ValueError("invalid nozzle footprint, angle, clearance or sample budget")


def _sampler_max_radius(sampler: RadialSurfaceSampler) -> float:
    # Triangles may have their maximum radius at ANY vertex, not merely
    # the first vertex in each triangle.
    vertices = (sampler.vertices, sampler.vertices + sampler.e0,
                sampler.vertices + sampler.e1)
    maximum = max((float(np.hypot(v[:,0],v[:,1]).max())
                   for v in vertices if len(v)), default=0.0)
    if maximum <= 0 or not math.isfinite(maximum):
        raise ValueError("sampler has no finite positive radial mesh extent")
    return maximum


def probe_nozzle_ribbon(
    sampler: RadialSurfaceSampler,
    positions: Sequence[tuple[float, float]],
    profile: NozzleTipProfile,
) -> dict:
    """Reject sampled contacts; an apparent pass is NEVER certified safe.

    Positions are continuous (theta-radians, axial-mm) centerline points.
    Target radii are sampled from exact source triangles, not interpolated
    from old constant-Z layers. A 2D regular grid covers the requested
    circular footprint at each sampled tool pose; unresolved neighbor
    rays remain explicitly unresolved rather than silently treated empty.
    """
    profile.validate()
    if not positions:
        raise ValueError("at least one candidate tool position is required")
    if not all(len(pt) == 2 and all(math.isfinite(v) for v in pt)
               for pt in positions):
        raise ValueError("candidate centerline positions must be finite")
    rmax = _sampler_max_radius(sampler)
    tips = []
    for index, position in enumerate(positions):
        hit = sampler.sample(*position)
        if hit is None:
            return {
                "status": "unresolved_centerline",
                "candidate_is_collision_certified": False,
                "first_invalid_position": index,
                "probes": 0,
            }
        if tips:
            previous = tips[-1]
            length = segment_surface_length(previous, hit)
            count = max(1, math.ceil(length / profile.path_step_mm))
            if count > 16384:
                return {"status": "budget_exceeded",
                        "candidate_is_collision_certified": False, "probes": 0}
            for step in range(1, count + 1):
                fraction = step / count
                theta = previous.theta + (hit.theta - previous.theta)*fraction
                axial = previous.axial + (hit.axial - previous.axial)*fraction
                middle = sampler.sample(theta, axial)
                if middle is None:
                    return {
                        "status": "unresolved_centerline",
                        "candidate_is_collision_certified": False,
                        "first_invalid_position": index,
                        "probes": 0,
                    }
                tips.append(middle)
        else:
            tips.append(hit)
    return _probe_sampled_tips(sampler, tips, profile, rmax)


def _probe_sampled_tips(sampler: RadialSurfaceSampler,
                        tips, profile: NozzleTipProfile, rmax: float) -> dict:
    tested = 0
    unknown = 0
    min_margin = float("inf")
    tan_cone = math.tan(math.radians(profile.cone_half_angle_deg))
    for pose_index, tip in enumerate(tips):
        # Any surface able to intrude above the current tip must have at
        # least the tip's radius. Such a point within the circular nozzle
        # footprint satisfies |delta_theta| <= asin(footprint/tip_radius).
        angular_extent = math.asin(min(1., profile.scan_radius_mm / tip.radius))
        angular_steps = max(1, math.ceil(angular_extent*rmax/profile.neighbor_step_mm))
        axial_steps = max(1, math.ceil(profile.scan_radius_mm/profile.neighbor_step_mm))
        count = (2*angular_steps+1)*(2*axial_steps+1)
        if tested + count > profile.max_total_probes:
            return {
                "status": "budget_exceeded",
                "candidate_is_collision_certified": False,
                "probes": tested,
                "path_samples": pose_index,
            }
        for j in range(-axial_steps,axial_steps+1):
            dx = profile.scan_radius_mm * j/axial_steps
            for i in range(-angular_steps,angular_steps+1):
                dtheta = angular_extent*i/angular_steps
                hit = sampler.sample(tip.theta + dtheta, tip.axial + dx)
                tested += 1
                if hit is None:
                    unknown += 1
                    continue
                lateral = math.hypot(dx,hit.radius*math.sin(dtheta))
                if lateral > profile.scan_radius_mm + 1e-9:
                    continue
                radial_above_tip = hit.radius*math.cos(dtheta) - tip.radius
                if lateral <= profile.flat_tip_radius_mm + 1e-9:
                    margin = profile.contact_allowance_mm - radial_above_tip
                else:
                    lower_head_surface = ((lateral-profile.flat_tip_radius_mm)
                                          / tan_cone)
                    margin = lower_head_surface - radial_above_tip - profile.safety_margin_mm
                min_margin = min(min_margin,margin)
                if margin < -1e-8:
                    return {
                        "status": "sampled_interference",
                        "candidate_is_collision_certified": False,
                        "probes": tested,
                        "path_samples": pose_index+1,
                        "unknown_samples": unknown,
                        "first_interference": {
                            "target_theta": tip.theta,
                            "target_axial_mm": tip.axial,
                            "target_radius_mm": tip.radius,
                            "neighbor_theta": hit.theta,
                            "neighbor_axial_mm": hit.axial,
                            "neighbor_radius_mm": hit.radius,
                            "lateral_offset_mm": lateral,
                            "radial_penetration_mm": -margin,
                        },
                    }
    return {
        "status": "unresolved_neighbors" if unknown else "no_sampled_interference",
        "candidate_is_collision_certified": False,
        "warning": "No continuous swept-tool, heater block, mandrel or printed-layer clearance proof.",
        "probes": tested,
        "path_samples": len(tips),
        "unknown_samples": unknown,
        "min_sampled_clearance_mm": min_margin if math.isfinite(min_margin) else None,
    }
