from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class SlicerConfig:
    layer_height: float = 0.5
    nozzle_width: float = 0.45
    wall_count: int = 2
    infill_spacing: float = 2.0
    support_spacing: float = 1.5
    theta_step_deg: float = 5.0
    axial_step: float = 1.0
    guard_columns: int = 2
    bed_radius: float = 0.0
    max_radius: float | None = None
    min_radius: float | None = None
    axial_margin: float = 0.0
    feed_mm_per_min: float = 1800.0
    travel_mm_per_min: float = 3000.0
    extrusion_per_mm: float = 0.035
    retract_mm: float = 0.0
    hotend_temp: int = 215
    min_path_length: float = 1e-4
    max_axial_segment: float = 1.5
    slice_workers: int = 1
    layer_backend: str = "raster"
    intersection_refine_factor: int = 0
    toolpath_workers: int = 1

    @property
    def theta_step(self) -> float:
        return math.radians(self.theta_step_deg)

    def validate(self) -> None:
        if self.layer_height <= 0:
            raise ValueError("layer_height must be positive")
        if self.nozzle_width <= 0:
            raise ValueError("nozzle_width must be positive")
        if self.wall_count < 1:
            raise ValueError("wall_count must be at least 1")
        if self.infill_spacing <= 0:
            raise ValueError("infill_spacing must be positive")
        if self.support_spacing <= 0:
            raise ValueError("support_spacing must be positive")
        if self.theta_step_deg <= 0 or self.theta_step_deg > 90:
            raise ValueError("theta_step_deg must be in (0, 90]")
        if self.axial_step <= 0:
            raise ValueError("axial_step must be positive")
        if self.guard_columns < 0:
            raise ValueError("guard_columns must be non-negative")
        if self.bed_radius < 0:
            raise ValueError("bed_radius must be non-negative")
        if self.max_axial_segment <= 0:
            raise ValueError("max_axial_segment must be positive")
        if self.slice_workers < 1:
            raise ValueError("slice_workers must be at least 1")
        if self.toolpath_workers < 1:
            raise ValueError("toolpath_workers must be at least 1")
        if self.layer_backend not in {"raster", "intersection"}:
            raise ValueError("layer_backend must be 'raster' or 'intersection'")
        if self.intersection_refine_factor < 0:
            raise ValueError("intersection_refine_factor must be non-negative")
        if self.max_radius is not None and self.max_radius <= 0:
            raise ValueError("max_radius must be positive when set")
        if self.min_radius is not None and self.min_radius < 0:
            raise ValueError("min_radius must be non-negative when set")
        if (
            self.min_radius is not None
            and self.max_radius is not None
            and self.min_radius >= self.max_radius
        ):
            raise ValueError("min_radius must be less than max_radius")
