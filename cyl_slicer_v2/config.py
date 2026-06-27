from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class SlicerConfig:
    layer_height: float = 0.5
    nozzle_width: float = 0.4
    filament_diameter: float = 1.75
    wall_count: int = 2
    top_layers: int = 3
    bottom_layers: int = 3
    infill_spacing: float = 2.0
    support_spacing: float = 1.5
    support_min_island_area: float = 0.05
    support_min_region_area: float = 1.0
    support_bearing_min_area: float = 1.0
    support_merge_gap: float = 0.0
    support_gap_layers: int = 1
    support_interface_layers: int = 2
    support_interface_spacing: float | None = None
    support_meander_max_angle_deg: float = 15.0
    raft_layers: int = 4
    raft_grid_spacing: float | None = None
    raft_axial_gap: float = 1.0
    theta_step_deg: float = 5.0
    axial_step: float = 1.0
    guard_columns: int = 2
    bed_radius: float = 0.0
    generate_gap_support: bool = True
    max_radius: float | None = None
    min_radius: float | None = None
    axial_margin: float = 0.0
    axial_min: float = 15.0
    feed_mm_per_min: float = 1800.0
    travel_mm_per_min: float = 3000.0
    extrusion_per_mm: float | None = None
    retract_mm: float = 1.0
    hotend_temp: int = 215
    min_path_length: float = 0.6
    max_axial_segment: float = 1.5
    slice_workers: int = 1
    layer_backend: str = "raster"
    intersection_refine_factor: int = 0
    toolpath_workers: int = 1

    def __post_init__(self) -> None:
        if self.extrusion_per_mm is not None or self.filament_diameter <= 0:
            return
        filament_area = math.pi * (self.filament_diameter / 2.0) ** 2
        object.__setattr__(
            self,
            "extrusion_per_mm",
            self.nozzle_width * self.layer_height / filament_area,
        )

    @property
    def theta_step(self) -> float:
        return math.radians(self.theta_step_deg)

    def validate(self) -> None:
        if self.layer_height <= 0:
            raise ValueError("layer_height must be positive")
        if self.nozzle_width <= 0:
            raise ValueError("nozzle_width must be positive")
        if self.filament_diameter <= 0:
            raise ValueError("filament_diameter must be positive")
        if self.extrusion_per_mm is None or self.extrusion_per_mm <= 0:
            raise ValueError("extrusion_per_mm must be positive")
        if self.wall_count < 1:
            raise ValueError("wall_count must be at least 1")
        if self.top_layers < 1:
            raise ValueError("top_layers must be at least 1")
        if self.bottom_layers < 1:
            raise ValueError("bottom_layers must be at least 1")
        if self.infill_spacing <= 0:
            raise ValueError("infill_spacing must be positive")
        if self.support_spacing <= 0:
            raise ValueError("support_spacing must be positive")
        if self.support_min_island_area < 0:
            raise ValueError("support_min_island_area must be non-negative")
        if self.support_min_region_area < 0:
            raise ValueError("support_min_region_area must be non-negative")
        if self.support_bearing_min_area < 0:
            raise ValueError("support_bearing_min_area must be non-negative")
        if self.support_merge_gap < 0:
            raise ValueError("support_merge_gap must be non-negative")
        if self.support_gap_layers < 0:
            raise ValueError("support_gap_layers must be non-negative")
        if self.support_interface_layers < 0:
            raise ValueError("support_interface_layers must be non-negative")
        if self.support_interface_spacing is not None and self.support_interface_spacing <= 0:
            raise ValueError("support_interface_spacing must be positive when set")
        if self.support_meander_max_angle_deg <= 0 or self.support_meander_max_angle_deg > 360:
            raise ValueError("support_meander_max_angle_deg must be in (0, 360]")
        if self.raft_layers < 0:
            raise ValueError("raft_layers must be non-negative")
        if self.raft_grid_spacing is not None and self.raft_grid_spacing <= 0:
            raise ValueError("raft_grid_spacing must be positive when set")
        if self.raft_axial_gap < 0:
            raise ValueError("raft_axial_gap must be non-negative")
        if self.theta_step_deg <= 0 or self.theta_step_deg > 90:
            raise ValueError("theta_step_deg must be in (0, 90]")
        if self.axial_step <= 0:
            raise ValueError("axial_step must be positive")
        if self.axial_min < 0:
            raise ValueError("axial_min must be non-negative")
        if self.guard_columns < 0:
            raise ValueError("guard_columns must be non-negative")
        if self.bed_radius < 0:
            raise ValueError("bed_radius must be non-negative")
        if self.max_axial_segment <= 0:
            raise ValueError("max_axial_segment must be positive")
        if self.min_path_length < 0:
            raise ValueError("min_path_length must be non-negative")
        if self.retract_mm < 0:
            raise ValueError("retract_mm must be non-negative")
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
