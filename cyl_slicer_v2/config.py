from __future__ import annotations

from dataclasses import dataclass
import math

from .calibration import RetractionSweep


@dataclass(frozen=True)
class SlicerConfig:
    layer_height: float = 0.16
    nozzle_width: float = 0.42  # Deposited track width; not physical nozzle orifice.
    nozzle_diameter: float | None = 0.4
    # "rounded" is the calibrated quality approximation. "rectangular"
    # remains available for old G-code/reproducibility profiles.
    bead_model: str = "rounded"
    filament_diameter: float = 1.75
    wall_count: int = 4
    wall_order: str = "outer-first"  # Inside-out for overhangs if selected.
    seam_policy: str = "nearest"  # nearest | aligned | corner, for closed loops only.
    seam_angle_deg: float = 180.0  # Cylinder phase target for aligned/corner.
    seam_exclude_deg: tuple[float, float] | None = None
    top_layers: int = 5
    bottom_layers: int = 5
    infill_spacing: float = 1.5
    sparse_infill_pattern: str = "alternate"  # optional alternating single-grid layer.
    solid_infill_pattern: str = "alternate"  # optional alternating exposed fill angle.
    monotonic_top_skin: bool = True  # One-way exposed top passes enabled for quality.
    # Fraction of line width that solid/sparse infill may enter the inner
    # perimeter bead, using the selected coupled bead cross-section model.
    infill_wall_overlap_fraction: float = 0.15
    support_spacing: float = 1.5
    support_print_feed_mm_per_min: float | None = None
    support_interface_print_feed_mm_per_min: float | None = None
    support_flow_ratio: float = 1.0
    support_interface_flow_ratio: float = 1.0
    support_min_island_area: float = 0.25
    support_min_region_area: float = 1.0
    support_bearing_min_area: float = 1.0
    support_merge_gap: float = 0.0
    support_gap_layers: int = 1
    support_interface_layers: int = 3
    support_interface_spacing: float | None = 0.4
    support_meander_max_angle_deg: float = 15.0
    support_cross_spacing: float | None = None
    raft_layers: int = 4
    raft_grid_spacing: float | None = None
    raft_axial_gap: float = 1.0
    theta_step_deg: float = 5.0
    axial_step: float = 1.0
    guard_columns: int = 2
    bed_radius: float = 30.0
    generate_gap_support: bool = True
    max_radius: float | None = None
    min_radius: float | None = None
    axial_margin: float = 0.0
    axial_min: float = 15.0
    feed_mm_per_min: float = 2700.0
    material_profile: str = "pla"  # Metadata; never guesses a material's K.
    linear_advance_k: float | None = None
    linear_advance_firmware: str | None = None  # Explicit firmware build declaration.
    travel_mm_per_min: float = 9000.0
    extrusion_per_mm: float | None = None
    retraction_sweep: RetractionSweep | None = None
    retract_speed_mm_per_min: float = 1200.0
    retract_mm: float = 1.4
    z_hop_height: float = 0.4
    z_hop_min_travel: float = 5.0
    z_hop_on_perimeter_crossing: bool = True
    route_travel_inside_island: bool = True
    route_max_extra_mm: float = 6.0  # Upper bound beyond direct travel length.
    route_max_factor: float = 1.65  # Multiplicative detour limit.
    route_max_vertices: int = 96  # Bound visibility-graph work.
    route_max_unretracted_mm: float = 5.0
    hotend_temp: int = 195
    fan_speed: int = 255
    fan_start_height: float = 1.0
    fan_ramp_height_mm: float = 1.0  # 0 restores immediate fan step.
    cooling_min_layer_time_s: float = 8.0  # 0 disables conservative slowdown.
    cooling_min_print_speed_mm_per_min: float = 300.0
    min_path_length: float = 0.1
    max_axial_segment: float = 1.5  # Secondary cap on curved-face intersection samples.
    intersection_curve_error_mm: float = 0.01  # Maximum physical unwrapped chord deviation.
    intersection_strict_recovery: bool = False  # Fail if shifted samples disagree.
    slice_workers: int = 1
    layer_backend: str = "intersection"
    intersection_refine_factor: int = 0
    toolpath_workers: int = 1

    def __post_init__(self) -> None:
        if self.extrusion_per_mm is not None or self.filament_diameter <= 0:
            return
        filament_area = math.pi * (self.filament_diameter / 2.0) ** 2
        object.__setattr__(
            self,
            "extrusion_per_mm",
            self.bead_area / filament_area,
        )

    @property
    def bead_area(self) -> float:
        """Deposited cross-sectional area (mm²) per mm of path.

        Rounded bead: a (w-h) by h rectangle and two half-circles of
        diameter h. Only valid when h <= w. This is the standard squashed
        strand approximation, not a literal cylindrical cross-section.
        """
        w, h = self.nozzle_width, self.layer_height
        if self.bead_model == "rounded":
            return h * (w - h) + math.pi * h * h / 4.0
        return w * h

    @property
    def dense_line_spacing(self) -> float:
        """Center-to-center pitch of adjacent dense beads, in physical mm.

        Matching pitch to A/h keeps the deposited volume per covered
        planar area consistent with bead_area. The outermost contour still
        sits w/2 inside the material boundary, independent of pitch.
        """
        return self.bead_area / self.layer_height

    def wall_coverage_depth(self, count: int) -> float:
        """Nominal inwards extent of count deposited perimeters (mm)."""
        if count < 1:
            raise ValueError("wall count must be at least one")
        return self.nozzle_width + (count - 1) * self.dense_line_spacing

    @property
    def theta_step(self) -> float:
        return math.radians(self.theta_step_deg)

    def validate(self) -> None:
        if not math.isfinite(self.layer_height) or self.layer_height <= 0:
            raise ValueError("layer_height must be finite and positive")
        if not math.isfinite(self.nozzle_width) or self.nozzle_width <= 0:
            raise ValueError("nozzle_width must be finite and positive")
        if self.nozzle_diameter is not None:
            if not math.isfinite(self.nozzle_diameter) or self.nozzle_diameter <= 0:
                raise ValueError("nozzle_diameter must be finite and positive")
            if self.layer_height > self.nozzle_diameter + 1e-9:
                raise ValueError("layer_height cannot exceed specified nozzle_diameter")
        if self.bead_model not in {"rectangular", "rounded"}:
            raise ValueError("bead_model must be 'rectangular' or 'rounded'")
        if self.bead_model == "rounded" and self.layer_height > self.nozzle_width:
            raise ValueError("rounded bead requires layer_height <= line width")
        if self.filament_diameter <= 0:
            raise ValueError("filament_diameter must be positive")
        if self.extrusion_per_mm is None or not math.isfinite(self.extrusion_per_mm) or self.extrusion_per_mm <= 0:
            raise ValueError("extrusion_per_mm must be finite and positive")
        if self.wall_count < 1:
            raise ValueError("wall_count must be at least 1")
        if self.wall_order not in {"outer-first", "inner-first"}:
            raise ValueError("wall_order must be 'outer-first' or 'inner-first'")
        if self.seam_policy not in {"nearest", "aligned", "corner"}:
            raise ValueError("seam_policy must be nearest, aligned or corner")
        if not math.isfinite(self.seam_angle_deg) or not 0 <= self.seam_angle_deg < 360:
            raise ValueError("seam_angle_deg must be in [0, 360)")
        if self.seam_exclude_deg is not None:
            pair = self.seam_exclude_deg
            if (len(pair) != 2 or any(not math.isfinite(x) or not 0 <= x < 360 for x in pair)
                    or pair[0] == pair[1]):
                raise ValueError("seam_exclude_deg must contain two distinct angles in [0, 360)")
        if self.top_layers < 1:
            raise ValueError("top_layers must be at least 1")
        if self.bottom_layers < 1:
            raise ValueError("bottom_layers must be at least 1")
        if self.infill_spacing <= 0:
            raise ValueError("infill_spacing must be positive")
        if self.sparse_infill_pattern not in {"both", "alternate"}:
            raise ValueError("sparse_infill_pattern must be 'both' or 'alternate'")
        if (self.sparse_infill_pattern == "alternate"
                and self.infill_spacing * 0.5 < self.dense_line_spacing - 1e-9):
            raise ValueError("alternate sparse infill requires infill_spacing >= 2 * dense_line_spacing")
        if self.solid_infill_pattern not in {"axial", "alternate"}:
            raise ValueError("solid_infill_pattern must be 'axial' or 'alternate'")
        if not isinstance(self.monotonic_top_skin, bool):
            raise ValueError("monotonic_top_skin must be bool")
        if not math.isfinite(self.infill_wall_overlap_fraction) or not 0 <= self.infill_wall_overlap_fraction <= 0.5:
            raise ValueError("infill_wall_overlap_fraction must be finite and within [0, 0.5]")
        if self.support_spacing <= 0:
            raise ValueError("support_spacing must be positive")
        for name in ("support_print_feed_mm_per_min", "support_interface_print_feed_mm_per_min"):
            speed = getattr(self, name)
            if speed is not None and (not math.isfinite(speed) or speed <= 0):
                raise ValueError(f"{name} must be finite and positive when set")
        for name in ("support_flow_ratio", "support_interface_flow_ratio"):
            flow = getattr(self, name)
            if not math.isfinite(flow) or not 0 < flow <= 2:
                raise ValueError(f"{name} must be finite and in (0, 2]")
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
        if self.support_cross_spacing is not None and self.support_cross_spacing <= 0:
            raise ValueError("support_cross_spacing must be positive when set")
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
        if not math.isfinite(self.intersection_curve_error_mm) or self.intersection_curve_error_mm <= 0:
            raise ValueError("intersection_curve_error_mm must be finite and positive")
        if self.min_path_length < 0:
            raise ValueError("min_path_length must be non-negative")
        if self.material_profile not in {"unprofiled", "pla", "petg", "abs", "tpu"}:
            raise ValueError("material_profile must be unprofiled, pla, petg, abs, or tpu")
        if self.linear_advance_firmware not in {None, "marlin-m900"}:
            raise ValueError("linear_advance_firmware must be None or marlin-m900")
        if self.linear_advance_k is not None:
            if not math.isfinite(self.linear_advance_k) or not (0 <= self.linear_advance_k <= 10):
                raise ValueError("linear_advance_k must be finite within [0, 10]")
            if self.linear_advance_firmware != "marlin-m900":
                raise ValueError("linear_advance_k requires explicit Marlin M900 firmware support")
        elif self.linear_advance_firmware is not None:
            raise ValueError("linear_advance_firmware requires an explicit K value")
        if not math.isfinite(self.feed_mm_per_min) or self.feed_mm_per_min <= 0:
            raise ValueError("feed_mm_per_min must be finite and positive")
        if not math.isfinite(self.travel_mm_per_min) or self.travel_mm_per_min <= 0:
            raise ValueError("travel_mm_per_min must be finite and positive")
        if not math.isfinite(self.retract_speed_mm_per_min) or self.retract_speed_mm_per_min <= 0:
            raise ValueError("retract_speed_mm_per_min must be finite and positive")
        if not math.isfinite(self.retract_mm) or self.retract_mm < 0:
            raise ValueError("retract_mm must be non-negative")
        if not isinstance(self.fan_speed, int) or not 0 <= self.fan_speed <= 255:
            raise ValueError("fan_speed must be an integer from 0 to 255")
        if not math.isfinite(self.fan_start_height) or self.fan_start_height < 0:
            raise ValueError("fan_start_height must be finite and non-negative")
        if not math.isfinite(self.fan_ramp_height_mm) or self.fan_ramp_height_mm < 0:
            raise ValueError("fan_ramp_height_mm must be finite and non-negative")
        if not math.isfinite(self.cooling_min_layer_time_s) or self.cooling_min_layer_time_s < 0:
            raise ValueError("cooling_min_layer_time_s must be finite and non-negative")
        if (not math.isfinite(self.cooling_min_print_speed_mm_per_min)
                or self.cooling_min_print_speed_mm_per_min <= 0):
            raise ValueError("cooling_min_print_speed_mm_per_min must be finite and positive")
        if (self.cooling_min_layer_time_s > 0
                and self.cooling_min_print_speed_mm_per_min > self.feed_mm_per_min):
            raise ValueError("cooling minimum print speed cannot exceed normal print speed")
        if self.z_hop_height < 0:
            raise ValueError("z_hop_height must be non-negative")
        if self.z_hop_min_travel < 0:
            raise ValueError("z_hop_min_travel must be non-negative")
        if not isinstance(self.route_travel_inside_island, bool):
            raise ValueError("route_travel_inside_island must be bool")
        if (not math.isfinite(self.route_max_extra_mm)
                or self.route_max_extra_mm < 0):
            raise ValueError("route_max_extra_mm must be finite and non-negative")
        if not math.isfinite(self.route_max_factor) or self.route_max_factor < 1:
            raise ValueError("route_max_factor must be finite and >= 1")
        if not isinstance(self.route_max_vertices, int) or not 4 <= self.route_max_vertices <= 512:
            raise ValueError("route_max_vertices must be between 4 and 512")
        if (not math.isfinite(self.route_max_unretracted_mm)
                or self.route_max_unretracted_mm < 0):
            raise ValueError("route_max_unretracted_mm must be finite and non-negative")
        if self.slice_workers < 1:
            raise ValueError("slice_workers must be at least 1")
        if self.toolpath_workers < 1:
            raise ValueError("toolpath_workers must be at least 1")
        if self.layer_backend not in {"raster", "intersection"}:
            raise ValueError("layer_backend must be 'raster' or 'intersection'")
        if not isinstance(self.intersection_strict_recovery, bool):
            raise ValueError("intersection_strict_recovery must be a bool")
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
