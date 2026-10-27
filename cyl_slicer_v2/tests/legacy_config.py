"""Explicit historical geometry fixture, NOT the application default.

The quality profile is the canonical SlicerConfig() default for production.
Tests written against pre-change dimensions use this dataclass subclass to
keep the same physical fixture geometry rather than silently weakening
the application's new defaults.
"""
from __future__ import annotations

from dataclasses import dataclass
from cyl_slicer_v2.config import SlicerConfig


@dataclass(frozen=True)
class LegacySlicerConfig(SlicerConfig):
    layer_height: float = 0.5
    nozzle_width: float = 0.4
    nozzle_diameter: float | None = None
    bead_model: str = "rectangular"
    wall_count: int = 2
    top_layers: int = 3
    bottom_layers: int = 3
    infill_spacing: float = 2.0
    sparse_infill_pattern: str = "both"
    solid_infill_pattern: str = "axial"
    monotonic_top_skin: bool = False
    gap_fill_enabled: bool = False
    support_min_island_area: float = 0.05
    support_interface_layers: int = 2
    support_interface_spacing: float | None = None
    raft_merge_gap_mm: float = 0.0
    bed_radius: float = 0.0
    feed_mm_per_min: float = 1800.0
    travel_mm_per_min: float = 3000.0
    retract_speed_mm_per_min: float = 1800.0
    retract_mm: float = 1.0
    hotend_temp: int = 215
    material_profile: str = "unprofiled"
    fan_start_height: float = 2.0
    fan_ramp_height_mm: float = 0.0
    cooling_min_layer_time_s: float = 0.0
    cooling_min_island_print_time_s: float = 0.0
    min_path_length: float = 0.6
    intersection_curve_error_mm: float = 0.025
    layer_backend: str = "raster"
