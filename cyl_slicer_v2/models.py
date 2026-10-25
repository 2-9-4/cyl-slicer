from __future__ import annotations

from dataclasses import dataclass, field

from shapely.geometry.base import BaseGeometry


@dataclass
class LayerSlice:
    index: int
    radius: float
    occupied: BaseGeometry
    top_surfaces: BaseGeometry
    bottom_surfaces: BaseGeometry
    needs_support: BaseGeometry
    support: BaseGeometry
    theta_window: tuple[float, float] = (0.0, 0.0)
    is_full_wrap: bool = False
    is_raft: bool = False
    raft_direction: str | None = None
    support_interface: BaseGeometry | None = None
    support_interface_direction: str | None = None
    support_interface_target: BaseGeometry | None = None


@dataclass
class Toolpath:
    role: str
    radius: float
    points: list[tuple[float, float]]
    closed: bool = False


@dataclass
class LayerToolpaths:
    layer_index: int
    radius: float
    paths_by_role: dict[str, list[Toolpath]] = field(default_factory=dict)


@dataclass
class GcodeState:
    theta_accum: float = 0.0
    theta_origin: float = 0.0
    axial: float = 0.0
    radial: float = 0.0
    extrusion: float = 0.0
    retracted: bool = False
    retracted_mm: float = 0.0
    retracted_speed_mm_per_min: float = 1800.0
    native_y_origin_degrees: float = 0.0
