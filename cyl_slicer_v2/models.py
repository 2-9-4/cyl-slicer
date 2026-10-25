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
    travel_region: BaseGeometry | None = None
    island_paths: list[dict[str, list[Toolpath]]] = field(default_factory=list)
    # radius above is the geometric slice SAMPLE plane, never the nozzle Z.
    # A physical deposited slab occupies [sample - thickness/2,
    # sample + thickness/2]; default calculation uses config.layer_height.
    nozzle_radius: float | None = None
    deposited_thickness: float | None = None

    def deposition_radius(self, configured_layer_height: float) -> float:
        if self.nozzle_radius is not None:
            return self.nozzle_radius
        height = self.deposited_thickness if self.deposited_thickness is not None else configured_layer_height
        return self.radius + height * 0.5


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
