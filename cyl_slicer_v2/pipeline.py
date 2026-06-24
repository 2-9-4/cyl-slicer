from __future__ import annotations

import trimesh

from .config import SlicerConfig
from .gcode import emit_gcode
from .layers import slice_mesh
from .occupancy_backends import OccupancyBackend
from .toolpaths import build_toolpaths


def slice_to_gcode(
    mesh: trimesh.Trimesh,
    config: SlicerConfig | None = None,
    max_layers: int | None = None,
    occupancy_backend: OccupancyBackend | None = None,
) -> str:
    config = config or SlicerConfig()
    layers = slice_mesh(mesh, config, max_layers=max_layers, occupancy_backend=occupancy_backend)
    layer_paths = build_toolpaths(layers, config)
    return emit_gcode(layer_paths, config)
