from __future__ import annotations

import trimesh

from .config import SlicerConfig
from .gcode import emit_gcode
from .layers import slice_mesh
from .toolpaths import build_toolpaths


def slice_to_gcode(mesh: trimesh.Trimesh, config: SlicerConfig | None = None, max_layers: int | None = None) -> str:
    config = config or SlicerConfig()
    layers = slice_mesh(mesh, config, max_layers=max_layers)
    layer_paths = build_toolpaths(layers, config)
    return emit_gcode(layer_paths, config)
