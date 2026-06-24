"""Segmented raster-based cylindrical slicer prototype."""

from .config import SlicerConfig
from .layers import slice_mesh
from .pipeline import slice_to_gcode

__all__ = ["SlicerConfig", "slice_mesh", "slice_to_gcode"]
