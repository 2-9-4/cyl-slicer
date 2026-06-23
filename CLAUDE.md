# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A **cylindrical FDM slicer** — it slices STL meshes radially (peeling concentric cylindrical shells outward) rather than the conventional planar layer-by-layer approach. The output is G-code for a printer that rotates a cylinder around the Y-axis and extrudes radially. There is also a companion C++ OpenGL viewer for the G-code.

## Running the slicer

```bash
# Activate the virtualenv first
source .venv/bin/activate

# Run the slicer (reads test.stl, writes output.gcode)
python slicer.py
```

Useful slicer options:
- `--max-layer N` — slice and emit only model layers `0..N`; layer `N` is treated as the outer/top limit.

Key toggle at the top of `slicer.py`:
- `plot_3d` / `plot_per_layer` — enable/disable matplotlib visualisation.

## Building the G-code viewer (C++)

```bash
cmake -S . -B build
cmake --build build
./build/gcode_viewer          # opens output.gcode by default
```

Depends on: OpenGL, GLFW3, GLEW (system packages). ImGui is fetched automatically via CMake FetchContent.

## Architecture

### Python slicer (`slicer.py` + `togcode.py`)

The pipeline runs entirely in `slicer.py` in two passes:

**Pass 1 — geometry / slicing loop** (iterates over radii from inner to outer):
1. Create a trimesh cylinder at `current_r` and intersect it with the loaded STL mesh.
2. Convert intersection vertices to cylindrical coords `(r, θ, z)`.
3. Identify boundary edges → reconstruct closed 2-D polygons in `(θ, z)` space.
4. Handle "full-circle" polygons (features that wrap all the way around) separately — they are paired up and unrolled into a `[0, 2π]` rectangle.
5. Build a containment hierarchy (`build_hierarchy`) to determine inside/outside levels.
6. Compute `external_peris` (union of all offset perimeter shapes, including ±2π shifted copies for seam handling) and detect `top_surfaces` / `bottom_surfaces` / `needs_support` by diffing adjacent layers.
   - `top_surfaces`: areas of layer `i-1` no longer covered by layer `i` (exposed outer face).
   - `bottom_surfaces`: areas of layer `i` not covered by layer `i-1` (exposed inner face / overhang floor). Layer 0 gets `bottom_surfaces = external_peris` (entire first shell).
   - `needs_support`: areas of layer `i` not covered by layer `i-1`'s `external_peris` — these are overhang regions requiring support material.

**Pass 1.5 — support propagation** (backward pass, outer → inner, after the slicing loop):
- Maintains `support_needed` — the union of all overhang regions not yet grounded by any inner model layer.
- At each layer going inward: subtract the layer's model (`external_peris`) from `support_needed` (model acts as support base), store the remainder as `layer_data[i]["support"]` (extra material this layer must print in empty space), then union in this layer's own `needs_support`.
- Gap layers: if the innermost model layer is above `bed_radius`, generates additional `gap_support_layers` (dicts with `radius` and `support`) printed as pure support shells from `bed_radius` up to the first model layer.

**Pass 2 — G-code generation loop** (iterates over `layer_data`):
1. Inset the external perimeter `wall_count` times for internal perimeters.
2. Compute `solid_infill` (radial lines through top/bottom surfaces), `support_infill` (radial lines through support regions), and `sparse_infill` (radial + axial grid lines for the interior).
3. Convert all geometry to `PathObject` instances via `shapely_to_paths()` in `togcode.py`.
4. Call `PathObject.to_gcode()` which outputs `G0`/`G1` commands where **X = axial position, Y = accumulated angle in degrees** (Y is never reset; Marlin always moves forward on the rotary axis).
5. Gap support layers are written first (before model layers), tagged `;LAYER:support_gap` / `;TYPE:SUPPORT`.

### Key coordinate system
Internally everything lives in `(θ, z)` space (theta in radians, z in mm). Before Shapely operations that depend on metric distances, the θ axis is scaled by `current_r` to convert arc-length to mm, then divided back afterwards (`do_transform` / `do_transform_inverse` closures scattered throughout `slicer.py`).

### Supporting modules
- `togcode.py` — `PathObject` class and `shapely_to_paths()`. Handles densification of arcs, duplicate-point removal, and seam-crossing in G-code (accumulated Y axis).
- `plotutil.py` — single helper `set_axes_equal` for 3D matplotlib plots.

### C++ viewer (`gcode_viewer.cpp`)
Standalone OpenGL app that reads `output.gcode`, parses `;TYPE:` / `;LAYER:` / `;Z:` comments, converts `(X, Y_deg)` back to cylindrical then Cartesian, and renders with per-section-type colours. Controls: left-drag = orbit, right/mid-drag = pan, scroll = zoom, R = reset, Q/Esc = quit.

## Slicing parameters (in `slicer.py`)

```python
slicing_params = {
    "layer_height": 0.2,       # radial step between shells (mm)
    "infill_spacing": 2,        # sparse infill grid spacing (mm)
    "nozzle_size": 0.4,         # extrusion width (mm)
    "wall_count": 3,            # number of perimeter shells
    "feed_mm_per_min": 6000,
    "travel_mm_per_min": 6000,
    "extrusion_per_mm": 0.05,   # E per mm of real travel
}
```

## Dependencies (Python)

Managed in `.venv/`. Key packages: `trimesh`, `numpy`, `scipy`, `shapely`, `matplotlib`, `plotly`, `vispy`, `alphashape`, `solid2`.
