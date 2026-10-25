# Cylindrical Slicer V2

This is a self-contained raster-based prototype. It does not import or mutate
the legacy `slicer.py`/`togcode.py` pipeline.

## Run

```bash
.venv/bin/python -m cyl_slicer_v2.cli --fixture tube --out /tmp/tube_v2.gcode --max-layers 2
.venv/bin/python -m cyl_slicer_v2.cli tube.stl --out /tmp/tube_stl_v2.gcode --max-layers 2
```

The default occupancy backend calls `trimesh.contains` directly. A faster
approximate voxel-grid backend can be enabled and cached:

```bash
.venv/bin/python -m cyl_slicer_v2.cli cubeandring.stl \
  --occupancy-backend grid \
  --grid-voxel-mm 1.0 \
  --grid-cache /tmp/cubeandring_v2_grid.npz \
  --out /tmp/cubeandring_grid_v2.gcode
```

Use `--rebuild-grid` to refresh an existing cache. Cached grids validate their
voxel size and mesh bounds before use.

There is also a narrow triangle/cylinder intersection layer backend. It bypasses
point occupancy sampling and extracts closed intersection loops directly:

```bash
.venv/bin/python -m cyl_slicer_v2.cli cubeandring.stl \
  --layer-backend intersection \
  --out /tmp/cubeandring_intersection_v2.gcode
```

Layer toolpath generation also runs in parallel in the CLI. By default it uses
the same worker count as `--slice-workers`; override it with
`--toolpath-workers N` when geometry slicing and toolpath generation need
different process counts.

When a model's innermost radius is above the printer build cylinder, v2 emits
support-only foundation layers with negative layer numbers before model layer 0.
Set the printer surface with `--build-cylinder-radius R` (`--bed-radius R` also
works), or pass `--no-gap-support` to disable those foundation layers.
The first four foundation layers default to a continuous raft envelope under the
part. Raft layers alternate one grid direction per layer, radial first, with
spacing defaulting to nozzle width. Change these with `--raft-layers N` and
`--raft-grid-spacing S`.

Extrusion is derived from line width, layer height, and filament cross-section:
`nozzle_width * layer_height / (pi * (filament_diameter / 2)^2)`. The default
filament diameter is 1.75 mm; override it with `--filament-diameter D`.

The intersection backend is experimental and intentionally small. Raster remains
the default. `--intersection-refine-factor N` optionally enables smooth topology
voting for mixed faces: sampled occupancy decides whether to keep or discard the
existing smooth polygonized face, but it does not replace boundaries with raster
cells. `0` disables it and is the default.

The emitted G-code preserves the current viewer convention:

- `X`: axial position in millimeters, shifted so the minimum emitted X is `--axial-min`
- `Y`: layer-local rotary angle in degrees after `G92 E0 Y0`
- `Z`: radial lift from `bed_radius`
- comments: `;LAYER`, `;Z`, `;TYPE`

## Retraction calibration

Both slicers support stepped retraction distance or speed sweeps over radial
height above the mandrel. See [the calibration guide](../calibration/README.md)
for the OpenSCAD coupon, commands, band mapping, and scoring procedure.
`--retract-speed` is in mm/s and defaults to 30; v2 no longer reuses travel
feedrate for E-only retraction moves. `--hotend-temp` sets a fixed temperature.

## Verify

```bash
.venv/bin/python -m unittest discover cyl_slicer_v2/tests
```

The tests are staged to match the implementation pipeline: coordinates/config,
fixtures, occupancy sampling, periodic geometry, layer classification, toolpath
generation, G-code emission, and CLI smoke tests.
