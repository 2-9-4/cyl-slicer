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
- `Y`: rotary degrees since the most recent Y reset (normally each layer)
- `Z`: radial lift from `bed_radius`
- comments: `;LAYER`, `;Z`, `;TYPE`

Straight moves are emitted as single extrusion blocks, including axial and rotary
moves. Redundant collinear mesh vertices are removed with a numerical tolerance
of 1e-8 mm; corners, reversals, and rotary winding are preserved. The internal
`max_axial_segment` setting controls intersection curve sampling only, rather than
limiting G-code block length.

Sparse radial infill can print in either rotary direction. Adjacent passes start
at the nearer endpoint, so they alternate positive and negative Y with a short
axial travel between rows. Reversed moves retain their intended angular span,
including across the periodic seam.

V2 targets the printer's Marlin 2.0.7 firmware with `POWER_LOSS_RECOVERY` enabled.
At startup and each layer, `M400`, `G92.9 Y<wrapped phase>`, and `G92 E0 Y0`
rebase the native rotary position by whole turns and reset logical Y without
rotating the motor. Retaining the angular phase avoids accumulating fractional
motor-step rounding at each layer's zero.
Ordinary `G92 Y0` alone only changes the workspace in this firmware and does not
bound its internal float32 planner position. See Marlin's
[G92 implementation](https://github.com/MarlinFirmware/Marlin/blob/2.0.7.2/Marlin/src/gcode/geometry/G92.cpp).

Additional Y-only resets happen inside large layers before the next command
would exceed ±3600 degrees. At that magnitude float32 spacing is approximately
0.000244 degrees, below the 0.001-degree output resolution. Moves too large to
fit after a reset are rejected before writing G-code. Mid-layer resets preserve
E, retraction state, and physical rotary winding. This requires `G92.9` support;
it is not a generic profile for firmware that lacks that Marlin feature.

`test_marlin_rotary_precision.py` replays emitted commands using Marlin's
float32 workspace math and signed 32-bit planner steps at the configured
52.7777778 steps/degree. It covers a print exceeding 80 million physical degrees,
large individual layers, reversals, small moves, existing workspace offsets,
and oversized-move rejection. It models these numeric operations rather than
executing a complete firmware image.

Z-hop defaults to 0.4 mm for travels of at least 5 mm or shorter travels crossing
a perimeter. Traveling along a perimeter (such as the 0.4 mm step between adjacent
solid infill lines) does not count as crossing it. Use `--z-hop-min-travel`,
`--z-hop-height`, and `--no-z-hop-perimeter-crossing` to adjust this policy.

## Retraction calibration

Part cooling uses standard Marlin `M106` commands: off at startup, full speed
once nominal radial Z reaches 2 mm above the mandrel, and off at shutdown.
Set `--fan-speed 0` to disable cooling, or choose PWM 0–255. Set
`--fan-start-height` to change the starting radial height; Z-hop does not trigger
cooling. For PLA, check that `M106 S255` actually starts the part fan on your printer.

At the end of a non-empty print, v2 raises radial Z by 10 mm from the final
position at 600 mm/min. It waits for the lift to finish (`M400`) before shutting
down the heater, fan, and motors. Re-slice existing G-code to include this lift.

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
