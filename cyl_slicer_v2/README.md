# Cylindrical Slicer V2

This is a self-contained raster-based prototype. It does not import or mutate
the legacy `slicer.py`/`togcode.py` pipeline.

## Machine and default profiles

The CLI retains the original `legacy-cli` defaults to avoid silently changing
existing prints: three walls, 30 mm build radius, and strong support filtering
and merging. `--profile library` instead selects the smaller support thresholds,
two walls, and zero build radius of bare `SlicerConfig()`. Explicit arguments
override either profile.

`--nozzle-diameter` specifies the physical nozzle orifice, while
`--line-width` specifies the deposited bead width. They are independent.
Supplying a nozzle diameter validates that layer height does not exceed the
physical orifice. Without a specified nozzle, hardware is marked `unknown`
in G-code and the old settings remain available for compatibility; this is
not evidence that the historic h0.5/w0.4 configuration is printable.
For a 0.4 mm nozzle, use `--nozzle-diameter 0.4 --layer-height 0.2` with a
separately calibrated line width.

## Perimeter wall order

The default `--wall-order outer-first` preserves the historical dimensional
outline order. For overhangs where inside walls should support the outer shell,
`--wall-order inner-first` prints every island’s internal perimeters before
its external perimeter, then its infill. This is per-island ordering: the
slicer still completes each island before traveling to another. A one-wall
model behaves identically under either setting. This does not change wall
geometry, extrusion flow, winding, or the selected seam.

## Optional firmware linear advance

The default G-code emits **no M900 commands**. The previously reviewed Marlin
source has `LIN_ADVANCE` disabled, so simply changing the slicer is not enough.
Verify that the **flashed** firmware actually supports `M900` before using
`--marlin-linear-advance --linear-advance-k K`. With both explicitly supplied,
the slicer issues `M900 K...` after heating, before printing, and `M900 K0`
at job end. It does not modify calculated extrusion distance, retranslations,
or the physical X/Y travel-speed correction.

`--material-profile pla|petg|abs|tpu|unprofiled` records the identity in
G-code; it deliberately assigns **no guessed** K factor, nozzle temperature,
fan override, or acceleration. Calibrate K separately for X-only, rotary-only,
and mixed axial/rotary moves at real print speeds on the actual firmware,
and confirm retract behavior. Extrusion-rate smoothing, if needed, requires
a separate volume-conserving motion planner and is not enabled here.

## Layer cooling controls (opt in)

Without options, the original fan step and full-speed motion are unchanged.
`--cooling-min-layer-time-s T` computes deposition-only time from all
cylindrical toolpath lengths in physical millimeters, then caps *that layer's*
print speed to approach at least T seconds. This is conservative: travel,
retraction and acceleration also consume time but are deliberately excluded
from the estimate. `--cooling-min-print-speed` (default 300 mm/min) bounds
the slowdown and takes precedence if T is unattainable; **no stationary nozzle
dwell is added**. G-code comments report chosen print speed and estimated
deposition time. Per-motion Marlin F conversion still uses the existing
physical surface-speed correction. Extrusion E remains strictly conserved.

`--fan-ramp-height-mm H` gradually ramps M106 PWM from zero at the already
configured fan-start height up to `--fan-speed` after H radial millimeters.
Zero retains the previous binary fan step. This is a radial-height profile,
**not** measured thermal cooling. Individual island revisit times, bridge
fan overrides, acceleration-aware timing and actual filament temperature
are still open validation/improvement items.

## Alternating infill (opt in)

The original sparse pattern crosshatches radial and axial passes **on every
layer**. `--sparse-infill-pattern alternate` deposits radial-only paths on
even model layers and axial-only on odd layers, reducing same-layer grid
crossings without assuming support that is not present. User-selected
`--infill-spacing` continues to control physical pitch.

`--solid-infill-pattern alternate` rotates dense skin/infill lines between
axial (even) and radial (odd) shells, retaining the unified bead-width/pitch
calibration. Default `axial` preserves current output. These settings operate
at the geometric toolpath generation stage, so periodic winding, wall order,
retraction and firmware speed correction are unchanged. Monotonic top-skin
ordering/seam tuning remains a separate task (#12/#13).

## Monotonic exposed top skin

`--monotonic-top-skin` separates the **exposed** upper-surface mask from
buried bottom solid infill, then prints its dense lines in a fixed physical
direction, progressing in stripe order. Unlike nearest-path or serpentine
printing, the nozzle never reverses during deposited top-skin stripes.
The travel back to the next stripe remains non-extruding and uses the same
island/seam-aware retract and Z-hop rules. Default behavior is unchanged.

Thin model regions marked both top and bottom use the top mask once; there
is no intentional double-extrusion. The same calibrated bead area and dense
pitch apply to both skin and buried fill. This is not yet a variable-width
skin/gap fill process and still needs surface/flow calibration coupons.

## Run

```bash
.venv/bin/python -m cyl_slicer_v2.cli --fixture tube --out /tmp/tube_v2.gcode --max-layers 2
.venv/bin/python -m cyl_slicer_v2.cli tube.stl --out /tmp/tube_stl_v2.gcode --max-layers 2
```

The physical deposition model is selected with `--bead-model`. The default
`rectangular` mode preserves the historical `width * layer_height` area and
full-width dense spacing. For calibrated printing with `layer_height <=
line_width`, use `--bead-model rounded --layer-height 0.2 --line-width 0.4`.
It models a flat-sided bead with round ends:

- area = `(w-h)*h + pi*h*h/4` mm² (extrusion mm filament/mm path
  is area divided by filament cross-sectional area);
- dense line and inner-wall center spacing = `area/h`;
- external perimeter centerline inset = `w/2`;
- inner fill region begins `w + (wall_count-1)*(area/h) - overlap` mm inside
  the original boundary.

The spacing and flow change together; switching only one creates under/overfill.
Sparse infill's explicit `--infill-spacing` still sets user-selected density.
Actual nozzle orifice diameter is not the same as deposited line width.
Existing slicing profiles retain the rectangular default until hardware
dimensions are known. G-code reports the chosen bead model and dense pitch.

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

Degenerate face/cylinder recovery checks nearby samples from both
sides of the requested radius, bounded by 10% of deposited layer thickness.
Strict mode additionally caps that shift to the intersection curve-error
budget. Compatibility mode explicitly warns if it exceeds that budget.
A two-sided match in island/hole topology
and physical boundary position is reported as verified. True radial
tangencies can cause topology to change exactly at the sample radius; the
default preserves bounded one-sided recovery for compatibility, but emits a
`RuntimeWarning` explicitly saying the substitution is unverified.
Use `--intersection-strict-recovery` to reject these ambiguous topology
transitions instead of printing them. Warnings include the radius shift and
boundary discrepancy; they are diagnostic output rather than changes to the
requested sample radius in G-code.

The intersection backend is experimental and intentionally small. Raster remains
the default. `--intersection-refine-factor N` optionally enables smooth topology
voting for mixed faces: sampled occupancy decides whether to keep or discard the
existing smooth polygonized face, but it does not replace boundaries with raster
cells. `0` disables it and is the default.

The emitted G-code preserves the current viewer convention:

- `X`: axial position in millimeters, shifted so the minimum emitted X is `--axial-min`
- `Y`: rotary degrees since the most recent Y reset (normally each layer)
- `Z`: radial nozzle clearance from `bed_radius`, at the OUTER surface
  of each deposited slab (first layer: Z = configured layer height for a
  sample at bed radius + half a layer). A slice's internal radius is its
  midpoint sample plane; `;SAMPLE_RADIUS` and `;DEPOSITED_THICKNESS` describe
  these distinct quantities in the output. Explicit machine Z calibration is
  independent of geometric slab sampling.
- comments: `;LAYER`, `;Z`, `;TYPE`

Straight moves are emitted as single extrusion blocks, including axial and rotary
moves. Redundant collinear mesh vertices are removed with a numerical tolerance
of 1e-8 mm; corners, reversals, and rotary winding are preserved. Intersection
face/cylinder curves are sampled to a physical unwrapped chord-error budget:
`--intersection-curve-error-mm` (default 0.025 mm), with
`max_axial_segment` an additional curved-face sample cap. Exact straight
surface paths do not need intermediate curve samples. These are sampling
controls, not a G-code block-length limit.

Sparse radial infill can print in either rotary direction. Adjacent passes start
at the nearer endpoint, so they alternate positive and negative Y with a short
axial travel between rows. Reversed moves retain their intended angular span,
including across the periodic seam.

Each model island is finished locally: external perimeter, internal perimeters,
then its solid or sparse infill. The next island is chosen from the current
nozzle position. Components connected across the rotary seam belong to one
island; support components are grouped separately so a support bridge does not
merge distinct model islands. Closed perimeters start at the nearest point on
their contour, including the middle of an edge, rather than their original
mesh vertex.

Short transfers between neighboring walls or fill passes stay primed when the entire
travel lies inside the layer's model or support footprint. The local distance
allowance is the larger pass spacing plus half a nozzle width. This also applies
across infill type changes; ending a type does not itself force retraction.
Longer travels and travel through gaps retain the ordinary retraction rule,
and layer ends still retract. Periodic copies of the footprint allow the same
decision across the rotary seam.

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

`--print-feed` and `--travel-feed` expose the commanded G-code F values (defaults
1800 and 3000). Travel includes radial Z hops; end parking stays at F600. These
are the existing feed settings, with no new radius-based scaling. X is millimeters
and Y is degrees, so F does not guarantee a particular physical nozzle speed on
rotary or mixed moves. Firmware acceleration and axis limits also constrain speed.
For stringing experiments, try faster travel while holding print feed fixed.

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
or
```bash
# Runs tests in parallel using all CPU cores, terminating any individual test that exceeds 10 seconds
PYTHONPATH=. pytest -n auto --timeout=30
```



The tests are staged to match the implementation pipeline: coordinates/config,
fixtures, occupancy sampling, periodic geometry, layer classification, toolpath
generation, G-code emission, and CLI smoke tests.
