# Cylindrical Slicer V2

This is a self-contained cylindrical slicing pipeline (triangle/cylinder intersection by default). It does not import or mutate
the legacy `slicer.py`/`togcode.py` pipeline.

## Quality defaults

`.venv/bin/python -m cyl_slicer_v2.cli part.stl --out part.gcode` now uses
the same quality-first defaults as `SlicerConfig()`. No quality flags are
required. These replace the old development/CLI defaults, so **reslice all old
G-code**; do not assume output from this version matches earlier jobs.

The default machine is a **30 mm radius mandrel, 0.4 mm physical nozzle,
1.75 mm PLA filament** at a 195°C setpoint. Geometry uses the
**triangle/cylinder intersection** backend with a 0.01 mm curve-error budget,
0.16 mm radial layers, 0.42 mm deposited line width and the rounded-bead
area/spacing model. Walls: four, outer-first. Top/bottom shells: five each.
Sparse grid: alternating directions at crosshatch-equivalent 1.5 mm spacing;
solid: alternating direction with meandering passes, including exposed top
surfaces. Small features are retained down to 0.1 mm and fine supports down to
0.25 mm². Support interface: three layers at 0.4 mm pitch and one-layer radial gap.

Printing targets 2700 mm/min (45 mm/s) physical print speed and
9000 mm/min (150 mm/s) surface travel, 8 s minimum estimated deposition time
per layer with 300 mm/min (5 mm/s) speed floor, full fan starting 1 mm above
mandrel with a further 1 mm fan ramp, 1.4 mm retraction at 20 mm/s. Actual
firmware acceleration and extrusion throughput may prevent reaching targets.

Use `--profile legacy-cli` to reproduce previous **CLI** defaults (0.5 mm
layers, rectangular flow, raster backend, 3 walls, etc.).
`--profile library` remains an alias for the new quality profile.
Explicit options override either profile; use `--monotonic-top-skin` to
enable one-way exposed top passes. Python can directly override fields:
`SlicerConfig(layer_backend="raster", monotonic_top_skin=True)`.

**Safety/calibration:** defaults represent the user's printer profile, not a
universal PLA recipe. Verify bed radius, nozzle, material, temperatures,
feed/acceleration limits, retraction and fan wiring before running new G-code.
Pressure advance (M900), strict intersection recovery, and fixed-position
seams remain opt-in because they require firmware calibration or depend on
the model. Support flow multipliers remain unity pending printed coupons.

## Experimental complete-terminal radial slabs

The optional `--radial-layer-mode complete-terminal` covers the requested mesh
radial range, including a fractional final layer ordinarily dropped by the
fixed-height scheduler. If the remainder is less than half the nominal layer
height, the final full slab and remainder are split into two equal-thickness
slabs, ensuring that the minimum printed slab is at least half-height. An
entire radial feature narrower than half a nominal layer is rejected rather
than silently depositing an unprintable sliver. The default remains
`--radial-layer-mode fixed` and preserves the previous sample locations.

Each slab is intersected at its own midpoint while the nozzle moves to its
outer radial face. Bead area, dense wall/skin spacing, variable-width paths,
and calibrated filament E scale consistently with that slab's actual height.
This is an **experimental prerequisite to adaptive layer selection**, not
slope-adaptive slicing yet. The approximate surface cross-section and contact
behavior still require real prints before promoting this mode to the default.
When a model floats above the mandrel, generated gap-support and raft
foundation slabs in this mode now partition the **actual free radial interval**
from the bed to the first occupied model slab's inner deposited face. They
rebalance sub-half-height remainders and carry individual physical thickness
through support toolpaths, nozzle positions, bead spacing and E. The foundation
ends at (but does not overlap) the model's inner face. An off-bed gap thinner
than half the nominal height is rejected rather than printed as an unsafe
sliver. The legacy fixed mode retains its original foundation scheduler.

This solves foundation overlap and terminal print-volume accounting, but
does not make arbitrary overhangs self-supporting. Independent model/support
contact gaps and material separation still require physical coupons, and
variable-height slope optimization remains experimental.

## Perimeter wall order

The default `--wall-order outer-first` preserves the historical dimensional
outline order. For overhangs where inside walls should support the outer shell,
`--wall-order inner-first` prints every island’s internal perimeters before
its external perimeter, then its infill. This is per-island ordering: the
slicer still completes each island before traveling to another. A one-wall
model behaves identically under either setting. This does not change wall
geometry, extrusion flow, winding, or the selected seam.

## Optional seam finish feed ramp

`--seam-end-ramp-mm MM --seam-end-ramp-fraction FRACTION` is an
**opt-in** closed-perimeter speed-shaping control for tuning seam starts/ends.
The default finish ramp length is zero (disabled), and ordinary open paths,
infill, support and bridges never receive a finish ramp. A short series of
Marlin mixed-axis speed-corrected G1 moves reduces commanded speed across
the last MM of each closed wall, with the final speed approaching FRACTION of
the normal role-limited feed. A configured start ramp may coexist with the
finish ramp; on very short loops they share the available length.

This does not coast, skip extrusion, offset E or change the selected seam,
wall winding, or variable-width bead area. Because changing feed may alter
nozzle pressure and thermal behavior, do not automatically enable it for
other machines. Measure axial, rotary and mixed wall loops and tune with
identical seam positions, material and temperature; pressure advance and
end speed shaping can interact. A slower finish is not necessarily a better
seam on every printer.

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
fan override, or acceleration.

### Optional M900 K radial-band calibration

After verifying that **the flashed firmware** enables Marlin `LIN_ADVANCE`,
you can slice a repetitive sharp-corner coupon with:

```bash
.venv/bin/python -m cyl_slicer_v2.cli calibration/retraction_posts.stl \
  --out /tmp/linear_advance_sweep.gcode --marlin-linear-advance \
  --calibration linear-advance-k --calibration-values 0,0.04,0.08,0.12,0.16 \
  --calibration-base-height 2 --calibration-band-height 3
```

The example coupon reuses the existing radial bands, but compare a dedicated
corner/acceleration coupon before committing a K value. The five K values are
**illustrative**, not material presets. A verified firmware capability flag is
mandatory, even for all-zero sweeps. Without `--linear-advance-k`, the
calibration sets explicit base K0, emits the next M900 only when the radial
band changes, then resets to K0 at the end of the print. An explicitly set
base K can also be supplied. Retraction sweeps remain separate.

Keep material, nozzle, temperature, path width, feed, cooling, retraction,
pressure transient layout and mandrel phase unchanged between compared
bands; print axial, rotary and mixed-axis features at realistic speeds.
M900 commands never adjust filament E, XY geometry or Marlin feed scaling.
Do **not** use this on firmware that lacks M900 support; the slicer cannot
inspect the controller or verify what was flashed. Calibrate K separately for X-only, rotary-only,
and mixed axial/rotary moves at real print speeds on the actual firmware,
and confirm retract behavior. Extrusion-rate smoothing, if needed, requires
a separate volume-conserving motion planner and is not enabled here.

## Layer cooling controls (opt in)

Without options, the original fan step and full-speed motion are unchanged.
`--cooling-min-layer-time-s T` computes deposition-only time from all
cylindrical toolpath lengths in physical millimeters, then caps *that layer's*
print speed to approach at least T seconds. This is conservative: travel,
retraction and acceleration also consume time but are deliberately excluded
from the estimate. Layer and connected-island timing use the **actual print-role speed
ceilings** (bridge, support and support interface) so slower material is
credited toward the cooling target instead of unnecessarily throttling all
other paths. `--cooling-min-print-speed` (default 300 mm/min) bounds
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

The optional `--sparse-infill-pattern both` crosshatches radial and axial passes **on every
layer** at the specified `--infill-spacing S`. In alternating mode the
slicer prints radial-only passes on even model layers and axial-only passes on
odd layers at *S/2*. This keeps approximate deposited line length (and thus
material volume) per layer comparable to the old two-direction crosshatch,
instead of silently halving infill density. Clipped-edge effects can cause
minor differences. Because the same physical material volume is concentrated
into one direction, alternation also avoids same-layer grid intersections.

A check rejects alternating sparse infill when S/2 would be smaller than the
calibrated dense-bead pitch; use solid infill instead of crossing the strands
at a physically overfull spacing. The configured S is the **crosshatch
equivalent spacing**, not the actual single-family pitch in alternate mode.

`--solid-infill-pattern alternate` rotates dense skin/infill lines between
axial (even) and radial (odd) shells, retaining the unified bead-width/pitch
calibration. Default `alternate` provides cross-layer reinforcement. These settings operate
at the geometric toolpath generation stage, so periodic winding, wall order,
retraction and firmware speed correction are unchanged. Monotonic top-skin
ordering/seam tuning remains a separate task (#12/#13).

## Monotonic exposed top skin

The opt-in `--monotonic-top-skin` separates the **exposed** upper-surface mask from
buried bottom solid infill, then prints its dense lines in a fixed physical
direction, progressing in stripe order. Unlike nearest-path or serpentine
printing, the nozzle never reverses during deposited top-skin stripes.
The travel back to the next stripe remains non-extruding and uses the same
island/seam-aware retract and Z-hop rules. Default top fill meanders instead.

Thin model regions marked both top and bottom use the top mask once; there
is no intentional double-extrusion. The same calibrated bead area and dense
pitch apply to both skin and buried fill. Use `--no-monotonic-top-skin` to disable it. This is not yet a variable-width
skin/gap fill process and still needs surface/flow calibration coupons.

## Predictable perimeter seams

For closed perimeters only, `--seam-policy nearest` preserves the existing
nearest-point entry. `--seam-policy aligned --seam-angle-deg 180` places the
start at the contour point nearest the selected cylinder phase (180° is the
default target); if the ray intersects an edge, it can split the edge without
changing its geometry or extrusion length. `--seam-policy corner` prefers
the sharpest actual polygon corner, with the target angle and travel distance
as tie-breakers. These policies are especially useful for stacking starts over
successive layers or hiding them at corners.

`--seam-exclude-deg START END` excludes start angles in the half-open cylinder
window [START,END) measured in degrees; START greater than END means the
window crosses the 0/360° rotary seam. When every candidate is excluded the
slicer raises an error instead of silently putting a seam in the forbidden
region. This applies to already closed perimeter loops only; open/wraparound
surface walls retain their original winding. No automatic coasting/restart
over-extrusion is enabled: calibrate pressure advance and retraction on the
machine first.

## Support and interface motion/flow tuning (opt in)

For a separately calibrated support recipe, set
`--support-print-feed` and `--support-interface-print-feed` in physical
mm/min. These are *maximum* role speeds, capped by any layer-wide cooling
slowdown. Per-motion rotary/axial Marlin F conversion remains in effect.

`--support-flow-ratio` and `--support-interface-flow-ratio` separately scale
the deposited filament on support and top support-interface toolpaths. Defaults
are 1.0, keeping the shared bead model and previous G-code unchanged. The
geometry, model walls, skin, line spacing, raft foundations, support gaps and
retract distance are unaffected. A flow ratio is **not** a substitute for
independently tuning interface density, contact gap, support temperature and
cooling on the actual printer. Avoid exceeding nozzle volumetric flow capacity.

## Bounded on-island travel routing

The quality profile now enables `--route-travel-inside-island`. When a
straight non-extruding move would cut across a model hole or leave a concave
island, the planner first tries following a short arc of its original boundary,
then a bounded local visibility graph in **physical** rotary/axial millimeters,
including periodic seam copies. Routing is allowed only when both endpoints lie
in the same island. Only vertices reachable within the travel distance budget
enter the graph, so distant outline detail does not prevent local routing.
The graph and returned boundary routes have at most 96 vertices, including the
start and destination. Every segment is checked against the original footprint.
A detour must be at
most 1.65 times the direct distance and at most 6 mm longer. Above those
limits, or with no verified inside route, the slicer retains its existing
straight travel/retract/Z-hop safety behavior. A successfully routed short
move may remain primed only if its full travel is ≤5 mm; longer detours
still retract. Straight support-to-support and support-interface connections
within one footprint can also remain primed up to that 5 mm limit. These are
deliberately bounded for performance and ooze.

Within each support island and role, open meanders may reverse while retaining
their segment widths and deposited filament. Ordering checks at most eight
nearby entries, preferring short verified travel over a nearer move across a
void. Model island order, wall winding, and configured perimeter seams remain
unchanged.

`--no-route-travel-inside-island` disables it, and
`--route-max-extra-mm`, `--route-max-factor`,
`--route-max-vertices`, and `--route-max-unretracted-mm` allow tuning.
`;ROUTED_TRAVEL:N` comments identify paths actually routed. Visibility
segments must remain inside the printed island footprint. The router does
not navigate disconnected islands, move at extrusion Z through empty air,
or change E, perimeter winding, or part geometry.

Retraction checks and emitted travel use the same physical island boundaries
and rotary coordinate frame, including after negative turns or rotary resets.
A 0.0000001 mm boundary tolerance handles floating-point drift along edges;
short connections must still fit inside a single island. Straight contained
travel bypasses visibility routing. Physical footprints are cached with a
bounded layer cache so repeated boundary passes do not rebuild the geometry.
Each cached physical island also retains at most 64 routing results, shared by
ordering, retraction checks, and emitted travel.

## Partial-wall bearing report (diagnostics only)

Use `--wall-support-report` to emit `;WALL_SUPPORT:...` comments beside
external/internal wall paths. For each printed wall it calculates the fraction
of the physical R×theta/axial perimeter centerline supported by the preceding
sampled layer (expanded by half a deposited line width), and the longest
uncovered continuous span in physical millimeters. A partly supported path is
marked `partial`; fully supported or entirely unsupported walls are reported
separately. It accounts for winding across the 0/360° seam, including an
uncovered span that straddles a closed loop's start/finish point.

This does **not** change walls, material, E, support generation, or the
`--bridge-mode` safety gate. A geometrically unsupported wall segment is not
automatically a safe gravity-independent bridge on a spinning cylinder. First
model layer diagnostics are omitted when the mandrel itself is the bearing.
The report is off by default to avoid additional geometry processing and
G-code noise. Print coupons at different rotary phases are still necessary
before exploiting bearing diagnostics for speed, bridging or support removal.

## Two-anchor bridge printing (experimental opt-in)

`--bridge-mode print` detects short, nearly straight **two-bearing**
unsupported corridors in the physical cylindrical distance metric. Both
anchor ends must meet the previous deposited layer; a one-ended cantilever,
a long span above `--bridge-max-span-mm` (default 8 mm), or a branched/irregular
gap retains normal support instead. Detection spans the periodic 0/360°
seam. Confirmed bridge infill is laid along the anchored direction at the
calibrated dense bead pitch, including a small `--bridge-anchor-overlap-mm`
(default 0.2) into solid bearing material, and excluded from other fill so
filament is not intentionally printed twice. Support propagation omits
only bridge regions; nearby overhangs retain their support.

`--bridge-print-feed` sets a speed limit in physical mm/min (default 900),
`--bridge-flow-ratio` scales bridge extrusion, and `--bridge-fan-pwm`
overrides the part fan for the bridge and restores it immediately afterward.
All are **opt-in via bridge mode**. The default `--bridge-mode off`
continues conventional support: bridging under a rotating mandrel depends
on gravity, radial phase and cooling hardware, not on horizontal-bed
heuristics. Verify bridge strips at relevant rotary phases before enabling
this in production. Partially supported walls, sag prediction, and variable
bridge bead width remain outside the conservative detector.

## Narrow-feature variable-width gap fill (quality default)

The slicer recovers conservative thin strips and tapered serif-like edges
whose fixed-width perimeter offsets disappear. It analyzes **unprinted
edge-band** material rather than intentionally sparse infill, measures
cross-section widths in physical R×theta/axial millimeters, and emits one
GAP_FILL centerline with a calibrated width and E for every segment. Default
minimum deposited gap width is 0.22 mm (at 0.16 mm layer height), bounded
above by the 0.42 mm nominal bead. Short, branched, sharply curved and
sub-minimum features are rejected rather than fabricated across empty space.
Gently curving narrow strips can now pass the local chord/inside-swept-footprint
checks even if their overall rotated bounding box is wider than a bead.
An area bound excludes broad plates and folded centerlines before sampling;
no extra plastic is generated across a separate island or sparse interior.

To reproduce the old slicer behavior use `--no-gap-fill-enabled` (also the
explicit legacy profile). Tune `--gap-fill-min-width-mm`,
`--gap-fill-sample-step-mm`, and `--gap-fill-min-length-mm`.
This is a staged conservative gap recovery, not a full Arachne/medial-axis
solver for complex curves or variable-width multiple walls.

## Physical support contact clearance (opt in)

The default `--support-gap-layers 1` remains exactly as before. To specify a
**physical radial** contact clearance instead, use
`--support-contact-gap-mm N`. The support generator now compares actual
radial *faces*: it permits a support slab only if its outer deposited face
is at least N mm below the originating model slab's inner face. This
uses the individual recorded slab thickness, including a variable-thickness
terminal layer, rather than converting N to nominal layer counts. The
printed gap is quantized upward by the real intervening slabs, never below
the requested clearance. At fixed 0.16 mm height, 0.17 mm still means two
empty layers. Passing 0 requests no empty radial slab; if unset, the old
integer `--support-gap-layers` setting controls.

Independent overhangs retain their own contact limits as the support
foundation propagates downward. This avoids prematurely printing one
island because a nearby overhang originates at a different radial height.
Clearance is radial-only; physical removability still requires print coupons.

`--support-lateral-gap-mm N` optionally removes ordinary support tracks
from an additional same-layer physical-mm buffer around model material,
measured in the unwrapped cylindrical R×theta/axial plane. It changes
printed support contact near model sidewalls, **not** the carried foundation
footprint on earlier layers; zero is the default and preserves geometry.
These are material/removability parameters, not universally safe automatic
quality improvements. Validate physical support removal and underside finish
before using non-default contact gaps. Fully adaptive radial sampling, variable-thickness gap-support/raft foundations,
bottom-interface patterns and printed underside contact validation remain
separate engineering tasks.

## Cooling per connected island (quality default)

The 8-second **layer** print-time target remains enabled. To avoid a tiny
feature being printed quickly just because a different large island fills
the layer, the quality profile additionally enables
`--cooling-min-island-print-time-s 2`. Each island's print-only physical
surface length and role-specific speed limits are used to set the
**fastest** permissible feed that reaches at least the requested deposition
time; other islands keep their existing layer-capped feed. A 300 mm/min
minimum speed still takes priority if the target cannot be reached. It never
adds stationary G4 commands over molten plastic and never modifies E,
geometry, bridge-flow multipliers or the temporary bridge-fan override.

`--cooling-min-island-print-time-s 0` disables per-island slowdown.
`--profile legacy-cli` also disables it. The `;COOLING_ISLAND:` G-code
annotations report conservative deposition-only duration estimates and
effective physical speed. They **do not guarantee** a between-layer
revisit/cooling interval: physical acceleration, motion blends, travel,
retractions, material behavior and physical tests remain important. Island
revisit scheduling and material-specific fan calibration remain open under
issue #9.

## Opt-in extrusion feed ramp (issue #8)

`--extrusion-feed-ramp-mm N` gradually raises commanded *physical speed*
from `--extrusion-feed-ramp-start-fraction 0.65` to the existing print
speed over the first N physical mm of each extrusion path. The slicer
inserts no more than four bounded ramp stages without changing original
path endpoints, winding, calibrated per-segment variable-width bead
areas, or total commanded filament E. It continues converting physical
feed to the Marlin X-mm/Y-degrees planner metric. Zero (default) leaves
all deposition moves unchanged; layer/island cooling and bridge/support
role caps still apply.

This is feed shaping, **not** pressure compensation or an E-volume ramp.
Firmware acceleration and junction blending may limit its impact.
Test X-only, rotary-only, mixed-direction and restart coupons before
enabling; Marlin M900 remains independently opt-in only when enabled in
the flashed firmware.

## Optional material fan ceilings

`--material-fan-policy preset` caps the user-requested fan PWM by material,
without increasing a user-specified value: PLA/unprofiled 255, PETG/TPU 180,
ABS 0. The fan start height and radial ramp still apply. This is an
**opt-in starting policy**, not a universally safe cooling recipe; enclosure
conditions, printer air ducts, part geometry, and filament grade must be
calibrated on the actual machine. The default `manual` policy preserves
existing cooling behavior and the PLA quality profile.

Example: `--material-profile petg --material-fan-policy preset`.
Bridge-specific fan overrides remain independent and are restored to the
material-capped baseline after each bridge. This change does not select
pressure advance, change the filament E length, or insert waiting dwell.

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

For the optional **raster** geometry backend, the mesh occupancy implementation
calls `trimesh.contains` directly. A faster
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
spacing defaulting to the selected bead model's dense line spacing. Change these
with `--raft-layers N` and `--raft-grid-spacing S`. Small angular gaps between
projected support strips merge into the raft when they measure at most 2 mm
at the innermost raft sample radius. This reduces fragmented first-layer
passes without filling larger angular cutouts or the configured axial relief
gap. Set `--raft-merge-gap-mm D` to change this distance, or `0` to retain all
angular gaps (`legacy-cli` defaults to `0`). Straight support and interface passes run
back and forth, keeping adjacent passes primed when the whole connecting travel
stays inside their footprint.

Above the raft, nearby support branches also merge into shared sparse bodies
instead of many separate columns with repeated meander walls. The quality
default `--support-merge-gap 2` closes gaps in the combined support footprint
in physical rotary/axial mm at each layer's radius. Set it to `0` to retain
separate branches, or choose another distance explicitly. Larger openings
remain open. The merge adds connecting material while retaining original
branches; it does not make the sparse support fill fully dense. Model/lateral
clearance and pending radial contact gaps are removed after merging. Dense
interface pads and their independent contact heights keep their original masks.

Extrusion is derived from the selected bead cross-section and filament area.
For the rounded default: `((w-h)*h + pi*h*h/4) /
(pi*(filament_diameter/2)^2)` filament-mm per toolpath-mm. The default
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

The intersection backend is the quality profile default. Raster remains
available with `--layer-backend raster` for debugging or compatibility. `--intersection-refine-factor N` optionally enables smooth topology
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

Sparse radial and solid infill can print in either rotary direction. Adjacent passes start
at the nearer endpoint, so they alternate positive and negative Y with a short
axial travel between rows. Reversed moves retain their intended angular span,
including across the periodic seam.
Exposed monotonic top skin retains its one-way passes.

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
At startup, `M400`, `G92.9 Y<wrapped phase>`, and `G92 E0 Y0` establish a
bounded native rotary coordinate. Ordinary layer changes use `G92 E0` only;
the Y frame remains continuous across layers. A Y reset is made only before a
move that would exceed the ±3600° precision budget. That reset uses `M400`,
`G92.9 Y<wrapped phase>`, and `G92 Y0` without rotating the motor. Retaining
the angular phase avoids accumulating fractional motor-step rounding.
Ordinary `G92 Y0` alone only changes the workspace in this firmware and does not
bound its internal float32 planner position. See Marlin's
[G92 implementation](https://github.com/MarlinFirmware/Marlin/blob/2.0.7.2/Marlin/src/gcode/geometry/G92.cpp).

Additional Y-only resets happen before any command that would exceed ±3600
degrees, whether within a layer or after a layer change. At that magnitude
float32 spacing is approximately 0.000244 degrees, below the 0.001-degree
output resolution. Moves too large to fit after a reset are rejected before
writing G-code. Mid-layer resets preserve
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

`--print-feed` and `--travel-feed` specify physical nozzle surface-speed
targets (quality defaults 2700 and 9000 mm/min). Marlin has mixed X-mm/Y-degree
axes; the slicer converts each move's F to account for its **physical** motion
length, but firmware acceleration and axis limits still constrain actual speed.
Radial Z hops and end parking use native mm/min.
For stringing experiments, try faster travel while holding print feed fixed.

Part cooling uses standard Marlin `M106` commands: off at startup, full speed
ramped over 1 mm of radial height starting 1 mm above the mandrel in the
quality profile, and off at shutdown.
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
