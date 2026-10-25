# Cylindrical Bowden calibration

The photograph shows threads between features. It is consistent with oozing on
travel, but cannot establish the cause by itself. Excess temperature, wet filament,
Bowden tube movement, and poor restart pressure can all complicate retraction tuning.
[Prusa's stringing guide](https://help.prusa3d.com/article/stringing-and-oozing_1805)
describes temperature, retraction, and filament moisture checks. Bowden retraction
depends on the particular hotend and tube; use the hotend manufacturer's permitted
distance and speed range rather than treating a generic Bowden number as a target.

## Coupon

`retraction_posts.scad` has a curved base seated on a 30 mm radius mandrel and
four separated radial posts. It needs no long bridges or floating cubes. Mesh Z
is the cylinder axis, matching both slicers' STL convention. Output machine Z is
radial lift above the mandrel.

- `layout="axial"`: two posts separated along the cylinder axis.
- `layout="rotary"`: two posts separated around its circumference.
- `layout="both"`: four posts, the default; exercises mixed travel too.

The base is 2 mm thick. Five 3 mm bands give a total radial height of 17 mm.
Post width is 6 mm, axial centre spacing 22 mm, and arc centre spacing at the
mandrel 20 mm. All dimensions are editable with OpenSCAD `-D` arguments. Set
`bed_radius` to the actual mandrel radius, and use the same radius in the slicer.
Optional grooves are off by default to keep band geometry uniform.

```bash
openscad -o /tmp/retraction_posts.stl calibration/retraction_posts.scad
# Example isolation test:
openscad -D 'layout="axial"' -o /tmp/retraction_axial.stl calibration/retraction_posts.scad
```

## Test sequence

1. Check that both Bowden tube ends stay seated during retract/unretract. Use dry
   filament, a clean nozzle, and verify ordinary extrusion and base adhesion.
2. Establish a sensible temperature within the filament maker's range. If needed,
   print separate identical short coupons at temperatures 5 degrees apart with
   fixed retraction. Do not sweep temperature and retraction together.
3. Sweep distance with a fixed speed. The example below uses **0, 1, 2, 3, 4 mm**
   at **30 mm/s**; these are experiment values, subject to the hotend's limits.
   Zero is a deliberately unretracted control. Before the first band, the ordinary
   `--retract-mm` setting applies (2 mm in legacy; 1 mm in v2).
4. Choose the lowest distance that removes substantial strings while retaining
   complete starts after travel. Then hold that distance and sweep speed.
5. Confirm the selected pair on a fresh coupon with calibration disabled, then on
   the original model. Compare axial and rotary layouts if one travel direction
   is much worse. Record actual travel duration: the same feed setting need not
   produce the same travel time on linear and rotary axes.

Keep flow, temperature, fan, travel speed, Z-hop, and path planning fixed within
each experiment. Do not add wiping/coasting/extra prime while isolating retraction.
If all bands look equally bad, revisit filament, temperature, tube play, and
whether the generated travels actually retract before increasing the range.

For severe stringing across the entire 0–4 mm sweep at 210 C, compare two short
fixed-setting coupons at 210 C and 200 C (if both are within the spool's range).
Use `--retract-mm 4 --retract-speed 30 --max-layers 35`, keep the other slicing
options, and omit all `--calibration...` arguments. This prints about 7 mm of
radial height including the 2 mm base. Change only `--hotend-temp` between the
two files. If lowering temperature does not help, check tube seating during
retraction and compare with known-dry PLA before extending the distance range.
Also verify part cooling physically: send `M106 S255` with the printer idle to
start the normal Marlin part fan; `M107` stops it. The old v2 test files did not
include a fan-on command. New v2 output explicitly starts with the fan off, enables
full cooling at radial Z=2 mm (after this coupon's base), and switches it off at
shutdown. Override with `--fan-speed 0..255` and `--fan-start-height MM`.
If the previous print ran without cooling, repeat the original 210 C sweep with
cooling first, changing only that variable before comparing temperatures.

New v2 output ends with a 10 mm radial Z lift at 600 mm/min and waits for the
move to finish before disabling motors. Existing files need to be re-sliced.

## Slice with the legacy slicer

Run from the repository root:

```bash
.venv/bin/python slicer.py /tmp/retraction_posts.stl --no-plot \
  --out /tmp/retraction_distance.gcode \
  --calibration retract-mm --calibration-values 0,1,2,3,4 \
  --calibration-base-height 2 --calibration-band-height 3 \
  --retract-speed 30 --hotend-temp 210
```

The legacy mandrel radius remains `printer_params["bed_radius"]` in `slicer.py`.
The example uses a fixed 210 C for PLA; adjust to the range on your spool.
Without `--hotend-temp`, the existing default is 215 C.
Legacy also writes its normal `stats.txt`; `--out` selects only the G-code path.

## Slice with v2

```bash
.venv/bin/python -m cyl_slicer_v2.cli /tmp/retraction_posts.stl \
  --layer-backend intersection --build-cylinder-radius 30 --layer-height 0.2 \
  --out /tmp/retraction_distance_v2.gcode \
  --calibration retract-mm --calibration-values 0,1,2,3,4 \
  --calibration-base-height 2 --calibration-band-height 3 \
  --retract-speed 30 --hotend-temp 210
```

V2's intersection backend is experimental; it is useful here to preserve narrow
post geometry without coarse raster sampling. Preview the output with the existing
viewer and your normal printer startup procedure before printing. These commands
retain the slicers' ordinary startup procedure and explicitly select 210 C.

For a speed sweep, use either slicer with:

```text
--calibration retract-speed --calibration-values 20,25,30,35,40 --retract-mm 2
```

Replace `2` with the distance selected in the first test. CLI speeds and sweep
speed values are **mm/s**. G-code E-only feedrates are **mm/min**. V2 now has a
separate retraction speed (default 30 mm/s) instead of reusing travel feedrate.

## Read the bands

| Radial height above mandrel | Distance in example |
| --- | --- |
| Below 2 mm | Ordinary base setting |
| 2 to below 5 mm | 0 mm |
| 5 to below 8 mm | 1 mm |
| 8 to below 11 mm | 2 mm |
| 11 to below 14 mm | 3 mm |
| 14 to 17 mm | 4 mm |

Band boundaries use nominal shell radius minus mandrel radius, never mesh axial Z
or temporary Z-hop. `;CALIBRATION:` comments report each shell's height, band, and
value. Header comments also record the base retraction settings, fixed temperature,
mandrel radius, and complete sweep schedule. The final value is held if the model extends beyond the requested bands;
match SCAD `base_height`, `band_height`, and `bands` to the sweep when changing them.
Shell spacing quantizes transitions; choose band heights divisible by layer height.
Foundation shells below the first band retain the normal setting.

Score the middle of each band, allowing a few transition shells to settle:

| Band/value | Fine hairs (0–3) | Thick strings/blobs (0–3) | Gaps at restarts (0–3) | Notes |
| --- | --- | --- | --- | --- |
| | | | | |

Zero means none, three means severe. Photograph against a dark background under
the same lighting. Reject settings that produce missing starts, clicking, or
grinding, even if hairs decrease. In v2, a pending retraction is repaid using its
original amount and speed after a layer reset, even when the new band is zero.

The angular footprint is fixed, so rotary gaps and post widths grow with radius.
Cooling and travel duration can therefore change through the test. Use the sweep
to shortlist values, then compare fixed-setting coupons of identical geometry and
height to confirm the result.

## Software verification

```bash
.venv/bin/python -m unittest discover cyl_slicer_v2/tests
```

Calibration tests cover base/band boundaries, invalid ranges, speed conversion,
exact retraction repayment across G92 resets and a zero-distance band, and unchanged
deposited filament. They do not establish physical printer performance.
