# Cylindrical FDM G-Code Format

## Coordinate System

The G-code uses a non-standard coordinate system adapted for cylindrical printing:

- **X axis**: Axial position (Z in cylindrical coords), in millimeters, shifted so the minimum emitted X is the configured output minimum
- **Y axis**: Rotation angle in **degrees**, relative to the most recent `G92 Y0`

The printer rotates around the Y-axis, extruding radially outward. `E` and `Y` are reset on each layer change with `G92 E0 Y0` to keep accumulating axes close to zero. The slicer tracks physical rotary position internally, so emitted Y values after the reset are local offsets from the current physical angle.

V2 targets Marlin with `POWER_LOSS_RECOVERY` enabled. Each rotary reset waits
with `M400`, rebases native Y by whole turns with `G92.9 Y<wrapped phase>`, then
sets logical Y to zero with ordinary `G92`. Retaining angular phase prevents
fractional-step drift across layer resets. Ordinary `G92 Y0` alone does not reset native coordinates
in this firmware. Extra Y-only resets within large layers bound native and
emitted Y to ±3600 degrees without resetting E or moving the rotary motor.
See the [v2 precision notes](cyl_slicer_v2/README.md).

## PathObject

Represents geometry as `(theta, x)` pairs internally:

- `theta`: angle in radians
- `x`: axial position (mm)
- `closed`: whether path forms a closed loop
- `kind`: "polygon_exterior", "polygon_hole", "line", "point"
- `radius`: cylinder radius for arc-length calculations
- `is_accumulated`: if True, thetas already in accumulated printer space

## G-Code Commands

### Travel (no extrusion)
```gcode
G0 X{x:.3f} Y{accum_degrees:.3f} F{travel_mm_per_min:.1f}
```

### Extrusion
```gcode
G1 X{x:.3f} Y{accum_degrees:.3f} E{e:.5f} F{feed_mm_per_min:.1f}
```

### Retraction Sequence
```gcode
;start of path
G1 E{current_e - retract_mm:.5f} F{retract_speed:.1f}  ; retract
G0 X{x:.3f} Y{accum_degrees:.3f} F{travel_mm_per_min:.1f}  ; travel
G1 E{current_e:.5f} F{retract_speed:.1f}  ; unretract
```

### Comments
```gcode
;start of path
;end of path
```

Calibration mode adds `;CALIBRATION:` comments after each model layer's radius
marker. Band selection uses nominal radius minus build-cylinder radius, ignoring
temporary Z-hop. Distance values are mm of filament; speed values in the comments
are mm/s, converted to mm/min in E-only G-code. See
[the calibration guide](calibration/README.md).

## Layer/Type Markers

Standard comments for layer identification:
```gcode
;LAYER:0
;LAYER:support_gap
;TYPE:WALL_OUTER
;TYPE:WALL_INNER
;TYPE:SKIN
;TYPE:INFILL
;TYPE:SUPPORT
```

## Key Transformations

1. **Theta unwrapping**: Canonical theta values (0–2π) are unwrapped from the nearest equivalent angle to avoid full-turn travel for tiny backtracks
2. **Seam handling**: Full-circle paths are rotated to start at current position via `rotate_full_circle_path()`
3. **Path rotation**: Closed perimeters are rotated to nearest point via `rotate_closed_path_to_nearest()`
4. **Local seam closure**: Segments crossing the periodic seam are handled by traveling to the equivalent endpoint and extruding the short chord

## Arc Length Calculation

Real arc length between two points:
```
dist = sqrt((x2 - x1)^2 + ((t2 - t1) * radius)^2)
```

Where `t2 - t1` uses accumulated (unwrapped) theta deltas to avoid ~2πr artifacts at seam crossings.

## Feed Rate Scaling

Feed rate is adjusted per-segment to maintain constant linear speed in parameter space:
```
f = feed_mm_per_min * fake_dist / real_dist
```
Where `fake_dist` uses accumulated angular deltas (in degrees).
