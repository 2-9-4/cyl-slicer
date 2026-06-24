# Cylindrical Slicer V2

This is a self-contained raster-based prototype. It does not import or mutate
the legacy `slicer.py`/`togcode.py` pipeline.

## Run

```bash
.venv/bin/python -m cyl_slicer_v2.cli --fixture tube --out /tmp/tube_v2.gcode --max-layers 2
.venv/bin/python -m cyl_slicer_v2.cli tube.stl --out /tmp/tube_stl_v2.gcode --max-layers 2
```

The emitted G-code preserves the current viewer convention:

- `X`: axial position in millimeters
- `Y`: accumulated rotary angle in degrees
- `Z`: radial lift from `bed_radius`
- comments: `;LAYER`, `;Z`, `;TYPE`

## Verify

```bash
.venv/bin/python -m unittest discover cyl_slicer_v2/tests
```

The tests are staged to match the implementation pipeline: coordinates/config,
fixtures, occupancy sampling, periodic geometry, layer classification, toolpath
generation, G-code emission, and CLI smoke tests.
