# Slicer development

- Treat SlicerConfig() as the canonical, quality-first default machine profile.
  The CLI must resolve unspecified settings from it; don't recreate
  divergent or outdated defaults in argparse. Keep explicit CLI overrides,
  quality/legacy profiles, and default-vs-explicit regression tests aligned.
- Prefer quality-improving algorithms as defaults when they preserve geometry
  and do not require calibration. Keep firmware commands, machine-specific
  tuning and potentially topology-changing experimental recovery opt-in.
- Work on the v2 slicer in `cyl_slicer_v2/` by default. The root-level
  `slicer.py` and `togcode.py` are legacy; do not change or use them to solve
  slicer issues unless the user explicitly requests legacy work.
- Run v2 with `.venv/bin/python -m cyl_slicer_v2.cli`.
- Verify slicer changes with `.venv/bin/python -m unittest discover cyl_slicer_v2/tests`.
- Machine axes are X = axial mm, Y = rotary degrees, Z = radial height above
  the build cylinder. Preserve rotary winding, geometry, and deposited filament
  when changing path planning or G-code output.
- Geometry sampling resolution must not force redundant G-code blocks on a
  straight move. Short travels along a perimeter are not perimeter crossings.
- Adjacent fill passes should stay primed for short travel contained in the
  model/support footprint. Check the whole connection, including across the seam;
  don't suppress retraction across holes or separate islands based on distance
  alone. A role/type change is not by itself a reason to retract.
- Finish each model island's walls and infill before moving to another island.
  Join periodic seam components for ordering; keep support components separate.
  Enter closed perimeters using the configured seam policy (nearest by default),
  preserving winding and deposited length. Keep intentional zero-distance calibration bands intact.
- Keep Marlin's native Y position bounded as well as emitted Y values. The target
  firmware enables `POWER_LOSS_RECOVERY`: use `M400`, `G92.9 Y<wrapped phase>`,
  then `G92 Y0` (or `G92 E0 Y0` at layer changes). Keep angular phase so resets
  don't accumulate fractional-step drift. Ordinary `G92 Y0` alone only shifts its
  workspace. Preserve the Marlin float32/planner regression tests and bidirectional
  infill when changing rotary coordinates.


test with ```PYTHONPATH=. .venv/bin/python3 -m pytest -n auto --timeout=30```