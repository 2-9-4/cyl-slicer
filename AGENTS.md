# Slicer development

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
