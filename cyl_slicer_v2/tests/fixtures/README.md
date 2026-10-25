# Regression geometry

- `cubeandring2.stl`: the existing intersection/layer-classification fixture.
- `cubeandring2_cutout.stl`: the ring and rotated cubes with the through cutout
  from `cubeandring2.scad`, saved for the first-layer regression at 0.16 mm.
- `cubeandring2_cubes.stl`: the rotated-cubes model without the ring base, saved
  from the repository version used by the intersection and support regressions.
  These layer-specific expectations must not depend on edits to the working
  `cubeandring2.stl` model.
- `cubeandring2_top_surfaces.wkt`: two valid polygons reduced from the exposed
  surface of the ring-base model at radius 54.634266 mm. A zero-width buffer of
  the MultiPolygon discards the large printable patch; buffering each component
  separately does not. Keep the coordinates at full precision.
