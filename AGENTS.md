# Slicer development

- Treat SlicerConfig() as the canonical, quality-first default machine profile.
  The CLI must resolve unspecified settings from it; don't recreate
  divergent or outdated defaults in argparse. Keep explicit CLI overrides,
  quality/legacy profiles, and default-vs-explicit regression tests aligned.
- Adaptive radial scheduling may only subdivide complete-terminal
  deposited slabs; preserve exact radial envelope, explicitly calibrated
  first slab, per-slab nozzle outer faces and bead-model E. Use a
  conservative projected facet-normal cusp bound across triangles active
  in physical radius. Fail when a requested cusp bound cannot meet the
  minimum printable thickness rather than silently creating slivers.
  Preserve default fixed schedules and complete-terminal foundations.
- Registered axial infill grids must use one 2π-periodic angular phase/step
  across radial layers and clipped islands. Fixed-outer preserves all tracks
  (can overfill small radii); nested dyadic insertion must preserve parent
  angular positions and never silently claim actual connected branching.
  Fixed-outer must reject gross small-radius overfill against its
  actual physical pitch and configurable density budget; do not silently
  permit extreme overlap. The nested mode can escape this by dropping
  child tracks. Keep both modes opt-in until print coupons validate.
- Narrow-feature fill must keep per-segment physical bead-width metadata and matched E; never lose widths through simplification, path reversal, or joining. Restrict thin-wall recovery to unprinted edge bands, not deliberately sparse interiors. A single 3- or 4-arm branch may be recovered only from a bounded, fully verified medial tree, with original polygon bead-sweep containment and nearly disjoint trunk/side-arm footprints. Each extra side arm must also avoid all previously accepted arms; reject multi-junction, holed, cyclic or ambiguous graph topology.
- Opt-in bottom support interfaces must seed only from exposed model
  outer radial faces with immediately adjacent support and propagate through
  consecutive support slabs. Keep top and bottom interface deposition disjoint
  and separate from sparse support; preserve wedges, periodic theta and role
  speed/flow/fan defaults. Independent bottom interface speed, flow and
  fan override must fall back to the top-interface policy when unspecified;
  preserve physical Marlin feed conversion, E exactly except intentional
  flow multipliers, and restore the layer fan after a role override.
  An optional support-bottom contact clearance must
  compare actual deposited radial faces (not rounded layer counts), carve only
  the identified contacting footprint, and never jump an existing gap.
- Cross-layer island cooling reordering must match islands by periodic
  physical position (not array indexes); avoid last-then-first consecutive
  visits only when an alternative exists. Preserve island-contiguous groups,
  safe retract/travel, all deposited geometry/E and one-direction start layers.
  Never claim exact elapsed revisit times from print-only estimates.
- Bridge classification must not promote entire irregular unsupported
  polygons based on a majority of two-ended samples. Check transverse
  physical-mm lanes for continuous spans and independent bearings at both
  ends; retain support in rejected/one-ended lanes, including across 2pi.
  A rectangular bridge may bypass stripe sampling only when both complete
  physical bearing lines are covered by the previous slab; never replace this
  with a few sampled points or a majority vote. Bridge print mode remains
  opt-in pending rotating-mandrel validation.
- Optional unsupported-wall feed limits must use the physical longest
  unsupported perimeter run from the previous radial slab; never use theta
  radians as millimeters or presume bearing on the first printed slab.
  Apply a per-path feed ceiling without modifying filament E/width or
  the bridge/support/fan policies. Leave unset until measured on hardware.
- Post-G-code cooling revisit reports must calculate a *lower bound* on
  actual elapsed time using emitted modal G0/G1 planner F, native rotary
  degree distances, standalone E-only retraction/priming, and non-moving
  G92/G92.9 coordinate rebases. Match periodic island positions and support
  classification. Do not certify a target from an optimistic or guessed
  acceleration model; `NOT_PROVEN` is not equivalent to a real-time violation.
  The report must leave every motion, E, fan and feed unchanged.
- For cooling, estimate each island's deposited path lengths in physical
  theta*radius/axial mm and respect role-specific maximum speeds and the
  machine's configured minimum speed. Cooling timing estimates must be
  conservative lower bounds (do not assume acceleration or travel time);
  do not insert stationary hot-nozzle dwells to meet impossible targets.
- Never change deposited filament E to smooth extrusion rate: when splitting an extrusion command, preserve physical endpoints, width-area metadata, and cumulative E. Hardware-dependent feed ramps remain opt-in.
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
- Bounded on-island travel routes must stay inside a single polygonal
  footprint in the physical cylindrical metric. Enforce finite graph/extra-
  distance budgets, correct periodic theta and Marlin rotary resets, and
  preserve safe retract/Z-hop fallback when routing cannot be verified.
- Adjacent fill passes should stay primed for short travel contained in the
  model/support footprint. Check the whole connection, including across the seam;
  don't suppress retraction across holes or separate islands based on distance
  alone. A role/type change is not by itself a reason to retract.
- Finish each model island's walls and infill before moving to another island.
  Opt-in island_print_order permutes infill/outer/inner within each island,
  preserving all emitted material, separate support groups, and the legacy
  wall_order behavior when unset.
  Join periodic seam components for ordering; keep support components separate.
  Enter closed perimeters using the configured seam policy (nearest by default),
  preserving winding and deposited length. Seam cuts on variable-width loops
  must split/rotate their segment_widths alongside edges, conserving E exactly.
  Physical axial targets and half-open exclusion windows must be tested as
  actual geometric edge intersections, not nearest sampled vertices. Hidden
  corner preference must classify concavity relative to the signed polygon
  orientation (reversed winding and 2pi translations invariant), keep
  exclusions strict, and fall back safely when a loop is convex.
  Keep intentional zero-distance calibration bands intact.
- Keep Marlin's native Y position bounded as well as emitted Y values. The target
  firmware enables `POWER_LOSS_RECOVERY`: use `M400`, `G92.9 Y<wrapped phase>`,
  then `G92 Y0` (or `G92 E0 Y0` at layer changes). Keep angular phase so resets
  don't accumulate fractional-step drift. Ordinary `G92 Y0` alone only shifts its
  workspace. Preserve the Marlin float32/planner regression tests and bidirectional
  infill when changing rotary coordinates.


test with ```PYTHONPATH=. .venv/bin/python3 -m pytest -n auto --timeout=30```