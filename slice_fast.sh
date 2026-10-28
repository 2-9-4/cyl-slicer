#!/usr/bin/env bash
# Fast PLA / 0.4 mm nozzle starting profile for the cylindrical v2 slicer.
# Usage: ./slice_fast.sh model.stl --out model.fast.gcode [slicer overrides...]
# Every CLI argument below can be overridden by passing it again after the script.
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

profile=(
    --profile quality
    --out output_fast.gcode
    # Coarser layers, ~0.9 mm walls and 0.84 mm top/bottom shells.
    --layer-height 0.28
    --line-width 0.48
    --nozzle-diameter 0.4
    --bead-model rounded
    --gap-fill-min-width-mm 0.28
    --wall-count 2
    --top-layers 3
    --bottom-layers 3
    # Alternate uses half this pitch per layer: 2.4 mm, ~17.5% infill.
    --infill-spacing 4.8
    --sparse-infill-pattern alternate
    --solid-infill-pattern alternate
    --no-monotonic-top-skin
    # Keep the foundation and dense contact pads; lighten only support bodies.
    --support-spacing 3.0
    --support-interface-layers 2
    --support-interface-spacing 0.42
    --raft-layers 2
    # Physical surface targets: 75 mm/s print, 200 mm/s travel.
    # The 6 mm³/s cap limits nominal beads to ~51 mm/s; thin beads may go faster.
    --print-feed 4500
    --travel-feed 12000
    --max-volumetric-flow 6
    --support-print-feed 3000
    --support-interface-print-feed 2100
    # Opt-in M204 P/T, in Marlin's native mixed-axis units/s², not surface mm/s².
    --print-accel 1200
    --travel-accel 2000
    --material-profile pla
    --hotend-temp 205
    # Keep the existing retraction calibration and short-island cooling safeguards.
    --retract-mm 1.4
    --retract-speed 20
    --fan-speed 255
    --cooling-min-layer-time-s 8
    --cooling-min-island-print-time-s 2
    --cooling-min-print-speed 300
)

# Keep input/output paths relative to the caller, even when invoked elsewhere.
export PYTHONPATH="${repo_root}${PYTHONPATH:+:${PYTHONPATH}}"
exec "${repo_root}/.venv/bin/python" -m cyl_slicer_v2.cli "${profile[@]}" "$@"
