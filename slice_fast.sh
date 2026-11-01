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
    --layer-height 0.12
    --first-layer-height 0.28
    --flow-percent 98
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
    --print-feed 4800
    --travel-feed 9600
    # Optional adhesion aid: pass --one-direction-first-layers 10 to avoid Y reversals.
    --support-print-feed 4800
    --support-interface-print-feed 4800
    # Physical mm/s², converted to M204 P/T for each layer radius and move direction.
    --print-accel 800
    --travel-accel 2000
    --material-profile pla
    --hotend-temp 200
    # Keep the existing retraction calibration and short-island cooling safeguards.
    --retract-mm 1.4
    --retract-speed 20
    --fan-speed 255
    --cooling-min-layer-time-s 8
    --cooling-min-island-print-time-s 2
    --cooling-min-print-speed 300
    --max-volumetric-flow 9 
    --axial-min 30 
    --support-gap-layers 2  
    --slow-start-layers 2
    --slow-start-print-feed 1200   
    --slow-start-travel-feed 4800 
    --one-direction-first-layers 0
    --support-wedge-count 12
)

# Keep input/output paths relative to the caller, even when invoked elsewhere.
export PYTHONPATH="${repo_root}${PYTHONPATH:+:${PYTHONPATH}}"
python_cmd="${repo_root}/.venv/bin/python"
# Developer/CI checkouts may have Python dependencies installed but no
# repository-local venv. Prefer the venv when available; otherwise run the
# environment-selected Python. Do not change the caller's working directory.
if [[ ! -x "${python_cmd}" ]]; then
    python_cmd="${PYTHON:-python3}"
fi
exec "${python_cmd}" -m cyl_slicer_v2.cli "${profile[@]}" "$@"
