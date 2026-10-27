"""Add a physical return to logical Y0 between layers of bounded v2 G-code.

This is an opt-in diagnostic pass for continuous-rotary programs. It keeps all
existing extrusion commands and their coordinates intact. The cylinder moves
to Y0 at a lifted, retracted nozzle; no coordinate reset is made at the layer
boundary. The input must already use one continuous absolute Y frame.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from .gcode import MAX_MARLIN_Y_DEGREES


def add_layer_unwinds(
    gcode: str,
    *,
    lift_mm: float = 0.8,
    rotary_feed_deg_per_min: float = 3600.0,
) -> str:
    if not math.isfinite(lift_mm) or lift_mm <= 0:
        raise ValueError("lift_mm must be finite and positive")
    if not math.isfinite(rotary_feed_deg_per_min) or rotary_feed_deg_per_min <= 0:
        raise ValueError("rotary_feed_deg_per_min must be finite and positive")

    lines = gcode.splitlines()
    if not lines or lines[0] != ";start of slice":
        raise ValueError("Expected v2 G-code starting with ';start of slice'")

    output = [lines[0], ";PHYSICAL_LAYER_UNWIND:Y0 between layers, no layer Y reset"]
    z = None
    y = 0.0
    retracted = False
    layers = 0
    unwind_count = 0
    absolute = False

    for line in lines[1:]:
        if line.startswith(";LAYER:"):
            if layers:
                if z is None or not retracted:
                    raise ValueError("Layer unwind requires a known Z and retracted filament")
                if abs(y) >= MAX_MARLIN_Y_DEGREES:
                    raise ValueError("Layer unwind requires bounded absolute Y coordinates")
                output.extend([
                    f";PHYSICAL_LAYER_UNWIND:{unwind_count + 1}",
                    f"G0 Z{z + lift_mm:.4f} F9000.0 ; lift before physical rotary unwind",
                    "M400 ; complete lift before rotary unwind",
                    f"G0 Y0.000 F{rotary_feed_deg_per_min:.1f} ; physically return to Y0",
                    "M400 ; complete physical rotary unwind",
                ])
                y = 0.0
                unwind_count += 1
            layers += 1

        words = line.split(";", 1)[0].split()
        if words:
            command = words[0]
            if command == "G90":
                absolute = True
            elif command == "G91":
                raise ValueError("Layer unwind requires absolute positioning throughout")
            if command in {"G0", "G1", "G92", "G92.9"}:
                values = {}
                for word in words[1:]:
                    if word[0] in "XYZE":
                        values[word[0]] = float(word[1:])
                if any(not math.isfinite(value) for value in values.values()):
                    raise ValueError("Non-finite machine coordinate in input G-code")
                if command in {"G92", "G92.9"} and layers and "Y" in values:
                    raise ValueError("Layer unwind requires no Y coordinate resets after startup")
                if command in {"G0", "G1"}:
                    if not absolute:
                        raise ValueError("Layer unwind requires G90 before motion")
                    z = values.get("Z", z)
                    y = values.get("Y", y)
                    if y is not None and abs(y) >= MAX_MARLIN_Y_DEGREES:
                        raise ValueError("Layer unwind requires bounded absolute Y coordinates")
                    if command == "G1" and "; retract" in line:
                        retracted = True
                    elif command == "G1" and "; unretract" in line:
                        retracted = False
        output.append(line)

    if layers < 2:
        raise ValueError("Layer unwind test needs at least two layers")
    return "\n".join(output) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Add physical Y0 returns between v2 G-code layers")
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--lift-mm", type=float, default=0.8)
    parser.add_argument("--rotary-feed-deg-per-min", type=float, default=3600.0)
    args = parser.parse_args(argv)
    if args.input.resolve() == args.out.resolve():
        parser.error("input and output must be different files")
    args.out.write_text(add_layer_unwinds(
        args.input.read_text(),
        lift_mm=args.lift_mm,
        rotary_feed_deg_per_min=args.rotary_feed_deg_per_min,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
