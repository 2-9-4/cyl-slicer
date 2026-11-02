"""Conservative revisit-time diagnostics from emitted Marlin motion commands.

This is a postprocessor for slicer-owned G-code, not a firmware simulator.
Native planner feed F is interpreted in the emitted X-mm/Y-degree/Z-mm units:
ideal constant-feed motion durations are *lower bounds* on actual runtime
where acceleration, junction slowdown and command processing add time.

Extrusion-only retract/prime moves are measured using native E-mm/F. G92,
G92.9 and M400 never become fabricated moving time. The parser expects
absolute XYZ/E as emitted by the slicer and rejects unsupported relative
coordinates instead of manufacturing an unsafe time certificate.
"""
from __future__ import annotations

import math
import re

TOKEN = re.compile(r"(?<![A-Za-z])([XYZEF])(-?(?:\\d+(?:\\.\\d*)?|\\.\\d+))")


def _coordinates(line: str) -> dict[str, float]:
    return {key: float(value) for key, value in TOKEN.findall(line.split(";", 1)[0])}


def annotate_revisit_timeline(code: str, *, match_mm: float, target_s: float = 0.0) -> str:
    """Annotate min elapsed between corresponding island deposits.

    Only positive E moves with nonzero geometric travel count as deposited
    material. The time is based on the actual emitted G-code, after speed
    ceilings, near-seam splits, travel routing and rotary coordinate resets.
    """
    if match_mm <= 0 or not math.isfinite(match_mm) or target_s < 0 or not math.isfinite(target_s):
        raise ValueError("Revisit matching radius must be positive and time target nonnegative")

    lines = code.splitlines()
    position: dict[str, float | None] = {"X": None, "Y": None, "Z": None, "E": 0.0}
    feed: float | None = None
    elapsed = 0.0
    layer: int | None = None
    radius: float | None = None
    group = None
    records: list[dict] = []
    annotations: dict[int, list[str]] = {}

    for index, line in enumerate(lines):
        if line.startswith(";LAYER:"):
            layer = int(line.split(":", 1)[1])
            group = None
            continue
        if line.startswith(";SAMPLE_RADIUS:"):
            radius = float(line.split(":", 1)[1])
            continue
        if line.startswith(";COOLING_REVISIT_GROUP:"):
            values = line.split()
            group = {
                "marker": index,
                "group": int(values[0].split(":", 1)[1]),
                "theta": float(values[1].split(":", 1)[1]) % (2 * math.pi),
                "axial": float(values[2].split(":", 1)[1]),
                "support": bool(int(values[3].split(":", 1)[1])),
                "layer": layer,
                "radius": radius,
                "first": None,
                "last": None,
            }
            records.append(group)
            continue
        command = line.split(";", 1)[0].strip()
        if not command:
            continue
        op = command.split(None, 1)[0]
        if op in {"G91", "M83"}:
            raise ValueError("Cannot certify G-code revisit timeline with relative coordinates")
        if op in {"G92", "G92.9"}:
            # G92.9 rebases native Marlin Y steps, then G92 rebases the
            # workspace. Neither reset moves the physical axis.
            for axis, value in _coordinates(command).items():
                if axis in position:
                    position[axis] = value
            continue
        if op not in {"G0", "G1", "G00", "G01"}:
            continue
        args = _coordinates(command)
        if "F" in args:
            feed = args["F"]
        before = dict(position)
        for axis in "XYZE":
            if axis in args:
                position[axis] = args[axis]
        moving = [
            position[axis] - before[axis]
            for axis in "XYZ"
            if position[axis] is not None and before[axis] is not None
        ]
        spatial = math.sqrt(sum(delta * delta for delta in moving))
        e_delta = (position["E"] - before["E"]
                   if position["E"] is not None and before["E"] is not None else 0)
        distance = spatial if spatial > 1e-10 else abs(e_delta)
        seconds = distance * 60.0 / feed if distance > 1e-10 and feed and feed > 0 else 0.0
        # These coordinates are native XYZ, including *degrees* for rotary Y.
        # Never scale Y by radius here: F was already converted by emitter.
        if group is not None and e_delta > 1e-9 and spatial > 1e-10:
            if group["first"] is None:
                group["first"] = elapsed
            group["last"] = elapsed + seconds
        elapsed += seconds

    previous: list[dict] = []
    for record in records:
        if record["first"] is None:
            continue
        matches = []
        for old in previous:
            if (old["layer"] is None or record["layer"] is None
                    or old["layer"] >= record["layer"]
                    or old["support"] != record["support"]
                    or record["radius"] is None):
                continue
            dt = math.remainder(record["theta"] - old["theta"], 2 * math.pi)
            distance = math.hypot(dt * record["radius"], record["axial"] - old["axial"])
            if distance <= match_mm + 1e-9:
                matches.append((old["layer"], -distance, old))
        if matches:
            matched = max(matches, key=lambda x: (x[0], x[1]))[2]
            minimum = max(0.0, record["first"] - matched["last"])
            result = [
                f";COOLING_REVISIT_LOWER_BOUND_S:{minimum:.4f} GROUP:{record['group']}"
                f" PRIOR_LAYER:{matched['layer']} PRIOR_GROUP:{matched['group']}"
            ]
            if target_s > 0:
                if minimum + 1e-9 >= target_s:
                    result.append(
                        f";COOLING_REVISIT_TARGET_CERTIFIED_LOWER_BOUND:"
                        f"TARGET_S:{target_s:g} MIN_S:{minimum:.4f}"
                    )
                else:
                    result.append(
                        f";COOLING_REVISIT_TARGET_NOT_PROVEN:"
                        f"TARGET_S:{target_s:g} MIN_S:{minimum:.4f}"
                    )
            annotations[record["marker"]] = result
        previous.append(record)
    out = []
    for index, line in enumerate(lines):
        out.append(line)
        out.extend(annotations.get(index, ()))
    return "\n".join(out) + ("\n" if code.endswith("\n") else "")
