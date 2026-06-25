#!/usr/bin/env python3
"""Measure whether cylindrical G-code for the rotated cube is actually planar.

The viewer maps each G-code point as:
    STL x = layer_radius * cos(Y degrees)
    STL y = layer_radius * sin(Y degrees)
    STL z = X

For a real cube, external perimeter samples should land close to one of the six
STL face planes after that conversion. This script samples the emitted rotary
centerlines and reports plane error, with periodic seam walls separated out.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
from statistics import median

import numpy as np
import trimesh


@dataclass(frozen=True)
class Segment:
    line_no: int
    layer: int
    role: str
    radius: float
    x0: float
    y0: float
    x1: float
    y1: float
    raw: str

    @property
    def length_mm(self) -> float:
        return math.hypot(self.x1 - self.x0, math.radians(self.y1 - self.y0) * self.radius)


@dataclass(frozen=True)
class Sample:
    distance: float
    segment: Segment
    frac: float
    point: np.ndarray


def parse_axis_words(line: str) -> dict[str, float]:
    words: dict[str, float] = {}
    for token in line.split()[1:]:
        if len(token) < 2:
            continue
        axis = token[0].upper()
        if axis not in {"X", "Y", "Z", "E", "F"}:
            continue
        try:
            words[axis] = float(token[1:])
        except ValueError:
            pass
    return words


def parse_extrusion_segments(gcode_path: Path, roles: set[str]) -> list[Segment]:
    segments: list[Segment] = []
    layer = -9999
    role = "UNKNOWN"
    radius: float | None = None
    cur_x: float | None = None
    cur_y: float | None = None
    cur_e: float | None = None

    for line_no, raw_line in enumerate(gcode_path.read_text().splitlines(), start=1):
        line = raw_line.split(";", 1)[0].strip()
        if raw_line.startswith(";LAYER:"):
            try:
                layer = int(raw_line.split(":", 1)[1])
            except ValueError:
                layer = -9999
            continue
        if raw_line.startswith(";Z:"):
            radius = float(raw_line.split(":", 1)[1])
            continue
        if raw_line.startswith(";TYPE:"):
            role = raw_line.split(":", 1)[1].strip()
            continue
        if not line:
            continue

        cmd = line.split(maxsplit=1)[0]
        if cmd == "G92":
            words = parse_axis_words(line)
            if "X" in words:
                cur_x = words["X"]
            if "Y" in words:
                cur_y = words["Y"]
            if "E" in words:
                cur_e = words["E"]
            continue
        if cmd not in {"G0", "G1"}:
            continue

        words = parse_axis_words(line)
        next_x = words.get("X", cur_x)
        next_y = words.get("Y", cur_y)
        next_e = words.get("E", cur_e)

        is_extrusion = (
            cmd == "G1"
            and "E" in words
            and cur_e is not None
            and next_e is not None
            and next_e > cur_e
        )
        if (
            is_extrusion
            and role in roles
            and radius is not None
            and cur_x is not None
            and cur_y is not None
            and next_x is not None
            and next_y is not None
        ):
            segments.append(
                Segment(
                    line_no=line_no,
                    layer=layer,
                    role=role,
                    radius=radius,
                    x0=cur_x,
                    y0=cur_y,
                    x1=next_x,
                    y1=next_y,
                    raw=raw_line,
                )
            )

        cur_x = next_x
        cur_y = next_y
        cur_e = next_e

    return segments


def stl_face_planes(stl_path: Path) -> list[tuple[np.ndarray, float]]:
    mesh = trimesh.load(stl_path, force="mesh")
    planes: list[tuple[np.ndarray, float]] = []
    for face, raw_normal in zip(mesh.faces, mesh.face_normals):
        normal = np.array(raw_normal, dtype=float, copy=True)
        norm = float(np.linalg.norm(normal))
        if norm <= 1e-12:
            continue
        normal /= norm
        offset = -float(np.dot(normal, mesh.vertices[face[0]]))
        if any(np.allclose(normal, n, atol=1e-5) and abs(offset - d) <= 1e-4 for n, d in planes):
            continue
        planes.append((normal, offset))
    return planes


def gcode_to_stl_point(radius: float, x_axis: float, y_degrees: float) -> np.ndarray:
    theta = math.radians(y_degrees % 360.0)
    return np.array(
        [
            radius * math.cos(theta),
            radius * math.sin(theta),
            x_axis,
        ],
        dtype=float,
    )


def distance_to_nearest_plane(point: np.ndarray, planes: list[tuple[np.ndarray, float]]) -> float:
    return min(abs(float(np.dot(normal, point) + offset)) for normal, offset in planes)


def is_periodic_seam_segment(segment: Segment, seam_window_deg: float) -> bool:
    for angle in (segment.y0 % 360.0, segment.y1 % 360.0):
        if min(angle, 360.0 - angle) > seam_window_deg:
            return False
    return True


def sample_segments(
    segments: list[Segment],
    planes: list[tuple[np.ndarray, float]],
    samples_per_segment: int,
) -> list[Sample]:
    samples: list[Sample] = []
    if samples_per_segment < 2:
        raise ValueError("samples_per_segment must be at least 2")
    for segment in segments:
        for index in range(samples_per_segment):
            frac = index / (samples_per_segment - 1)
            x_axis = segment.x0 + (segment.x1 - segment.x0) * frac
            y_degrees = segment.y0 + (segment.y1 - segment.y0) * frac
            point = gcode_to_stl_point(segment.radius, x_axis, y_degrees)
            samples.append(
                Sample(
                    distance=distance_to_nearest_plane(point, planes),
                    segment=segment,
                    frac=frac,
                    point=point,
                )
            )
    return samples


def percentile(values: list[float], pct: float) -> float:
    if not values:
        return float("nan")
    return float(np.percentile(np.asarray(values, dtype=float), pct))


def print_stats(label: str, samples: list[Sample], flat_tolerance: float) -> None:
    distances = [sample.distance for sample in samples]
    if not distances:
        print(f"{label}: no samples")
        return
    over = sum(value > flat_tolerance for value in distances)
    print(
        f"{label}: {len(samples)} samples, "
        f"median={median(distances):.4f} mm, "
        f"p95={percentile(distances, 95):.4f} mm, "
        f"p99={percentile(distances, 99):.4f} mm, "
        f"max={max(distances):.4f} mm, "
        f">{flat_tolerance:g} mm={over}"
    )


def print_worst(label: str, samples: list[Sample], limit: int) -> None:
    if not samples or limit <= 0:
        return
    print(f"\nWorst {label} samples:")
    for sample in sorted(samples, key=lambda item: item.distance, reverse=True)[:limit]:
        segment = sample.segment
        angle = (segment.y0 + (segment.y1 - segment.y0) * sample.frac) % 360.0
        x_axis = segment.x0 + (segment.x1 - segment.x0) * sample.frac
        px, py, pz = sample.point
        print(
            f"  {sample.distance:8.4f} mm  "
            f"line={segment.line_no:<6d} layer={segment.layer:<5d} "
            f"role={segment.role:<20s} frac={sample.frac:.2f} "
            f"r={segment.radius:.4f} theta={angle:.3f} X={x_axis:.3f} "
            f"stl=({px:.3f}, {py:.3f}, {pz:.3f})"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check whether rotated-cube cylindrical G-code lies on the STL cube planes."
    )
    parser.add_argument("gcode", nargs="?", default="output_v2.gcode", help="G-code file to inspect.")
    parser.add_argument("stl", nargs="?", default="single_rotated_cube.stl", help="Planar STL to compare against.")
    parser.add_argument(
        "--role",
        action="append",
        default=None,
        help="G-code ;TYPE role to inspect. May be repeated. Default: EXTERNAL_PERIMETER.",
    )
    parser.add_argument(
        "--samples-per-segment",
        type=int,
        default=5,
        help="Samples along each emitted segment, including endpoints.",
    )
    parser.add_argument(
        "--seam-window-deg",
        type=float,
        default=0.01,
        help="Canonical theta window used to separate periodic seam segments.",
    )
    parser.add_argument(
        "--flat-tolerance",
        type=float,
        default=0.05,
        help="Plane-distance tolerance used in the summary counts.",
    )
    parser.add_argument("--worst", type=int, default=8, help="Number of worst samples to print per category.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    gcode_path = Path(args.gcode)
    stl_path = Path(args.stl)
    roles = set(args.role or ["EXTERNAL_PERIMETER"])

    planes = stl_face_planes(stl_path)
    if not planes:
        raise SystemExit(f"No face planes found in {stl_path}")

    segments = parse_extrusion_segments(gcode_path, roles)
    if not segments:
        raise SystemExit(f"No extrusion segments found in {gcode_path} for roles {sorted(roles)}")

    samples = sample_segments(segments, planes, args.samples_per_segment)
    seam_segments = [segment for segment in segments if is_periodic_seam_segment(segment, args.seam_window_deg)]
    seam_ids = {id(segment) for segment in seam_segments}
    seam_samples = [sample for sample in samples if id(sample.segment) in seam_ids]
    nonseam_samples = [sample for sample in samples if id(sample.segment) not in seam_ids]

    total_length = sum(segment.length_mm for segment in segments)
    seam_length = sum(segment.length_mm for segment in seam_segments)

    print(f"Compared {gcode_path} against {stl_path}")
    print(f"STL face planes: {len(planes)}")
    print(f"Roles: {', '.join(sorted(roles))}")
    print(f"Extrusion segments: {len(segments)}")
    print(
        f"Periodic seam segments: {len(seam_segments)} "
        f"({seam_length:.2f} mm of {total_length:.2f} mm, {100.0 * seam_length / total_length:.1f}%)"
    )
    print()
    print_stats("All inspected samples", samples, args.flat_tolerance)
    print_stats("Non-seam samples", nonseam_samples, args.flat_tolerance)
    print_stats("Periodic seam samples", seam_samples, args.flat_tolerance)
    print_worst("non-seam", nonseam_samples, args.worst)
    print_worst("periodic seam", seam_samples, args.worst)

    nonseam_max = max((sample.distance for sample in nonseam_samples), default=0.0)
    seam_max = max((sample.distance for sample in seam_samples), default=0.0)
    print("\nVerdict:")
    if nonseam_max > args.flat_tolerance:
        print(
            "  The visible curvature is present in the emitted G-code centerlines: "
            f"non-seam samples are up to {nonseam_max:.4f} mm from every STL cube face plane."
        )
    else:
        print("  Non-seam external samples are planar within the requested tolerance.")
    if seam_segments and seam_max > args.flat_tolerance:
        print(
            "  The G-code also extrudes periodic seam wall segments at theta=0/360; "
            f"those reach {seam_max:.4f} mm from the real cube face planes and are not a viewer artifact."
        )
    elif seam_segments:
        print("  Periodic seam samples are also planar within the requested tolerance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
