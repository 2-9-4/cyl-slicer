"""Inspect candidate nonplanar surface ribbons without emitting machine moves.

Example:
python -m cyl_slicer_v2.nonplanar_probe part.stl \
    --point 20,15 --point 30,25 --max-segment-mm 0.4

Coordinates are theta degrees and axial machine-X millimeters.
JSON results are candidates ONLY, not collision-certified printable paths.
"""
from __future__ import annotations

import argparse
import json
import math
from typing import Sequence

import trimesh

from .mesh_io import load_mesh
from .nonplanar_surface import RadialSurfaceSampler, sample_surface_ribbon
from .nonplanar_surface import segment_surface_length


def degree_point(value: str) -> tuple[float, float]:
    try:
        fields = value.split(",")
        if len(fields) != 2:
            raise ValueError
        deg, axial = map(float, fields)
        if not math.isfinite(deg) or not math.isfinite(axial):
            raise ValueError
        return math.radians(deg), axial
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--point must be finite THETA_DEGREES,AXIAL_MM") from exc


def analyze_ribbon(mesh: trimesh.Trimesh,
                   points: Sequence[tuple[float, float]], *,
                   max_segment_mm: float = 0.5,
                   max_chord_error_mm: float = 0.015,
                   max_slope: float = 0.3,
                   max_excursion_mm: float = 1.0,
                   max_points: int = 4096) -> dict:
    """Produce a numerical report; never claim swept-tool clearance."""
    try:
        sampler = RadialSurfaceSampler(mesh)
        ribbon = sample_surface_ribbon(
            sampler, points,
            max_segment_mm=max_segment_mm,
            max_radial_chord_error_mm=max_chord_error_mm,
            max_radial_slope=max_slope,
            max_total_radial_excursion_mm=max_excursion_mm,
            max_points=max_points,
        )
    except ValueError as exc:
        return {
            "status": "rejected",
            "reason": str(exc),
            "nonplanar_extrusion_enabled": False,
        }
    radii = [point.radius for point in ribbon]
    lengths = [segment_surface_length(a, b) for a, b in zip(ribbon, ribbon[1:])]
    return {
        "status": "geometry_candidate_only",
        "nonplanar_extrusion_enabled": False,
        "warning": "No printhead swept-collision, bead ownership, or actual deposition validation.",
        "sample_count": len(ribbon),
        "min_radius_mm": min(radii),
        "max_radius_mm": max(radii),
        "radial_excursion_mm": max(radii) - min(radii),
        "path_length_3d_mm": sum(lengths),
        "max_segment_3d_mm": max(lengths) if lengths else 0.0,
        "min_outward_radial_normal": min(point.outward_normal_radial for point in ribbon),
        "samples": [
            {"theta_deg": math.degrees(point.theta),
             "axial_mm": point.axial, "radius_mm": point.radius}
            for point in ribbon
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Inspect cylindrical conformal top-surface candidates; NEVER emits G-code."
    )
    parser.add_argument("mesh", help="Input STL or supported mesh")
    parser.add_argument("--point", action="append", type=degree_point, required=True,
                        help="Repeat THETA_DEGREES,AXIAL_MM to define a continuous ribbon.")
    parser.add_argument("--max-segment-mm", type=float, default=0.5)
    parser.add_argument("--max-chord-error-mm", type=float, default=0.015)
    parser.add_argument("--max-slope", type=float, default=0.3)
    parser.add_argument("--max-excursion-mm", type=float, default=1.0)
    parser.add_argument("--max-points", type=int, default=4096)
    args = parser.parse_args(argv)
    if len(args.point) < 2:
        parser.error("at least two --point positions are required")
    result = analyze_ribbon(
        load_mesh(args.mesh), args.point,
        max_segment_mm=args.max_segment_mm,
        max_chord_error_mm=args.max_chord_error_mm,
        max_slope=args.max_slope,
        max_excursion_mm=args.max_excursion_mm,
        max_points=args.max_points,
    )
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["status"] == "geometry_candidate_only" else 2


if __name__ == "__main__":
    raise SystemExit(main())
