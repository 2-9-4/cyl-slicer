"""Inspect candidate nonplanar surfaces without generating printable G-code.

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
                   max_points: int = 4096,
                   motion_bead_area_mm2: float | None = None,
                   motion_feed_mm_per_min: float = 2700,
                   motion_filament_diameter_mm: float = 1.75,
                   motion_max_e_error_fraction: float = 0.01,
                   motion_max_segments: int = 8192,
                   motion_max_volumetric_flow_mm3_s: float | None = None,
                   motion_xyz_commands_preview: bool = False,
                   motion_bed_radius_mm: float = 30.0) -> dict:
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
    result = {
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
    if motion_bead_area_mm2 is not None:
        from .nonplanar_plan import plan_conformal_motion
        try:
            plan = plan_conformal_motion(
                ribbon,
                physical_feed_mm_per_min=motion_feed_mm_per_min,
                bead_area_mm2=motion_bead_area_mm2,
                filament_diameter_mm=motion_filament_diameter_mm,
                max_local_volume_error_fraction=motion_max_e_error_fraction,
                max_volumetric_flow_mm3_s=motion_max_volumetric_flow_mm3_s,
                max_segments=motion_max_segments,
            )
        except ValueError as exc:
            return {
                "status": "rejected",
                "reason": str(exc),
                "nonplanar_extrusion_enabled": False,
            }
        if motion_xyz_commands_preview:
            from .config import SlicerConfig
            from .models import Toolpath
            from .nonplanar_emission import preview_conformal_toolpath_xyz
            config = SlicerConfig(
                bed_radius=motion_bed_radius_mm,
                filament_diameter=motion_filament_diameter_mm,
                feed_mm_per_min=motion_feed_mm_per_min,
                max_volumetric_flow_mm3_s=(
                    motion_max_volumetric_flow_mm3_s
                    if motion_max_volumetric_flow_mm3_s is not None else 0.0
                ),
            )
            candidate = Toolpath(
                "TOP_SKIN", ribbon[0].radius,
                [(point.theta, point.axial) for point in ribbon],
                point_radii=[point.radius for point in ribbon],
                segment_areas_mm2=[motion_bead_area_mm2]*(len(ribbon)-1),
            )
            try:
                commands = preview_conformal_toolpath_xyz(
                    candidate, config,
                    max_local_e_error_fraction=motion_max_e_error_fraction,
                    max_segments=motion_max_segments,
                )
            except ValueError as exc:
                return {
                    "status": "rejected",
                    "reason": str(exc),
                    "nonplanar_extrusion_enabled": False,
                }
            result["xyz_commands_preview"] = {
                "status": commands.status,
                "nonplanar_extrusion_enabled": False,
                "warning": "UNCERTIFIED MOTION SNIPPET: no homing, machine origin, initial travel, shell ownership, deposited-layer/mandrel/nozzle collision, retraction or finish. NEVER run as printer G-code.",
                "command_count": commands.segment_count,
                "source_segment_count": commands.source_segment_count,
                "physical_length_mm": commands.exact_physical_length_mm,
                "ideal_filament_length_mm": commands.ideal_filament_length_mm,
                "emitted_filament_length_mm": commands.emitted_filament_length_mm,
                "max_local_e_error_fraction": commands.maximum_local_e_error_fraction,
                "commands": list(commands.commands),
            }
        result["motion_plan"] = {
            "status": "analysis_only_no_gcode",
            "warning": "Machine-linear conformal E/feed analysis only: NO approved nozzle clearance, support or shell replacement.",
            "native_marlin_y_units": "degrees",
            "segment_count": len(plan.segments),
            "original_segment_count": plan.input_segment_count,
            "split_count": plan.split_count,
            "physical_length_mm": plan.physical_length_mm,
            "filament_length_mm": plan.filament_length_mm,
            "max_local_e_error_fraction": plan.largest_volume_error_fraction,
            "segments": [
                {
                    "start": {"theta_deg": math.degrees(seg.start[0]),
                              "axial_mm": seg.start[1],
                              "radius_mm": seg.start[2]},
                    "end": {"theta_deg": math.degrees(seg.end[0]),
                            "axial_mm": seg.end[1],
                            "radius_mm": seg.end[2]},
                    "physical_length_mm": seg.metrics.length_mm,
                    "marlin_native_length_mm": seg.metrics.marlin_native_length,
                    "marlin_feed_mm_per_min": seg.metrics.marlin_feed_mm_per_min,
                    "cumulative_filament_mm": seg.cumulative_e_end_mm,
                    "max_local_e_error_fraction": seg.metrics.max_local_volume_error_fraction,
                }
                for seg in plan.segments
            ],
        }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Inspect cylindrical conformal top-surface candidates; NEVER emits G-code."
    )
    parser.add_argument("mesh", help="Input STL or supported mesh")
    parser.add_argument("--point", action="append", type=degree_point,
                        help="Repeat THETA_DEGREES,AXIAL_MM to define a continuous ribbon.")
    parser.add_argument("--scan", action="store_true",
                        help="Discover connected candidate surface patches automatically (no G-code).")
    parser.add_argument("--theta-range-deg", nargs=2, type=float, metavar=("START", "END"),
                        default=(0,360), help="Unwrapped angular scan interval, up to 360 degrees.")
    parser.add_argument("--axial-range-mm", nargs=2, type=float, metavar=("START", "END"),
                        help="Axial scan interval in physical mm; defaults to source mesh bounds.")
    parser.add_argument("--theta-cells", type=int, default=72)
    parser.add_argument("--axial-cells", type=int, default=32)
    parser.add_argument("--max-cell-excursion-mm", type=float, default=0.5,
                        help="Maximum variation of surface radius within one scan cell.")
    parser.add_argument("--sharp-dihedral-deg", type=float, default=55.,
                        help="Reject scan cells touched by mesh edges with sharper dihedral.")
    parser.add_argument("--min-patch-cells", type=int, default=4)
    parser.add_argument("--max-scan-cells", type=int, default=12000)
    parser.add_argument("--raster-preview", action="store_true",
                        help="Suggest surface-following axial finish lines for accepted scan cells; JSON only.")
    parser.add_argument("--raster-pitch-mm", type=float, default=.42)
    parser.add_argument("--raster-edge-inset-mm", type=float, default=.21)
    parser.add_argument("--raster-min-path-mm", type=float, default=.5)
    parser.add_argument("--raster-max-paths", type=int, default=2048)
    parser.add_argument("--raster-max-points", type=int, default=100000)
    parser.add_argument("--raster-include-points", action="store_true")
    parser.add_argument("--max-segment-mm", type=float, default=0.5)
    parser.add_argument("--max-chord-error-mm", type=float, default=0.015)
    parser.add_argument("--max-slope", type=float, default=0.3)
    parser.add_argument("--max-excursion-mm", type=float, default=1.0)
    parser.add_argument("--max-points", type=int, default=4096)
    parser.add_argument("--motion-plan", action="store_true",
                        help="Also analyze XYZ segment F/E (JSON only; never a printable job).")
    parser.add_argument("--xyz-motion-preview", action="store_true",
                        help="Include unexecutable G1 XYZ/E/F motion snippets in JSON for inspection (NOT printable G-code).")
    parser.add_argument("--motion-bed-radius-mm", type=float, default=30.0,
                        help="Build cylinder radius for non-printable XYZ preview Z coordinates.")
    parser.add_argument("--bead-area-mm2", type=float,
                        help="Required with --motion-plan; measured bead cross-section in mm².")
    parser.add_argument("--physical-feed-mm-per-min", type=float, default=2700)
    parser.add_argument("--filament-diameter-mm", type=float, default=1.75)
    parser.add_argument("--motion-max-e-error-fraction", type=float, default=0.01)
    parser.add_argument("--motion-max-segments", type=int, default=8192)
    parser.add_argument("--motion-volumetric-limit-mm3-s", type=float)
    args = parser.parse_args(argv)
    if args.scan and args.point:
        parser.error("--scan and --point cannot be used together")
    if args.raster_preview and not args.scan:
        parser.error("--raster-preview requires --scan")
    if args.raster_include_points and not args.raster_preview:
        parser.error("--raster-include-points requires --raster-preview")
    if args.scan and args.motion_plan:
        parser.error("--motion-plan requires a specified --point ribbon, not --scan")
    if args.motion_plan and args.bead_area_mm2 is None:
        parser.error("--motion-plan requires explicit --bead-area-mm2")
    if args.xyz_motion_preview and not args.motion_plan:
        parser.error("--xyz-motion-preview requires --motion-plan and --bead-area-mm2")
    if args.bead_area_mm2 is not None and not args.motion_plan:
        parser.error("--bead-area-mm2 requires --motion-plan")
    if not args.scan and (not args.point or len(args.point) < 2):
        parser.error("provide at least two --point positions or use --scan")
    mesh = load_mesh(args.mesh)
    if args.scan:
        from .nonplanar_patches import scan_surface_patches
        angular = [math.radians(v) for v in args.theta_range_deg]
        axial = (args.axial_range_mm if args.axial_range_mm is not None
                 else (float(mesh.vertices[:,2].min()), float(mesh.vertices[:,2].max())))
        try:
            result = scan_surface_patches(
                mesh, theta_min=angular[0], theta_max=angular[1],
                axial_min=axial[0], axial_max=axial[1],
                theta_cells=args.theta_cells, axial_cells=args.axial_cells,
                max_slope=args.max_slope,
                max_chord_error_mm=args.max_chord_error_mm,
                max_cell_excursion_mm=args.max_cell_excursion_mm,
                sharp_dihedral_deg=args.sharp_dihedral_deg,
                max_cells=args.max_scan_cells,
                min_patch_cells=args.min_patch_cells,
            )
            if args.raster_preview:
                from .nonplanar_raster import propose_axial_conformal_rasters
                result["raster_preview"] = propose_axial_conformal_rasters(
                    mesh, result,
                    nominal_surface_pitch_mm=args.raster_pitch_mm,
                    allowed_surface_slope=args.max_slope,
                    edge_inset_mm=args.raster_edge_inset_mm,
                    min_path_length_mm=args.raster_min_path_mm,
                    max_segment_mm=args.max_segment_mm,
                    max_radial_chord_error_mm=args.max_chord_error_mm,
                    max_radial_excursion_mm=args.max_excursion_mm,
                    max_paths=args.raster_max_paths,
                    max_total_points=args.raster_max_points,
                    include_points=args.raster_include_points,
                )
        except ValueError as exc:
            parser.error(str(exc))
    else:
        result = analyze_ribbon(
            mesh, args.point, max_segment_mm=args.max_segment_mm,
            max_chord_error_mm=args.max_chord_error_mm,
            max_slope=args.max_slope,
            max_excursion_mm=args.max_excursion_mm,
            max_points=args.max_points,
            motion_bead_area_mm2=(args.bead_area_mm2 if args.motion_plan else None),
            motion_feed_mm_per_min=args.physical_feed_mm_per_min,
            motion_filament_diameter_mm=args.filament_diameter_mm,
            motion_max_e_error_fraction=args.motion_max_e_error_fraction,
            motion_max_segments=args.motion_max_segments,
            motion_max_volumetric_flow_mm3_s=args.motion_volumetric_limit_mm3_s,
            motion_xyz_commands_preview=args.xyz_motion_preview,
            motion_bed_radius_mm=args.motion_bed_radius_mm,
        )
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["status"] in {"geometry_candidate_only", "geometry_candidates_only"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
