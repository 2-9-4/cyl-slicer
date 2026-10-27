from __future__ import annotations

import argparse
import os
from pathlib import Path
import time

from .calibration import add_calibration_arguments, sweep_from_args
from .config import SlicerConfig
from .fixtures import FIXTURES, make_fixture
from .gcode import emit_gcode
from .layers import slice_mesh
from .mesh_io import load_mesh
from .occupancy_backends import build_occupancy_backend
from .toolpaths import build_toolpaths


# Snapshot of the pre-quality CLI defaults for reproduction of old slices.
# New production defaults live in SlicerConfig, so Python and CLI agree.
CLI_PROFILE_OVERRIDES = {
    "layer_backend": "raster",
    "layer_height": 0.5,
    "nozzle_width": 0.4,
    "nozzle_diameter": None,
    "bead_model": "rectangular",
    "wall_count": 3,
    "top_layers": 3,
    "bottom_layers": 3,
    "infill_spacing": 2.0,
    "sparse_infill_pattern": "both",
    "solid_infill_pattern": "axial",
    "monotonic_top_skin": False,
    "gap_fill_enabled": False,
    "min_path_length": 0.6,
    "intersection_curve_error_mm": 0.025,
    "bed_radius": 30.0,
    "support_min_island_area": 10.0,
    "support_min_region_area": 10.0,
    "support_bearing_min_area": 25.0,
    "support_merge_gap": 10.0,
    "raft_merge_gap_mm": 0.0,
    "support_interface_layers": 2,
    "support_interface_spacing": None,
    "feed_mm_per_min": 1800.0,
    "travel_mm_per_min": 3000.0,
    "cooling_min_layer_time_s": 0.0,
    "cooling_min_island_print_time_s": 0.0,
    "cooling_min_print_speed_mm_per_min": 300.0,
    "fan_speed": 255,
    "fan_start_height": 2.0,
    "fan_ramp_height_mm": 0.0,
    "hotend_temp": 215,
    "material_profile": "unprofiled",
    "retract_speed_mm_per_min": 1800.0,
    "retract_mm": 1.0,
}
_CLI_ALIASES = {
    "nozzle_width": "line_width",
    "feed_mm_per_min": "print_feed",
    "travel_mm_per_min": "travel_feed",
    "retract_speed_mm_per_min": "retract_speed",
    "bridge_print_feed_mm_per_min": "bridge_print_feed",
    "cooling_min_print_speed_mm_per_min": "cooling_min_print_speed",
}


def _profile_value(args, key: str):
    # A default of None in argparse distinguishes explicit flags from the
    # active profile. Physical retraction speed is configured in mm/min but
    # the CLI exposes mm/s.
    option = _CLI_ALIASES.get(key, key)
    specified = getattr(args, option)
    if specified is not None:
        return specified * 60.0 if key == "retract_speed_mm_per_min" else specified
    if args.profile == "legacy-cli" and key in CLI_PROFILE_OVERRIDES:
        return CLI_PROFILE_OVERRIDES[key]
    return getattr(SlicerConfig(), key)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Cylindrical slicer v2 (quality profile by default)")
    parser.add_argument("input", nargs="?", help="STL/mesh input path")
    parser.add_argument("--fixture", choices=sorted(FIXTURES), help="Use a generated fixture mesh")
    parser.add_argument("--out", default="output_v2.gcode", help="G-code output path")
    parser.add_argument("--max-layers", type=int, default=None)
    parser.add_argument("--profile", choices=("quality", "legacy-cli", "library"), default="quality",
                        help="Quality (default, same as SlicerConfig); library is its alias; legacy-cli reproduces previous CLI defaults.")
    parser.add_argument("--layer-height", type=float, default=None)
    parser.add_argument("--radial-layer-mode", choices=("fixed", "complete-terminal"), default=None,
                        help="Fixed midpoint slabs (default), or experimental complete-terminal slabs with bounded partial thickness.")
    parser.add_argument("--filament-diameter", type=float, default=1.75, help="Filament diameter in mm.")
    parser.add_argument("--nozzle-diameter", type=float, default=None,
                        help="Physical nozzle orifice diameter in mm; independent of line width.")
    parser.add_argument("--theta-step-deg", type=float, default=5.0)
    parser.add_argument("--axial-step", type=float, default=1.0)
    parser.add_argument("--axial-min", type=float, default=15.0, help="Minimum emitted X position in mm.")
    parser.add_argument("--wall-count", type=int, default=None)
    parser.add_argument("--wall-order", choices=("outer-first", "inner-first"), default="outer-first",
                        help="Per-island printing order; inner-first prints supporting walls before visible shell.")
    parser.add_argument("--seam-policy", choices=("nearest", "aligned", "corner"), default="nearest",
                        help="Closed-wall start selection; nearest preserves current toolpath behavior.")
    parser.add_argument("--seam-angle-deg", type=float, default=180.0,
                        help="Cylinder angular phase target for aligned seams and corner tie-breaking.")
    parser.add_argument("--seam-exclude-deg", type=float, nargs=2, default=None, metavar=("START", "END"),
                        help="Exclude closed-loop starts in this angular window (wraps through zero if START>END).")
    parser.add_argument("--seam-end-ramp-mm", type=float, default=None,
                        help="Opt-in physical-mm feed reduction over the end of closed perimeter loops.")
    parser.add_argument("--seam-end-ramp-fraction", type=float, default=None,
                        help="Target final feed fraction at the seam; no filament is removed.")
    parser.add_argument("--wall-support-report", action=argparse.BooleanOptionalAction,
                        default=None, help="Annotate per-wall support fraction and longest unsupported span in physical mm.")
    parser.add_argument("--top-layers", type=int, default=None, help="Number of solid shells on the outer surface.")
    parser.add_argument("--bottom-layers", type=int, default=None, help="Number of solid shells on the inner surface.")
    parser.add_argument("--infill-spacing", type=float, default=None)
    parser.add_argument("--sparse-infill-pattern", choices=("both", "alternate"), default=None,
                        help="Legacy crosshatched grid on every layer or one alternating grid axis per layer.")
    parser.add_argument("--solid-infill-pattern", choices=("axial", "alternate"), default=None,
                        help="Legacy axial solid lines or alternate radial/axial lines by layer.")
    parser.add_argument("--monotonic-top-skin", action=argparse.BooleanOptionalAction,
                        default=None,
                        help="Opt-in one-way exposed top-surface passes; default top fill meanders.")
    parser.add_argument("--gap-fill-enabled", action=argparse.BooleanOptionalAction,
                        default=None, help="Recover thin features with variable-width centerlines (quality default on).")
    parser.add_argument("--gap-fill-min-width-mm", type=float, default=None)
    parser.add_argument("--gap-fill-sample-step-mm", type=float, default=None)
    parser.add_argument("--gap-fill-min-length-mm", type=float, default=None)
    parser.add_argument("--line-width", type=float, default=None,
                        help="Nominal deposited bead width in mm (not physical nozzle diameter).")
    parser.add_argument("--bead-model", choices=("rectangular", "rounded"), default=None,
                        help="Coupled extrusion/spacing model; rounded requires layer height <= line width.")
    parser.add_argument("--infill-wall-overlap-fraction", type=float, default=0.15,
                        help="Fraction of bead width infill may extend into the inner wall footprint.")
    parser.add_argument("--min-path-length", type=float, default=None,
                        help="Reject entire unprintably short features, not corners within a path.")
    add_calibration_arguments(parser)
    parser.add_argument("--retract-speed", type=float, default=None, help="Retraction speed in mm/s.")
    parser.add_argument("--hotend-temp", type=int, default=None)
    parser.add_argument("--material-profile", choices=("unprofiled", "pla", "petg", "abs", "tpu"),
                        default=None, help="Material identity metadata, not an assumed pressure-advance K.")
    parser.add_argument("--marlin-linear-advance", action="store_true",
                        help="Explicitly declare flashed firmware supports Marlin LIN_ADVANCE M900.")
    parser.add_argument("--linear-advance-k", type=float, default=None,
                        help="Calibrated Marlin K, 0..10; requires --marlin-linear-advance.")
    parser.add_argument("--extrusion-feed-ramp-mm", type=float, default=None,
                        help="Optional physical-mm path-start commanded feed ramp (default disabled).")
    parser.add_argument("--extrusion-feed-ramp-start-fraction", type=float, default=None,
                        help="Start/feed fraction for opt-in ramp, in (0,1], default 0.65.")
    parser.add_argument("--print-feed", type=float, default=None,
                        help="Target physical surface printing speed (mm/min); Marlin F varies per move.")
    parser.add_argument("--travel-feed", type=float, default=None,
                        help="Target physical XY travel mm/min; radial Z moves retain native mm/min.")
    parser.add_argument("--material-fan-policy", choices=("manual", "preset"), default=None,
                        help="Manual fan PWM (default), or opt-in material cooling preset with a conservative cap.")
    parser.add_argument("--fan-speed", type=int, default=None,
                        help="Part-cooling PWM 0–255; 255 in quality profile.")
    parser.add_argument("--fan-start-height", type=float, default=None,
                        help="Radial Z height above the build cylinder where cooling starts (quality: 1 mm).")
    parser.add_argument("--fan-ramp-height-mm", type=float, default=None,
                        help="Ramp fan PWM from 0 to --fan-speed over this radial height; 0 is an immediate step.")
    parser.add_argument("--cooling-min-layer-time-s", type=float, default=None,
                        help="Conservatively slow extrusions to this minimum estimated printing time per layer; 0 disables.")
    parser.add_argument("--cooling-min-island-print-time-s", type=float, default=None,
                        help="Target minimum print duration per connected island (default 2 s; 0 disables, no dwell).")
    parser.add_argument("--cooling-min-print-speed", type=float, default=None,
                        help="Physical mm/min floor on cooling slowdown (default 300); never adds dwell.")
    parser.add_argument("--retract-mm", type=float, default=None, help="Filament retraction distance in mm between paths.")
    parser.add_argument("--z-hop-height", type=float, default=0.4, help="Radial Z-hop lift in mm for risky travel moves.")
    parser.add_argument("--z-hop-min-travel", type=float, default=5.0, help="Minimum travel distance in mm that triggers Z-hop.")
    parser.add_argument("--route-travel-inside-island", action=argparse.BooleanOptionalAction,
                        default=None, help="Try bounded detours around voids within a connected island (default on).")
    parser.add_argument("--route-max-extra-mm", type=float, default=None,
                        help="Maximum physical mm of extra detour beyond direct travel (default 6).")
    parser.add_argument("--route-max-factor", type=float, default=None,
                        help="Maximum routed/direct distance factor (default 1.65).")
    parser.add_argument("--route-max-vertices", type=int, default=None,
                        help="Maximum visibility graph nodes; above cap fall back safely (default 96).")
    parser.add_argument("--route-max-unretracted-mm", type=float, default=None,
                        help="Safe on-island route length that may remain primed (default 5 mm).")
    parser.add_argument(
        "--no-z-hop-perimeter-crossing",
        action="store_false",
        dest="z_hop_on_perimeter_crossing",
        help="Disable Z-hop for travel moves that cross layer perimeters.",
    )
    parser.add_argument("--bridge-mode", choices=("off", "print"), default=None,
                        help="Opt-in short two-bearing spans on a rotating mandrel (requires hardware calibration).")
    parser.add_argument("--bridge-max-span-mm", type=float, default=None)
    parser.add_argument("--bridge-anchor-overlap-mm", type=float, default=None)
    parser.add_argument("--bridge-print-feed", type=float, default=None,
                        help="Bridge speed ceiling in physical mm/min.")
    parser.add_argument("--bridge-flow-ratio", type=float, default=None,
                        help="Bridge extrusion multiplier relative to calibrated bead.")
    parser.add_argument("--bridge-fan-pwm", type=int, default=None,
                        help="Temporary M106 PWM for bridges; layer fan is restored afterward.")
    parser.add_argument("--support-spacing", type=float, default=1.5)
    parser.add_argument("--support-print-feed", type=float, default=None,
                        help="Physical mm/min max print speed for support; default standard speed.")
    parser.add_argument("--support-interface-print-feed", type=float, default=None,
                        help="Physical mm/min max print speed for support interface; default standard speed.")
    parser.add_argument("--support-flow-ratio", type=float, default=1.0,
                        help="Support extrusion multiplier relative to calibrated bead (default 1.0).")
    parser.add_argument("--support-interface-flow-ratio", type=float, default=1.0,
                        help="Support-interface extrusion multiplier (default 1.0).")
    parser.add_argument(
        "--raft-layers",
        type=int,
        default=None,
        help="Number of bed-side gap support layers printed as a continuous raft.",
    )
    parser.add_argument(
        "--raft-grid-spacing",
        type=float,
        default=None,
        help="Raft grid spacing in mm; defaults to nozzle width.",
    )
    parser.add_argument(
        "--raft-axial-gap",
        type=float,
        default=None,
        help="Axial relief gap kept in the bottom raft footprint in mm.",
    )
    parser.add_argument(
        "--raft-merge-gap-mm", type=float, default=None,
        help="Merge raft angular gaps up to this physical mm at the innermost raft radius (default 2; 0 disables).",
    )
    parser.add_argument(
        "--support-min-island-area",
        type=float,
        default=None,
        help="Drop support fragments smaller than this metric area in mm^2.",
    )
    parser.add_argument(
        "--support-min-region-area",
        type=float,
        default=None,
        help="Inflate surviving support fragments below this metric area in mm^2.",
    )
    parser.add_argument(
        "--support-bearing-min-area",
        type=float,
        default=None,
        help="Minimum contact area in mm^2 before lower geometry counts as load-bearing.",
    )
    parser.add_argument(
        "--support-merge-gap",
        type=float,
        default=None,
        help="Merge combined support bodies across physical gaps up to this mm (default 2; 0 disables); preserve model/contact clearance.",
    )
    parser.add_argument("--support-contact-gap-mm", type=float, default=None,
                        help="Physical radial contact clearance mm, overrides --support-gap-layers.")
    parser.add_argument("--support-lateral-gap-mm", type=float, default=None,
                        help="Physical same-layer model/support clearance in R*theta/axial mm.")
    parser.add_argument(
        "--support-gap-layers",
        type=int,
        default=1,
        help="Empty radial layers between top support interface and supported model surface.",
    )
    parser.add_argument(
        "--support-interface-layers",
        type=int,
        default=None,
        help="Dense top support interface layers below the support gap.",
    )
    parser.add_argument(
        "--support-interface-spacing",
        type=float,
        default=None,
        help="Top support interface line spacing in mm; defaults to nozzle width.",
    )
    parser.add_argument(
        "--support-meander-max-angle-deg",
        type=float,
        default=15.0,
        help="Maximum angular width of local support meander groups in degrees.",
    )
    parser.add_argument(
        "--support-cross-spacing",
        type=float,
        default=None,
        help="Coarse axial tie spacing for sparse supports in mm; defaults to 4x support spacing.",
    )
    parser.add_argument(
        "--build-cylinder-radius",
        "--bed-radius",
        dest="bed_radius",
        type=float,
        default=None,
        help="Build cylinder radius (quality/legacy-cli: 30 mm).",
    )
    parser.add_argument(
        "--no-gap-support",
        action="store_false",
        dest="generate_gap_support",
        help="Disable support-only layers from the build cylinder up to the first model layer.",
    )
    parser.add_argument(
        "--slice-workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Parallel worker processes for geometry slicing (default: CPU count).",
    )
    parser.add_argument(
        "--toolpath-workers",
        type=int,
        default=None,
        help="Parallel worker processes for toolpath generation (default: --slice-workers).",
    )
    parser.add_argument(
        "--layer-backend",
        choices=("raster", "intersection"),
        default=None,
        help="Layer geometry backend: raster occupancy sampling or narrow triangle/cylinder intersection.",
    )
    parser.add_argument(
        "--intersection-curve-error-mm",
        type=float,
        default=None,
        help="Maximum physical mm deviation of sampled face/cylinder curves from chord.",
    )
    parser.add_argument("--intersection-strict-recovery", action="store_true",
                        help="Reject ambiguous radial topology transitions instead of warning and using bounded one-sided recovery.")
    parser.add_argument(
        "--intersection-refine-factor",
        type=int,
        default=0,
        help="Optional smooth topology voting density for mixed intersection faces; 0 disables it.",
    )
    parser.add_argument(
        "--occupancy-backend",
        choices=("mesh", "grid"),
        default="mesh",
        help="Point occupancy backend: exact mesh.contains or cached voxel grid.",
    )
    parser.add_argument("--grid-voxel-mm", type=float, default=1.0, help="Voxel size for --occupancy-backend grid.")
    parser.add_argument("--grid-cache", default=None, help="Optional .npz cache path for grid backend.")
    parser.add_argument("--rebuild-grid", action="store_true", help="Rebuild --grid-cache even when it exists.")
    return parser


def main(argv: list[str] | None = None) -> int:
    total_start = time.perf_counter()
    parser = build_parser()
    args = parser.parse_args(argv)
    sweep = sweep_from_args(args, parser)
    if sweep is not None and sweep.parameter == "linear-advance-k" and not args.marlin_linear_advance:
        parser.error("linear-advance-k calibration requires --marlin-linear-advance "
                     "and verified flashed firmware M900 support")

    timings: list[tuple[str, float]] = []

    def timed(label: str, fn):
        start = time.perf_counter()
        result = fn()
        timings.append((label, time.perf_counter() - start))
        return result

    if args.fixture:
        mesh = timed("build fixture mesh", lambda: make_fixture(args.fixture))
    elif args.input:
        mesh = timed("load mesh", lambda: load_mesh(args.input))
    else:
        raise SystemExit("provide an input mesh path or --fixture")

    config = SlicerConfig(
        layer_height=_profile_value(args, "layer_height"),
        radial_layer_mode=_profile_value(args, "radial_layer_mode"),
        filament_diameter=args.filament_diameter,
        nozzle_diameter=_profile_value(args, "nozzle_diameter"),
        theta_step_deg=args.theta_step_deg,
        axial_step=args.axial_step,
        axial_min=args.axial_min,
        wall_count=_profile_value(args, "wall_count"),
        wall_order=args.wall_order,
        seam_policy=args.seam_policy,
        seam_angle_deg=args.seam_angle_deg,
        seam_exclude_deg=(tuple(args.seam_exclude_deg) if args.seam_exclude_deg is not None else None),
        seam_end_ramp_mm=_profile_value(args, "seam_end_ramp_mm"),
        seam_end_ramp_fraction=_profile_value(args, "seam_end_ramp_fraction"),
        wall_support_report=_profile_value(args, "wall_support_report"),
        top_layers=_profile_value(args, "top_layers"),
        bottom_layers=_profile_value(args, "bottom_layers"),
        infill_spacing=_profile_value(args, "infill_spacing"),
        sparse_infill_pattern=_profile_value(args, "sparse_infill_pattern"),
        solid_infill_pattern=_profile_value(args, "solid_infill_pattern"),
        monotonic_top_skin=_profile_value(args, "monotonic_top_skin"),
        gap_fill_enabled=_profile_value(args, "gap_fill_enabled"),
        gap_fill_min_width_mm=_profile_value(args, "gap_fill_min_width_mm"),
        gap_fill_sample_step_mm=_profile_value(args, "gap_fill_sample_step_mm"),
        gap_fill_min_length_mm=_profile_value(args, "gap_fill_min_length_mm"),
        nozzle_width=_profile_value(args, "nozzle_width"),
        bead_model=_profile_value(args, "bead_model"),
        infill_wall_overlap_fraction=args.infill_wall_overlap_fraction,
        min_path_length=_profile_value(args, "min_path_length"),
        retract_mm=_profile_value(args, "retract_mm"),
        retract_speed_mm_per_min=_profile_value(args, "retract_speed_mm_per_min"),
        retraction_sweep=sweep,
        hotend_temp=_profile_value(args, "hotend_temp"),
        feed_mm_per_min=_profile_value(args, "feed_mm_per_min"),
        material_profile=_profile_value(args, "material_profile"),
        # For a K sweep, a known no-advance base (K0) is explicitly set.
        # Every M900 remains gated by the verified firmware capability flag.
        linear_advance_k=(0.0 if sweep is not None
                          and sweep.parameter == "linear-advance-k"
                          and args.linear_advance_k is None
                          else args.linear_advance_k),
        extrusion_feed_ramp_mm=_profile_value(args, "extrusion_feed_ramp_mm"),
        extrusion_feed_ramp_start_fraction=_profile_value(args, "extrusion_feed_ramp_start_fraction"),
        linear_advance_firmware=("marlin-m900" if args.marlin_linear_advance else None),
        travel_mm_per_min=_profile_value(args, "travel_mm_per_min"),
        fan_speed=_profile_value(args, "fan_speed"),
        material_fan_policy=_profile_value(args, "material_fan_policy"),
        fan_start_height=_profile_value(args, "fan_start_height"),
        fan_ramp_height_mm=_profile_value(args, "fan_ramp_height_mm"),
        cooling_min_layer_time_s=_profile_value(args, "cooling_min_layer_time_s"),
        cooling_min_island_print_time_s=_profile_value(args, "cooling_min_island_print_time_s"),
        cooling_min_print_speed_mm_per_min=_profile_value(args, "cooling_min_print_speed_mm_per_min"),
        z_hop_height=args.z_hop_height,
        z_hop_min_travel=args.z_hop_min_travel,
        route_travel_inside_island=_profile_value(args, "route_travel_inside_island"),
        route_max_extra_mm=_profile_value(args, "route_max_extra_mm"),
        route_max_factor=_profile_value(args, "route_max_factor"),
        route_max_vertices=_profile_value(args, "route_max_vertices"),
        route_max_unretracted_mm=_profile_value(args, "route_max_unretracted_mm"),
        z_hop_on_perimeter_crossing=args.z_hop_on_perimeter_crossing,
        bridge_mode=_profile_value(args, "bridge_mode"),
        bridge_max_span_mm=_profile_value(args, "bridge_max_span_mm"),
        bridge_anchor_overlap_mm=_profile_value(args, "bridge_anchor_overlap_mm"),
        bridge_print_feed_mm_per_min=_profile_value(args, "bridge_print_feed_mm_per_min"),
        bridge_flow_ratio=_profile_value(args, "bridge_flow_ratio"),
        bridge_fan_pwm=_profile_value(args, "bridge_fan_pwm"),
        support_spacing=args.support_spacing,
        support_print_feed_mm_per_min=args.support_print_feed,
        support_interface_print_feed_mm_per_min=args.support_interface_print_feed,
        support_flow_ratio=args.support_flow_ratio,
        support_interface_flow_ratio=args.support_interface_flow_ratio,
        support_min_island_area=_profile_value(args, "support_min_island_area"),
        support_min_region_area=_profile_value(args, "support_min_region_area"),
        support_bearing_min_area=_profile_value(args, "support_bearing_min_area"),
        support_merge_gap=_profile_value(args, "support_merge_gap"),
        support_gap_layers=_profile_value(args, "support_gap_layers"),
        support_contact_gap_mm=_profile_value(args, "support_contact_gap_mm"),
        support_lateral_gap_mm=_profile_value(args, "support_lateral_gap_mm"),
        support_interface_layers=_profile_value(args, "support_interface_layers"),
        support_interface_spacing=_profile_value(args, "support_interface_spacing"),
        support_meander_max_angle_deg=args.support_meander_max_angle_deg,
        support_cross_spacing=args.support_cross_spacing,
        raft_layers=_profile_value(args, "raft_layers"),
        raft_grid_spacing=_profile_value(args, "raft_grid_spacing"),
        raft_axial_gap=_profile_value(args, "raft_axial_gap"),
        raft_merge_gap_mm=_profile_value(args, "raft_merge_gap_mm"),
        bed_radius=_profile_value(args, "bed_radius"),
        generate_gap_support=args.generate_gap_support,
        slice_workers=args.slice_workers,
        toolpath_workers=args.toolpath_workers if args.toolpath_workers is not None else args.slice_workers,
        layer_backend=_profile_value(args, "layer_backend"),
        intersection_refine_factor=args.intersection_refine_factor,
        intersection_curve_error_mm=_profile_value(args, "intersection_curve_error_mm"),
        intersection_strict_recovery=args.intersection_strict_recovery,
    )
    try:
        config.validate()
    except ValueError as exc:
        parser.error(str(exc))
    backend = None
    if config.layer_backend == "raster":
        backend = timed(
            "occupancy backend",
            lambda: build_occupancy_backend(
                mesh,
                kind=args.occupancy_backend,
                grid_voxel_mm=args.grid_voxel_mm,
                grid_cache=args.grid_cache,
                rebuild_grid=args.rebuild_grid,
            ),
        )
    layers = timed("slice geometry", lambda: slice_mesh(mesh, config, max_layers=args.max_layers, occupancy_backend=backend))
    layer_paths = timed("generate toolpaths", lambda: build_toolpaths(layers, config))
    gcode = timed("emit gcode", lambda: emit_gcode(layer_paths, config))
    timed("write output", lambda: Path(args.out).write_text(gcode))
    timings.append(("total", time.perf_counter() - total_start))

    print("Timings:")
    for label, seconds in timings:
        print(f"  {label:<18} {seconds:8.3f}s")
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
