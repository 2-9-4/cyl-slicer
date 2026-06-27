from __future__ import annotations

import argparse
import os
from pathlib import Path
import time

from .config import SlicerConfig
from .fixtures import FIXTURES, make_fixture
from .gcode import emit_gcode
from .layers import slice_mesh
from .mesh_io import load_mesh
from .occupancy_backends import build_occupancy_backend
from .toolpaths import build_toolpaths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Raster-based cylindrical slicer v2")
    parser.add_argument("input", nargs="?", help="STL/mesh input path")
    parser.add_argument("--fixture", choices=sorted(FIXTURES), help="Use a generated fixture mesh")
    parser.add_argument("--out", default="output_v2.gcode", help="G-code output path")
    parser.add_argument("--max-layers", type=int, default=None)
    parser.add_argument("--layer-height", type=float, default=0.5)
    parser.add_argument("--filament-diameter", type=float, default=1.75, help="Filament diameter in mm.")
    parser.add_argument("--theta-step-deg", type=float, default=5.0)
    parser.add_argument("--axial-step", type=float, default=1.0)
    parser.add_argument("--axial-min", type=float, default=15.0, help="Minimum emitted X position in mm.")
    parser.add_argument("--wall-count", type=int, default=3)
    parser.add_argument("--top-layers", type=int, default=3, help="Number of solid shells on the outer surface.")
    parser.add_argument("--bottom-layers", type=int, default=3, help="Number of solid shells on the inner surface.")
    parser.add_argument("--infill-spacing", type=float, default=2.0)
    parser.add_argument("--min-path-length", type=float, default=0.6, help="Drop printable path segments shorter than this mm length.")
    parser.add_argument("--retract-mm", type=float, default=1.0, help="Filament retraction distance in mm between paths.")
    parser.add_argument("--support-spacing", type=float, default=1.5)
    parser.add_argument(
        "--raft-layers",
        type=int,
        default=4,
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
        default=1.0,
        help="Axial relief gap kept in the bottom raft footprint in mm.",
    )
    parser.add_argument(
        "--support-min-island-area",
        type=float,
        default=3,
        help="Drop support fragments smaller than this metric area in mm^2.",
    )
    parser.add_argument(
        "--support-min-region-area",
        type=float,
        default=3.0,
        help="Inflate surviving support fragments below this metric area in mm^2.",
    )
    parser.add_argument(
        "--support-bearing-min-area",
        type=float,
        default=3.0,
        help="Minimum contact area in mm^2 before lower geometry counts as load-bearing.",
    )
    parser.add_argument(
        "--support-merge-gap",
        type=float,
        default=3.0,
        help="Merge support branches closer than this metric gap in mm.",
    )
    parser.add_argument(
        "--support-gap-layers",
        type=int,
        default=1,
        help="Empty radial layers between top support interface and supported model surface.",
    )
    parser.add_argument(
        "--support-interface-layers",
        type=int,
        default=2,
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
        "--build-cylinder-radius",
        "--bed-radius",
        dest="bed_radius",
        type=float,
        default=30.0,
        help="Radius of the printer build cylinder; support foundations start from this radius.",
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
        default="raster",
        help="Layer geometry backend: raster occupancy sampling or narrow triangle/cylinder intersection.",
    )
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
    args = build_parser().parse_args(argv)

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
        layer_height=args.layer_height,
        filament_diameter=args.filament_diameter,
        theta_step_deg=args.theta_step_deg,
        axial_step=args.axial_step,
        axial_min=args.axial_min,
        wall_count=args.wall_count,
        top_layers=args.top_layers,
        bottom_layers=args.bottom_layers,
        infill_spacing=args.infill_spacing,
        min_path_length=args.min_path_length,
        retract_mm=args.retract_mm,
        support_spacing=args.support_spacing,
        support_min_island_area=args.support_min_island_area,
        support_min_region_area=args.support_min_region_area,
        support_bearing_min_area=args.support_bearing_min_area,
        support_merge_gap=args.support_merge_gap,
        support_gap_layers=args.support_gap_layers,
        support_interface_layers=args.support_interface_layers,
        support_interface_spacing=args.support_interface_spacing,
        support_meander_max_angle_deg=args.support_meander_max_angle_deg,
        raft_layers=args.raft_layers,
        raft_grid_spacing=args.raft_grid_spacing,
        raft_axial_gap=args.raft_axial_gap,
        bed_radius=args.bed_radius,
        generate_gap_support=args.generate_gap_support,
        slice_workers=args.slice_workers,
        toolpath_workers=args.toolpath_workers if args.toolpath_workers is not None else args.slice_workers,
        layer_backend=args.layer_backend,
        intersection_refine_factor=args.intersection_refine_factor,
    )
    backend = None
    if args.layer_backend == "raster":
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
