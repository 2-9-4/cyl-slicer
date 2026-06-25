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
    parser.add_argument("--theta-step-deg", type=float, default=5.0)
    parser.add_argument("--axial-step", type=float, default=1.0)
    parser.add_argument("--wall-count", type=int, default=2)
    parser.add_argument("--infill-spacing", type=float, default=2.0)
    parser.add_argument("--support-spacing", type=float, default=1.5)
    parser.add_argument(
        "--slice-workers",
        type=int,
        default=os.cpu_count() or 1,
        help="Parallel worker processes for geometry slicing (default: CPU count).",
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
        theta_step_deg=args.theta_step_deg,
        axial_step=args.axial_step,
        wall_count=args.wall_count,
        infill_spacing=args.infill_spacing,
        support_spacing=args.support_spacing,
        slice_workers=args.slice_workers,
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
