from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np
import trimesh


class OccupancyBackend(Protocol):
    def contains(self, points: np.ndarray) -> np.ndarray:
        ...


@dataclass
class MeshContainsBackend:
    mesh: trimesh.Trimesh

    def contains(self, points: np.ndarray) -> np.ndarray:
        return self.mesh.contains(points)


@dataclass
class GridOccupancyBackend:
    origin: np.ndarray
    voxel_mm: float
    occupied: np.ndarray
    mesh_bounds: np.ndarray | None = None

    def contains(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        idx = np.rint((points - self.origin) / self.voxel_mm).astype(np.int64)
        valid = (
            (idx[:, 0] >= 0)
            & (idx[:, 0] < self.occupied.shape[0])
            & (idx[:, 1] >= 0)
            & (idx[:, 1] < self.occupied.shape[1])
            & (idx[:, 2] >= 0)
            & (idx[:, 2] < self.occupied.shape[2])
        )
        result = np.zeros(len(points), dtype=bool)
        good = idx[valid]
        result[valid] = self.occupied[good[:, 0], good[:, 1], good[:, 2]]
        return result

    def save(self, path: str | Path) -> None:
        np.savez_compressed(
            path,
            origin=self.origin,
            voxel_mm=np.array([self.voxel_mm], dtype=float),
            occupied=self.occupied,
            mesh_bounds=np.asarray(self.mesh_bounds, dtype=float) if self.mesh_bounds is not None else np.array([], dtype=float),
        )

    @classmethod
    def load(cls, path: str | Path) -> "GridOccupancyBackend":
        data = np.load(path)
        mesh_bounds = np.asarray(data["mesh_bounds"], dtype=float) if "mesh_bounds" in data else np.array([], dtype=float)
        return cls(
            origin=np.asarray(data["origin"], dtype=float),
            voxel_mm=float(data["voxel_mm"][0]),
            occupied=np.asarray(data["occupied"], dtype=bool),
            mesh_bounds=mesh_bounds.reshape((2, 3)) if mesh_bounds.size else None,
        )

    @classmethod
    def build(
        cls,
        mesh: trimesh.Trimesh,
        voxel_mm: float,
        *,
        padding_voxels: int = 2,
        chunk_points: int = 250_000,
    ) -> "GridOccupancyBackend":
        if voxel_mm <= 0:
            raise ValueError("voxel_mm must be positive")
        bounds = np.asarray(mesh.bounds, dtype=float)
        pad = voxel_mm * padding_voxels
        lo = bounds[0] - pad
        hi = bounds[1] + pad
        axes = [
            np.arange(lo[axis], hi[axis] + voxel_mm * 0.5, voxel_mm, dtype=float)
            for axis in range(3)
        ]
        occupied = np.zeros((len(axes[0]), len(axes[1]), len(axes[2])), dtype=bool)
        yz_count = len(axes[1]) * len(axes[2])
        x_chunk = max(1, chunk_points // max(1, yz_count))
        yy, zz = np.meshgrid(axes[1], axes[2], indexing="ij")
        yz = np.column_stack([yy.ravel(), zz.ravel()])

        for start in range(0, len(axes[0]), x_chunk):
            stop = min(start + x_chunk, len(axes[0]))
            xs = axes[0][start:stop]
            points = np.empty((len(xs) * yz_count, 3), dtype=float)
            points[:, 0] = np.repeat(xs, yz_count)
            points[:, 1:] = np.tile(yz, (len(xs), 1))
            inside = mesh.contains(points).reshape((len(xs), len(axes[1]), len(axes[2])))
            occupied[start:stop, :, :] = inside

        return cls(
            origin=np.array([axes[0][0], axes[1][0], axes[2][0]], dtype=float),
            voxel_mm=voxel_mm,
            occupied=occupied,
            mesh_bounds=bounds,
        )


def build_occupancy_backend(
    mesh: trimesh.Trimesh,
    *,
    kind: str = "mesh",
    grid_voxel_mm: float = 1.0,
    grid_cache: str | Path | None = None,
    rebuild_grid: bool = False,
) -> OccupancyBackend:
    if kind == "mesh":
        return MeshContainsBackend(mesh)
    if kind != "grid":
        raise ValueError(f"unknown occupancy backend {kind!r}")

    cache_path = Path(grid_cache) if grid_cache else None
    if cache_path is not None and cache_path.exists() and not rebuild_grid:
        backend = GridOccupancyBackend.load(cache_path)
        if abs(backend.voxel_mm - grid_voxel_mm) > 1e-9:
            raise ValueError(
                f"grid cache voxel size {backend.voxel_mm} does not match requested {grid_voxel_mm}; use --rebuild-grid"
            )
        if backend.mesh_bounds is not None and not np.allclose(backend.mesh_bounds, mesh.bounds, atol=1e-9):
            raise ValueError("grid cache mesh bounds do not match input mesh; use --rebuild-grid")
        return backend

    backend = GridOccupancyBackend.build(mesh, grid_voxel_mm)
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        backend.save(cache_path)
    return backend
