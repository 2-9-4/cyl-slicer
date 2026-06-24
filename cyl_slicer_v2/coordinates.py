from __future__ import annotations

import math
from typing import Iterable

import numpy as np

TWO_PI = 2.0 * math.pi


def wrap_theta(theta: float) -> float:
    return theta % TWO_PI


def cylindrical_to_cartesian(theta: np.ndarray, axial: np.ndarray, radius: float) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    axial = np.asarray(axial, dtype=float)
    return np.stack(
        [radius * np.cos(theta), radius * np.sin(theta), axial],
        axis=-1,
    )


def cartesian_to_cylindrical(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    radius = np.hypot(points[:, 0], points[:, 1])
    theta = np.mod(np.arctan2(points[:, 1], points[:, 0]), TWO_PI)
    axial = points[:, 2]
    return np.stack([radius, theta, axial], axis=-1)


def cylindrical_distance(p0: tuple[float, float], p1: tuple[float, float], radius: float) -> float:
    t0, x0 = p0
    t1, x1 = p1
    return math.hypot((t1 - t0) * radius, x1 - x0)


def _forward_delta(theta: float, reference: float, tol: float = 1e-9) -> float:
    delta = (theta - reference) % TWO_PI
    if delta <= tol or abs(delta - TWO_PI) <= tol:
        return 0.0
    return delta


def unwrap_forward(thetas: Iterable[float], start: float | None = None, tol: float = 1e-9) -> list[float]:
    values = [float(t) for t in thetas]
    if not values:
        return []
    if start is None:
        out = [values[0]]
    else:
        first = start + _forward_delta(values[0], start, tol)
        out = [first]
    previous_raw = values[0]
    for theta in values[1:]:
        raw_delta = theta - previous_raw
        if raw_delta >= -tol:
            delta = max(0.0, raw_delta)
        else:
            delta = _forward_delta(theta, previous_raw, tol)
        out.append(out[-1] + delta)
        previous_raw = theta
    return out


def unwrap_path_from_current(thetas: Iterable[float], current_accum: float) -> list[float]:
    return unwrap_forward(thetas, current_accum)
