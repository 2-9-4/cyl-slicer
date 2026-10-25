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


def coalesce_collinear_points(
    points: list[tuple[float, float]],
    radius: float,
    tol: float = 1e-8,
) -> list[tuple[float, float]]:
    """Remove redundant vertices in continuous (theta, axial) machine space.

    The tolerance is in physical mm and only absorbs numerical noise. Corners,
    reversals, endpoints, and unwrapped rotary winding must survive.
    """
    out: list[tuple[float, float]] = []
    for point in points:
        if out and cylindrical_distance(out[-1], point, radius) <= 1e-12:
            continue
        while len(out) >= 2:
            a, b = out[-2:]
            ab = ((b[0] - a[0]) * radius, b[1] - a[1])
            bc = ((point[0] - b[0]) * radius, point[1] - b[1])
            ac_length = cylindrical_distance(a, point, radius)
            cross = ab[0] * bc[1] - ab[1] * bc[0]
            dot = ab[0] * bc[0] + ab[1] * bc[1]
            if dot <= 0 or abs(cross) > tol * ac_length:
                break
            out.pop()
        out.append(point)
    return out


def _forward_delta(theta: float, reference: float, tol: float = 1e-9) -> float:
    delta = (theta - reference) % TWO_PI
    if delta <= tol or abs(delta - TWO_PI) <= tol:
        return 0.0
    return delta


def _nearest_delta(theta: float, reference: float, tol: float = 1e-9) -> float:
    delta = (theta - reference + math.pi) % TWO_PI - math.pi
    if abs(delta) <= tol:
        return 0.0
    if abs(delta + math.pi) <= tol and theta > reference:
        return math.pi
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
    values = [float(t) for t in thetas]
    if not values:
        return []
    out = [current_accum + _nearest_delta(values[0], current_accum)]
    previous_raw = values[0]
    for theta in values[1:]:
        raw_delta = theta - previous_raw
        if raw_delta >= -1e-9:
            delta = max(0.0, raw_delta)
        else:
            delta = _forward_delta(theta, previous_raw)
        out.append(out[-1] + delta)
        previous_raw = theta
    return out
