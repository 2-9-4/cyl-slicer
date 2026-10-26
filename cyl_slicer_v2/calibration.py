"""Radial-height calibration bands, independent of temporary travel Z-hop."""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class RetractionSweep:
    parameter: str
    values: tuple[float, ...]
    base_height: float = 2.0
    band_height: float = 3.0

    def __post_init__(self):
        if self.parameter not in {"retract-mm", "retract-speed", "linear-advance-k"}:
            raise ValueError("calibration must be retract-mm, retract-speed or linear-advance-k")
        if not self.values or any(not math.isfinite(v) or v < 0 for v in self.values):
            raise ValueError("calibration values must be finite, non-negative numbers")
        if self.parameter == "retract-speed" and any(v <= 0 for v in self.values):
            raise ValueError("retraction speeds must be positive (mm/s)")
        if self.parameter == "linear-advance-k" and any(v > 10 for v in self.values):
            raise ValueError("Marlin linear-advance calibration K must be in [0, 10]")
        if not math.isfinite(self.base_height) or self.base_height < 0:
            raise ValueError("calibration base height must be finite and non-negative")
        if not math.isfinite(self.band_height) or self.band_height <= 0:
            raise ValueError("calibration band height must be finite and positive")

    def band_at(self, height: float) -> int | None:
        relative = height - self.base_height
        if relative < -1e-8:
            return None
        # Radial shell arithmetic can put exact boundaries a few ULPs below an integer.
        return min(len(self.values) - 1, max(0, math.floor(relative / self.band_height + 1e-8)))

    def overrides(self, height: float) -> dict[str, float]:
        band = self.band_at(height)
        if band is None:
            return {}
        key = {
            "retract-mm": "retract_mm",
            "retract-speed": "retract_speed_mm_per_min",
            "linear-advance-k": "linear_advance_k",
        }[self.parameter]
        return {key: self.values[band] * (60 if self.parameter == "retract-speed" else 1)}

    def comment(self, height: float) -> str:
        band = self.band_at(height)
        if band is None:
            return f";CALIBRATION:base height={height:.4f}\n"
        unit = {"retract-mm": "mm", "retract-speed": "mm/s",
                "linear-advance-k": "K"}[self.parameter]
        return (f";CALIBRATION:{self.parameter} band={band} "
                f"value={self.values[band]:g} unit={unit} height={height:.4f}\n")

    def header(self) -> str:
        values = ",".join(f"{value:g}" for value in self.values)
        return (f";CALIBRATION_SCHEDULE:{self.parameter} values={values} "
                f"base_height={self.base_height:g} band_height={self.band_height:g}\n")


def add_calibration_arguments(parser):
    parser.add_argument("--calibration", choices=("retract-mm", "retract-speed", "linear-advance-k"))
    parser.add_argument("--calibration-values", help="Comma-separated radial band values; retract speed in mm/s, M900 K unitless.")
    parser.add_argument("--calibration-base-height", type=float, default=2.0,
                        help="Radial height above build cylinder where the first band starts (mm).")
    parser.add_argument("--calibration-band-height", type=float, default=3.0,
                        help="Radial thickness of each band (mm); final value is held above the last band.")


def sweep_from_args(args, parser) -> RetractionSweep | None:
    if not args.calibration and not args.calibration_values:
        return None
    if not args.calibration or not args.calibration_values:
        parser.error("--calibration and --calibration-values must be supplied together")
    try:
        return RetractionSweep(args.calibration,
                               tuple(float(v) for v in args.calibration_values.split(",")),
                               args.calibration_base_height, args.calibration_band_height)
    except ValueError as exc:
        parser.error(str(exc))
