"""Radial-height calibration bands, independent of temporary travel Z-hop."""
from __future__ import annotations

from dataclasses import dataclass
import math


CALIBRATION_PARAMETERS = ("retract-mm", "retract-speed", "linear-advance-k",
                          "print-speed", "print-accel", "speed-accel")


@dataclass(frozen=True)
class CalibrationSweep:
    parameter: str
    values: tuple[float, ...] | tuple[tuple[float, float], ...]
    base_height: float = 2.0
    band_height: float = 3.0

    def __post_init__(self):
        if self.parameter not in CALIBRATION_PARAMETERS:
            raise ValueError("unknown calibration parameter")
        if self.parameter == "speed-accel":
            if any(not isinstance(v, tuple) or len(v) != 2 for v in self.values):
                raise ValueError("speed-accel values must be speed:acceleration pairs (mm/s:mm/s²)")
            numbers = tuple(n for pair in self.values for n in pair)
        else:
            numbers = self.values
        if not numbers or any(not math.isfinite(v) or v < 0 for v in numbers):
            raise ValueError("calibration values must be finite, non-negative numbers")
        if self.parameter in {"retract-speed", "print-speed", "print-accel", "speed-accel"}:
            if any(v <= 0 for v in numbers):
                raise ValueError("calibration speeds and accelerations must be positive")
        speeds = (tuple(v[0] for v in self.values) if self.parameter == "speed-accel"
                  else self.values if self.parameter in {"retract-speed", "print-speed"} else ())
        if any(not math.isfinite(v * 60) for v in speeds):
            raise ValueError("calibration speed is too large to convert to mm/min")
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
        if self.parameter == "speed-accel":
            speed, accel = self.values[band]
            return {"feed_mm_per_min": speed * 60, "print_accel": accel}
        key = {
            "retract-mm": "retract_mm",
            "retract-speed": "retract_speed_mm_per_min",
            "linear-advance-k": "linear_advance_k",
            "print-speed": "feed_mm_per_min",
            "print-accel": "print_accel",
        }[self.parameter]
        return {key: self.values[band] * (60 if self.parameter in {"retract-speed", "print-speed"} else 1)}

    def _value_text(self, band: int) -> str:
        value = self.values[band]
        return f"{value[0]:g}:{value[1]:g}" if isinstance(value, tuple) else f"{value:g}"

    def comment(self, height: float) -> str:
        band = self.band_at(height)
        if band is None:
            return f";CALIBRATION:base height={height:.4f}\n"
        unit = {"retract-mm": "mm", "retract-speed": "mm/s",
                "linear-advance-k": "K", "print-speed": "mm/s",
                "print-accel": "mm/s²", "speed-accel": "mm/s:mm/s²"}[self.parameter]
        return (f";CALIBRATION:{self.parameter} band={band} "
                f"value={self._value_text(band)} unit={unit} height={height:.4f}\n")

    def header(self) -> str:
        values = ",".join(self._value_text(i) for i in range(len(self.values)))
        return (f";CALIBRATION_SCHEDULE:{self.parameter} values={values} "
                f"base_height={self.base_height:g} band_height={self.band_height:g}\n")


# Preserve the existing Python API and SlicerConfig.retraction_sweep field.
RetractionSweep = CalibrationSweep


def add_calibration_arguments(parser):
    parser.add_argument("--calibration", choices=CALIBRATION_PARAMETERS)
    parser.add_argument("--calibration-values", help="Comma-separated band values: speeds in mm/s, accelerations in mm/s², K unitless; speed-accel uses speed:acceleration pairs.")
    parser.add_argument("--calibration-base-height", type=float, default=2.0,
                        help="Radial height above build cylinder where the first band starts (mm).")
    parser.add_argument("--calibration-band-height", type=float, default=3.0,
                        help="Radial thickness of each band (mm); final value is held above the last band.")


def sweep_from_args(args, parser) -> CalibrationSweep | None:
    if not args.calibration and not args.calibration_values:
        return None
    if not args.calibration or not args.calibration_values:
        parser.error("--calibration and --calibration-values must be supplied together")
    try:
        values = tuple(tuple(float(n) for n in v.split(":"))
                       if args.calibration == "speed-accel" else float(v)
                       for v in args.calibration_values.split(","))
        return CalibrationSweep(args.calibration, values,
                                args.calibration_base_height, args.calibration_band_height)
    except ValueError as exc:
        parser.error(str(exc))
