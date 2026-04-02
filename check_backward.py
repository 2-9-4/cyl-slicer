"""
check_backward.py — detect "thread jump" artifacts in cylindrical slicer G-code.

A thread jump is an extrusion segment that leaps axially (large ΔX) while barely
moving angularly (small Δθ canonical), connecting two different helix threads at
the seam instead of staying on the same thread.

Usage:
    python check_backward.py [output.gcode] [--axial-thresh MM] [--angle-thresh DEG]

Outputs:
  - Per-flag: layer, type, approximate canonical θ, ΔX, Δθ, coordinates
  - Summary count and worst offenders
  - Optional matplotlib plot of flagged segments in (θ, X) space
"""

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

TWO_PI = 2 * math.pi
DEG2RAD = math.pi / 180.0


# ---------------------------------------------------------------------------
# G-code parser
# ---------------------------------------------------------------------------

@dataclass
class Segment:
    """One G1 extrusion segment."""
    layer: int
    kind: str                  # TYPE tag active at this segment
    # Accumulated Y (degrees) — as in the file
    y_start_deg: float
    y_end_deg: float
    x_start: float
    x_end: float

    @property
    def dy_deg(self) -> float:
        return self.y_end_deg - self.y_start_deg

    @property
    def dx(self) -> float:
        return self.x_end - self.x_start

    @property
    def theta_start_deg(self) -> float:
        """Canonical angle [0, 360) at segment start."""
        return self.y_start_deg % 360.0

    @property
    def theta_end_deg(self) -> float:
        """Canonical angle [0, 360) at segment end."""
        return self.y_end_deg % 360.0

    def near_seam(self, window_deg: float = 30.0) -> bool:
        """True if either endpoint is within window_deg of the seam (θ=0/360)."""
        for t in (self.theta_start_deg, self.theta_end_deg):
            if t <= window_deg or t >= 360.0 - window_deg:
                return True
        return False


def parse_gcode(path: Path) -> list[Segment]:
    segments: list[Segment] = []
    current_layer = -9999
    current_kind = "UNKNOWN"
    cur_x: Optional[float] = None
    cur_y: Optional[float] = None  # accumulated degrees

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(";LAYER:"):
                try:
                    current_layer = int(line.split(":")[1])
                except ValueError:
                    pass
                continue

            if line.startswith(";TYPE:"):
                current_kind = line.split(":", 1)[1].strip()
                continue

            if not (line.startswith("G0") or line.startswith("G1")):
                continue

            # Parse X and Y tokens
            tokens = line.split()
            cmd = tokens[0]
            new_x, new_y = cur_x, cur_y
            has_e = False

            for tok in tokens[1:]:
                if tok.startswith("X"):
                    try:
                        new_x = float(tok[1:])
                    except ValueError:
                        pass
                elif tok.startswith("Y"):
                    try:
                        new_y = float(tok[1:])
                    except ValueError:
                        pass
                elif tok.startswith("E"):
                    has_e = True

            # Only record G1 moves with extrusion where both coords are known
            if cmd == "G1" and has_e and cur_x is not None and cur_y is not None:
                if new_x is not None and new_y is not None:
                    if new_x != cur_x or new_y != cur_y:
                        segments.append(Segment(
                            layer=current_layer,
                            kind=current_kind,
                            y_start_deg=cur_y,
                            y_end_deg=new_y,
                            x_start=cur_x,
                            x_end=new_x,
                        ))

            if new_x is not None:
                cur_x = new_x
            if new_y is not None:
                cur_y = new_y

    return segments


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

@dataclass
class Flag:
    seg: Segment
    reason: str


def detect_thread_jumps(
    segments: list[Segment],
    axial_thresh: float = 2.0,
    angle_thresh: float = 30.0,
    seam_window: float = 45.0,
    all_angles: bool = False,
) -> list[Flag]:
    """
    Flag segments that look like thread jumps.

    Criteria:
      - |ΔX| > axial_thresh  (significant axial leap)
      - AND  |Δθ_canonical| < angle_thresh  (barely moved angularly within the
        segment itself, where Δθ is the UNWRAPPED delta — the actual motion)
      - AND (by default) the segment is near the seam (θ ≈ 0/360)

    For a legitimate perimeter, each segment moves the nozzle around the
    cylinder; ΔX per-segment stays small (≤ wall extrusion width, ~0.4 mm).
    A thread-jump segment moves purely axially by ~pitch while staying at the
    same angle.
    """
    flags: list[Flag] = []
    for seg in segments:
        ax = abs(seg.dx)
        if ax < axial_thresh:
            continue

        # Use the raw accumulated delta — no unwrap needed since Y only increases
        ang = abs(seg.dy_deg)
        if ang >= angle_thresh:
            continue

        if not all_angles and not seg.near_seam(seam_window):
            continue
        
        reason = (
            f"|ΔX|={ax:.2f} mm  |Δθ|={ang:.2f}°  "
            f"seam={'yes' if seg.near_seam(seam_window) else 'no'}"
        )
        flags.append(Flag(seg=seg, reason=reason))

    return flags


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report(flags: list[Flag], top_n: int = 20) -> None:
    if not flags:
        print("No thread jumps detected.")
        return

    print(f"\n{'='*70}")
    print(f"  THREAD JUMP REPORT — {len(flags)} segments flagged")
    print(f"{'='*70}\n")

    # Sort by |ΔX| descending
    flags_sorted = sorted(flags, key=lambda f: abs(f.seg.dx), reverse=True)

    for i, flag in enumerate(flags_sorted[:top_n]):
        s = flag.seg
        print(
            f"  [{i+1:3d}] layer={s.layer:5d}  type={s.kind:<20s}\n"
            f"        X: {s.x_start:.3f} → {s.x_end:.3f}  (ΔX={s.dx:+.3f} mm)\n"
            f"        Y_accum: {s.y_start_deg:.2f}° → {s.y_end_deg:.2f}°"
            f"  θ_canon: {s.theta_start_deg:.1f}° → {s.theta_end_deg:.1f}°\n"
            f"        {flag.reason}\n"
        )

    if len(flags) > top_n:
        print(f"  ... and {len(flags) - top_n} more (sorted by |ΔX|)\n")

    # Layer histogram
    from collections import Counter
    by_layer = Counter(f.seg.layer for f in flags)
    print(f"\n  Layers with most flags (top 10):")
    for layer, count in by_layer.most_common(10):
        bar = "#" * min(count, 40)
        print(f"    layer {layer:5d}: {bar} ({count})")

    # Type breakdown
    by_type = Counter(f.seg.kind for f in flags)
    print(f"\n  By type:")
    for kind, count in by_type.most_common():
        print(f"    {kind:<25s} {count}")

    print()


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_flags(
    all_segments: list[Segment],
    flags: list[Flag],
    layer: Optional[int] = None,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        print("matplotlib not available — skipping plot")
        return

    flagged_ids = {id(f.seg): i for i, f in enumerate(flags)}

    # Distinct color per flagged segment using a qualitative colormap
    n_flags = len(flags)
    cmap = cm._colormaps["tab20" if n_flags <= 20 else "hsv"].resampled(max(n_flags, 1))
    flag_colors = [cmap(i / max(n_flags - 1, 1)) for i in range(n_flags)]

    fig, ax = plt.subplots(figsize=(14, 6))

    for seg in all_segments:
        if layer is not None and seg.layer != layer:
            continue
        fid = flagged_ids.get(id(seg))
        if fid is not None:
            color = flag_colors[fid]
            alpha, lw = 1.0, 2.5
        else:
            color, alpha, lw = "steelblue", 0.15, 0.5

        t0 = seg.theta_start_deg
        t1 = t0 + seg.dy_deg
        ax.plot([t0, t1], [seg.x_start, seg.x_end],
                color=color, alpha=alpha, linewidth=lw)

        # Label at segment midpoint
        if fid is not None:
            tmid = (t0 + t1) / 2
            xmid = (seg.x_start + seg.x_end) / 2
            ax.annotate(
                str(fid + 1),
                xy=(tmid, xmid),
                fontsize=7,
                fontweight="bold",
                color=color,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=color,
                          lw=0.8, alpha=0.85),
            )

    ax.axvline(0,   color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(360, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Canonical θ (degrees)")
    ax.set_ylabel("X / axial (mm)")
    title = "Thread-jump detection"
    if layer is not None:
        title += f"  — layer {layer}"
    ax.set_title(title)

    from matplotlib.lines import Line2D
    legend = [Line2D([0], [0], color="steelblue", lw=1, label="Normal")]
    for i, f in enumerate(flags):
        legend.append(Line2D([0], [0], color=flag_colors[i], lw=2,
                             label=f"[{i+1}] {f.seg.kind} ΔX={f.seg.dx:+.1f}"))
    ax.legend(handles=legend, fontsize=7, loc="upper right",
              ncol=max(1, n_flags // 15 + 1))
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Detect thread-jump artifacts in cylindrical slicer G-code."
    )
    parser.add_argument(
        "gcode", nargs="?", default="output.gcode",
        help="Path to G-code file (default: output.gcode)",
    )
    parser.add_argument(
        "--axial-thresh", type=float, default=10.0, metavar="MM",
        help="Min |ΔX| (mm) to flag a segment (default: 2.0)",
    )
    parser.add_argument(
        "--angle-thresh", type=float, default=30.0, metavar="DEG",
        help="Max |Δθ| (deg) within a segment to flag it (default: 30.0)",
    )
    parser.add_argument(
        "--seam-window", type=float, default=45.0, metavar="DEG",
        help="Angular window around seam (θ=0/360) to restrict search (default: 45.0)",
    )
    parser.add_argument(
        "--all-angles", action="store_true",
        help="Check all angles, not just near the seam",
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Show matplotlib plot of flagged segments",
    )
    parser.add_argument(
        "--plot-layer", type=int, default=None, metavar="N",
        help="Restrict plot to a single layer number",
    )
    parser.add_argument(
        "--top", type=int, default=20, metavar="N",
        help="Show top N worst offenders (default: 20)",
    )
    args = parser.parse_args()

    gcode_path = Path(args.gcode)
    if not gcode_path.exists():
        print(f"Error: {gcode_path} not found", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing {gcode_path} …", end=" ", flush=True)
    segments = parse_gcode(gcode_path)
    print(f"{len(segments)} extrusion segments found.")
    flags = detect_thread_jumps(
        segments,
        axial_thresh=args.axial_thresh,
        angle_thresh=args.angle_thresh,
        seam_window=args.seam_window,
        all_angles=args.all_angles,
    )

    report(flags, top_n=args.top)

    if args.plot:
        plot_flags(segments, flags, layer=args.plot_layer)


if __name__ == "__main__":
    main()
