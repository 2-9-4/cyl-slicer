#!/usr/bin/env python3
"""Check for theta gaps in solid infill G-code for a given layer.

Usage: python check_infill_gaps.py [gcode_file] [layer]
Defaults: output.gcode, layer 25
"""

import re, sys, math
import numpy as np
import matplotlib.pyplot as plt

gcode_file = sys.argv[1] if len(sys.argv) > 1 else "output.gcode"
target_layer = int(sys.argv[2]) if len(sys.argv) > 2 else 25

lines = open(gcode_file).readlines()

# --- Parse G-code into per-layer, per-type segments ---
current_layer = None
current_type = None
prev_x, prev_y = None, None

# List of (theta_canon_deg, x_start, x_end) for solid infill in target layer
infill_segs = []

for line in lines:
    l = line.strip()
    if l.startswith(';LAYER:'):
        try:
            current_layer = int(l.split(':')[1])
        except ValueError:
            pass
        current_type = None
    elif l.startswith(';TYPE:'):
        current_type = l.split(':')[1]

    if current_layer != target_layer:
        # Still track position
        if l.startswith('G0') or l.startswith('G1'):
            xm = re.search(r'X([-\d.]+)', l)
            ym = re.search(r'Y([-\d.]+)', l)
            if xm: prev_x = float(xm.group(1))
            if ym: prev_y = float(ym.group(1))
        continue

    if not (l.startswith('G0') or l.startswith('G1')):
        continue

    xm = re.search(r'X([-\d.]+)', l)
    ym = re.search(r'Y([-\d.]+)', l)
    x = float(xm.group(1)) if xm else prev_x
    y = float(ym.group(1)) if ym else prev_y

    if current_type == 'SOLID_INFILL' and 'E' in l and l.startswith('G1'):
        if prev_x is not None and prev_y is not None and x is not None and y is not None:
            theta_canon = y % 360.0  # degrees
            infill_segs.append((theta_canon, prev_x, x))

    if x is not None: prev_x = x
    if y is not None: prev_y = y

print(f"Layer {target_layer} SOLID_INFILL: {len(infill_segs)} extrusion segments")

if not infill_segs:
    print("No solid infill found for this layer.")
    sys.exit(0)

# Each axial infill line is a single G1 segment (constant Y, varying X).
# Filter to axial lines: dY ≈ 0 (same theta), significant dX.
# The canonical theta of each line is its Y mod 360.
axial = [(theta, min(x0,x1), max(x0,x1)) for theta, x0, x1 in infill_segs
         if abs(x1 - x0) > 0.1]  # non-trivial X movement

print(f"Axial infill segments (|ΔX|>0.1): {len(axial)}")

if not axial:
    print("No axial segments found.")
    sys.exit(0)

thetas = np.array(sorted(set(round(t, 4) for t, _, _ in axial)))
print(f"Unique canonical thetas: {len(thetas)}")

gaps = np.diff(thetas)
wrap_gap = (360.0 - thetas[-1]) + thetas[0]
all_gaps = np.append(gaps, wrap_gap)

median_step = float(np.median(gaps))
print(f"Median step: {median_step:.4f}°")

threshold = 1.1 * median_step
large_gaps = [(float(thetas[i]), float(thetas[(i+1) % len(thetas)]), float(g))
              for i, g in enumerate(all_gaps) if g > threshold]

print(f"\nGap histogram (in units of median step):")
buckets = {}
for g in all_gaps:
    n = round(g / median_step)
    buckets[n] = buckets.get(n, 0) + 1
for n in sorted(buckets):
    bar = '#' * min(buckets[n], 60)
    print(f"  {n:2d}x ({n*median_step:.3f}°): {buckets[n]:4d}  {bar}")

if large_gaps:
    print(f"\nGaps > {threshold:.3f}° (>{1.6:.1f}x step): {len(large_gaps)}")
    for t0, t1, g in sorted(large_gaps, key=lambda x: -x[2])[:30]:
        print(f"  θ {t0:.3f}° → {t1:.3f}°  gap={g:.3f}°  ({g/median_step:.1f}x step)")
else:
    print(f"\nNo gaps > 1.6x step. Coverage looks uniform.")

plt.plot(all_gaps,'o-')
plt.show()