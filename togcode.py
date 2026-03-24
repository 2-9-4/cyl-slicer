from typing import List, Tuple
from shapely.geometry import Polygon, LineString, MultiPolygon, MultiLineString, GeometryCollection, Point
from matplotlib import pyplot
import shapely
import numpy as np
import plotly.graph_objects as go
from vispy import scene
from vispy.scene import Line
from scipy.spatial import cKDTree
Coord = Tuple[float, float]
def cylinder_to_cartesian(theta, z, r):
    """
    Convert cylindrical coords (theta, z) with radius r
    to Cartesian (x, y, z).
    
    theta: array-like, angle in radians
    z: array-like, height
    r: float or array-like, cylinder radius
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return np.stack((x, y, z), axis=-1)
def unwrap_angle(theta, prev_theta):
    """Adjust theta so it's continuous with prev_theta."""
    d = theta - prev_theta
    # shift d into [-pi, pi)
    d = (d + np.pi) % (2*np.pi) - np.pi
    return prev_theta + d

def rotate_full_circle_path(path: 'PathObject', current_pos: tuple) -> 'PathObject':
    """Rotate a full-circle (2π-spanning, constant-z) open path to start at the
    current printer position, avoiding a travel back to θ=0.

    current_pos: (accumulated_theta_radians, x_mm) — the printer's current position.

    Each point's canonical theta is mapped to its forward distance from the
    current canonical angle, then added to the accumulated printer position.
    The result spans [current_pos[0], current_pos[0] + 2π] in accumulated space
    and has is_accumulated=True so to_gcode skips the unwrap step.

    Returns the original path unchanged if it isn't a full-circle radial line.
    """
    if not path.points or len(path.points) < 2:
        return path

    pts = np.asarray(path.points, dtype=float)
    TWO_PI = 2 * np.pi

    # Must span ~2π in theta and be constant z
    theta_span = pts[-1, 0] - pts[0, 0]
    if not np.isclose(theta_span, TWO_PI, atol=0.1):
        return path
    if np.ptp(pts[:, 1]) > 0.01:
        return path

    accum_start = current_pos[0]              # printer's accumulated theta
    canon_start = accum_start % TWO_PI        # canonical [0, 2π)

    # Forward distance of each original point from canon_start around the circle.
    # Points exactly at canon_start get relative=0; all others are in (0, 2π].
    relative = (pts[:, 0] - canon_start + TWO_PI) % TWO_PI

    # Sort by forward distance so the path runs canon_start → canon_start+2π.
    order = np.argsort(relative, kind='stable')
    new_pts = list(zip(accum_start + relative[order], pts[order, 1]))

    # Guarantee exact start and end in accumulated space (argsort alone leaves
    # both the original 0 and 2π endpoints collapsed to the same relative value).
    z0 = float(pts[0, 1])
    if not np.isclose(new_pts[0][0], accum_start, atol=1e-9):
        new_pts.insert(0, (accum_start, z0))
    if not np.isclose(new_pts[-1][0], accum_start + TWO_PI, atol=1e-9):
        new_pts.append((accum_start + TWO_PI, z0))

    result = PathObject(new_pts, path.closed, path.kind, path.radius)
    result.is_accumulated = True
    return result

def rotate_closed_path_to_nearest(path: 'PathObject', current_pos: tuple) -> 'PathObject':
    """Rotate a closed path to start at the point nearest to current position.

    For closed paths (like perimeters), we can start at any point along the path.
    This finds the point closest to the current printer position and rotates the
    path to start there, minimizing travel distance.

    path: PathObject (closed or open perimeter path)
    current_pos: (accumulated_theta_radians, x_mm) - printer's current position

    Returns a new PathObject with rotated points (in ACCUMULATED theta space,
    near current_pos[0]), or the original if not suitable.
    """
    if not path.points or len(path.points) < 2:
        return path

    # Check if this is a perimeter path (either closed or spanning ~2π)
    pts = np.asarray(path.points, dtype=float)
    TWO_PI = 2 * np.pi
    
    theta_span = np.max(pts[:, 0]) - np.min(pts[:, 0])
    is_perimeter = path.closed or theta_span >= TWO_PI * 0.9
    
    if not is_perimeter:
        return path  # Not a perimeter, don't rotate

    current_theta, current_x = current_pos

    # Find the nearest point on the path to current position.
    # For each edge, compute the closest point on that edge segment.
    min_dist = float('inf')
    best_t = 0.0  # Interpolation parameter along best edge
    best_edge_start = 0
    best_closest_t_accum = 0.0  # In accumulated space
    best_closest_x = 0.0

    for i in range(len(pts) - 1):
        t1_orig, x1 = pts[i, 0], pts[i, 1]
        t2_orig, x2 = pts[i + 1, 0], pts[i + 1, 1]
        
        # Compute edge span from original (possibly unwrapped) theta values
        edge_dt_orig = t2_orig - t1_orig
        
        # Convert to canonical for wrap calculation
        t1_canon = t1_orig % TWO_PI
        
        # Find the wrap of t1 that's closest to current_theta
        t1_accum = t1_canon + round((current_theta - t1_canon) / TWO_PI) * TWO_PI
        
        # t2_accum = t1_accum + edge_dt_orig (preserving the original edge span)
        t2_accum = t1_accum + edge_dt_orig
        
        # Ensure the edge contains current_theta by adjusting the wrap if needed
        if edge_dt_orig > 0:
            # Edge goes from t1 to t2 (increasing theta)
            if current_theta < t1_accum:
                t1_accum -= TWO_PI
                t2_accum -= TWO_PI
            elif current_theta > t2_accum:
                t1_accum += TWO_PI
                t2_accum += TWO_PI
        else:
            # Edge goes from t1 to t2 (decreasing theta)
            if current_theta > t1_accum:
                t1_accum += TWO_PI
                t2_accum += TWO_PI
            elif current_theta < t2_accum:
                t1_accum -= TWO_PI
                t2_accum -= TWO_PI
        
        # Vector from (t1_accum, x1) to (t2_accum, x2)
        edge_dt = t2_accum - t1_accum
        edge_dx = x2 - x1
        edge_len_sq = edge_dt * edge_dt * path.radius * path.radius + edge_dx * edge_dx
        
        # Vector from (t1_accum, x1) to current position
        to_curr_dt = current_theta - t1_accum
        to_curr_dx = current_x - x1
        
        # Project current position onto edge (clamped to [0, 1])
        if edge_len_sq > 1e-10:
            t_param = (to_curr_dt * edge_dt * path.radius * path.radius + to_curr_dx * edge_dx) / edge_len_sq
            t_param = max(0, min(1, t_param))
        else:
            t_param = 0
        
        # Closest point on this edge (in accumulated space)
        closest_t_accum = t1_accum + t_param * edge_dt
        closest_x = x1 + t_param * edge_dx
        
        # Distance from current position to closest point
        dtheta = closest_t_accum - current_theta
        dx = closest_x - current_x
        dist_sq = dx * dx + (dtheta * path.radius) ** 2
        
        if dist_sq < min_dist:
            min_dist = dist_sq
            best_t = t_param
            best_edge_start = i
            best_closest_t_accum = closest_t_accum
            best_closest_x = closest_x

    # Insert a new vertex at the closest point if it's not already a vertex
    # Then rotate to start at that point
    if best_t > 0.01 and best_t < 0.99:
        # Insert new vertex at closest point
        new_pts = []
        for i in range(len(pts) - 1):
            # Convert to accumulated space
            t_accum = (pts[i, 0] % TWO_PI) + round((current_theta - (pts[i, 0] % TWO_PI)) / TWO_PI) * TWO_PI
            new_pts.append([t_accum, pts[i, 1]])
            if i == best_edge_start:
                # Insert the new point
                new_pts.append([best_closest_t_accum, best_closest_x])
        # Don't close the loop - preserve the original span
        new_pts.append([new_pts[0][0] + (pts[-1, 0] - pts[0, 0]), new_pts[0][1]])
        pts = np.array(new_pts)
        best_idx = best_edge_start + 1  # The new point is at this index
    else:
        # Closest point is near a vertex, convert all to accumulated space and use that vertex
        new_pts = []
        for i in range(len(pts)):
            t_accum = (pts[i, 0] % TWO_PI) + round((current_theta - (pts[i, 0] % TWO_PI)) / TWO_PI) * TWO_PI
            new_pts.append([t_accum, pts[i, 1]])
        # Preserve the original span
        new_pts.append([new_pts[0][0] + (pts[-1, 0] - pts[0, 0]), new_pts[0][1]])
        pts = np.array(new_pts)

        best_idx = best_edge_start
        if best_t >= 0.5:
            best_idx = best_edge_start + 1
            if best_idx >= len(pts) - 1:
                best_idx = 0

    if best_idx == 0:
        return path  # Already starts at nearest point

    # Build rotated path: start at best_idx, continue to end, then wrap to beginning
    # Output in accumulated theta space (near current_theta)
    rotated_pts = []
    orig_span = pts[-1, 0] - pts[0, 0]  # Original theta span (should be ~2π for perimeters)
    
    # Start at best_idx
    start_theta_accum = pts[best_idx, 0]
    rotated_pts.append([start_theta_accum, pts[best_idx, 1]])
    
    # Go from best_idx+1 to end (excluding the duplicate end point if present)
    last_theta = start_theta_accum
    for i in range(best_idx + 1, len(pts)):
        curr_theta = pts[i, 0]
        # Choose the representation closest to last_theta
        best_candidate = curr_theta
        best_diff = abs(curr_theta - last_theta)
        for wrap in [TWO_PI, -TWO_PI]:
            candidate = curr_theta + wrap
            diff = abs(candidate - last_theta)
            if diff < best_diff:
                best_candidate = candidate
                best_diff = diff
        rotated_pts.append([best_candidate, pts[i, 1]])
        last_theta = best_candidate
    
    # Continue from start to best_idx-1 (wrapping around)
    for i in range(best_idx):
        curr_theta = pts[i, 0]
        # Choose the representation closest to last_theta
        best_candidate = curr_theta
        best_diff = abs(curr_theta - last_theta)
        for wrap in [TWO_PI, -TWO_PI]:
            candidate = curr_theta + wrap
            diff = abs(candidate - last_theta)
            if diff < best_diff:
                best_candidate = candidate
                best_diff = diff
        rotated_pts.append([best_candidate, pts[i, 1]])
        last_theta = best_candidate
    
    # Ensure path ends at start_theta + original_span
    if rotated_pts:
        rotated_pts[-1] = [rotated_pts[0][0] + orig_span, rotated_pts[-1][1]]

    result = PathObject(rotated_pts, path.closed, path.kind, path.radius)
    result.is_accumulated = True
    return result


class PathObject:
    def __init__(self, points: List[Coord], closed: bool, kind: str, radius: float):
        self.points = points  # list of (x,y)
        self.closed = closed
        self.kind = kind      # "polygon_exterior", "polygon_hole", "line", "point"
        self.radius = radius
        self.is_accumulated = False  # True → thetas already in accumulated printer space

    def __str__(self):
        return self.to_gcode()
    
    def dist(self,x2,x1,t2,t1):
        return ((x2 - x1)**2 + ((t2 - t1) * self.radius)**2)**0.5
    
    # def to_gcode(self,
    #              slicing_params,current_pos,current_e) -> str:
    #     """Turn this path into G-code commands (basic FDM style)."""
    #     if not self.points:
    #         return "",current_pos,current_e

    #     if not self.closed and self.dist(current_pos[1],self.points[-1][1],current_pos[0],self.points[-1][0]) < self.dist(current_pos[1],self.points[0][1],current_pos[0],self.points[0][0]):
    #         self.points = self.points[::-1]
    #     g = []
    #     # travel to start
    #     t, x = self.points[0]
    #     t_offset = (current_pos[0] - t + np.pi) % (2*np.pi) - np.pi
    #     g.append(f"G0 X{x:.3f} Y{(t+t_offset)*180/np.pi:.3f} F{slicing_params["travel_mm_per_min"]:.1f}")
    #     e = current_e
    #     last = (t, x)
    #     def seg_move(pt1, pt2):
    #         nonlocal e
    #         t1, x1 = pt1
    #         t2, x2 = pt2
    #         # distance in parameter space (theta treated linear)
    #         fake_dist = ((x2 - x1)**2 + ((t2 - t1)*180/np.pi)**2)**0.5
    #         # distance in real space (theta scaled by radius)
    #         real_dist = self.dist(x2,x1,t2,t1)
    #         e += real_dist * slicing_params["extrusion_per_mm"]
    #         assert(real_dist > 0)
    #         fval = slicing_params["feed_mm_per_min"] * fake_dist / real_dist
    #         return f"G1 X{x2:.3f} Y{(t2+t_offset)*180/np.pi:.3f} E{e:.5f} F{fval:.1f}"

    #     for pt in self.points[1:]:
    #         g.append(seg_move(last, pt))
    #         last = pt

    #     if self.closed and self.points[-1] != self.points[0]:
    #         g.append(seg_move(last, self.points[0]))
    #         last = self.points[0]
    #     return "\n".join(g), last, e
    def to_gcode(self, slicing_params, current_pos, current_e):
        """Turn this path into G-code commands (basic FDM style)."""
        if not self.points:
            return "", current_pos, current_e
        
        pts = np.asarray(self.points, dtype=float)  # shape (N, 2), [theta, x]
        
        # Filter duplicate points (vectorized) then remove collinear interior points.
        if len(pts) > 1:
            dx   = pts[1:, 1] - pts[:-1, 1]
            dt   = (pts[1:, 0] - pts[:-1, 0]) * self.radius
            keep = np.concatenate([[True], np.hypot(dx, dt) > 5e-3])
            pts  = pts[keep]
            if len(pts) < 2:
                return "", current_pos, current_e

        if len(pts) > 2:
            pts = _reduce_collinear(pts, self.radius)

        # current_pos[0] is the printer's ACCUMULATED Y in radians (never resets to
        # canonical [0,2π) — it only increases).  This lets Marlin always move forward.
        # The viewer recovers canonical position with Y_gcode % 360.

        is_accumulated = self.is_accumulated  # set by rotate_*_path helpers

        # For non-accumulated paths: pick the closer endpoint to minimise travel.
        # Skip for accumulated paths — their start was already chosen by the
        # rotation helper and reversing would send the printer 2π away.
        if not is_accumulated:
            current_canonical = current_pos[0] % (2 * np.pi)
            if (not self.closed and len(pts) >= 2 and
                _travel_dist(self.radius, current_pos[1], current_canonical,
                             pts[-1,1], pts[-1,0] % (2*np.pi)) <
                _travel_dist(self.radius, current_pos[1], current_canonical,
                             pts[0,1],  pts[0,0]  % (2*np.pi))):
                pts = pts[::-1]

        # --- Travel (G0) to path start ---
        if is_accumulated:
            # Thetas already in accumulated printer space — use directly.
            t0_accum = pts[0, 0]
        else:
            # Canonical space: advance by shortest arc from current accumulated Y.
            t0_canonical = pts[0, 0]
            accum_y      = current_pos[0]
            dt           = ((t0_canonical - accum_y) + np.pi) % (2 * np.pi) - np.pi
            t0_accum     = accum_y + dt
        x0           = pts[0, 1]
        _retract_mm    = slicing_params.get("retract_mm", 0.0)
        _retract_speed = slicing_params.get("retract_speed_mm_per_min", 1800)
        _min_retract   = slicing_params.get("min_retract_travel_mm", 0.0)
        travel_dist    = self.dist(current_pos[1], x0, current_pos[0], t0_accum)
        if _retract_mm > 0 and travel_dist >= _min_retract:
            g = [";start of path\n",
                 f"G1 E{current_e - _retract_mm:.5f} F{_retract_speed:.1f}",
                 f"G0 X{x0:.3f} Y{t0_accum * 180/np.pi:.3f} F{slicing_params['travel_mm_per_min']:.1f}",
                 f"G1 E{current_e:.5f} F{_retract_speed:.1f}"]
        else:
            g = [";start of path\n",
                 f"G0 X{x0:.3f} Y{t0_accum * 180/np.pi:.3f} F{slicing_params['travel_mm_per_min']:.1f}"]

        if len(pts) == 1:
            g.append(";end of path\n")
            return "\n".join(g) + "\n", (t0_accum, x0), current_e

        # --- Unwrap canonical thetas within the path ---
        # np.unwrap removes >π jumps (handles seam crossings inside a path),
        # then shift so the sequence starts at t0_accum.
        if is_accumulated:
            # Thetas are already continuous in accumulated space — use directly.
            accum_thetas = pts[:, 0]
        else:
            # Canonical space: unwrap then shift to start at t0_accum.
            # Improved unwrap: np.unwrap collapses an exact ±2π step (full circle)
            # to 0.  Detect those steps and restore them so full-circle perimeters
            # aren't lost.
            canon_thetas = pts[:, 0]
            _diffs  = np.diff(canon_thetas)
            _wdiffs = (_diffs + np.pi) % (2 * np.pi) - np.pi
            _fc     = (np.abs(_diffs) > np.pi) & (np.abs(_wdiffs) < 1e-9)
            _wdiffs[_fc] = _diffs[_fc]
            unwrapped    = np.concatenate([[canon_thetas[0]], canon_thetas[0] + np.cumsum(_wdiffs)])
            accum_thetas = unwrapped + (t0_accum - unwrapped[0])

        p1 = pts[:-1]
        p2 = pts[1:]

        # Real distances: use accum_thetas (already seam-unwrapped) so a segment
        # whose canonical endpoints straddle the 0→2π boundary gets the correct
        # arc length rather than the full-circumference (~2πr) "long-way-around"
        # distance.  Using raw pts[:,0] here was the cause of ~6 E-units extruded
        # over a 0.02° move and F=0.4 mm/min on seam-crossing support lines.
        try:
            real_dist = self.dist(p2[:,1], p1[:,1], accum_thetas[1:], accum_thetas[:-1])
        except (TypeError, AttributeError):
            real_dist = np.array([self.dist(p2[i,1], p1[i,1], accum_thetas[i+1], accum_thetas[i])
                                  for i in range(len(p1))])

        # fake_dist uses ACCUMULATED deltas so seam-crossing segments get the
        # correct (small) angular delta for Marlin feed-rate scaling, not the
        # large canonical wrap-around value.
        dt_accum  = np.diff(accum_thetas) * 180 / np.pi   # degrees, accumulated
        dx        = p2[:,1] - p1[:,1]
        fake_dist = np.sqrt(dx*dx + dt_accum*dt_accum)

        safe_real_dist = np.where(real_dist > 1e-10, real_dist, 1e-10)
        e_incr  = real_dist * slicing_params["extrusion_per_mm"] * slicing_params.get("flow_rate", 1.0)
        e_vals  = np.cumsum(e_incr) + current_e
        fvals   = slicing_params["feed_mm_per_min"] * fake_dist / safe_real_dist

        # Guard: warn on any segment with absurd extrusion or crawling feedrate —
        # both are symptoms of a seam-coordinate bug where real_dist is computed
        # as ~2πr instead of the true tiny arc length.
        _exp_e_per_mm = slicing_params["extrusion_per_mm"] * slicing_params.get("flow_rate", 1.0)
        _bad_e = (real_dist > 1e-6) & (e_incr / np.maximum(real_dist, 1e-12) > _exp_e_per_mm * 10)
        _bad_f = (safe_real_dist > 1e-3) & (fvals < slicing_params.get("min_feed_mm_per_min", 10.0))
        if np.any(_bad_e) | np.any(_bad_f):
            import sys
            print(
                f"WARNING to_gcode: degenerate segment — "
                f"e/mm max={( e_incr / np.maximum(real_dist, 1e-12)).max():.4f} "
                f"(expected ≤{_exp_e_per_mm * 10:.4f}), "
                f"feed min={fvals.min():.1f} mm/min",
                file=sys.stderr,
            )

        angles_deg = accum_thetas[1:] * 180 / np.pi   # accumulated degrees for G1

        g.extend(
            f"G1 X{x:.3f} Y{a:.3f} E{e:.5f} F{f:.1f}"
            for (_, x), a, e, f in zip(p2, angles_deg, e_vals, fvals)
        )

        # Handle closed path closing segment
        if self.closed and len(pts) > 1:
            t1_canon, x1 = pts[-1]
            t2_canon, x2 = pts[0]
            if not np.allclose([t1_canon, x1], [t2_canon, x2], atol=1e-10):
                real_dist_close = self.dist(x2, x1, accum_thetas[0], accum_thetas[-1])
                if real_dist_close > 1e-10:
                    current_e = e_vals[-1] + real_dist_close * slicing_params["extrusion_per_mm"] * slicing_params.get("flow_rate", 1.0)
                    dt_close  = (accum_thetas[0] - accum_thetas[-1]) * 180 / np.pi
                    fval_close = (slicing_params["feed_mm_per_min"]
                                  * np.sqrt((x2-x1)**2 + dt_close**2) / real_dist_close)
                    g.append(f"G1 X{x2:.3f} Y{accum_thetas[0]*180/np.pi:.3f} E{current_e:.5f} F{fval_close:.1f}")
                    last = (accum_thetas[0], x2)
                else:
                    last = (accum_thetas[-1], x1)
                    current_e = e_vals[-1] if len(e_vals) else current_e
            else:
                last = (accum_thetas[-1], x1)
                current_e = e_vals[-1] if len(e_vals) else current_e
        else:
            last = (accum_thetas[-1], pts[-1, 1])
            current_e = e_vals[-1] if len(e_vals) else current_e
        g.append(";end of path\n")
        return "\n".join(g) + "\n", last, current_e
    
    def plot(self,color,ax): 
        if not len(self.points):
            return
        pts_3d = self.get_pts_3d()
        ax.plot(pts_3d[:,0],pts_3d[:,1],pts_3d[:,2],color)
    
    def plotly(self, color="blue"):
        if not len(self.points):
            return None
        pts_3d = self.get_pts_3d()
        return go.Scatter3d(
            x=pts_3d[:, 0],
            y=pts_3d[:, 1],
            z=pts_3d[:, 2],
            mode="lines",
            line=dict(color=color, width=4)
        )
    def plot_vis(self, color, scene):
        if not len(self.points):
            return None
        return Line(pos=self.get_pts_3d(), color='red', width=2, parent=scene)
    def get_pts_3d(self):
        if not len(self.points):
            return None
        pt = np.array(self.points)
        pts_3d = cylinder_to_cartesian(pt[:, 0], pt[:, 1], self.radius)
        return pts_3d
        
def _travel_dist(radius, x1, t1_canon, x2, t2_canon):
    """Shortest-arc travel distance between two canonical (θ, x) positions."""
    dt = ((t2_canon - t1_canon) + np.pi) % (2 * np.pi) - np.pi
    return ((x2 - x1) ** 2 + (dt * radius) ** 2) ** 0.5


def optimize_path_order(paths: List[PathObject], current_pos) -> List[PathObject]:
    """Greedy nearest-neighbor reordering to minimize travel distance within a section.

    current_pos is (accumulated_theta_radians, x_mm) — the printer's position
    before this section starts.  Open paths may be traversed in either
    direction; closed paths are left in their original orientation.

    Uses a cKDTree for O(n log n) average complexity instead of O(n²).
    Theta is scaled by radius to convert arc-length to mm for distance queries.
    Three theta wraps (±2π) are included so seam-crossing paths are found correctly.
    """
    if len(paths) <= 1:
        return paths

    n = len(paths)
    TWO_PI = 2 * np.pi
    canon_theta = current_pos[0] % TWO_PI
    cur_x = current_pos[1]
    r = next((p.radius for p in paths if p.points), 1.0)

    # Build KD-tree over scaled (theta*r, x) coordinates.
    # Each path contributes up to 6 points: start+end × 3 wraps.
    kd_pts = []
    kd_path_idx = []  # path index for each kd point
    kd_is_end = []    # True if this kd point represents the end endpoint

    for i, path in enumerate(paths):
        if not path.points:
            continue
        pts = np.asarray(path.points)
        s_tr = (pts[0, 0] % TWO_PI) * r
        s_x  = pts[0, 1]
        e_tr = (pts[-1, 0] % TWO_PI) * r
        e_x  = pts[-1, 1]
        for wrap_r in (-TWO_PI * r, 0.0, TWO_PI * r):
            kd_pts.append([s_tr + wrap_r, s_x])
            kd_path_idx.append(i)
            kd_is_end.append(False)
            if not path.closed:
                kd_pts.append([e_tr + wrap_r, e_x])
                kd_path_idx.append(i)
                kd_is_end.append(True)

    if not kd_pts:
        return paths

    kd_pts_arr = np.array(kd_pts)
    kd_path_idx = np.array(kd_path_idx)
    kd_is_end   = np.array(kd_is_end, dtype=bool)
    tree = cKDTree(kd_pts_arr)

    visited = np.zeros(n, dtype=bool)
    result  = []
    cur_pt  = np.array([canon_theta * r, cur_x])

    for step in range(n):
        # Query with progressively larger k until we find an unvisited path.
        # k=6 covers one path in the best case; multiply up through 24, 96, then
        # fall back to a full linear scan only as a last resort.
        chosen_pi  = -1
        chosen_rev = False
        for k_try in (6, 24, 96, len(kd_pts_arr)):
            k = min(k_try, len(kd_pts_arr))
            _, idxs = tree.query(cur_pt, k=k)
            if np.isscalar(idxs):
                idxs = (idxs,)
            for ki in idxs:
                pi = kd_path_idx[ki]
                if not visited[pi]:
                    chosen_pi  = pi
                    chosen_rev = bool(kd_is_end[ki])
                    break
            if chosen_pi >= 0:
                break

        if chosen_pi < 0:
            # Should never happen, but guard anyway
            for i in range(n):
                if not visited[i]:
                    result.append(paths[i])
                    visited[i] = True
            break

        visited[chosen_pi] = True
        result.append(paths[chosen_pi])

        # Advance to the endpoint we'll finish at after to_gcode reversal logic.
        pts = np.asarray(paths[chosen_pi].points)
        end_pt = pts[0] if chosen_rev else pts[-1]
        cur_pt = np.array([(end_pt[0] % TWO_PI) * r, end_pt[1]])

    return result


def _reduce_collinear(pts: np.ndarray, radius: float, tol: float = 0.005) -> np.ndarray:
    """Remove points collinear (within tol mm) with their neighbours.

    Runs multiple passes until stable so that consecutive collinear runs
    (e.g. densified straight lines) are fully collapsed.
    pts: (N, 2) array of (theta, x) — theta in radians.
    """
    while len(pts) > 2:
        scaled = np.column_stack([pts[:, 0] * radius, pts[:, 1]])
        a = scaled[:-2]
        b = scaled[1:-1]
        c = scaled[2:]
        ab = b - a
        ac = c - a
        ac_len = np.hypot(ac[:, 0], ac[:, 1])
        safe_len = np.where(ac_len > 1e-10, ac_len, 1e-10)
        perp = np.abs(ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]) / safe_len
        keep_mid = (perp >= tol) | (ac_len <= 1e-10)
        keep = np.empty(len(pts), dtype=bool)
        keep[0] = True
        keep[-1] = True
        keep[1:-1] = keep_mid
        if keep.all():
            break
        pts = pts[keep]
    return pts


def shapely_to_paths(geom, radius, max_seg_len=None) -> List[PathObject]:
    """Convert a shapely geometry to a list of PathObjects.

    max_seg_len: if set (mm), long segments are densified into chunks of that
    length.  Defaults to None (no densification) because the printer's
    coordinate space (X=axial, Y=angle) is already linear in (theta, z), so
    intermediate points on a straight segment add no accuracy.
    """

    def densify(coords, closed):
        if not coords:
            return []
        coords = list(coords)
        if max_seg_len is None:
            out = coords[:]
            if closed and out[-1] != out[0]:
                out.append(out[0])
            return out
        out = []
        for i in range(len(coords) - 1):
            a = [coords[i][0] * radius,  coords[i][1]]
            b = [coords[i+1][0] * radius, coords[i+1][1]]
            dx, dy = b[0] - a[0], b[1] - a[1]
            dist = (dx*dx + dy*dy)**0.5
            steps = max(1, int(dist / max_seg_len))
            for s in range(steps):
                t = s / steps
                out.append(((a[0] + t*dx) / radius, a[1] + t*dy))
        out.append(coords[-1])
        if closed and out[-1] != out[0]:
            out.append(out[0])
        return out

    paths = []
    if isinstance(geom, Polygon):
        paths.append(PathObject(densify(list(geom.exterior.coords), True), True, "polygon_exterior", radius))
        for ring in geom.interiors:
            paths.append(PathObject(densify(list(ring.coords), True), True, "polygon_hole", radius))
    elif isinstance(geom, MultiPolygon):
        for poly in geom.geoms:
            paths.extend(shapely_to_paths(poly, radius, max_seg_len))
    elif isinstance(geom, LineString):
        paths.append(PathObject(densify(list(geom.coords), geom.is_ring), geom.is_ring, "line", radius))
    elif isinstance(geom, MultiLineString):
        for line in geom.geoms:
            paths.extend(shapely_to_paths(line, radius, max_seg_len))
    elif isinstance(geom, GeometryCollection):
        for g in geom.geoms:
            paths.extend(shapely_to_paths(g, radius, max_seg_len))
    elif isinstance(geom, List):
        for g in geom:
            paths.extend(shapely_to_paths(g, radius, max_seg_len))
    elif isinstance(geom, Point):
        print("got point, ignoring")
        # assert(False)
    elif geom is None:
        pass
    else:
        print(type(geom))
        assert(False)
    return paths