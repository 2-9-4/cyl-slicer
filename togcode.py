from typing import List, Tuple
from shapely.geometry import Polygon, LineString, MultiPolygon, MultiLineString, GeometryCollection, Point
from matplotlib import pyplot
import shapely
import numpy as np
import plotly.graph_objects as go
from vispy import scene
from vispy.scene import Line
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
class PathObject:
    def __init__(self, points: List[Coord], closed: bool, kind: str, radius: float):
        self.points = points  # list of (x,y)
        self.closed = closed
        self.kind = kind      # "polygon_exterior", "polygon_hole", "line", "point"
        self.radius = radius

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
        
        # Filter out duplicate points first
        if len(pts) > 1:
            # Calculate real distances between consecutive points
            p1 = pts[:-1]
            p2 = pts[1:]
            real_dists = np.array([self.dist(p2[i,1], p1[i,1], p2[i,0], p1[i,0]) 
                                for i in range(len(p1))])
            
            # Keep points where distance > 0 (plus always keep first point)
            keep_mask = np.concatenate([[True], real_dists > 1e-10])  # Small epsilon for floating point
            pts = pts[keep_mask]
            
            # If we filtered everything except first point, return early
            if len(pts) < 2:
                return "", current_pos, current_e
        
        # current_pos[0] is the printer's ACCUMULATED Y in radians (never resets to
        # canonical [0,2π) — it only increases).  This lets Marlin always move forward.
        # The viewer recovers canonical position with Y_gcode % 360.

        # For the closest-endpoint check use canonical current angle.
        current_canonical = current_pos[0] % (2 * np.pi)
        if (not self.closed and len(pts) >= 2 and
            self.dist(current_pos[1], pts[-1,1], current_canonical, pts[-1,0]) <
            self.dist(current_pos[1], pts[0,1], current_canonical, pts[0,0])):
            pts = pts[::-1]

        # --- Travel (G0) to path start ---
        # Advance by the shortest FORWARD arc from current accumulated Y to t0.
        t0_canonical = pts[0, 0]            # canonical radians [0, 2π)
        accum_y      = current_pos[0]       # printer's accumulated Y (radians)
        dt           = ((t0_canonical - accum_y) + np.pi) % (2 * np.pi) - np.pi  # shortest-path, [-π, π)
        t0_accum     = accum_y + dt
        x0           = pts[0, 1]
        _retract_mm    = slicing_params.get("retract_mm", 0.0)
        _retract_speed = slicing_params.get("retract_speed_mm_per_min", 1800)
        if _retract_mm > 0:
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
        canon_thetas  = pts[:, 0]
        unwrapped     = np.unwrap(canon_thetas)                       # still starts near t0_canonical
        accum_thetas  = unwrapped + (t0_accum - unwrapped[0])        # shift to match t0_accum

        p1 = pts[:-1]
        p2 = pts[1:]

        # Real distances use canonical thetas (correct arc length regardless of accumulation)
        try:
            real_dist = self.dist(p2[:,1], p1[:,1], p2[:,0], p1[:,0])
        except (TypeError, AttributeError):
            real_dist = np.array([self.dist(p2[i,1], p1[i,1], p2[i,0], p1[i,0])
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

        angles_deg = accum_thetas[1:] * 180 / np.pi   # accumulated degrees for G1

        gcode_lines = []
        for (t, x), e_val, fval, angle_deg in zip(p2, e_vals, fvals, angles_deg):
            gcode_lines.append(f"G1 X{x:.3f} Y{angle_deg:.3f} E{e_val:.5f} F{fval:.1f}")
        g.extend(gcode_lines)

        # Handle closed path closing segment
        if self.closed and len(pts) > 1:
            t1_canon, x1 = pts[-1]
            t2_canon, x2 = pts[0]
            if not np.allclose([t1_canon, x1], [t2_canon, x2], atol=1e-10):
                real_dist_close = self.dist(x2, x1, t2_canon, t1_canon)
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
    """
    if len(paths) <= 1:
        return paths

    remaining = list(range(len(paths)))
    result = []
    canon_theta = current_pos[0] % (2 * np.pi)
    x = current_pos[1]

    while remaining:
        best_idx = None
        best_dist = float("inf")
        best_reversed = False

        for idx in remaining:
            path = paths[idx]
            if not path.points:
                continue
            pts = np.asarray(path.points)
            s = pts[0]   # (theta, x)
            e = pts[-1]  # (theta, x)
            r = path.radius

            d_start = _travel_dist(r, x, canon_theta, s[1], s[0] % (2 * np.pi))
            if path.closed:
                d, rev = d_start, False
            else:
                d_end = _travel_dist(r, x, canon_theta, e[1], e[0] % (2 * np.pi))
                if d_end < d_start:
                    d, rev = d_end, True
                else:
                    d, rev = d_start, False

            if d < best_dist:
                best_dist = d
                best_idx = idx
                best_reversed = rev

        if best_idx is None:
            for idx in remaining:
                result.append(paths[idx])
            break

        remaining.remove(best_idx)
        result.append(paths[best_idx])

        # Advance position to the endpoint we'll finish at.
        pts = np.asarray(paths[best_idx].points)
        end_pt = pts[0] if best_reversed else pts[-1]
        canon_theta = end_pt[0] % (2 * np.pi)
        x = end_pt[1]

    return result


def shapely_to_paths(geom, radius, max_seg_len=0.5) -> List[PathObject]:
    """Convert a shapely geometry to a list of PathObjects (densified)."""

    def densify(coords, closed):
        if(len(coords) == 0):
            return []
        out = []
        for i in range(len(coords)-1):
            a, b = list(coords[i]), list(coords[i+1])
            a[0] *= radius
            b[0] *= radius
            dx, dy = b[0]-a[0], b[1]-a[1]
            dist = (dx*dx + dy*dy)**0.5
            steps = max(1, int(dist/max_seg_len))
            for s in range(steps):
                t = s/steps
                out.append(((a[0]+t*dx)/radius, a[1]+t*dy))
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