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
        
        # If open path, ensure we start closest to current_pos
        if (not self.closed and len(pts) >= 2 and
            self.dist(current_pos[1], pts[-1,1], current_pos[0], pts[-1,0]) <
            self.dist(current_pos[1], pts[0,1], current_pos[0], pts[0,0])):
            pts = pts[::-1]

        # Travel move to start - use continuous theta for Marlin compatibility
        t0, x0 = pts[0]
        # Calculate offset to make theta continuous from current position
        t_offset = (current_pos[0] - t0 + np.pi) % (2*np.pi) - np.pi
        # Output Y as continuous angle in degrees (can exceed 360)
        t0_cont = t0 + t_offset
        t0_deg = t0_cont * 180 / np.pi

        g = [";start of path\n",f"G0 X{x0:.3f} Y{t0_deg:.3f} F{slicing_params['travel_mm_per_min']:.1f}"]
        
        if len(pts) == 1:
            g.append(";end of path\n")
            return "\n".join(g) + "\n", tuple(pts[0]), current_e
        
        # Vectorized segment calculations
        p1 = pts[:-1]
        p2 = pts[1:]
        
        # Calculate all real distances in one vectorized call if possible
        # (This assumes self.dist can be vectorized - if not, keep the loop)
        try:
            # Try vectorized version first
            real_dist = self.dist(p2[:,1], p1[:,1], p2[:,0], p1[:,0])
        except (TypeError, AttributeError):
            # Fall back to loop if self.dist doesn't support vectorization
            real_dist = np.array([self.dist(p2[i,1], p1[i,1], p2[i,0], p1[i,0])
                                for i in range(len(p1))])
        
        # Vectorized calculations
        dt = (p2[:,0] - p1[:,0]) * 180/np.pi
        dx = p2[:,1] - p1[:,1]
        fake_dist = np.sqrt(dx*dx + dt*dt)
        
        # Avoid division by zero
        safe_real_dist = np.where(real_dist > 1e-10, real_dist, 1e-10)
        
        # Extrusion and feed calculations
        e_incr = real_dist * slicing_params["extrusion_per_mm"]
        e_vals = np.cumsum(e_incr) + current_e
        fvals = slicing_params["feed_mm_per_min"] * fake_dist / safe_real_dist
        
        # Pre-calculate angle conversions - continuous theta for Marlin
        angles_deg = (p2[:,0] + t_offset) * 180 / np.pi

        # Build G-code strings efficiently
        gcode_lines = []
        for i, ((t, x), e_val, fval, angle_deg) in enumerate(zip(p2, e_vals, fvals, angles_deg)):
            gcode_lines.append(f"G1 X{x:.3f} Y{angle_deg:.3f} E{e_val:.5f} F{fval:.1f}")

        g.extend(gcode_lines)

        # Handle closed path
        if self.closed and len(pts) > 1:
            t1, x1 = pts[-1]
            t2, x2 = pts[0]

            # Only add closing segment if endpoints are different
            if not np.allclose([t1, x1], [t2, x2], atol=1e-10):
                real_dist_close = self.dist(x2, x1, t2, t1)
                if real_dist_close > 1e-10:  # Only if there's actual distance
                    current_e = e_vals[-1] + real_dist_close * slicing_params["extrusion_per_mm"]
                    fake_dist_close = np.sqrt((x2-x1)**2 + ((t2-t1)*180/np.pi)**2)
                    fval_close = slicing_params["feed_mm_per_min"] * fake_dist_close / real_dist_close
                    t2_cont = t2 + t_offset
                    g.append(f"G1 X{x2:.3f} Y{t2_cont*180/np.pi:.3f} E{current_e:.5f} F{fval_close:.1f}")
                    last = (t2, x2)
                else:
                    last = (t1, x1)
                    current_e = e_vals[-1] if len(e_vals) else current_e
            else:
                last = (t1, x1)
                current_e = e_vals[-1] if len(e_vals) else current_e
        else:
            last = tuple(pts[-1])
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