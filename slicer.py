import sys
import trimesh
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, box
import numpy as np
import plotly.graph_objects as go
import shapely
from shapely.plotting import plot_polygon, plot_line
from collections import defaultdict, deque
import time
import math
from vispy.scene import visuals
from togcode import shapely_to_paths, optimize_path_order, rotate_full_circle_path, rotate_closed_path_to_nearest
from concurrent.futures import ThreadPoolExecutor
from plotutil import set_axes_equal

from vispy import scene
from vispy.scene import Line
from vispy import app, gloo
import solid2
def find_polygons(edges, num_points):
    # Build undirected edge map
    adjacency = defaultdict(list)
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)

    # Traverse to extract polygons
    visited_edges = set()
    polygons = []

    def edge_key(a, b):
        return tuple(sorted((a, b)))

    for start in range(num_points):
        neighbors = adjacency[start]
        for neighbor in neighbors:
            ek = edge_key(start, neighbor)
            if ek in visited_edges:
                continue

            # Begin a new polygon
            poly = [start]
            curr = neighbor
            came_from = ek  # track the specific edge we arrived on, not just prev vertex
            visited_edges.add(ek)

            while curr != start:
                poly.append(curr)
                for next_node in adjacency[curr]:
                    ek2 = edge_key(curr, next_node)
                    if ek2 == came_from or ek2 in visited_edges:
                        continue
                    visited_edges.add(ek2)
                    came_from = ek2
                    curr = next_node
                    break
                else:
                    break  # degenerate topology — abandon this polygon

            if curr == start:  # only keep properly closed loops
                polygons.append(poly)

    return polygons

def shift_polygon(poly, dx):
    return shapely.transform(poly,lambda x: x + [dx,0])

def build_hierarchy(polygons, period=2 * np.pi):
    from shapely import STRtree
    hierarchy = [[] for _ in polygons]
    if len(polygons) < 2:
        return hierarchy
    tree = STRtree(polygons)
    for i, outer_poly in enumerate(polygons):
        seen = set()
        for shift in [0, -period, period]:
            shifted = shift_polygon(outer_poly, shift) if shift != 0 else outer_poly
            for j in tree.query(shifted, predicate='contains'):
                if j != i and j not in seen:
                    seen.add(j)
                    hierarchy[i].append(j)
    return hierarchy

def unwrap_polygon(polygon, period=2*np.pi):
    coords = np.array(polygon)
    theta = coords[:, 0]
    
    # Unwrap theta to remove jumps across 0/2pi
    theta_unwrapped = np.unwrap(theta, discont=np.pi)  # default jump is pi
    coords[:, 0] = theta_unwrapped
    coords -= [np.floor(max(0.0, np.min(coords[:,0]))/(2*np.pi))*2*np.pi,0]
    return coords

def find_boundary_edges(triangles):
    """
    Given a (n, 3) array of triangle vertex indices,
    return a list of boundary edges (edges that appear only once).
    """
    edge_list = []

    for tri in triangles:
        # Generate edges with sorted endpoints to avoid direction mismatches
        edges = [tuple(sorted((tri[i], tri[(i+1) % 3]))) for i in range(3)]
        edge_list.extend(edges)

    # Count edge occurrences
    edge_counts = Counter(edge_list)

    # Keep only edges that occur once (i.e., boundary edges)
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    return boundary_edges

def cartesian_to_cylindrical(points):
    """
    Convert Nx3 numpy array of Cartesian (x, y, z) points 
    to cylindrical coordinates (r, theta, z)
    
    theta is in radians, in the range [-pi, pi]
    """
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    
    return np.stack((r, theta, z), axis=-1)

def find_connected_shapes(shapes, triangles):
    # Build a map from vertex index to shape index
    vertex_to_shape = {}
    for shape_idx, shape in enumerate(shapes):
        for v in shape:
            vertex_to_shape[v] = shape_idx

    # Build graph of shape connections
    shape_graph = defaultdict(set)
    for tri in triangles:
        # Get the set of unique shape indices this triangle touches
        shape_indices = set(vertex_to_shape[v] for v in tri if v in vertex_to_shape)
        # Connect all pairs of shapes in this triangle
        for s1 in shape_indices:
            for s2 in shape_indices:
                if s1 != s2:
                    shape_graph[s1].add(s2)
                    shape_graph[s2].add(s1)

    # Find connected components of the shape graph
    visited = set()
    connected_components = []

    def dfs(node, component):
        visited.add(node)
        component.append(node)
        for neighbor in shape_graph[node]:
            if neighbor not in visited:
                dfs(neighbor, component)

    for shape_idx in range(len(shapes)):
        if shape_idx not in visited:
            component = []
            dfs(shape_idx, component)
            connected_components.append(component)

    return connected_components

def is_full_circle(polygon, tol=0.1):
    def total_theta_sweep(polygon, period=2*np.pi):
        thetas = np.array(polygon)[:, 0]
        diffs = np.diff(thetas, append=thetas[0])  # close the loop
        diffs = (diffs + np.pi) % (2 * np.pi) - np.pi  # wrap to [-π, π]
        return np.abs(diffs.sum())
    return np.isclose(total_theta_sweep(polygon), 2*np.pi, atol=tol)


def cylindrical_cross_section(mesh, radius):
    """
    Compute the cross-section of mesh at cylindrical radius using directed marching triangles.

    Uses face winding order to orient each crossing segment (inside→outside = FROM,
    outside→inside = TO).  After building the directed graph, T-junction gaps introduced
    by subdivide_to_size (boundary edges with no manifold partner) are stitched by
    nearest-neighbour matching of open chain endpoints.

    Returns (pts_cartesian, edges) where pts_cartesian is (N,3) Cartesian points
    at r≈radius and edges is a list of (i,j) index pairs forming closed loops.
    """
    from scipy.spatial import cKDTree

    verts = mesh.vertices
    r_vals = np.sqrt(verts[:, 0]**2 + verts[:, 1]**2)

    pts = []
    edge_cache = {}

    def get_crossing(vi, vj):
        key = (min(vi, vj), max(vi, vj))
        if key not in edge_cache:
            ri, rj = r_vals[vi], r_vals[vj]
            t = (radius - ri) / (rj - ri)
            edge_cache[key] = len(pts)
            pts.append(verts[vi] + t * (verts[vj] - verts[vi]))
        return edge_cache[key]

    # next_of[from_idx] = to_idx encodes the directed polygon graph.
    # Each node has at most one outgoing edge; ambiguous nodes (mesh self-intersections)
    # are marked None and skipped during loop extraction.
    next_of = {}

    for v0, v1, v2 in mesh.faces:
        r0, r1, r2 = r_vals[v0], r_vals[v1], r_vals[v2]
        in0, in1, in2 = r0 < radius, r1 < radius, r2 < radius
        if (in0 == in1 == in2):
            continue  # all inside or all outside

        from_idx = to_idx = None
        for a, b in ((v0, v1), (v1, v2), (v2, v0)):
            a_in = r_vals[a] < radius
            b_in = r_vals[b] < radius
            if a_in and not b_in:      # inside → outside: start of contour segment
                from_idx = get_crossing(a, b)
            elif not a_in and b_in:    # outside → inside: end of contour segment
                to_idx = get_crossing(a, b)

        if from_idx is None or to_idx is None or from_idx == to_idx:
            continue

        if from_idx in next_of:
            next_of[from_idx] = None   # ambiguous: two segments leave same point
        else:
            next_of[from_idx] = to_idx

    # Stitch open chain endpoints caused by T-junctions from subdivide_to_size.
    # T-junctions create boundary edges whose crossing points have no manifold partner:
    #   open tail  = from_idx that is never a to_idx (chain source, no incoming edge)
    #   open head  = to_idx that is never a from_idx (chain sink, no outgoing edge)
    # Nearest-neighbour matching in 3-D re-connects these broken links.
    if pts:
        pts_arr = np.array(pts)
        valid_keys = {k for k, v in next_of.items() if v is not None}
        valid_vals = {v for v in next_of.values() if v is not None}
        open_tails = list(valid_keys - valid_vals)   # no incoming edge
        open_heads = list(valid_vals - valid_keys)   # no outgoing edge

        # Iterative stitching: greedy nearest-neighbour can fail to match a pair when
        # two tails tie for the same head (one grabs it, leaving the other's true partner
        # unclaimed).  Re-running with the updated open set resolves these cascades.
        stitch_tol = radius * 0.5   # generous — only exact T-junction pairs will be close
        prev_n_open = float('inf')
        while open_tails and open_heads and len(open_tails) < prev_n_open:
            prev_n_open = len(open_tails)
            tail_coords = pts_arr[open_tails]
            head_coords = pts_arr[open_heads]
            tree = cKDTree(head_coords)
            dists, idxs = tree.query(tail_coords, distance_upper_bound=stitch_tol)
            stitched_heads = set()
            for tail, dist, hi in zip(open_tails, dists, idxs):
                if dist < stitch_tol and hi < len(open_heads):
                    head = open_heads[hi]
                    if head not in next_of and head not in stitched_heads:
                        next_of[head] = next_of[tail]
                        next_of[tail] = head
                        stitched_heads.add(head)
            # Recompute open endpoints for the next iteration
            valid_keys = {k for k, v in next_of.items() if v is not None}
            valid_vals = {v for v in next_of.values() if v is not None}
            open_tails = list(valid_keys - valid_vals)
            open_heads = list(valid_vals - valid_keys)
    else:
        pts_arr = np.empty((0, 3))

    # Trace closed loops in the functional directed graph.
    state = {}   # 0 unvisited | 1 in-progress | 2 done
    closed_edges = []

    for start in next_of:
        if state.get(start, 0) != 0:
            continue
        chain = []
        pos = {}
        curr = start
        while (next_of.get(curr) is not None
               and state.get(curr, 0) == 0
               and curr not in pos):
            state[curr] = 1
            pos[curr] = len(chain)
            chain.append(curr)
            curr = next_of[curr]

        if curr in pos:  # closed loop — curr is the cycle entry point
            cycle = chain[pos[curr]:]
            for i in range(len(cycle)):
                closed_edges.append((cycle[i], cycle[(i + 1) % len(cycle)]))

        for node in chain:
            state[node] = 2

    return pts_arr, closed_edges


def fast_extract_faces(cut_mesh, mask):
    mask = np.asarray(mask)
    mask_set = set(mask)

    # Step 1: Get all faces involving any of the vertices in the mask
    vertex_faces = cut_mesh.vertex_faces[mask]
    vertex_faces = vertex_faces[vertex_faces != -1]

    # Step 2: Remove duplicate face indices
    face_indices = np.unique(vertex_faces)

    # Step 3: Build a mask for which of those faces have ALL vertices in `mask`
    face_verts = cut_mesh.faces[face_indices]
    keep_mask = np.all(np.isin(face_verts, mask), axis=1)

    # Step 4: Extract the actual face data
    triangles = face_verts[keep_mask]

    return triangles

def trim_polygon_left_right(poly: Polygon, percent: float = 0.05) -> Polygon:
    # Get bounding box
    minx, miny, maxx, maxy = poly.bounds
    
    # Width of the polygon
    width = maxx - minx
    margin = percent * width
    
    # New bounding box excluding left and right margins
    trimmed_box = box(minx + margin, miny, maxx - margin, maxy)
    
    # Intersect polygon with trimmed box
    trimmed_poly = poly.intersection(trimmed_box)
    
    return trimmed_poly
def fix_mesh(mesh1):
    mesh1.fix_normals()
    mesh1.remove_duplicate_faces()
    mesh1.remove_degenerate_faces()
    mesh1.remove_infinite_values()
    mesh1.remove_unreferenced_vertices()
    mesh1.fill_holes()  
    mesh1.merge_vertices()
    return mesh1.process(validate=True)

def _collect_lines(geom, result, min_length=0.0):
    """Append all non-empty LineStrings from any geometry type into result."""
    if geom is None or geom.is_empty:
        return
    if isinstance(geom, shapely.LineString):
        if geom.length >= min_length:
            result.append(geom)
    elif hasattr(geom, 'geoms'):
        for g in geom.geoms:
            _collect_lines(g, result, min_length)


def _snap_radial_to_full_circle(seg, period=2 * np.pi, tol=0.02):
    """Snap a constant-z radial segment to [0, period] if it nearly spans a full circle.

    Root cause: when a support/infill polygon has a tiny gap at the θ=0/2π seam
    (a numerical artefact from fold_to_canonical or ±2π union ops), the intersection
    with the horizontal line [[0,z],[2π,z]] produces a single LineString whose
    endpoints are just inside [0, 2π] (e.g. θ=0.000175 → θ=6.28299).
    stitch_radial_seam ignores this (tol=1e-9), so PathObject stores it with
    pts = [(6.28299,z),(0.000175,z)] after to_gcode reversal.  np.unwrap then maps
    the −6.28 rad diff to the shortest arc (+0.000364 rad), producing a 0.019° move
    that triggers the seam-cross real_dist bug.  Snapping to [0, 2π] instead lets
    rotate_full_circle_path handle it correctly as a genuine full-circle arc.

    tol: tolerated gap in radians (default 0.02 rad ≈ 1.1°).
    """
    if not isinstance(seg, shapely.geometry.LineString) or seg.is_empty:
        return seg
    coords = np.array(seg.coords)
    if len(coords) < 2 or np.ptp(coords[:, 1]) > 1e-4:  # not a constant-z segment
        return seg
    if coords[:, 0].max() - coords[:, 0].min() < period - tol:
        return seg
    return shapely.geometry.LineString([[0.0, coords[0, 1]], [period, coords[0, 1]]])

def clean_support_region(geom_theta_z, r_base, min_area_mm2, merge_gap_mm):
    """
    1. Inflate sub-polygons below min_area_mm2 to at least that area.
    2. Merge patches closer than merge_gap_mm by buffering then eroding.
    All size parameters are in real mm / mm².  Geometry in (θ,z) space.
    """
    if geom_theta_z is None or geom_theta_z.is_empty:
        return geom_theta_z

    # Work in metric space: scale θ axis by r_base so distances are in mm
    metric = shapely.affinity.scale(geom_theta_z, xfact=r_base, yfact=1.0, origin=(0, 0))

    if min_area_mm2 > 0:
        parts = list(metric.geoms) if hasattr(metric, 'geoms') else [metric]
        inflated = []
        for part in parts:
            if part.is_empty:
                continue
            if part.area < min_area_mm2:
                r_target = (min_area_mm2 / np.pi) ** 0.5
                r_current = (part.area / np.pi) ** 0.5
                part = part.buffer(r_target - r_current)
            inflated.append(part)
        if not inflated:
            return shapely.GeometryCollection()
        metric = shapely.unary_union(inflated)

    if merge_gap_mm > 0:
        metric = metric.buffer(merge_gap_mm / 2).buffer(-(merge_gap_mm / 2))

    return shapely.affinity.scale(metric, xfact=1.0 / r_base, yfact=1.0, origin=(0, 0))



def stitch_radial_seam(segments, period=2*np.pi, tol=1e-9):
    """
    Merge radial (constant-z) infill segments split at the theta=0 seam.

    A radial line [0, 2pi] intersecting a seam-crossing polygon returns two
    pieces: one ending at theta~2pi and one starting at theta~0 (physically
    adjacent). Shift the theta=0 piece by +period so endpoints align, then
    linemerge into one continuous LineString.
    """
    lines_in = [g for g in segments
                if isinstance(g, shapely.geometry.LineString) and not g.is_empty]
    if len(lines_in) <= 1:
        return lines_in

    from shapely.ops import linemerge as _lm
    left, right, other = [], [], []
    for ln in lines_in:
        coords = np.array(ln.coords)
        if np.min(coords[:, 0]) < tol:
            shifted = coords.copy()
            shifted[:, 0] += period
            left.append(shapely.geometry.LineString(shifted))
        elif np.max(coords[:, 0]) > period - tol:
            right.append(ln)
        else:
            other.append(ln)

    if not left:
        return lines_in

    to_merge = [l for l in right + left if not l.is_empty]
    merged = _lm(shapely.geometry.MultiLineString(to_merge))
    if isinstance(merged, shapely.geometry.LineString):
        return [merged] + other
    return list(merged.geoms) + other


def stitch_axial_gap(segments, max_gap):
    """
    Merge axial (constant-θ) LineString segments whose z-gap is ≤ max_gap.

    Mesh-intersection and buffer artifacts occasionally leave a tiny notch in the
    interior polygon at a specific z, splitting one axial line into two close
    segments.  Bridging gaps up to max_gap (typically nozzle_size) eliminates the
    resulting travel move / visible gap without merging intentional voids.
    """
    lines = [g for g in segments
             if isinstance(g, shapely.geometry.LineString) and not g.is_empty]
    if len(lines) <= 1:
        return lines

    def z_range(seg):
        coords = np.array(seg.coords)
        return coords[:, 1].min(), coords[:, 1].max()

    lines.sort(key=lambda l: z_range(l)[0])

    result = [lines[0]]
    for seg in lines[1:]:
        prev_z_min, prev_z_max = z_range(result[-1])
        seg_z_min, seg_z_max = z_range(seg)
        if seg_z_min - prev_z_max <= max_gap:
            theta = np.array(result[-1].coords)[0, 0]
            result[-1] = shapely.LineString(
                [[theta, prev_z_min], [theta, max(prev_z_max, seg_z_max)]]
            )
        else:
            result.append(seg)
    return result


def fold_to_canonical(geom, clip_box, period=2 * np.pi):
    """Fold (θ,z)-space geometry into [0, 2π] by unioning ±period shifted copies, then clipping."""
    if geom is None or geom.is_empty:
        return geom
    geom = shapely.make_valid(geom)
    geom = shapely.set_precision(geom, 1e-8)
    folded = shapely.unary_union([
        geom,
        shift_polygon(geom, period),
        shift_polygon(geom, -period),
    ])
    return shapely.set_precision(folded, 1e-8).intersection(clip_box)


def _make_sine_wave_ls(theta_center, z_min, z_max, amplitude_rad, period_mm, pts_per_period=24):
    """
    Build a Shapely LineString for a single sine wave in (θ, z) space.

    The wave runs along z with:
        θ(z) = theta_center + amplitude_rad * sin(2π * z / period_mm)
    """
    n_pts = max(4, int(np.ceil((z_max - z_min) / period_mm * pts_per_period)) + 1)
    zs = np.linspace(z_min, z_max, n_pts)
    thetas = theta_center + amplitude_rad * np.sin(2 * np.pi * zs / period_mm)
    return shapely.LineString(np.column_stack([thetas, zs]))


def do_support_infill(support_geom, slicing_params, current_r, direction="body"):
    """
    Generate support infill in (θ,z) space.

    direction="radial"  — radial lines only (constant z) at support_interface_spacing.
    direction="axial"   — axial lines only (constant θ) at support_interface_spacing.
    direction="body_radial" — radial lines only (constant z) at support_body_spacing.
    direction="body_axial"  — axial lines only (constant θ) at support_body_spacing.
    direction="body"    — full body pattern (sine or grid) at support_body_spacing.

    Interface layers spread the grid over multiple layers by alternating direction
    based on depth from the nearest model layer (depth 1 → radial, depth 2 → axial, …).

    Returns a list of Shapely LineStrings.
    """
    if support_geom is None or support_geom.is_empty:
        return []

    result = []
    bounds = shapely.bounds(support_geom)
    z_min, z_max = bounds[1], bounds[3]

    # Anchor θ positions to r_base so lines stay at the same physical angles on every layer,
    # matching the strategy used by do_sparse_infill.
    r_base = slicing_params["r_base"]

    # Support strips are ~layer_height thick in z, so their infill intersections are
    # naturally short.  Only filter true zero-length degenerate geometry.
    min_length = 1e-6

    if direction in ("radial", "body_radial"):
        spacing = slicing_params["support_interface_spacing"] if direction == "radial" else slicing_params["support_body_spacing"]
        # Radial lines only (constant z, varying θ) — stitch seam-split arcs
        z_lo = np.floor(z_min / spacing) * spacing
        for z in np.arange(z_lo, z_max + spacing, spacing):
            ls = shapely.LineString([[0, z], [2 * np.pi, z]])
            inter = shapely.intersection(ls, support_geom)
            if hasattr(inter, 'geoms'):
                for seg in stitch_radial_seam(inter.geoms):
                    _collect_lines(_snap_radial_to_full_circle(seg), result, min_length)
            else:
                _collect_lines(_snap_radial_to_full_circle(inter), result, min_length)

    elif direction in ("axial", "body_axial"):
        spacing = slicing_params["support_interface_spacing"] if direction == "axial" else slicing_params["support_body_spacing"]
        _ts = spacing / r_base
        _n = int(np.floor(2 * np.pi / _ts))
        thetas = np.arange(_n) * _ts + _ts / 2   # offset by half-step to avoid θ=0 seam
        # Axial lines only (constant θ, varying z)
        for theta in thetas:
            ls = shapely.LineString([[theta, z_min - spacing], [theta, z_max + spacing]])
            _collect_lines(shapely.intersection(ls, support_geom), result, min_length)

    else:  # "body"
        spacing = slicing_params["support_body_spacing"]
        _ts = spacing / r_base
        n_body = int(np.floor(2 * np.pi / _ts))
        thetas = np.arange(n_body) * _ts + _ts / 2   # offset by half-step to avoid θ=0 seam

        mode = slicing_params.get("support_infill_mode", "grid")

        if mode == "sine":
            amplitude_rad = slicing_params["support_sine_amplitude_mm"] / current_r
            period_mm = slicing_params["support_sine_period_mm"]
            # Extend z slightly so waves fully cover the geometry bounds
            z_lo = z_min - period_mm
            z_hi = z_max + period_mm
            for theta in thetas:
                # Only add a seam-wrap copy when this specific wave can actually
                # cross the 0/2π boundary. Unconditionally emitting ±2π copies
                # triples support density on canonical gap layers because most
                # shifted waves still intersect the same support footprint.
                shifts = [0.0]
                if theta - amplitude_rad < 0.0:
                    shifts.append(2 * np.pi)
                if theta + amplitude_rad > 2 * np.pi:
                    shifts.append(-2 * np.pi)
                for shift in shifts:
                    ls = _make_sine_wave_ls(theta + shift, z_lo, z_hi, amplitude_rad, period_mm)
                    _collect_lines(shapely.intersection(ls, support_geom), result, min_length)
        else:
            # Axial lines (constant θ, varying z)
            for theta in thetas:
                ls = shapely.LineString([[theta, z_min - spacing], [theta, z_max + spacing]])
                _collect_lines(shapely.intersection(ls, support_geom), result, min_length)

            # Radial lines (constant z, varying θ) — stitch seam-split arcs
            z_lo = np.floor(z_min / spacing) * spacing
            for z in np.arange(z_lo, z_max + spacing, spacing):
                ls = shapely.LineString([[0, z], [2 * np.pi, z]])
                inter = shapely.intersection(ls, support_geom)
                if hasattr(inter, 'geoms'):
                    for seg in stitch_radial_seam(inter.geoms):
                        _collect_lines(_snap_radial_to_full_circle(seg), result, min_length)
                else:
                    _collect_lines(_snap_radial_to_full_circle(inter), result, min_length)

    return result

def do_solid_infill(peri, slicing_params, current_r, plot_infill=False, ax=None, direction="axial"):
    """Generate solid infill lines clipped to *peri* (in raw (θ,z) space).

    direction="axial"  — constant-θ lines (runs along the cylinder axis)
    direction="radial" — constant-z lines (arcs around the cylinder)
    """
    if peri.is_empty:
        return
    def do_transform(x):
        return shapely.transform(x,lambda p: p * [current_r,1])
    def do_transform_inverse(x):
        return shapely.transform(x,lambda p: p * [1/current_r,1])
    infill = []
    bounds = shapely.bounds(peri)
    theta_min, z_min, theta_max, z_max = bounds
    multi_turn = (theta_max - theta_min) > (2 * np.pi + 1e-6)

    def _collect(inter):
        geoms = None
        if type(inter) == shapely.MultiLineString:
            geoms = inter.geoms
        elif type(inter) == shapely.geometry.collection.GeometryCollection:
            geoms = inter.geoms
        if geoms is not None:
            for seg in geoms:
                if plot_infill:
                    plot_line(do_transform(seg), ax=ax, color="greenyellow", add_points=False)
                infill.append(seg)
        else:
            if type(inter) == shapely.geometry.multipoint.MultiPoint:
                print("got MultiPoint in solid infill")
            else:
                if plot_infill:
                    plot_line(do_transform(inter), ax=ax, color="greenyellow", add_points=False)
                infill.append(inter)

    if direction == "axial":
        # Constant-θ lines — runs along the cylinder axis (z direction)
        _solid_step = slicing_params["nozzle_size"] / current_r
        _n_solid = int(np.floor(2 * np.pi / _solid_step))
        # Round step so lines divide evenly into 2π (keeps grid aligned across layers).
        _solid_step = 2 * np.pi / _n_solid
        theta_origin = _solid_step / 2

        # Single-turn seam-copy geometry should still be folded onto one canonical
        # turn so we do not emit stacked duplicate lines at equivalent angles.
        if multi_turn:
            _peri_scan = peri
            _k_start = int(np.ceil((theta_min - theta_origin) / _solid_step))
            _k_end = int(np.floor((theta_max - theta_origin) / _solid_step))
            _all_thetas = theta_origin + np.arange(_k_start, _k_end + 1) * _solid_step
        else:
            _TWO_PI = 2 * np.pi
            _n_lo = int(np.floor(theta_min / _TWO_PI))
            _n_hi = int(np.ceil(theta_max / _TWO_PI))
            _pieces = []
            _clip_box = shapely.box(0, z_min, _TWO_PI, z_max)
            for _n in range(_n_lo, _n_hi + 1):
                _shifted = shapely.transform(peri, lambda p: p + [-_n * _TWO_PI, 0])
                _clipped = _shifted.intersection(_clip_box)
                if not _clipped.is_empty:
                    _pieces.append(_clipped)
            _peri_scan = shapely.unary_union(_pieces) if _pieces else peri
            _all_thetas = np.arange(_n_solid) * _solid_step + theta_origin

        for theta in _all_thetas:
            ls = shapely.LineString([[theta, z_min], [theta, z_max]])
            _collect(shapely.intersection(ls, _peri_scan))
    else:
        # Constant-z lines — arcs around the cylinder (radial / circumferential)
        z_step = slicing_params["nozzle_size"]   # mm; z axis is already in mm
        z_lo = np.floor(z_min / z_step) * z_step
        for z in np.arange(z_lo, z_max + z_step, z_step):
            ls = shapely.LineString([[theta_min, z], [theta_max, z]])
            _collect(shapely.intersection(ls, peri))

    return infill

def do_sparse_infill(innermost_peri, slicing_params, current_r, plot_infill=False, ax=None):
    """
    Generate sparse grid infill for the interior of a cylindrical shell layer.

    Works in (θ, z) space (theta in radians, z in mm). The grid has two families
    of lines:
      - "Radial" lines  → constant-z horizontal lines in (θ,z) space, i.e. the
                          printer rotates while extruding — these become arcs on
                          the cylinder surface.
      - "Axial" lines   → constant-θ vertical lines in (θ,z) space, i.e. the
                          printhead moves along the cylinder axis while extruding.

    Axial lines are anchored to a fixed θ-grid derived from r_base (the innermost
    shell). The first n_base lines land at the same θ positions on every layer.
    Extra lines beyond n_base are "spawned" fresh when the radius grows enough to
    fit them; they are arc-aligned but have no counterpart on inner layers.

    Returns three lists of Shapely LineString segments (all in raw (θ,z) space):
      fully_radial_infill — radial lines that span the full interior in one shot
                            (no multi-part intersection; single connected segment)
      radial_infill       — radial line segments clipped to the interior polygon
                            (multi-part intersections produce multiple entries)
      axial_infill        — axial line segments clipped to the interior polygon
    """
    if innermost_peri.is_empty:
        return [], [], []

    # Scale θ axis by current_r to convert arc-length → mm for plotting only.
    def do_transform(x):
        return shapely.transform(x, lambda p: p * [current_r, 1])
    def do_transform_inverse(x):
        return shapely.transform(x, lambda p: p * [1/current_r, 1])

    fully_radial_infill = []  # single-segment radial lines (full crossing)
    radial_infill = []        # multi-segment radial lines (partial crossings)
    axial_infill = []         # axial (constant-θ) line segments

    # ── Snap the z bounding-box to the infill grid spacing ──────────────────
    # Floor the bottom z-bound and ceil the top so grid lines always start/end
    # outside the polygon, guaranteeing full coverage without partial first/last gaps.
    lowest = shapely.bounds(innermost_peri)[1]          # min z of the polygon
    lowest /= slicing_params["infill_spacing"]
    lowest = np.floor(lowest) * slicing_params["infill_spacing"]   # snap down

    highest = shapely.bounds(innermost_peri)[3]         # max z of the polygon
    highest /= slicing_params["infill_spacing"]
    highest = np.ceil(highest) * slicing_params["infill_spacing"]  # snap up

    theta_min = shapely.bounds(innermost_peri)[0]
    theta_max = shapely.bounds(innermost_peri)[2]
    multi_turn = (theta_max - theta_min) > (2 * np.pi + 1e-6)

    # ── Pass 1: Radial infill (constant-z lines, θ = 0 → 2π) ───────────────
    # Sweep horizontal lines from lowest to highest z, spaced by infill_spacing.
    for current_line_v in np.arange(lowest, highest + 1, slicing_params["infill_spacing"]):
        # Build a full-width horizontal line at this z value across the whole
        # occupied theta span, not just one canonical turn.
        lp_x, lp_y = [theta_min, theta_max], [current_line_v, current_line_v]
        ls = shapely.LineString(np.array([lp_x, lp_y]).swapaxes(0, 1))

        # Clip the line to the interior polygon.
        inter = shapely.intersection(ls, innermost_peri)

        geoms = None
        if type(inter) == shapely.MultiLineString:
            # The line crossed the interior in multiple disjoint segments —
            # collect each segment individually.
            geoms = inter.geoms
        elif type(inter) == shapely.geometry.collection.GeometryCollection:
            # Mixed geometry result (can include points + lines); iterate all.
            geoms = inter.geoms

        if geoms is not None:
            stitched = list(geoms) if multi_turn else stitch_radial_seam(geoms)
            if len(stitched) == 1:
                # Seam-crossing arc stitched into one segment → fully radial.
                seg = stitched[0] if multi_turn else _snap_radial_to_full_circle(stitched[0])
                if plot_infill:
                    plot_line(do_transform(seg), ax=ax, color="red", add_points=False)
                fully_radial_infill.append(seg)
            else:
                for line_segment in stitched:
                    if not multi_turn:
                        line_segment = _snap_radial_to_full_circle(line_segment)
                    if plot_infill:
                        plot_line(do_transform(line_segment), ax=ax, color="orange", add_points=False)
                    radial_infill.append(line_segment)
        else:
            if type(inter) == shapely.geometry.multipoint.MultiPoint:
                # Degenerate: line only grazed polygon vertices — skip, just log.
                print("got MultiPoint")
                shapely.plotting.plot_points(do_transform(inter), color="red")
            else:
                # Single continuous segment spanning the full interior crossing.
                if not multi_turn:
                    inter = _snap_radial_to_full_circle(inter)
                if plot_infill:
                    plot_line(do_transform(inter), ax=ax, color="red", add_points=False)
                fully_radial_infill.append(inter)

    # ── Pass 2: Axial infill (constant-θ lines, z = lowest → highest) ───────
    # Strategy: anchor the entire grid at r_base (the innermost shell radius,
    # fixed for the whole print). This gives a constant θ-step that never
    # changes, so the first n_base lines land at exactly the same θ positions
    # on every layer. As the radius grows, floor(2π * current_r / infill_spacing)
    # eventually exceeds n_base — any extra lines are "spawned" fresh on that
    # layer at their natural positions (arc-aligned but not previously visited).
    # Fixed axial grid anchored at r_base — same lines on every layer, no spawning.
    # Spawning extra lines as radius grows creates unsupported strands (spaghetti):
    # each spawned midpoint is ~half the infill arc-gap from its nearest inner support.
    # 30° slope limit: supportable only if arctan(half_gap / layer_height) ≤ 30°,
    # i.e. half_gap ≤ layer_height·tan(30°) ≈ 0.115mm — far below typical arc-gaps.
    r_base = slicing_params["r_base"]
    fixed_step = slicing_params["infill_spacing"] / r_base     # never changes
    # Offset by half a step so no line falls exactly on the θ=0/2π polygon boundary.
    # After the morphological open in arc-length space the polygon's left edge retracts
    # slightly inward from θ=0, so a line at exactly 0 would miss the polygon.
    theta_origin = fixed_step / 2
    k_start = int(np.ceil((theta_min - theta_origin) / fixed_step))
    k_end = int(np.floor((theta_max - theta_origin) / fixed_step))
    all_thetas = theta_origin + np.arange(k_start, k_end + 1) * fixed_step

    for current_line_h in all_thetas:
        # Build a vertical line in (θ,z) space at this fixed θ value.
        lp_x, lp_y = [current_line_h, current_line_h], [lowest, highest]
        ls = shapely.LineString(np.array([lp_x, lp_y]).swapaxes(0, 1))

        # Clip the line to the interior polygon.
        inter = shapely.intersection(ls, innermost_peri)

        geoms = None
        if type(inter) == shapely.MultiLineString:
            geoms = inter.geoms
        elif type(inter) == shapely.geometry.collection.GeometryCollection:
            geoms = inter.geoms

        if geoms is not None:
            # Multiple disjoint segments — stitch tiny z-gaps before appending.
            stitched = stitch_axial_gap(geoms, max_gap=slicing_params["nozzle_size"])
            for line_segment in stitched:
                if plot_infill:
                    plot_line(do_transform(line_segment), ax=ax, color="darkviolet", add_points=False)
                axial_infill.append(line_segment)
        else:
            # Single segment — straightforward full crossing at this θ.
            if plot_infill:
                plot_line(do_transform(inter), ax=ax, color="darkviolet", add_points=False)
            axial_infill.append(inter)

    return fully_radial_infill, radial_infill, axial_infill
from shapely.geometry import Polygon, LineString, box, MultiLineString
from shapely.ops import linemerge

def _stitch_lines(lines, lower=0.0, period=2*math.pi, tol=1e-9):
    """
    Given a list of LineStrings produced by clipping a ring at [lower, lower+period],
    shift any 'left-boundary' pieces (those whose minimum x is within tol of lower)
    by +period so their endpoints align with the 'right-boundary' pieces, then
    linemerge the result.  Only shifts when multiple pieces exist (i.e. the ring
    actually crossed the seam) to avoid wrongly shifting un-split geometry.
    """
    lines = [l for l in lines if not l.is_empty]
    if not lines:
        return []
    if len(lines) > 1:
        wrapped = []
        for line in lines:
            coords = np.array(line.coords)
            if np.min(coords[:, 0]) < lower + tol:
                coords[:, 0] += period
            wrapped.append(LineString(coords))
    else:
        wrapped = lines

    wrapped = [l for l in wrapped if not l.is_empty]
    if not wrapped:
        return []
    merged = linemerge(MultiLineString(wrapped))
    if isinstance(merged, LineString):
        return [merged]
    return [l for l in merged.geoms if isinstance(l, LineString)]

def _split_theta_monotone_runs(line, tol=1e-9):
    """Split a LineString into maximal runs with constant theta-direction sign."""
    if line is None or line.is_empty:
        return []

    coords = np.asarray(line.coords, dtype=float)
    if len(coords) < 2:
        return []

    # Drop consecutive duplicate vertices from clipping / merge artifacts.
    step = np.abs(np.diff(coords, axis=0))
    keep = np.concatenate([[True], np.any(step > tol, axis=1)])
    coords = coords[keep]
    if len(coords) < 2:
        return []

    def _edge_sign(dtheta):
        if dtheta > tol:
            return 1
        if dtheta < -tol:
            return -1
        return 0

    runs = []

    def _append_run(run_coords, sign):
        arr = np.asarray(run_coords, dtype=float)
        if len(arr) < 2:
            return
        if sign < 0:
            arr = arr[::-1]
        ls = LineString(arr)
        if not ls.is_empty and ls.length > 0:
            runs.append(ls)

    run_start = 0
    run_sign = _edge_sign(coords[1, 0] - coords[0, 0])
    for i in range(1, len(coords) - 1):
        next_sign = _edge_sign(coords[i + 1, 0] - coords[i, 0])
        if next_sign != run_sign:
            _append_run(coords[run_start:i + 1], run_sign)
            run_start = i
            run_sign = next_sign
    _append_run(coords[run_start:], run_sign)
    return runs


def _restore_theta_walk(coords, theta_tol=1e-9, z_tol=1e-6):
    """Rebuild a continuous theta walk without collapsing real multi-turn edges.

    `np.unwrap` treats any >pi jump as a seam wrap, but on multi-turn helices a
    single boundary edge can legitimately span far more than one turn while also
    advancing in Z. Preserve those large pitch-carrying jumps; only shorten
    wrap-only jumps to their nearest seam-equivalent.
    """
    if len(coords) < 2:
        return coords

    out = np.asarray(coords, dtype=float).copy()
    diffs = np.diff(out[:, 0])
    dz = np.diff(out[:, 1])
    wrapped = (diffs + np.pi) % (2 * np.pi) - np.pi

    preserve_raw = (np.abs(diffs) > np.pi + theta_tol) & (np.abs(dz) > z_tol)
    step = np.where(preserve_raw, diffs, wrapped)
    out[:, 0] = np.concatenate([[out[0, 0]], out[0, 0] + np.cumsum(step)])
    return out

def _same_cyl_point(p0, p1, period=2 * math.pi, tol=1e-6):
    """Return True when two (theta, z) points are the same modulo one turn."""
    if abs(p0[1] - p1[1]) > tol:
        return False
    dtheta = ((p0[0] - p1[0]) + period / 2) % period - period / 2
    return abs(dtheta) <= tol

def _order_ring_segments(lines, period=2 * math.pi, tol=1e-6):
    """Orient seam-split ring segments so adjacent paths chain in ring order.

    `linemerge` / `unary_union` preserve geometry but not traversal order. For
    helical perimeters that means valid segments can come back reversed or
    permuted, producing start/end jumps in G-code. Reconstruct one walk around
    the ring by matching segment endpoints modulo one full turn, then unwrap
    each oriented segment so theta increases monotonically.
    """
    segs = []
    for line in lines:
        coords = np.asarray(line.coords, dtype=float)
        if len(coords) >= 2:
            if coords[-1, 0] < coords[0, 0] - tol:
                coords = coords[::-1]
            segs.append(coords)
    if len(segs) <= 1:
        return [LineString(segs[0])] if segs else []

    used = [False] * len(segs)

    start_idx = min(range(len(segs)), key=lambda i: np.min(segs[i][:, 0]))
    current = segs[start_idx].copy()
    ordered = [current]
    used[start_idx] = True
    current_end = current[-1].copy()

    while not all(used):
        found = False
        for i, coords in enumerate(segs):
            if used[i]:
                continue
            for reverse in (False, True):
                cand = coords[::-1] if reverse else coords
                if cand[-1, 0] < cand[0, 0] - tol:
                    cand = cand[::-1]
                if not _same_cyl_point(cand[0], current_end, period=period, tol=tol):
                    continue
                shift = round((current_end[0] - cand[0, 0]) / period) * period
                cand = cand.copy()
                cand[:, 0] += shift
                cand[0, 0] = current_end[0]
                ordered.append(cand)
                used[i] = True
                current_end = cand[-1].copy()
                found = True
                break
            if found:
                break
        if not found:
            # Fallback: start a new chain for any remaining unmatched fragment.
            # This keeps output deterministic instead of silently dropping data.
            next_idx = next(i for i, flag in enumerate(used) if not flag)
            current = segs[next_idx].copy()
            ordered.append(current)
            used[next_idx] = True
            current_end = current[-1].copy()

    return [LineString(coords) for coords in ordered if len(coords) >= 2]

def _ring_to_wrapped_lines(ring, lower=0.0, upper=2*math.pi, tol=1e-9):
    """Convert one polygon ring into printable lines in unwrapped theta space.

    For seam-crossing / multi-turn rings, clip the boundary line itself against the
    canonical strip instead of clipping the polygon.  That preserves the real edge
    geometry and avoids creating synthetic seam/endcap edges that later need to be
    guessed away with geometry-specific heuristics.
    """
    coords = np.asarray(ring.coords, dtype=float)
    if len(coords) < 2:
        return []

    coords = _restore_theta_walk(coords)
    theta_span = coords[:, 0].max() - coords[:, 0].min()
    if theta_span < 1e-3:
        return []

    period = upper - lower
    ring_ls = LineString(coords)

    # Simple case: ring already fits in one canonical turn, so preserve it as-is.
    if coords[:, 0].min() >= lower - tol and coords[:, 0].max() <= upper + tol and theta_span < period - tol:
        return [ring_ls] if not ring_ls.is_empty else []

    clip_box = box(lower, -1e15, upper, 1e15)
    k_min = int(math.ceil((lower - coords[:, 0].max()) / period))
    k_max = int(math.floor((upper - coords[:, 0].min()) / period))

    pieces = []
    for k in range(k_min, k_max + 1):
        shifted = shift_polygon(ring_ls, k * period) if k else ring_ls
        clipped = shifted.intersection(clip_box)
        if clipped.is_empty:
            continue
        if k:
            clipped = shift_polygon(clipped, -k * period)
        _collect_lines(clipped, pieces, min_length=0.0)

    if not pieces:
        return []

    merged = linemerge(shapely.unary_union(MultiLineString(pieces)))
    merged_lines = []
    _collect_lines(merged, merged_lines, min_length=0.0)

    raw_runs = []
    for line in merged_lines:
        line_coords = np.asarray(line.coords, dtype=float)
        line_span = line_coords[:, 0].max() - line_coords[:, 0].min()
        if line.is_ring and line_span < period - tol:
            raw_runs.append(line)
        else:
            raw_runs.extend(_split_theta_monotone_runs(line, tol=tol))
    return _order_ring_segments(raw_runs, period=period, tol=tol)

def polygon_to_wrapped_lines(poly, lower=0.0, upper=2*math.pi):
    """
    Convert a polygon (or collection) to LineStrings for gcode/plotting.

    external_peris stores each physical polygon THREE times (canonical + ±2π
    copies).  We keep only the canonical copy by checking that the centroid
    theta is within [lower, upper).  Rings that fit inside one canonical turn
    are returned unchanged; rings that cross the seam or extend past one turn
    are clipped as boundary lines, shifted back into unwrapped theta space, and
    split only where their theta direction actually changes.
    """
    if poly is None or (hasattr(poly, 'is_empty') and poly.is_empty):
        return []

    # Recurse into collections / lists
    if isinstance(poly, list):
        ret = []
        for each in poly:
            ret += polygon_to_wrapped_lines(each, lower, upper)
        return ret

    if isinstance(poly, (shapely.MultiPolygon, shapely.GeometryCollection)):
        ret = []
        for geom in poly.geoms:
            ret += polygon_to_wrapped_lines(geom, lower, upper)
        return ret

    if isinstance(poly, (LineString, shapely.MultiLineString)):
        # Already a line (e.g. infill passed through here) — flatten and return as-is
        ret = []
        _collect_lines(poly, ret, min_length=0.0)
        return ret

    if not isinstance(poly, Polygon):
        return []

    # Skip ±2π duplicate copies — only keep the canonical polygon
    if not (lower <= poly.centroid.x < upper):
        return []

    result = []
    for ring in [poly.exterior] + list(poly.interiors):
        result.extend(_ring_to_wrapped_lines(ring, lower=lower, upper=upper))
    return result


def keep_canonical_polygons(geom, lower=0.0, upper=2 * math.pi):
    """Keep only polygon components whose centroid lies in the canonical turn.

    `external_peris` intentionally carries seam-neighbor copies at ±2π so support
    and boundary extraction can see across the wrap. Print-phase model geometry
    should use only the canonical physical copy, otherwise buffers/differences
    triple the interior area and later canonical clipping chops real multi-turn
    features back down to one turn.
    """
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()

    if isinstance(geom, Polygon):
        return geom if lower <= geom.centroid.x < upper else shapely.GeometryCollection()

    if isinstance(geom, shapely.MultiPolygon):
        polys = [g for g in geom.geoms if isinstance(g, Polygon) and lower <= g.centroid.x < upper]
        if not polys:
            return shapely.GeometryCollection()
        if len(polys) == 1:
            return polys[0]
        return shapely.MultiPolygon(polys)

    if isinstance(geom, shapely.GeometryCollection):
        kept = [keep_canonical_polygons(g, lower=lower, upper=upper) for g in geom.geoms]
        kept = [g for g in kept if g is not None and not g.is_empty]
        if not kept:
            return shapely.GeometryCollection()
        return shapely.unary_union(kept)

    return geom
    
def sum_poly(p):
    if isinstance(p, Polygon):
        return p
    if len(p) == 0:
        return MultiLineString()
    f = p[0]
    for e in p[1:]:
        f = f.union(e)
    return f


def densify_cylindrical(coords, max_angle_step=np.pi / 36):
    """
    Insert intermediate points along each segment so that arcs follow the
    cylinder surface instead of cutting through it as a chord.
    max_angle_step: maximum angular gap between consecutive plotted points (radians).
    """
    if len(coords) == 0:
        return coords
    result = [coords[0]]
    for i in range(len(coords) - 1):
        t0, z0 = coords[i]
        t1, z1 = coords[i + 1]
        n = max(2, int(np.ceil(abs(t1 - t0) / max_angle_step)) + 1)
        thetas = np.linspace(t0, t1, n)
        zs = np.linspace(z0, z1, n)
        result.extend(zip(thetas[1:], zs[1:]))
    return np.array(result)


def remove_collinear_points_angle(coords, angle_tolerance=1e-6):
    """
    Remove collinear points using angle change method (vectorized).
    """
    coords = np.asarray(coords)
    if len(coords) <= 2:
        return coords.tolist()

    v = np.diff(coords, axis=0)                          # (n-1, 2)
    norms = np.linalg.norm(v, axis=1)                    # (n-1,)
    safe = norms > 1e-10
    v_norm = v / np.where(norms[:, None] > 1e-10, norms[:, None], 1.0)

    dots = np.sum(v_norm[:-1] * v_norm[1:], axis=1)     # (n-2,)
    dots = np.clip(dots, -1.0, 1.0)
    angle_changes = np.arccos(np.abs(dots))

    keep_middle = (angle_changes > angle_tolerance) & safe[:-1] & safe[1:]
    keep = np.concatenate([[True], keep_middle, [True]])
    return coords[keep].tolist()


slicing_params = {
    "layer_height": 0.2,
    "infill_spacing" : 3,
    "nozzle_size": 0.4,
    "wall_count": 3,
    "feed_mm_per_min" : 3600,
    "travel_mm_per_min" : 4800,
    "support_infill_mode": "sine",   # "grid" or "sine"
    "support_sine_amplitude_mm": 2.0, # mm — sine wave half-amplitude (grid mode ignores this)
    "support_sine_period_mm": 8.0,   # mm — sine wave period along z (grid mode ignores this)
    "support_interface_spacing": 2,  # mm — line spacing for interface layers
    "support_interface_layers": 2,   # number of layers at each end of a support column that use alternating interface pattern
    "support_gap_layers": 1,         # number of empty (air-gap) layers between the top of support and the part bottom
    "support_body_spacing": 2.5,     # mm — line spacing for support interior (arc-length mm)
    "support_speed_factor": 1.5,     # multiply feed/travel speed for support sections (faster = less sag on bridges)
    "support_flow_rate": 0.8,        # flow multiplier for support (reduce to shrink blobs on contact)
    "infill_speed_factor": 1.3,      # multiply feed/travel speed for sparse infill (reduces bridging blobs)
    "infill_flow_rate": 0.9,         # flow multiplier for sparse infill
    "support_min_area": 1.0,         # mm² — inflate overhang patches below this area to at least this size
    "support_merge_gap": 15.0,       # mm — merge support branches closer than this into a shared trunk
    "support_taper_per_layer_mm": -0.1,  # mm — how much the support footprint shrinks per layer going inward
    "support_xy_distance_mm": 0.3,   # mm — lateral gap between support and model walls
    "retract_mm": 2.0,               # mm — retraction distance on travel moves
    "retract_speed_mm_per_min": 1800, # mm/min — retraction speed (30 mm/s)
    "min_retract_travel_mm": 2.0,    # mm — skip retraction for travel moves shorter than this
    "flow_rate": 1.0,                # multiplier on extrusion_per_mm (tune for over/underextrusion)
    "first_layer_speed_factor": 0.5, # fraction of feed/travel speeds for the first model layer
    "fan_start_layer": 2,            # model layer index at which part cooling fan turns on
    "top_layers": 3,                 # number of solid shells on the outer (top) surface
    "bottom_layers": 3,              # number of solid shells on the inner (bottom) surface
}

printer_params = {
    "bed_radius": 30.0,    # mm — mandrel/cylinder surface radius; no paths printed below this
    "max_radius": 120.0,    # mm — maximum printable radius; slicing stops here
    "axial_length": 165.0, # mm — printable length along the cylinder axis
    "filament_diameter": 1.75, # mm — filament diameter (1.75 or 2.85)
    "hotend_temp": 215,    # °C — PLA
    "min_x_pos": 0.0,    # mm — X coordinate where the printable bed begins (prime line here)
    "min_x_margin": 5.0,  # mm — gap between bed start and model start
}
# Derived: how many mm of filament per mm of real travel
slicing_params["extrusion_per_mm"] = (
    slicing_params["nozzle_size"] * slicing_params["layer_height"]
    / (np.pi * (printer_params["filament_diameter"] / 2) ** 2)
)
max_real_layers = None  # Set to an int to cap layers for testing, None for full print
mesh = trimesh.load_mesh(sys.argv[1] if len(sys.argv) > 1 else "helix.stl")
print(f"Loaded mesh: {len(mesh.faces)} faces, watertight={mesh.is_watertight}")
# Remove duplicate faces before anything else — Blender STL exports each face twice,
# which makes every edge appear 4× instead of 2×, causing the directed cross-section
# algorithm to mark every crossing point as ambiguous and produce zero geometry.
_orig_face_count = len(mesh.faces)
mesh.update_faces(mesh.unique_faces())
mesh.remove_unreferenced_vertices()
if len(mesh.faces) < _orig_face_count:
    print(f"  removed {_orig_face_count - len(mesh.faces)} duplicate faces ({_orig_face_count} → {len(mesh.faces)})")
if not mesh.is_watertight:
    # The directed cross-section algorithm handles open boundaries without repair.
    # Just ensure normals are consistently oriented so FROM/TO direction is reliable.
    trimesh.repair.fix_normals(mesh)
    print(f"  normals fixed (non-watertight mesh — open chains will be silently dropped)")
pts = cartesian_to_cylindrical(mesh.vertices)
# Shift mesh so its axial start aligns with min_x_pos + min_x_margin
_z_shift = printer_params["min_x_pos"] + printer_params["min_x_margin"] - np.min(pts[:, 2])
mesh.apply_translation([0, 0, _z_shift])

# Subdivide large faces so that cylindrical cross-sections have enough sample points
# along long edges (avoids jitter from straight-line approximation of curved arcs).
# IMPORTANT: trimesh.remesh.subdivide_to_size is NOT manifold-preserving — it splits
# each face's edges independently, creating T-junctions where one face is split but its
# neighbour still has the original unsplit edge.  T-junctions cause the directed
# cross-section algorithm to produce two outgoing crossing points (both tails) instead
# of a tail+head pair, so stitching can't close them, and entire polygon loops are
# silently discarded.  Loop subdivision (trimesh.remesh.subdivide) splits all edges of
# each face and maintains manifold topology, at the cost of a few extra passes.
_max_edge = slicing_params["layer_height"] * 20  # subdivide any edge longer than 5 layers
_orig_faces = len(mesh.faces)
_sub_verts, _sub_faces = mesh.vertices, mesh.faces
for _ in range(8):
    _cur_max = trimesh.Trimesh(_sub_verts, _sub_faces, process=False).edges_unique_length.max()
    if _cur_max <= _max_edge:
        break
    _sub_verts, _sub_faces = trimesh.remesh.subdivide(_sub_verts, _sub_faces)
mesh = trimesh.Trimesh(_sub_verts, _sub_faces, process=False)
print(f"Subdivided mesh: {_orig_faces} → {len(mesh.faces)} faces (max edge {_max_edge:.2f} mm)")
pts = cartesian_to_cylindrical(mesh.vertices)
max_z = np.max(pts[:,2])
min_z = np.min(pts[:,2])
part_height = (max_z-min_z)

# Z bounds of the printable volume in the same coordinate space as the mesh
z_clip_min = min_z
z_clip_max = min_z + printer_params["axial_length"]

# Only clip the model in z here. A fixed theta window truncates valid multi-turn
# geometry once the unwrapped perimeter exceeds that hardcoded span.
axial_clip_box = box(-1e15, z_clip_min, 1e15, z_clip_max)
# Canonical one-rotation clip box for support geometry
theta_clip_box = box(0, z_clip_min, 2 * np.pi, z_clip_max)

# for current_r in np.arange(max_r-layer_height*35,max_r,layer_height):
layer_data = {}
sliced_layer_count = 0
# Store the innermost shell radius once so do_sparse_infill can anchor its
# axial-line grid at a fixed θ-step that never drifts across layers.
r_start = max(np.min(pts[:,0]) + slicing_params["layer_height"], printer_params["bed_radius"])
r_end   = min(np.max(pts[:,0]) - 1e-3, printer_params["max_radius"])
print(f"Mesh radial range: {np.min(pts[:,0]):.2f} – {np.max(pts[:,0]):.2f} mm")
print(f"Slicing range:     {r_start:.2f} – {r_end:.2f} mm  ({max(0, (r_end-r_start)/slicing_params['layer_height']):.0f} layers)")
if r_start >= r_end:
    raise ValueError(f"bed_radius ({printer_params['bed_radius']} mm) is at or beyond the mesh max radius ({np.max(pts[:,0]):.2f} mm) — nothing to slice.")
slicing_params["r_base"] = r_start
_total_layers = max(0, int((r_end - r_start) / slicing_params["layer_height"]))
print(f"\nPass 1 — slicing geometry ({_total_layers} layers, r={r_start:.2f}–{r_end:.2f} mm)")
for layer_i,current_r in enumerate( np.arange(r_start, r_end, slicing_params["layer_height"])):
    if max_real_layers is not None and sliced_layer_count >= max_real_layers:
        print(f"\nSliced {sliced_layer_count} layers, stopping.")
        break

    times = [(time.perf_counter(), "")]
    print(f"  [{sliced_layer_count+1}/{_total_layers}] r={current_r:.3f} mm", end="\r", flush=True)
    # Offset by a tiny fraction of layer_height so we never slice exactly through
    # mesh boundary vertices (e.g. the outer cone bottom edge at r=42.5mm).
    # The 2e-5 mm shift is negligible for geometry but prevents tangled zero-length
    # segments that arise when mesh vertices straddle the exact slice radius.
    _slice_r = current_r + slicing_params["layer_height"] * 1e-4
    cross_pts_cart, outer_edges = cylindrical_cross_section(mesh, _slice_r)
    if len(cross_pts_cart) < 3:
        print(f"empty layer for {layer_i}")
        continue

    layer_data[layer_i] = {"radius": current_r}
    sliced_layer_count += 1

    times.append((time.perf_counter(), "intersection"))
    layer_pts = cartesian_to_cylindrical(cross_pts_cart)
    triangles = []
    times.append((time.perf_counter(), "mask generation"))

    # try:
    poly_idxs = find_polygons(outer_edges, len(layer_pts))
    # except Exception as e:
    #     print(e)
    #     cut_mesh.show()
    #     print(triangles.shape)
    #     for tri in triangles:
    #         if 1914 not in tri:
    #             continue
    #         print(tri)
    #         tri2 = tri.tolist()+[tri[0]]
    #         plt.plot(layer_pts[tri2,1],layer_pts[tri2,2],"o-b", markersize=3,alpha=0.2)
    #     plt.plot(layer_pts[1914,1],layer_pts[1914,2],"o", markersize=15)
    #     plt.show()
    #     raise e
    poly_pts_raw = [unwrap_polygon(layer_pts[each,1:]).reshape(-1, 2) for each in poly_idxs]
    poly_pts = [p for p in poly_pts_raw if len(p) >= 3]
    # for poly in poly_pts:
    #     plt.plot(poly[:,0],poly[:,1],"-")
    # plt.show()
    full_circles = [is_full_circle(poly) for poly in poly_pts]
    # triangles is always [] with the cylindrical_cross_section approach, so
    # find_connected_shapes can't pair full-circle loops via shared faces.
    # Instead, sort by average z and pair consecutive loops (bottom/top of each
    # solid region).
    fc_original_indices = [i for i, fc in enumerate(full_circles) if fc]
    if len(fc_original_indices) >= 2:
        fc_z_mean = [np.mean(np.array(poly_pts[i])[:, 1]) for i in fc_original_indices]
        sorted_order = np.argsort(fc_z_mean)
        sorted_fc = [fc_original_indices[j] for j in sorted_order]
        full_circle_pairs = [np.array([sorted_fc[k], sorted_fc[k+1]]) for k in range(0, len(sorted_fc) - 1, 2)]
    else:
        full_circle_pairs = []
    full_circle_pairs = np.array(full_circle_pairs) if full_circle_pairs else np.empty((0, 2), dtype=int)
    polys = [Polygon(p) for i,p in enumerate(poly_pts) if not full_circles[i]]
    mapping = []
    for pair in full_circle_pairs:
        p1 = poly_pts[pair[0]]
        p2 = poly_pts[pair[1]]
        # split_idx1 = np.argmax(np.abs(np.diff(p1[:,0])))
        # split_idx2 = np.argmax(np.abs(np.diff(p2[:,0])))
        # p1 = np.roll(p1,-split_idx1)
        # p2 = np.roll(p2,-split_idx2)
        p1 -= [np.round(np.min(p1[:,0])/(2*np.pi))*2*np.pi,0]
        p2 -= [np.round(np.min(p2[:,0])/(2*np.pi))*2*np.pi,0]
        if(p1[0][0]-p1[-1][0] >0):
            p1 = p1[::-1]
        if(p2[0][0]-p2[-1][0] >0):
            p2 = p2[::-1]
        # p1 = p1[np.argsort(p1[:,0])]
        # p2 = p2[np.argsort(p2[:,0])]
        final_poly = np.concatenate([p1-[2*np.pi,0],p1,p1+[2*np.pi,0],(p2+[2*np.pi,0])[::-1,:],p2[::-1,:],(p2-[2*np.pi,0])[::-1,:]])
        # plt.plot(final_poly[:,0],final_poly[:,1],"o-")
        # plt.show()
        final_poly = shapely.make_valid(Polygon(final_poly))
        final_poly = trim_polygon_left_right(final_poly,0.1)
        # make_valid + trim can return a MultiPolygon when the concatenated
        # 6π-wide ring has self-intersections.  Flatten into individual
        # Polygons so the containment hierarchy loop (which only handles
        # Polygon instances) processes each component.
        fp_parts = list(final_poly.geoms) if hasattr(final_poly, "geoms") else [final_poly]
        for part in fp_parts:
            if not isinstance(part, Polygon) or part.is_empty:
                continue
            mapping += [len(polys)]
            polys.append(part)
    # for poly in poly_pts:
    #     if is_full_circle(poly):
    #         plt.plot(poly[:,0],poly[:,1])
    # plt.show()
    hier = build_hierarchy(polys)
    times.append((time.perf_counter(),"build_hierarchy"))

    flat_hier = [each2 for each in hier for each2 in each]
    levels = [0]*len(polys)
    _levels_i,_levels_count = np.unique(flat_hier,return_counts=True)
    for i,count in zip(_levels_i,_levels_count):
        levels[i] = count
    external_perimeters = []
    times.append((time.perf_counter(),"count levels"))


    plot_perimeters = False
    def do_transform(x):
        return shapely.transform(x, lambda p: p * [current_r, 1])
    def do_transform_inverse(x):
        return shapely.transform(x, lambda p: p * [1/current_r, 1])

    def _shift_periods_to_canonical_centroid(geom, lower=0.0, upper=2 * np.pi):
        """Return an integer turn shift that places the centroid in [lower, upper)."""
        if geom is None or geom.is_empty:
            return 0
        period = upper - lower
        centroid_x = geom.centroid.x
        k = math.floor((lower - centroid_x) / period)
        while centroid_x + k * period < lower:
            k += 1
        while centroid_x + k * period >= upper:
            k -= 1
        return k

    for poly_i,sub_polys in enumerate(hier):
        if levels[poly_i] % 2 == 1:
            continue
        if not isinstance(polys[poly_i], Polygon) or polys[poly_i].is_empty:
            continue
        sub_poly_list = [np.array(polys[sub_poly].exterior.coords).astype(float) for sub_poly in sub_polys if levels[sub_poly] == levels[poly_i]+1 and isinstance(polys[sub_poly], Polygon) and not polys[sub_poly].is_empty]
        sub_poly_list = [[each - [2*np.pi,0],each,each + [2*np.pi,0]] for each in sub_poly_list]
        sub_poly_list = [remove_collinear_points_angle(each2*np.array([current_r,1]),0.01)*np.array([1/current_r,1]) for each in sub_poly_list for each2 in each]
        sub_poly_list = [h for h in sub_poly_list if len(h) >= 4]
        poly_exterior = remove_collinear_points_angle(polys[poly_i].exterior.coords*np.array([current_r,1]),0.01)*np.array([1/current_r,1])
        if len(poly_exterior) < 4:
            continue
        if not len(sub_poly_list):
            poly_with_holes = Polygon(poly_exterior)
        else:
            poly_with_holes = Polygon(poly_exterior,holes = sub_poly_list)

        if plot_perimeters:
            plot_polygon(do_transform(poly_with_holes),add_points=False,facecolor=(0, 0, 0, 0))

        poly_with_holes = do_transform_inverse(do_transform(poly_with_holes).buffer(-slicing_params["nozzle_size"]/2,join_style="mitre"))
        poly_with_holes = poly_with_holes.intersection(axial_clip_box)
        if poly_with_holes.is_empty:
            continue

        # Normalize each physical perimeter into a canonical turn before adding
        # seam-neighbor copies. This keeps one stored copy selectable by the
        # centroid-based canonical filter even when the raw geometry spans many
        # full turns.
        _k_canon = _shift_periods_to_canonical_centroid(poly_with_holes)
        if _k_canon:
            poly_with_holes = shift_polygon(poly_with_holes, _k_canon * 2 * np.pi)

        external_perimeters.append(poly_with_holes)
        if poly_i not in mapping:
            external_perimeters.append(shift_polygon(poly_with_holes,2*np.pi))
            external_perimeters.append(shift_polygon(poly_with_holes,-2*np.pi))
        if plot_perimeters:
            plot_polygon(do_transform(external_perimeters[-1]),add_points=False,color="green",facecolor=(0, 0, 0, 0))
    if plot_perimeters:
        plt.show()
    times.append((time.perf_counter(),"form poly with holes"))
    
    final = shapely.unary_union(external_perimeters)
    layer_data[layer_i]["external_peris"] = final
    times.append((time.perf_counter(),"full union"))

    if layer_i == 0:
        layer_data[layer_i]["bottom_surfaces"] = layer_data[layer_i]["external_peris"]
        continue

    e = do_transform(final)

    e -= do_transform(layer_data[layer_i-1]["external_peris"]).buffer(slicing_params["nozzle_size"]/2)
    e -= do_transform(shift_polygon(layer_data[layer_i-1]["external_peris"],-np.pi*2)).buffer(slicing_params["nozzle_size"]/2)

    if not e.is_empty:
        ns = do_transform_inverse(e)
        ns = clean_support_region(
            ns, current_r,
            min_area_mm2=slicing_params["support_min_area"],
            merge_gap_mm=slicing_params["support_merge_gap"],
        )
        # Fold shifted fragments into the canonical [0, 2π] range and remove
        # areas already occupied by the model at this layer.
        ns = fold_to_canonical(ns, theta_clip_box)
        if ns.area > 1e-6:
            layer_data[layer_i]["needs_support"] = ns
    times.append((time.perf_counter(),"support check"))

    e = do_transform(layer_data[layer_i-1]["external_peris"])

    e = e.buffer(-slicing_params["nozzle_size"],join_style="mitre")# shrink the area of the top layer so we don't flag slicing_params["wall_count"] as needing top layers, TODO tweak this metric slicing_params["nozzle_size"]*slicing_params["wall_count"] is probably too little overlap between slicing_params["wall_count"]
    e = do_transform_inverse(e)

    e -= final
    e -= shift_polygon(final,-np.pi*2)
    times.append((time.perf_counter(),"top surface check"))

    layer_data[layer_i-1]["top_surfaces"] = e

    # Bottom surfaces of current layer: areas not covered by the layer below
    e = do_transform(final)
    e = e.buffer(-slicing_params["nozzle_size"], join_style="mitre")
    e = do_transform_inverse(e)
    e -= layer_data[layer_i-1]["external_peris"]
    e -= shift_polygon(layer_data[layer_i-1]["external_peris"], np.pi*2)
    layer_data[layer_i]["bottom_surfaces"] = e
    times.append((time.perf_counter(), "bottom surface check"))
    # for i in range(1,len(times)):
    #     print(f"{times[i][1]} took {times[i][0]-times[i-1][0]:.2f} seconds")




layer_data[list(layer_data)[-1]]["top_surfaces"] = layer_data[list(layer_data)[-1]]["external_peris"]

# ── Propagate top/bottom solid surfaces inward/outward ──────────────────
# top_layers: the outermost shell already has top_surfaces set; extend it
# to the next (top_layers-1) inner shells so they also print solid there.
# bottom_layers: same idea for inner-surface solid shells, propagating outward.
# Snapshot original source layers BEFORE modifying any, so propagated geometry
# doesn't cascade further on subsequent loop iterations.
_sorted_keys = sorted(layer_data.keys())
_empty_geom = shapely.GeometryCollection()
_top_sources  = [(pos, key, layer_data[key]["top_surfaces"])
                 for pos, key in enumerate(_sorted_keys)
                 if not layer_data[key].get("top_surfaces", _empty_geom).is_empty]
_bot_sources  = [(pos, key, layer_data[key]["bottom_surfaces"])
                 for pos, key in enumerate(_sorted_keys)
                 if not layer_data[key].get("bottom_surfaces", _empty_geom).is_empty]

for _pos, _key, _ts in _top_sources:
    for _depth in range(1, slicing_params["top_layers"]):
        _inner_pos = _pos - _depth
        if _inner_pos < 0:
            break
        _ik = _sorted_keys[_inner_pos]
        _prop = _ts.intersection(layer_data[_ik]["external_peris"])
        if _prop.is_empty:
            break
        _existing = layer_data[_ik].get("top_surfaces", _empty_geom)
        layer_data[_ik]["top_surfaces"] = _existing.union(_prop)

for _pos, _key, _bs in _bot_sources:
    for _depth in range(1, slicing_params["bottom_layers"]):
        _outer_pos = _pos + _depth
        if _outer_pos >= len(_sorted_keys):
            break
        _ok = _sorted_keys[_outer_pos]
        _prop = _bs.intersection(layer_data[_ok]["external_peris"])
        if _prop.is_empty:
            break
        _existing = layer_data[_ok].get("bottom_surfaces", _empty_geom)
        layer_data[_ok]["bottom_surfaces"] = _existing.union(_prop)

print(f"\nPass 1 complete — {sliced_layer_count} layers sliced.")

# ═══════════════════════════════════════════════════════════════════════
# Support propagation: backward pass (outer → inner)
# ═══════════════════════════════════════════════════════════════════════
sorted_layer_keys = sorted(layer_data.keys())
print(f"Pass 1.5 — support propagation (backward pass over {len(sorted_layer_keys)} layers)...")

for i in sorted_layer_keys:
    layer_data[i]["support"] = shapely.GeometryCollection()

support_needed = shapely.GeometryCollection()   # (θ,z) regions needing a support column below

taper_mm   = slicing_params.get("support_taper_per_layer_mm", 0.5)
xy_dist_mm = slicing_params.get("support_xy_distance_mm", 0.3)
merge_gap  = slicing_params.get("support_merge_gap", 0)

for i in reversed(sorted_layer_keys):
    current_ext = layer_data[i]["external_peris"]
    current_r   = layer_data[i]["radius"]

    # Ground: subtract model, expanded by xy gap so support doesn't press against walls.
    # xy_dist is a lateral (arc-length) gap only — clamp to the model's original Z bounds
    # so the buffer doesn't also push support away from the cone's top/bottom surfaces in
    # Z, which would create a spurious vertical gap at the narrow end of tapered geometry.
    if xy_dist_mm > 0 and not current_ext.is_empty:
        metric_ext     = shapely.affinity.scale(current_ext, xfact=current_r, yfact=1.0, origin=(0, 0))
        metric_exp     = metric_ext.buffer(xy_dist_mm)
        _zb            = metric_ext.bounds          # (xmin, z_min, xmax, z_max)
        metric_exp     = metric_exp.intersection(shapely.geometry.box(-1e9, _zb[1], 1e9, _zb[3]))
        model_expanded = shapely.affinity.scale(metric_exp, xfact=1.0 / current_r, yfact=1.0, origin=(0, 0))
        support_needed = support_needed.difference(model_expanded)
    else:
        support_needed = support_needed.difference(current_ext)

    # Record what this layer must print as support.
    # Store with ±2π shifted copies (like external_peris) so seam-crossing
    # regions remain as one continuous polygon rather than two split pieces.
    if not support_needed.is_empty:
        clipped = support_needed.intersection(theta_clip_box)
        layer_data[i]["support"] = shapely.unary_union([
            clipped,
            shift_polygon(clipped, 2*np.pi),
            shift_polygon(clipped, -2*np.pi),
        ])

    # Add this layer's overhangs so inner layers know to support them
    if "needs_support" in layer_data[i]:
        support_needed = support_needed.union(layer_data[i]["needs_support"])
        support_needed = fold_to_canonical(support_needed, theta_clip_box)

    # Organic taper: erode the support footprint so columns narrow toward the bed
    if not support_needed.is_empty :
        metric = shapely.affinity.scale(support_needed, xfact=current_r, yfact=1.0, origin=(0, 0))
        metric = metric.buffer(-taper_mm)
        if metric.is_empty:
            support_needed = shapely.GeometryCollection()
        else:
            support_needed = shapely.affinity.scale(metric, xfact=1.0 / current_r, yfact=1.0, origin=(0, 0))
            support_needed = fold_to_canonical(support_needed, theta_clip_box)

    # Merge nearby branches that have tapered into proximity
    if merge_gap > 0 and not support_needed.is_empty:
        support_needed = clean_support_region(
            support_needed, current_r,
            min_area_mm2=0,
            merge_gap_mm=merge_gap,
        )

print("Support propagation complete.")

# ── Gap layers: model floating above the bed ────────────────────────────
# If the innermost model layer is above bed_radius, print support shells
# from bed_radius up to the first model layer.
gap_support_layers = []
first_r = layer_data[sorted_layer_keys[0]]["radius"]
if first_r > printer_params["bed_radius"]:
    print(f"Gap support: model starts at r={first_r:.2f} mm (bed r={printer_params['bed_radius']:.2f} mm) — generating gap layers...")
    # Union the accumulated ungrounded support with the full first-layer footprint
    # (the entire first layer floats, so all of it needs a support base)
    gap_footprint = support_needed.union(
        layer_data[sorted_layer_keys[0]]["external_peris"]
    ).intersection(theta_clip_box)
    if not gap_footprint.is_empty:
        gap_radii = list(np.arange(printer_params["bed_radius"],
                                   first_r,
                                   slicing_params["layer_height"]))
        n_gap = len(gap_radii)
        for k, r in enumerate(gap_radii):
            # Organic flare: layers nearer the bed are wider (distance from top increases)
            layers_from_top = n_gap - 1 - k
            if taper_mm > 0 and layers_from_top > 0:
                m = shapely.affinity.scale(gap_footprint, xfact=r, yfact=1.0, origin=(0, 0))
                m = m.buffer(taper_mm * layers_from_top)
                tapered = shapely.affinity.scale(m, xfact=1.0 / r, yfact=1.0, origin=(0, 0))
                tapered = tapered.intersection(theta_clip_box)
            else:
                tapered = gap_footprint
            gap_support_layers.append({"radius": r, "support": tapered})

with open("stats.txt", "w") as _sf:
    _sf.write(
        "# Column definitions (areas are in (θ,z) space — radians×mm; multiply by radius for metric mm²)\n"
        "#   layer     — layer index; 0 = innermost model layer, counts outward\n"
        "#   radius    — nozzle distance from rotation axis in mm (the printer's Z axis)\n"
        "#   ext_area  — cross-section area of the model at this radius; grows as model widens outward\n"
        "#   ns_area   — needs_support area: portion of this layer not covered by the layer below (overhang);\n"
        "#               non-zero means the model grew here in a way the inner layer can't support\n"
        "#   sup_area  — support geometry area assigned to print at this radius: union of all outer-layer\n"
        "#               overhang footprints not yet grounded by the model; what actually gets deposited\n"
        "#   sup_paths — number of support infill line segments generated from sup_area;\n"
        "#               should be non-zero wherever sup_area > 0 (if zero, support is missing)\n"
        "#   solid_ln  — solid infill line count (dense radial sweeps for top/bottom surfaces)\n"
        "#   sparse_ln — sparse infill line count (radial + axial + fully-radial interior fill)\n"
        "#\n"
        "# Key diagnostic: ns_area > 0 at layer N → expect sup_area > 0 at all layers 0..N-1\n"
        "#                 sup_area > 0 but sup_paths = 0 → support geometry exists but no infill generated\n"
        "\n"
    )
    _sf.write(f"mesh: {mesh.metadata.get('file_path', 'unknown')}  vertices={len(mesh.vertices)}\n")
    _sf.write(f"bed_radius={printer_params['bed_radius']} mm  max_radius={printer_params['max_radius']} mm\n")
    _sf.write(f"layer_height={slicing_params['layer_height']} mm  nozzle={slicing_params['nozzle_size']} mm  walls={slicing_params['wall_count']}\n")
    _sf.write(f"first_r={first_r:.3f} mm  r_start={slicing_params['r_base']:.3f} mm\n")
    _sf.write(f"model_layers={len(sorted_layer_keys)}  gap_layers={len(gap_support_layers)}\n")
    if gap_support_layers:
        _sf.write(f"gap_r: {gap_support_layers[0]['radius']:.3f} → {gap_support_layers[-1]['radius']:.3f} mm"
                  f"  footprint_area(θz)={gap_support_layers[0]['support'].area:.6f}\n")
    _sf.write("\n")
    _sf.write(f"{'layer':>6}  {'radius':>8}  {'ext_area':>10}  {'ns_area':>10}  {'sup_area':>10}"
              f"  {'sup_paths':>9}  {'solid_ln':>8}  {'sparse_ln':>9}\n")
    _sf.write("-" * 80 + "\n")

plot_3d = False       # live combined 3D figure (non-blocking)
plot_per_layer = False  # show blocking per-layer 2D+3D window; set False for full gcode runs
ax_3d = None

def on_key(event):
    """Ctrl-C or Q in any plot window exits the program cleanly."""
    if event.key in ('ctrl+c', 'q', 'escape'):
        plt.close('all')
        sys.exit(0)

if plot_3d:
    fig_combined = plt.figure(figsize=(8, 8))
    ax_3d = fig_combined.add_subplot(1, 1, 1, projection='3d')
    fig_combined.canvas.mpl_connect('key_press_event', on_key)
    plt.ion()

# output_scad = solid2.cube()
print(f"\nPass 2 — generating G-code ({len(gap_support_layers)} gap support + {len(sorted_layer_keys)} model layers)...")
gcodes = []
last_pos = [0,0]
last_e = 0
real_layer_count = 0
if max_real_layers is not None:
    print(f"\n=== CAPPED AT {max_real_layers} LAYERS FOR TESTING ===\n")
_bed_r = printer_params["bed_radius"]
_prime_e = 2 * np.pi * _bed_r * slicing_params["extrusion_per_mm"] * slicing_params["flow_rate"]
_first_feed = slicing_params["feed_mm_per_min"] * slicing_params["first_layer_speed_factor"]
_min_x = printer_params["min_x_pos"]
_start_gcode = f""";start of slice
;LAYER_HEIGHT:{slicing_params['layer_height']}
;LINE_WIDTH:{slicing_params['nozzle_size']}
M82 ; absolute extrusion mode
G28 ; home all axes
M104 S{printer_params['hotend_temp']} ; start hotend heating
M109 S{printer_params['hotend_temp']} ; wait for hotend
G92 E0 ; reset extruder
M107 ; fan off
G0 Z0.000 F600 ; move to print radius (Z=0 at bed surface)
;LAYER:-9999
;Z:{_bed_r:.4f}
;TYPE:SUPPORT
G0 X{_min_x:.3f} Y0.000 F{slicing_params['travel_mm_per_min']:.1f} ; move to prime position
G1 Y360.000 E{_prime_e:.5f} F{_first_feed:.1f} ; prime rotation (one full circle)
G92 E0 ; reset extruder after prime
"""
with open("output.gcode","w") as f:
    f.write(_start_gcode)

# Write gap support layers (printed before the model, from bed outward).
# Use negative layer numbers (-N … -1) so the viewer layer slider includes them.
if gap_support_layers:
    print(f"  Writing {len(gap_support_layers)} gap support layer(s)...")
first_model_ext = fold_to_canonical(
    layer_data[sorted_layer_keys[0]]["external_peris"], theta_clip_box
) if sorted_layer_keys else shapely.GeometryCollection()

for gap_layer_idx, gap_layer in enumerate(gap_support_layers):
    layer_num = gap_layer_idx - len(gap_support_layers)   # e.g. -5, -4, …, -1
    r = gap_layer["radius"]
    n_gap = len(gap_support_layers)
    gap_sup = gap_layer["support"]
    is_first = gap_layer_idx == 0
    is_last = gap_layer_idx == n_gap - 1

    interface_n = slicing_params.get("support_interface_layers", 2)
    # depth_from_bottom: 1-based distance from bed side (innermost = 1, printed first)
    # depth_from_top:    1-based distance from model side (outermost = 1, printed last)
    depth_from_bottom = gap_layer_idx + 1
    depth_from_top    = n_gap - gap_layer_idx

    # Gap support layers are the bed-to-model foundation pillar — the air-gap feature
    # does not apply here; every shell is printed.
    if depth_from_top <= interface_n and depth_from_top <= depth_from_bottom:
        # Top interface: depth 1 is outermost (printed last) — flip parity
        direction = "axial" if depth_from_top % 2 == 1 else "radial"
        support_infill = do_support_infill(gap_sup, slicing_params, r, direction=direction)
    elif depth_from_bottom <= interface_n:
        # Bottom interface: depth 1 is innermost (printed first) — radial first
        direction = "radial" if depth_from_bottom % 2 == 1 else "axial"
        support_infill = do_support_infill(gap_sup, slicing_params, r, direction=direction)
    else:
        # Gap-foundation shells already stack across many radii; printing the full
        # 2D body pattern on every shell over-densifies support dramatically.
        # Use body-spacing radial sweeps here so negative layers honor the support
        # spacing setting instead of fragmenting into hundreds of short paths.
        support_infill = do_support_infill(gap_sup, slicing_params, r, direction="body_radial")
    with open("output.gcode", "a") as f:
        f.write(f";LAYER:{layer_num}\n")
        f.write(f";Z:{r:.4f}\n")
        f.write(f"G0 Z{r - _bed_r:.4f}\n")
        _y_accum_deg  = math.degrees(last_pos[0])
        _y_canon_deg  = _y_accum_deg % 360.0
        f.write(f"G92 Y{_y_canon_deg:.6f} ; reset Y accumulator to canonical\n")
        last_pos = [math.radians(_y_canon_deg), last_pos[1]]
        gcode_paths = [p for p in shapely_to_paths(support_infill, r) if p is not None]
        gcode_paths = optimize_path_order(gcode_paths, last_pos)
        for path in gcode_paths:
            f.write(";TYPE:SUPPORT\n")
            g, last_pos, last_e = path.to_gcode(slicing_params, last_pos, last_e)
            f.write(g)

_sorted_model_keys = sorted(layer_data.keys())
for _model_idx, i in enumerate(_sorted_model_keys):
    real_layer_count += 1
    if max_real_layers is not None and real_layer_count > max_real_layers:
        print(f"\nReached {max_real_layers} layer cap, exiting.")
        break

    total = len(layer_data)
    print(f">>> Layer {real_layer_count}/{total} (layer_data key: {i})", end="\r",flush=True)

    # Per-layer figure: 2D on left, 3D on right (only when plot_per_layer)
    ax_2d_layer = None
    ax_3d_layer = None
    if plot_per_layer:
        fig_layer = plt.figure(figsize=(14, 6))
        ax_2d_layer = fig_layer.add_subplot(1, 2, 1)
        ax_3d_layer = fig_layer.add_subplot(1, 2, 2, projection='3d')
        fig_layer.canvas.mpl_connect('key_press_event', on_key)

    current_r = layer_data[i]["radius"]
    def do_transform(x):
        return shapely.transform(x,lambda p: p * [current_r,1])
    def do_transform_inverse(x):
        return shapely.transform(x,lambda p: p * [1/current_r,1])
    current_ext_print = keep_canonical_polygons(layer_data[i]["external_peris"])
    current_top_print = keep_canonical_polygons(layer_data[i].get("top_surfaces", shapely.GeometryCollection()))
    current_bottom_print = keep_canonical_polygons(layer_data[i].get("bottom_surfaces", shapely.GeometryCollection()))
    if plot_per_layer:
        plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(current_ext_print))),ax=ax_2d_layer,add_points=False,color="blue")
    final_perimeter = do_transform(current_ext_print)


    layer_data[i]["internal_peris"] = []
    for wall_number in range(1,slicing_params["wall_count"]):
        internal_perimeter = do_transform(current_ext_print).buffer(-slicing_params["nozzle_size"]*(wall_number),join_style="mitre")
        if plot_per_layer:
            plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(do_transform_inverse(internal_perimeter)))),ax=ax_2d_layer,add_points=False,color="green")
        layer_data[i]["internal_peris"].append(do_transform_inverse(internal_perimeter))
        final_perimeter = internal_perimeter



    final_top_surface = do_transform(current_top_print).buffer(slicing_params["nozzle_size"]*slicing_params["wall_count"],join_style="mitre").intersection(final_perimeter)
    if plot_per_layer:
        plot_polygon(final_top_surface,ax=ax_2d_layer,add_points=False,color="indigo")
    _solid_direction = "axial"
    layer_data[i]["solid_infill"] = do_solid_infill(do_transform_inverse(final_top_surface),slicing_params,current_r,plot_per_layer,ax=ax_2d_layer,direction=_solid_direction) or []

    final_bottom_surface = do_transform(current_bottom_print).buffer(slicing_params["nozzle_size"]*slicing_params["wall_count"],join_style="mitre").intersection(final_perimeter)
    if plot_per_layer:
        plot_polygon(final_bottom_surface,ax=ax_2d_layer,add_points=False,color="cyan")
    bottom_solid_infill = do_solid_infill(do_transform_inverse(final_bottom_surface),slicing_params,current_r,plot_per_layer,ax=ax_2d_layer,direction=_solid_direction) or []
    layer_data[i]["solid_infill"] += bottom_solid_infill

    # Support infill: interface layers (near model) alternate radial/axial direction across
    # support_interface_layers shells; body layers use the body pattern.
    support_geom = layer_data[i].get("support", shapely.GeometryCollection())
    if not support_geom.is_empty:
        interface_n = slicing_params.get("support_interface_layers", 2)
        gap_n = slicing_params.get("support_gap_layers", 0)

        # Check distance to upper model boundary (i+d) and lower model boundary (i-d) separately.
        # For bottom interface, depth 1 is the innermost layer (printed first) → radial.
        # For top interface, depth 1 is the outermost layer (printed last) → parity is flipped
        # so the first-printed layer is still radial.
        def _iface_depth(delta_sign):
            for d in range(1, interface_n + 1):
                ext = layer_data.get(i + delta_sign * d, {}).get("external_peris", shapely.GeometryCollection())
                if not ext.is_empty:
                    candidate = support_geom.intersection(fold_to_canonical(ext, theta_clip_box))
                    if not candidate.is_empty:
                        return d
            return None

        top_depth = _iface_depth(+1)   # distance to model above (outward)
        bot_depth = _iface_depth(-1)   # distance to model below (inward)

        support_clipped = support_geom.intersection(theta_clip_box)

        # Gap detection: this layer is in the air-gap zone if there is no further support
        # gap_n shells outward — i.e. we are at the top of the support column here.
        # This is independent of how close the model is, so it works correctly for sloped
        # surfaces (like a cone) where the model is always at d=1 above every support layer.
        is_gap_layer = gap_n > 0 and layer_data.get(i + gap_n, {}).get("support", shapely.GeometryCollection()).is_empty

        if is_gap_layer:
            # Air-gap zone: leave these layers empty so the part separates cleanly
            layer_data[i]["support_infill"] = []
        else:
            if top_depth is not None and (bot_depth is None or top_depth <= bot_depth):
                # Top interface: depth 1 is outermost (printed last), so flip parity
                # → odd depth = axial, even depth = radial (first-printed layer = radial)
                direction = "axial" if top_depth % 2 == 1 else "radial"
            elif bot_depth is not None:
                # Bottom interface: depth 1 is innermost (printed first) → radial
                direction = "radial" if bot_depth % 2 == 1 else "axial"
            else:
                direction = "body"

            layer_data[i]["support_infill"] = do_support_infill(support_clipped, slicing_params, current_r, direction=direction)

        if plot_per_layer:
            plot_polygon(do_transform(support_geom), ax=ax_2d_layer, add_points=False, color="magenta")
    else:
        layer_data[i]["support_infill"] = []

    if "needs_support" in layer_data[i] and plot_per_layer:
        plot_polygon(do_transform(layer_data[i]["needs_support"]),ax=ax_2d_layer,add_points=False,color="red")

    # plt.savefig(f"frames/{i}.png",dpi=600)

    # Inset by half a nozzle so sparse infill lines land just inside the innermost
    # perimeter without overlapping it (its center line is at the boundary of final_perimeter).
    _infill_boundary = final_perimeter.buffer(-slicing_params["nozzle_size"] * 0.5, join_style="mitre")
    _diff = _infill_boundary - final_top_surface - final_bottom_surface
    if not _diff.is_empty:
        # Morphological opening: removes any remaining floating-point slivers
        # narrower than half a nozzle width.
        _half_noz = slicing_params["nozzle_size"] * 0.5
        _diff = _diff.buffer(-_half_noz, join_style="mitre").buffer(_half_noz, join_style="mitre")
    shape_to_infill = do_transform_inverse(_diff)

    fully_radial_infill, radial_infill, axial_infill = do_sparse_infill(shape_to_infill,slicing_params,current_r,plot_per_layer,ax=ax_2d_layer)
    layer_data[i]["fully_radial_infill"] = fully_radial_infill
    layer_data[i]["radial_infill"] = radial_infill
    layer_data[i]["axial_infill"] = axial_infill

    # 3D plotting helpers
    def plot_lines_3d(shapely_obj, ax, color, linewidth, label=None):
        """Polygon/MultiPolygon: clip to [0,2pi], densify arcs, convert cylindrical -> Cartesian."""
        lines = polygon_to_wrapped_lines(shapely_obj)
        for k, line in enumerate(lines):
            coords = densify_cylindrical(np.array(line.coords))
            x = current_r * np.cos(coords[:, 0])
            y = current_r * np.sin(coords[:, 0])
            ax.plot(x, y, coords[:, 1], color=color, linewidth=linewidth,
                    label=label if k == 0 else None)

    def plot_linelist_3d(line_list, ax, color, linewidth, label=None):
        """List of LineStrings already in (theta, z) space — densify arcs before plotting."""
        if not line_list:
            return
        for k, line in enumerate(line_list):
            if line is None or line.is_empty:
                continue
            coords = densify_cylindrical(np.array(line.coords))
            x = current_r * np.cos(coords[:, 0])
            y = current_r * np.sin(coords[:, 0])
            ax.plot(x, y, coords[:, 1], color=color, linewidth=linewidth,
                    label=label if k == 0 else None)

    def plot_line_endpoints_3d(line_list, ax, marker_size=3, start_color='lime', end_color='red'):
        """Plot spheres at start and end of each LineString in 3D space."""
        if not line_list:
            return
        for line in line_list:
            if line is None or line.is_empty:
                continue
            coords = densify_cylindrical(np.array(line.coords))
            if len(coords) < 2:
                continue
            # Convert to Cartesian
            x = current_r * np.cos(coords[:, 0])
            y = current_r * np.sin(coords[:, 0])
            z = coords[:, 1]
            # Start point (green sphere)
            ax.scatter(x[0], y[0], z[0], c=start_color, s=marker_size**2, alpha=0.7)
            # End point (red sphere)
            ax.scatter(x[-1], y[-1], z[-1], c=end_color, s=marker_size**2, alpha=0.7)

    def set_equal_aspect_3d(ax):
        limits = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
        spans = limits[:, 1] - limits[:, 0]
        spans = np.where(spans > 0, spans, 1.0)
        ax.set_box_aspect(spans / spans.max())

    top_surf = layer_data[i].get("top_surfaces")

    # --- Per-layer 3D (only when plot_per_layer) ---
    if plot_3d and plot_per_layer:
        plot_lines_3d(layer_data[i]["external_peris"], ax_3d_layer, 'blue', 2, label='Ext Peri')
        for j, int_peri in enumerate(layer_data[i].get("internal_peris", [])):
            plot_lines_3d(int_peri, ax_3d_layer, 'green', 1.5, label='Int Peri' if j == 0 else None)
        if top_surf is not None and not top_surf.is_empty:
            plot_lines_3d(top_surf, ax_3d_layer, 'red', 1.5, label='Top Surf')
        plot_linelist_3d(layer_data[i].get("solid_infill", []),        ax_3d_layer, 'greenyellow', 0.8, label='Solid Infill')
        plot_linelist_3d(layer_data[i].get("radial_infill", []),       ax_3d_layer, 'orange',      0.8, label='Radial Infill')
        plot_linelist_3d(layer_data[i].get("axial_infill", []),        ax_3d_layer, 'darkviolet',  0.8, label='Axial Infill')
        plot_linelist_3d(layer_data[i].get("fully_radial_infill", []), ax_3d_layer, 'orangered',   0.8, label='Full Radial')
        # Plot 3D endpoint markers
        plot_line_endpoints_3d(layer_data[i].get("solid_infill", []), ax_3d_layer, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("radial_infill", []), ax_3d_layer, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("axial_infill", []), ax_3d_layer, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("fully_radial_infill", []), ax_3d_layer, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("support_infill", []), ax_3d_layer, marker_size=3, start_color='cyan', end_color='magenta')
        set_equal_aspect_3d(ax_3d_layer)
        ax_3d_layer.set_xlabel('X'); ax_3d_layer.set_ylabel('Y'); ax_3d_layer.set_zlabel('Z')
        ax_3d_layer.set_title(f'Layer {i} (r={current_r:.2f})')
        ax_3d_layer.legend()

    # --- Combined 3D (always when plot_3d, non-blocking live update) ---
    if plot_3d:
        plot_lines_3d(layer_data[i]["external_peris"], ax_3d, 'blue', 2, label=f'Ext Peri L{i}')
        for j, int_peri in enumerate(layer_data[i].get("internal_peris", [])):
            plot_lines_3d(int_peri, ax_3d, 'green', 1.5, label=f'Int Peri L{i}' if j == 0 else None)
        if top_surf is not None and not top_surf.is_empty:
            plot_lines_3d(top_surf, ax_3d, 'red', 1.5, label=f'Top Surf L{i}')
        plot_linelist_3d(layer_data[i].get("solid_infill", []),        ax_3d, 'greenyellow', 0.8, label=f'Solid Infill L{i}')
        plot_linelist_3d(layer_data[i].get("radial_infill", []),       ax_3d, 'orange',      0.8, label=f'Radial Infill L{i}')
        plot_linelist_3d(layer_data[i].get("axial_infill", []),        ax_3d, 'darkviolet',  0.8, label=f'Axial Infill L{i}')
        plot_linelist_3d(layer_data[i].get("fully_radial_infill", []), ax_3d, 'orangered',   0.8, label=f'Full Radial L{i}')
        # Plot 3D endpoint markers
        plot_line_endpoints_3d(layer_data[i].get("solid_infill", []), ax_3d, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("radial_infill", []), ax_3d, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("axial_infill", []), ax_3d, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("fully_radial_infill", []), ax_3d, marker_size=4, start_color='lime', end_color='red')
        plot_line_endpoints_3d(layer_data[i].get("support_infill", []), ax_3d, marker_size=3, start_color='cyan', end_color='magenta')
        set_equal_aspect_3d(ax_3d)
        ax_3d.set_title(f'Up to layer {i} (r={current_r:.2f})')
        fig_combined.canvas.draw()
        plt.pause(0.001)

    # --- Finalize and show per-layer figure (blocks until closed) ---
    if plot_per_layer:
        ax_2d_layer.set_title(f'Layer {i}')
        ax_2d_layer.set_aspect('equal')
        fig_layer.tight_layout()
        fig_layer.canvas.draw()
        print(f"Close layer {i} window to continue...")
        while plt.fignum_exists(fig_layer.number):
            plt.pause(0.1)

    # G-code generation
    traces = []
    sections = [
        ("external_peris", "EXTERNAL_PERIMETER"),
        ("internal_peris", "INTERNAL_PERIMETER"),
        ("solid_infill", "SOLID_INFILL"),
        ("support_infill", "SUPPORT"),
        ("fully_radial_infill", "INFILL_RADIAL_FULL"),
        ("radial_infill", "INFILL_RADIAL"),
        ("axial_infill", "INFILL_AXIAL")
    ]
    # Temporarily scale speeds for the first model layer
    _speed_factor = slicing_params["first_layer_speed_factor"] if i == 0 else 1.0
    _layer_params = slicing_params.copy()
    _layer_params["feed_mm_per_min"]   *= _speed_factor
    _layer_params["travel_mm_per_min"] *= _speed_factor

    with open("output.gcode","a") as f:
        f.write(f";LAYER:{i}\n")
        f.write(f";Z:{current_r}\n")
        if i == slicing_params["fan_start_layer"]:
            f.write("M106 S255 ; part cooling fan on\n")
        f.write(f"G0 Z{current_r - _bed_r:.4f}\n")
        # Reset accumulated Y to its canonical [0, 360) equivalent at each layer
        # boundary.  At high layer counts last_pos[0] grows to millions of degrees,
        # causing float64 precision loss that manifests as infill gaps and rendering
        # loops in the viewer.  G92 Y<canon> tells the printer "this position is now
        # <canon>" without moving, and we update last_pos to match.
        _y_accum_deg  = math.degrees(last_pos[0])
        _y_canon_deg  = _y_accum_deg % 360.0
        f.write(f"G92 Y{_y_canon_deg:.6f} ; reset Y accumulator to canonical\n")
        last_pos = [math.radians(_y_canon_deg), last_pos[1]]

        # Precompute shapely_to_paths + optimize_path_order for all sections in
        # parallel. All sections get the same last_pos snapshot (start of layer);
        # actual last_pos is then updated sequentially below during to_gcode.
        _snap_pos = last_pos
        def _prep_section(subsection):
            if subsection in ("external_peris", "internal_peris"):
                raw = shapely_to_paths(polygon_to_wrapped_lines(layer_data[i][subsection]), current_r)
                return [p for p in raw if p is not None]
            else:
                raw = shapely_to_paths(layer_data[i].get(subsection, []), current_r)
                paths = [p for p in raw if p is not None]
                return optimize_path_order(paths, _snap_pos)

    
        with ThreadPoolExecutor() as _ex:
            _futs = {s: _ex.submit(_prep_section, s) for s, _ in sections}
            _prepped = {s: _futs[s].result() for s, _ in sections}

        _support_params = _layer_params.copy()
        _support_params["feed_mm_per_min"]   *= slicing_params.get("support_speed_factor", 1.0)
        _support_params["travel_mm_per_min"] *= slicing_params.get("support_speed_factor", 1.0)
        _support_params["flow_rate"]          = slicing_params.get("support_flow_rate", _layer_params.get("flow_rate", 1.0))

        _infill_params = _layer_params.copy()
        _infill_params["feed_mm_per_min"]   *= slicing_params.get("infill_speed_factor", 1.0)
        _infill_params["travel_mm_per_min"] *= slicing_params.get("infill_speed_factor", 1.0)
        _infill_params["flow_rate"]          = slicing_params.get("infill_flow_rate", _layer_params.get("flow_rate", 1.0))

        _infill_sections = {"fully_radial_infill", "radial_infill", "axial_infill"}

        for subsection, section_name in sections:
            if subsection == "support_infill":
                _path_params = _support_params
            elif subsection in _infill_sections:
                _path_params = _infill_params
            else:
                _path_params = _layer_params
            for path in _prepped[subsection]:
                f.write(f";TYPE:{section_name}\n")
                # For fully_radial_infill (full circles), rotate to start at current theta
                # to avoid unnecessary travel back to theta=0
                if subsection == "fully_radial_infill":
                    path = rotate_full_circle_path(path, last_pos)
                # For perimeters (external/internal), rotate to start at nearest point
                # Check by subsection name since paths from polygon_to_wrapped_lines are open LineStrings
                elif subsection in ("external_peris", "internal_peris"):
                    path = rotate_closed_path_to_nearest(path, last_pos)
                g,last_pos,last_e = path.to_gcode(_path_params,last_pos,last_e)
                f.write(g)

    _ns  = layer_data[i].get("needs_support")
    _sup = layer_data[i].get("support", shapely.GeometryCollection())
    _solid_n  = len(layer_data[i].get("solid_infill", []))
    _sparse_n = (len(layer_data[i].get("radial_infill", [])) +
                 len(layer_data[i].get("axial_infill", [])) +
                 len(layer_data[i].get("fully_radial_infill", [])))
    _sup_n = len(layer_data[i].get("support_infill", []))
    with open("stats.txt", "a") as _sf:
        _sf.write(f"{i:>6}  {current_r:>8.3f}"
                  f"  {layer_data[i]['external_peris'].area:>10.5f}"
                  f"  {(_ns.area if _ns else 0.0):>10.5f}"
                  f"  {_sup.area:>10.5f}"
                  f"  {_sup_n:>9}"
                  f"  {_solid_n:>8}"
                  f"  {_sparse_n:>9}\n")

    # traces = [trace for trace in traces if trace is not None]
    # traces = np.vstack(traces)
    # for each in traces:
    #     output_scad += solid2.translate(each)(solid2.sphere())
    # print(traces.shape)
    # scatter = visuals.Markers()
    # scatter.set_data(
    #     traces,             # Nx3 array of positions
    #     face_color="red",    # can be a single color or Nx4 array
    #     size=2               # size in pixels
    # )
    # scatter.parent = view.scene  # attach to the view


    # fig = go.Figure(data=traces)
    # fig.show()
# set_axes_equal(ax)
    # Helper functions for endpoint markers
    def plot_line_endpoints_2d(line_list, ax, marker_size=8, color='red'):
        """Plot circles at start (green) and end (red) of each LineString."""
        if not line_list:
            return
        for line in line_list:
            if line is None or line.is_empty:
                continue
            coords = np.array(line.coords)
            if len(coords) < 2:
                continue
            # Transform to metric space for plotting
            coords_metric = coords * np.array([current_r, 1])
            # Start point (green circle)
            ax.plot(coords_metric[0, 0], coords_metric[0, 1], 'go', 
                    markersize=marker_size, alpha=0.7)
            # End point (red circle)
            ax.plot(coords_metric[-1, 0], coords_metric[-1, 1], 'ro', 
                    markersize=marker_size, alpha=0.7)

    # Plot endpoints for all infill types in 2D
    if plot_per_layer:
        plot_line_endpoints_2d(layer_data[i].get("solid_infill", []), ax_2d_layer, marker_size=6, color='greenyellow')
        plot_line_endpoints_2d(layer_data[i].get("radial_infill", []), ax_2d_layer, marker_size=6, color='orange')
        plot_line_endpoints_2d(layer_data[i].get("axial_infill", []), ax_2d_layer, marker_size=6, color='darkviolet')
        plot_line_endpoints_2d(layer_data[i].get("fully_radial_infill", []), ax_2d_layer, marker_size=6, color='orangered')
        plot_line_endpoints_2d(layer_data[i].get("support_infill", []), ax_2d_layer, marker_size=5, color='magenta')

    # plt.clf()

with open("output.gcode", "a") as f:
    f.write(
        "\nM104 S0 ; hotend off\n"
        "M106 S0 ; fan off\n"
        "G91 ; relative mode\n"
        "G0 Z5 F600 ; retract radially away from print\n"
        "G90 ; absolute mode\n"
        "M84 ; disable steppers\n"
        ";end of slice\n"
    )
print(f"\nDone. Wrote {real_layer_count} model layer(s) to output.gcode.")

if plot_3d and ax_3d:
    ax_3d.set_xlabel('X')
    ax_3d.set_ylabel('Y')
    ax_3d.set_zlabel('Z (height)')
    ax_3d.set_title('All Layers')
    ax_3d.legend()
    plt.ioff()
    fig_combined.tight_layout()
    plt.show()

# mesh.show()
# output_scad.save_as_scad()
