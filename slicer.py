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
from togcode import shapely_to_paths, optimize_path_order
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
    coords -= [np.floor(np.min(coords[:,0])/(2*np.pi))*2*np.pi,0]
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
    Compute the cross-section of mesh at cylindrical radius using marching triangles.
    Works on any mesh regardless of manifold status — no boolean ops needed.
    Returns (pts_cartesian, edges) where pts_cartesian is (N,3) Cartesian points
    at r≈radius and edges is a list of (i,j) index pairs forming closed loops.
    """
    verts = mesh.vertices
    r_vals = np.sqrt(verts[:, 0]**2 + verts[:, 1]**2)

    pts = []
    edges = []
    edge_cache = {}

    def get_crossing(vi, vj):
        key = (min(vi, vj), max(vi, vj))
        if key not in edge_cache:
            ri, rj = r_vals[vi], r_vals[vj]
            t = (radius - ri) / (rj - ri)
            edge_cache[key] = len(pts)
            pts.append(verts[vi] + t * (verts[vj] - verts[vi]))
        return edge_cache[key]

    for v0, v1, v2 in mesh.faces:
        inside = [r_vals[v] < radius for v in (v0, v1, v2)]
        if all(inside) or not any(inside):
            continue
        crossing = []
        for a, b in [(v0, v1), (v1, v2), (v2, v0)]:
            if (r_vals[a] < radius) != (r_vals[b] < radius):
                crossing.append(get_crossing(a, b))
        if len(crossing) == 2 and crossing[0] != crossing[1]:
            edges.append(tuple(crossing))

    # Prune dangling edges caused by open (non-watertight) mesh boundaries.
    # Boundary mesh edges belong to only one face; their crossing points end up
    # with odd degree (degree=1) in the polygon graph, breaking find_polygons.
    # Iteratively remove every edge that touches an odd-degree vertex until all
    # remaining vertices have even degree (i.e. form closed loops).
    while True:
        degree = defaultdict(int)
        for a, b in edges:
            degree[a] += 1
            degree[b] += 1
        odd = {v for v, d in degree.items() if d % 2 != 0}
        if not odd:
            break
        edges = [(a, b) for a, b in edges if a not in odd and b not in odd]

    return (np.array(pts) if pts else np.empty((0, 3))), edges


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


def do_support_infill(support_geom, slicing_params, current_r, is_interface=False):
    """
    Generate support infill in (θ,z) space.

    Interface layers (top/bottom of a support column):
        axial + radial grid at support_interface_spacing.
    Body layers, mode "grid":
        axial + radial grid at support_body_spacing.
    Body layers, mode "sine":
        parallel sine waves along z, spaced support_body_spacing apart in arc-length.
        Amplitude and period controlled by support_sine_amplitude_mm / support_sine_period_mm.

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

    if is_interface:
        spacing = slicing_params["support_interface_spacing"]
        _ts = spacing / r_base
        thetas = np.arange(int(np.floor(2 * np.pi / _ts))) * _ts

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
                    _collect_lines(seg, result, min_length)
            else:
                _collect_lines(inter, result, min_length)
    else:
        spacing = slicing_params["support_body_spacing"]
        _ts = spacing / r_base
        thetas = np.arange(int(np.floor(2 * np.pi / _ts))) * _ts

        mode = slicing_params.get("support_infill_mode", "grid")

        if mode == "sine":
            amplitude_rad = slicing_params["support_sine_amplitude_mm"] / current_r
            period_mm = slicing_params["support_sine_period_mm"]
            # Extend z slightly so waves fully cover the geometry bounds
            z_lo = z_min - period_mm
            z_hi = z_max + period_mm
            for theta in thetas:
                # Generate ±2π shifted copies so waves that cross the 0/2π seam
                # are captured correctly — same seam-handling strategy as the rest
                # of the codebase.
                for shift in (0.0, -2 * np.pi, 2 * np.pi):
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
                        _collect_lines(seg, result, min_length)
                else:
                    _collect_lines(inter, result, min_length)

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
        # linspace distributes the leftover evenly instead of leaving a single
        # wide gap at the theta=0/2π seam.
        for theta in np.linspace(0, 2 * np.pi, _n_solid, endpoint=False):
            ls = shapely.LineString([[theta, bounds[1]], [theta, bounds[3]]])
            _collect(shapely.intersection(ls, peri))
    else:
        # Constant-z lines — arcs around the cylinder (radial / circumferential)
        z_step = slicing_params["nozzle_size"]   # mm; z axis is already in mm
        z_lo = np.floor(bounds[1] / z_step) * z_step
        for z in np.arange(z_lo, bounds[3] + z_step, z_step):
            ls = shapely.LineString([[0, z], [2 * np.pi, z]])
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

    # ── Pass 1: Radial infill (constant-z lines, θ = 0 → 2π) ───────────────
    # Sweep horizontal lines from lowest to highest z, spaced by infill_spacing.
    for current_line_v in np.arange(lowest, highest + 1, slicing_params["infill_spacing"]):
        # Build a full-width horizontal line at this z value across the whole θ range.
        lp_x, lp_y = [0, 2*np.pi], [current_line_v, current_line_v]
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
            stitched = stitch_radial_seam(geoms)
            if len(stitched) == 1:
                # Seam-crossing arc stitched into one segment → fully radial.
                if plot_infill:
                    plot_line(do_transform(stitched[0]), ax=ax, color="red", add_points=False)
                fully_radial_infill.append(stitched[0])
            else:
                for line_segment in stitched:
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
    n_base = int(np.floor(2 * np.pi / fixed_step))             # same on every layer
    all_thetas = np.linspace(0, 2 * np.pi, n_base, endpoint=False)

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
            # Multiple disjoint segments (line crossed a hollow region mid-way).
            for line_segment in geoms:
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

def polygon_to_wrapped_lines(poly, lower=0.0, upper=2*math.pi):
    """
    Convert a polygon (or collection) to LineStrings for gcode/plotting.

    external_peris stores each physical polygon THREE times (canonical + ±2π
    copies).  We keep only the canonical copy by checking that the centroid
    theta is within [lower, upper).  The ring coords are then unwrapped with
    np.unwrap so any seam-crossing polygon is returned as one continuous
    LineString — no clipping, no stitching.
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
        # Already a line (e.g. infill passed through here) — return as-is
        return [poly] if not poly.is_empty else []

    if not isinstance(poly, Polygon):
        return []

    # Skip ±2π duplicate copies — only keep the canonical polygon
    if not (lower <= poly.centroid.x < upper):
        return []

    result = []
    for ring in [poly.exterior] + list(poly.interiors):
        coords = np.array(ring.coords)
        raw_span = coords[:, 0].max() - coords[:, 0].min()
        if raw_span > 2 * math.pi - 0.5:
            # Full-circle ring: spans more than one period.  Clip to [0, 2π],
            # then emit only the arc segments (varying θ) as separate open
            # LineStrings — skip the constant-θ seam edges so no extra material
            # is deposited at the seam.  to_gcode handles the ±2π step in each
            # arc via its improved unwrap logic.
            canonical_clip = shapely.geometry.box(0, -1e15, 2 * math.pi, 1e15)
            clipped = poly.intersection(canonical_clip)
            if clipped.is_empty or not isinstance(clipped, Polygon):
                continue
            rc = np.array(clipped.exterior.coords)
            for k in range(len(rc) - 1):
                p0, p1 = rc[k], rc[k + 1]
                dtheta = p1[0] - p0[0]
                if abs(dtheta) > 0.1:          # arc segment (not a seam line)
                    seg = [p1, p0] if dtheta < 0 else [p0, p1]  # always θ-increasing
                    result.append(LineString(seg))
            continue
        else:
            coords[:, 0] = np.unwrap(coords[:, 0])
        # Skip degenerate rings whose θ-extent is negligibly thin (e.g. a
        # near-tangent cylindrical cross-section of a flat face).  Such rings
        # produce a back-and-forth line in G-code rather than a closed loop.
        # Threshold: 1e-3 rad ≈ 0.06° — well below any printable feature but
        # above floating-point noise.
        if coords[:, 0].max() - coords[:, 0].min() < 1e-3:
            continue
        ls = LineString(coords)
        if not ls.is_empty:
            result.append(ls)
    return result
    
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
    "feed_mm_per_min" : 2400,
    "travel_mm_per_min" : 2400,
    "support_infill_mode": "sine",   # "grid" or "sine"
    "support_sine_amplitude_mm": 2.0, # mm — sine wave half-amplitude (grid mode ignores this)
    "support_sine_period_mm": 8.0,   # mm — sine wave period along z (grid mode ignores this)
    "support_interface_spacing": 1,  # mm — grid spacing for top/bottom surfaces of support
    "support_body_spacing": 5,       # mm — line spacing for support interior (arc-length mm)
    "support_min_area": 1.0,         # mm² — inflate overhang patches below this area to at least this size
    "support_merge_gap": 15.0,       # mm — merge support branches closer than this into a shared trunk
    "support_taper_per_layer_mm": -0.1,  # mm — how much the support footprint shrinks per layer going inward
    "support_xy_distance_mm": 0.3,   # mm — lateral gap between support and model walls
    "retract_mm": 2.0,               # mm — retraction distance on travel moves
    "retract_speed_mm_per_min": 1800, # mm/min — retraction speed (30 mm/s)
    "min_retract_travel_mm": 6.0,    # mm — skip retraction for travel moves shorter than this
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
    "hotend_temp": 205,    # °C — PLA
    "min_x_pos": 0.0,    # mm — X coordinate where the printable bed begins (prime line here)
    "min_x_margin": 5.0,  # mm — gap between bed start and model start
}
# Derived: how many mm of filament per mm of real travel
slicing_params["extrusion_per_mm"] = (
    slicing_params["nozzle_size"] * slicing_params["layer_height"]
    / (np.pi * (printer_params["filament_diameter"] / 2) ** 2)
)
max_real_layers = None  # Set to an int to cap layers for testing, None for full print
mesh = trimesh.load_mesh(sys.argv[1] if len(sys.argv) > 1 else "offset_cube.stl")
if not mesh.is_watertight:
    print(f"Mesh is not watertight — attempting repair...")
    trimesh.repair.fill_holes(mesh)
    trimesh.repair.fix_normals(mesh)
    print(f"  after repair: watertight={mesh.is_watertight}")
pts = cartesian_to_cylindrical(mesh.vertices)
# Shift mesh so its axial start aligns with min_x_pos + min_x_margin
_z_shift = printer_params["min_x_pos"] + printer_params["min_x_margin"] - np.min(pts[:, 2])
mesh.apply_translation([0, 0, _z_shift])

# Subdivide large faces so that cylindrical cross-sections have enough sample points
# along long edges (avoids jitter from straight-line approximation of curved arcs).
_max_edge = slicing_params["layer_height"] * 20  # subdivide any edge longer than 5 layers
_orig_faces = len(mesh.faces)
mesh = trimesh.Trimesh(*trimesh.remesh.subdivide_to_size(mesh.vertices, mesh.faces, _max_edge), process=False)
print(f"Subdivided mesh: {_orig_faces} → {len(mesh.faces)} faces (max edge {_max_edge:.2f} mm)")
pts = cartesian_to_cylindrical(mesh.vertices)
max_z = np.max(pts[:,2])
min_z = np.min(pts[:,2])
part_height = (max_z-min_z)

# Z bounds of the printable volume in the same coordinate space as the mesh
z_clip_min = min_z
z_clip_max = min_z + printer_params["axial_length"]

# layer (θ,z)-space clipping box — wide enough to cover all θ shifts (±2π copies)
axial_clip_box = box(-2*np.pi, z_clip_min, 4*np.pi, z_clip_max)
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
    cross_pts_cart, outer_edges = cylindrical_cross_section(mesh, current_r)
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
    poly_pts = [unwrap_polygon(layer_pts[each,1:]).squeeze() for each in poly_idxs]
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
        mapping += [len(polys)]

        final_poly = Polygon(final_poly)
        final_poly = trim_polygon_left_right(final_poly,0.1)
        polys.append(final_poly)
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

    for poly_i,sub_polys in enumerate(hier):
        if levels[poly_i] % 2 == 1:
            continue
        sub_poly_list = [np.array(polys[sub_poly].exterior.coords).astype(float) for sub_poly in sub_polys if levels[sub_poly] == levels[poly_i]+1]
        sub_poly_list = [[each - [2*np.pi,0],each,each + [2*np.pi,0]] for each in sub_poly_list]
        sub_poly_list = [remove_collinear_points_angle(each2*np.array([current_r,1]),0.01)*np.array([1/current_r,1]) for each in sub_poly_list for each2 in each]
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
        if not ns.is_empty:
            model_mask = fold_to_canonical(final, theta_clip_box)
            if not model_mask.is_empty:
                ns = ns.difference(model_mask)
        if not ns.is_empty:
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

    # Ground: subtract model, expanded by xy gap so support doesn't press against walls
    if xy_dist_mm > 0 and not current_ext.is_empty:
        metric_ext     = shapely.affinity.scale(current_ext, xfact=current_r, yfact=1.0, origin=(0, 0))
        metric_exp     = metric_ext.buffer(xy_dist_mm)
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

    if n_gap <= 2:
        # Very few gap layers — treat everything as interface
        support_infill = do_support_infill(gap_sup, slicing_params, r, is_interface=True)
    elif is_last:
        # Top of gap column: only the part under the model is interface
        top_iface = gap_sup.intersection(first_model_ext) if not first_model_ext.is_empty else shapely.GeometryCollection()
        body = gap_sup.difference(top_iface) if not top_iface.is_empty else gap_sup
        support_infill = (
            do_support_infill(top_iface, slicing_params, r, is_interface=True) +
            do_support_infill(body, slicing_params, r, is_interface=False)
        )
    elif is_first:
        support_infill = do_support_infill(gap_sup, slicing_params, r, is_interface=True)
    else:
        support_infill = do_support_infill(gap_sup, slicing_params, r, is_interface=False)
    with open("output.gcode", "a") as f:
        f.write(f";LAYER:{layer_num}\n")
        f.write(f";Z:{r:.4f}\n")
        f.write(f"G0 Z{r - _bed_r:.4f}\n")
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
    if plot_per_layer:
        plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(layer_data[i]["external_peris"]))),ax=ax_2d_layer,add_points=False,color="blue")
    final_perimeter = do_transform(layer_data[i]["external_peris"])
    layer_data[i]["internal_peris"] = []
    for wall_number in range(1,slicing_params["wall_count"]):
        internal_perimeter = do_transform(layer_data[i]["external_peris"]).buffer(-slicing_params["nozzle_size"]*(wall_number),join_style="mitre")
        if plot_per_layer:
            plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(do_transform_inverse(internal_perimeter)))),ax=ax_2d_layer,add_points=False,color="green")
        layer_data[i]["internal_peris"].append(do_transform_inverse(internal_perimeter))
        final_perimeter = internal_perimeter



    final_top_surface = do_transform(layer_data[i]["top_surfaces"]).buffer(slicing_params["nozzle_size"]*slicing_params["wall_count"],join_style="mitre").intersection(final_perimeter)
    # Clip to canonical x=[0,2π·r] to remove thin-strip artifacts where the buffer
    # bleeds into the ±2π copies of final_perimeter, which would create double infill lines near the seam.
    _canonical_x_clip = shapely.geometry.box(0, z_clip_min - 10, 2 * np.pi * current_r, z_clip_max + 10)
    final_top_surface = final_top_surface.intersection(_canonical_x_clip)
    if plot_per_layer:
        plot_polygon(final_top_surface,ax=ax_2d_layer,add_points=False,color="indigo")
    _solid_direction = "axial" if _model_idx % 2 == 0 else "radial"
    layer_data[i]["solid_infill"] = do_solid_infill(do_transform_inverse(final_top_surface),slicing_params,current_r,plot_per_layer,ax=ax_2d_layer,direction=_solid_direction) or []

    final_bottom_surface = do_transform(layer_data[i].get("bottom_surfaces", shapely.GeometryCollection())).buffer(slicing_params["nozzle_size"]*slicing_params["wall_count"],join_style="mitre").intersection(final_perimeter)
    final_bottom_surface = final_bottom_surface.intersection(_canonical_x_clip)
    if plot_per_layer:
        plot_polygon(final_bottom_surface,ax=ax_2d_layer,add_points=False,color="cyan")
    bottom_solid_infill = do_solid_infill(do_transform_inverse(final_bottom_surface),slicing_params,current_r,plot_per_layer,ax=ax_2d_layer,direction=_solid_direction) or []
    layer_data[i]["solid_infill"] += bottom_solid_infill

    # Support infill: interface (top/bottom of column) gets a grid; body gets axial lines only.
    support_geom = layer_data[i].get("support", shapely.GeometryCollection())
    if not support_geom.is_empty:
        # Top interface: support here that directly contacts model on the next (outer) layer
        next_ext = layer_data.get(i + 1, {}).get("external_peris", shapely.GeometryCollection())
        if not next_ext.is_empty:
            top_interface = support_geom.intersection(fold_to_canonical(next_ext, theta_clip_box))
        else:
            top_interface = shapely.GeometryCollection()
        # Bottom interface: support here that directly contacts model on the previous (inner) layer
        prev_ext = layer_data.get(i - 1, {}).get("external_peris", shapely.GeometryCollection())
        if not prev_ext.is_empty:
            bottom_interface = support_geom.intersection(fold_to_canonical(prev_ext, theta_clip_box))
        else:
            bottom_interface = shapely.GeometryCollection()

        interface = top_interface.union(bottom_interface).intersection(theta_clip_box)
        body = support_geom.difference(interface).intersection(theta_clip_box)

        layer_data[i]["support_infill"] = (
            do_support_infill(interface, slicing_params, current_r, is_interface=True) +
            do_support_infill(body,      slicing_params, current_r, is_interface=False)
        )
        if plot_per_layer:
            plot_polygon(do_transform(support_geom), ax=ax_2d_layer, add_points=False, color="magenta")
    else:
        layer_data[i]["support_infill"] = []

    if "needs_support" in layer_data[i] and plot_per_layer:
        plot_polygon(do_transform(layer_data[i]["needs_support"]),ax=ax_2d_layer,add_points=False,color="red")

    # plt.savefig(f"frames/{i}.png",dpi=600)

    _diff = final_perimeter - final_top_surface - final_bottom_surface
    # Clip to the canonical [0, 2π·r] arc-length window.  external_peris carries
    # ±2π shifted copies so final_perimeter spans [−2π·r, 4π·r], but
    # final_top_surface is already clipped to [0, 2π·r].  Without this step the
    # difference leaves two strips (one at x<0, one at x>2π·r) whose boundaries
    # sit exactly at θ=0 and θ=2π — causing the first sparse-infill axial scan
    # to graze them and emit a stray line at the seam.
    _diff = _diff.intersection(_canonical_x_clip)
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

        # Precompute shapely_to_paths + optimize_path_order for all sections in
        # parallel. All sections get the same last_pos snapshot (start of layer);
        # actual last_pos is then updated sequentially below during to_gcode.
        _snap_pos = last_pos
        def _prep_section(subsection):
            if subsection in ("external_peris", "internal_peris"):
                raw = shapely_to_paths(polygon_to_wrapped_lines(layer_data[i][subsection]), current_r)
            else:
                raw = shapely_to_paths(layer_data[i].get(subsection, []), current_r)
            paths = [p for p in raw if p is not None]
            return optimize_path_order(paths, _snap_pos)

    
        with ThreadPoolExecutor() as _ex:
            _futs = {s: _ex.submit(_prep_section, s) for s, _ in sections}
            _prepped = {s: _futs[s].result() for s, _ in sections}

        for subsection, section_name in sections:
            for path in _prepped[subsection]:
                f.write(f";TYPE:{section_name}\n")
                g,last_pos,last_e = path.to_gcode(_layer_params,last_pos,last_e)
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
