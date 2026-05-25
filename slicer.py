import sys
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
import argparse

def parse_layer_spec(spec):
    """Parse comma-separated layer ids/ranges like '0,5,10-12'."""
    layers = set()
    if not spec:
        return layers
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part[1:]:
            start_s, end_s = part.split("-", 1)
            start = int(start_s)
            end = int(end_s)
            step = 1 if end >= start else -1
            layers.update(range(start, end + step, step))
        else:
            layers.add(int(part))
    return layers


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", nargs="?", default="cubeandring.stl")
    parser.add_argument(
        "--debug-layers",
        default="",
        help="Comma-separated model layer keys/ranges to export as PNGs, e.g. '0,12,20-25'.",
    )
    parser.add_argument(
        "--debug-dir",
        default="frames",
        help="Directory for --debug-layers PNG output.",
    )
    parser.add_argument(
        "--debug-plot-3d",
        action="store_true",
        help="Include a 3D subplot in exported per-layer debug images.",
    )
    parser.add_argument(
        "--show-debug-layers",
        action="store_true",
        help="Show selected debug layer windows interactively after saving.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable interactive plot windows. Debug PNGs are still written.",
    )
    parser.add_argument(
        "--max-layer",
        type=int,
        default=None,
        help=(
            "Highest 0-based model layer key to slice and emit. "
            "Layers above this are ignored, so this layer becomes the outer/top limit."
        ),
    )
    parsed = parser.parse_args()
    if parsed.max_layer is not None and parsed.max_layer < 0:
        parser.error("--max-layer must be >= 0")
    return parsed


args = parse_args()
debug_layers = parse_layer_spec(args.debug_layers)
defer_single_debug_layer_plots = len(debug_layers) == 1
if debug_layers:
    os.makedirs(args.debug_dir, exist_ok=True)

import gc
from dataclasses import dataclass
import trimesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from collections import Counter
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString, MultiLineString, box
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

def save_raw_stl_layer_plot(layer_i, current_r, slice_r, poly_pts_raw):
    """Save STL cross-section loops before slicer polygon repair/buffering."""
    if layer_i not in debug_layers:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    for poly in poly_pts_raw:
        coords = np.asarray(poly, dtype=float)
        if len(coords) < 2:
            continue
        metric = coords * np.array([current_r, 1])
        ax.plot(metric[:, 0], metric[:, 1], "-o", markersize=2, linewidth=1)
    ax.set_title(f"Raw STL cross-section layer {layer_i} (r={current_r:.3f}, sample={slice_r:.3f})")
    ax.set_xlabel("arc length theta*r (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    debug_path = os.path.join(args.debug_dir, f"layer_{layer_i:04d}_raw_stl.png")
    fig.savefig(debug_path, dpi=200)
    if args.no_plot:
        plt.close(fig)
    print(f"\nSaved raw STL debug plot: {debug_path}")


def _plot_theta_z_loop(ax, coords, current_r, style="-o", color=None, label=None, alpha=1.0):
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 2:
        return
    metric = coords * np.array([current_r, 1])
    ax.plot(metric[:, 0], metric[:, 1], style, color=color, label=label,
            markersize=2, linewidth=1, alpha=alpha)


def _plot_theta_z_polygon(ax, poly, current_r, facecolor, edgecolor, label=None, alpha=0.35):
    if poly is None or poly.is_empty or not isinstance(poly, Polygon):
        return
    if not getattr(ax, "_theta_z_seam_lines_drawn", False):
        for theta in (0.0, 2 * np.pi):
            ax.axvline(theta * current_r, color="black", linewidth=0.9, alpha=0.85, zorder=5)
        ax._theta_z_seam_lines_drawn = True
    exterior = np.asarray(poly.exterior.coords, dtype=float) * np.array([current_r, 1])
    ax.scatter(exterior[:, 0], exterior[:, 1])
    ax.fill(exterior[:, 0], exterior[:, 1], facecolor=facecolor, edgecolor=edgecolor,
            alpha=alpha, linewidth=1.2, label=label)
    for interior in poly.interiors:
        hole = np.asarray(interior.coords, dtype=float) * np.array([current_r, 1])
        ax.fill(hole[:, 0], hole[:, 1], facecolor="white", edgecolor="#e63946",
                alpha=1.0, linewidth=1.0)


def _finish_debug_stage_plot(fig, ax, layer_i, current_r, filename_suffix, title):
    ax.set_title(f"{title} layer {layer_i} (r={current_r:.3f})")
    ax.set_xlabel("arc length theta*r (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_aspect("equal", adjustable="box")
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="best")
    fig.tight_layout()
    debug_path = os.path.join(args.debug_dir, f"layer_{layer_i:04d}_{filename_suffix}.png")
    fig.savefig(debug_path, dpi=200)
    if args.no_plot:
        plt.close(fig)
    print(f"\nSaved {title.lower()} debug plot: {debug_path}")


def save_classified_loops_layer_plot(layer_i, current_r, poly_pts, fc_original_indices, simple_fc_indices, complex_fc_loops):
    """Save loops after filtering and full-circle classification."""
    if layer_i not in debug_layers:
        return
    simple_set = set(simple_fc_indices)
    complex_ids = {id(loop) for loop in complex_fc_loops}
    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, poly in enumerate(poly_pts):
        if idx in simple_set:
            _plot_theta_z_loop(ax, poly, current_r, color="#ffb703", label="simple full-circle")
        elif id(poly) in complex_ids:
            _plot_theta_z_loop(ax, poly, current_r, color="#d00000", label="complex full-circle")
        elif idx in fc_original_indices:
            _plot_theta_z_loop(ax, poly, current_r, color="#9d4edd", label="other full-circle")
        else:
            _plot_theta_z_loop(ax, poly, current_r, color="#1d3557", label="ordinary loop")
    _finish_debug_stage_plot(fig, ax, layer_i, current_r, "classified_loops", "Classified STL loops")


def save_reconstructed_candidates_layer_plot(layer_i, current_r, base_polys, complex_parts, paired_parts):
    """Save polygon candidates after full-circle reconstruction/pairing, before hierarchy."""
    if layer_i not in debug_layers:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    for poly in base_polys:
        _plot_theta_z_polygon(ax, poly, current_r, "#a8dadc", "#1d3557", label="ordinary polygon", alpha=0.30)
    for poly in complex_parts:
        _plot_theta_z_polygon(ax, poly, current_r, "#ffadad", "#d00000", label="complex reconstruction", alpha=0.35)
        _plot_theta_z_polygon(ax, shapely.affinity.translate(poly,2*np.pi,0,0), current_r, "#ffadad", "#d00000", label="complex reconstruction", alpha=0.35)
        _plot_theta_z_polygon(ax, shapely.affinity.translate(poly,-2*np.pi,0,0), current_r, "#ffadad", "#d00000", label="complex reconstruction", alpha=0.35)
    for poly in paired_parts:
        _plot_theta_z_polygon(ax, poly, current_r, "#ffd166", "#bc6c25", label="simple full-circle pair", alpha=0.35)
    _finish_debug_stage_plot(fig, ax, layer_i, current_r, "reconstructed_candidates", "Reconstructed polygon candidates")


def save_hierarchy_layer_plot(layer_i, current_r, polys, hierarchy, levels):
    """Save polygon hierarchy levels before holes are punched."""
    if layer_i not in debug_layers:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ["#8ecae6", "#ffb703", "#90be6d", "#f28482", "#b8c0ff"]
    for idx, poly in enumerate(polys):
        if not isinstance(poly, Polygon) or poly.is_empty:
            continue
        level = levels[idx] if idx < len(levels) else 0
        _plot_theta_z_polygon(
            ax,
            poly,
            current_r,
            colors[level % len(colors)],
            "#1d3557" if level % 2 == 0 else "#e63946",
            label=f"level {level}",
            alpha=0.25 if level % 2 == 0 else 0.18,
        )
        centroid = poly.centroid
        ax.text(centroid.x * current_r, centroid.y, str(idx), fontsize=8, ha="center", va="center")
        for child_idx in hierarchy[idx]:
            if child_idx >= len(polys) or not isinstance(polys[child_idx], Polygon) or polys[child_idx].is_empty:
                continue
            child_centroid = polys[child_idx].centroid
            ax.plot(
                [centroid.x * current_r, child_centroid.x * current_r],
                [centroid.y, child_centroid.y],
                color="#6c757d",
                linewidth=0.7,
                alpha=0.6,
            )
    _finish_debug_stage_plot(fig, ax, layer_i, current_r, "hierarchy_levels", "Hierarchy levels before holes")


def save_paired_full_circle_traces_layer_plot(layer_i, current_r, pairs, leftovers=None, period=2 * np.pi):
    """Save normalized full-span boundary traces before polygon construction."""
    if layer_i not in debug_layers:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ["#d00000", "#1d3557", "#2a9d8f", "#bc6c25", "#6a4c93"]
    copy_styles = [
        (-1, "#8ecae6", 0.35, "-2pi copy"),
        (0, None, 0.9, "base"),
        (1, "#ffb703", 0.35, "+2pi copy"),
    ]

    def _plot_trace(trace, style, base_color, label_prefix, annotate=None, linewidth=1.2):
        coords = np.asarray(trace.coords, dtype=float)
        if len(coords) < 2:
            return
        for turn, copy_color, alpha, suffix in copy_styles:
            shifted = coords.copy()
            shifted[:, 0] += turn * period
            metric = shifted * np.array([current_r, 1])
            color = base_color if copy_color is None else copy_color
            ax.plot(
                metric[:, 0],
                metric[:, 1],
                style,
                color=color,
                linewidth=linewidth if turn == 0 else max(0.8, linewidth * 0.75),
                alpha=alpha,
                label=f"{label_prefix} base" if turn == 0 else suffix,
            )
            if turn == 0 and annotate:
                mid = metric[len(metric) // 2]
                ax.text(mid[0], mid[1], annotate, fontsize=7, ha="center", va="center")

    for pair_i, (lower, upper) in enumerate(pairs or []):
        color = colors[pair_i % len(colors)]
        for role, trace, style in (
            ("lower", lower, "-"),
            ("upper", upper, "--"),
        ):
            _plot_trace(
                trace,
                style,
                color,
                f"pair {pair_i} {role}",
                annotate=f"{pair_i}:{role[0]}",
            )
    for leftover_i, trace in enumerate(leftovers or []):
        _plot_trace(
            trace,
            ":",
            "#6c757d",
            "unpaired trace",
            annotate=f"u{leftover_i}",
            linewidth=1.0,
        )
    _finish_debug_stage_plot(
        fig,
        ax,
        layer_i,
        current_r,
        "complex_paired_traces",
        "Complex full-circle paired traces",
    )


def save_stitched_holes_layer_plot(layer_i, current_r, polygons):
    """Save polygons after exterior stitching and hole assignment."""
    if layer_i not in debug_layers or not polygons:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    for geom in polygons:
        geom = polygonal_only(geom)
        if geom is None or geom.is_empty:
            continue
        parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        for poly in parts:
            if not isinstance(poly, Polygon) or poly.is_empty:
                continue
            exterior = np.asarray(poly.exterior.coords, dtype=float) * np.array([current_r, 1])
            ax.fill(exterior[:, 0], exterior[:, 1], facecolor="#8ecae6", edgecolor="#1d3557", alpha=0.45, linewidth=1.2)
            for interior in poly.interiors:
                hole = np.asarray(interior.coords, dtype=float) * np.array([current_r, 1])
                ax.fill(hole[:, 0], hole[:, 1], facecolor="white", edgecolor="#e63946", alpha=1.0, linewidth=1.0)
    ax.set_title(f"Stitched polygons with holes layer {layer_i} (r={current_r:.3f})")
    ax.set_xlabel("arc length theta*r (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    debug_path = os.path.join(args.debug_dir, f"layer_{layer_i:04d}_stitched_holes.png")
    fig.savefig(debug_path, dpi=200)
    if args.no_plot:
        plt.close(fig)
    elif defer_single_debug_layer_plots:
        fig.canvas.draw_idle()
    else:
        plt.show()
    print(f"\nSaved stitched/holes debug plot: {debug_path}")


def extract_extrusion_segments_from_gcode(gcode, start_e, section_name):
    """Extract actual emitted extrusion moves from a generated G-code chunk."""
    segments = []
    last_xy = None
    current_e = float(start_e)
    for raw_line in gcode.splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        if not parts or parts[0] not in {"G0", "G1"}:
            continue

        values = {}
        for token in parts[1:]:
            key = token[:1]
            if key not in {"X", "Y", "E"} or len(token) < 2:
                continue
            try:
                values[key] = float(token[1:])
            except ValueError:
                continue

        if "X" in values and "Y" in values:
            new_xy = (math.radians(values["Y"]), values["X"])
            if (
                parts[0] == "G1"
                and "E" in values
                and last_xy is not None
                and values["E"] > current_e + 1e-9
            ):
                segments.append((section_name, np.asarray([last_xy, new_xy], dtype=float)))
            last_xy = new_xy
        if "E" in values:
            current_e = values["E"]
    return segments


def save_generated_lines_layer_plot(layer_i, current_r, emitted_segments, period=2 * np.pi):
    """Save the actual emitted G-code extrusion moves for a selected model layer."""
    if layer_i not in debug_layers:
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = {
        "EXTERNAL_PERIMETER": "#1d4ed8",
        "INTERNAL_PERIMETER": "#16a34a",
        "SOLID_INFILL": "#84cc16",
        "SUPPORT": "#d946ef",
        "INFILL_RADIAL_FULL": "#dc2626",
        "INFILL_RADIAL": "#f97316",
        "INFILL_AXIAL": "#7c3aed",
    }
    for section_name in colors:
        section_segments = []
        for seg_section, coords in emitted_segments:
            if seg_section != section_name:
                continue
            coords = np.asarray(coords, dtype=float).copy()
            if coords.shape != (2, 2):
                continue
            # G-code Y is accumulated. Fold only by whole turns for a layer-local
            # physical view while preserving each segment's actual emitted delta.
            turn_offset = math.floor(float(np.min(coords[:, 0])) / period) * period
            coords[:, 0] -= turn_offset
            section_segments.append(coords * np.array([current_r, 1.0]))
        if not section_segments:
            continue
        collection = LineCollection(
            section_segments,
            colors=colors[section_name],
            linewidths=0.75 if "INFILL" in section_name or section_name == "SUPPORT" else 1.2,
            alpha=0.85,
            label=section_name,
        )
        ax.add_collection(collection)

    if emitted_segments:
        ax.autoscale()
    else:
        ax.text(0.5, 0.5, "no extrusion moves", ha="center", va="center", transform=ax.transAxes)
    _finish_debug_stage_plot(
        fig,
        ax,
        layer_i,
        current_r,
        "generated_lines",
        "Generated G-code extrusion lines",
    )
    if args.no_plot:
        plt.close(fig)
    else:
        plt.show()

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


def periodic_copies(geom, period=2 * np.pi, turns=1):
    """Return geom plus integer-turn theta copies covering neighboring seams."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()
    geom = polygonal_only(shapely.make_valid(geom))
    if geom.is_empty:
        return shapely.GeometryCollection()
    copies = [
        shift_polygon(geom, k * period) if k else geom
        for k in range(-turns, turns + 1)
    ]
    try:
        return shapely.unary_union(copies)
    except shapely.errors.GEOSException:
        return shapely.GeometryCollection(copies)


def periodic_copies_for_bounds(geom, target_geom, period=2 * np.pi):
    """Return enough periodic copies of geom to cover target_geom's theta span."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()
    if target_geom is None or target_geom.is_empty:
        return polygonal_only(shapely.make_valid(geom))
    geom = polygonal_only(shapely.make_valid(geom))
    if geom.is_empty:
        return shapely.GeometryCollection()
    geom_min, _, geom_max, _ = geom.bounds
    target_min, _, target_max, _ = target_geom.bounds
    k_min = int(math.floor((target_min - geom_max) / period)) - 1
    k_max = int(math.ceil((target_max - geom_min) / period)) + 1
    copies = [shift_polygon(geom, k * period) if k else geom for k in range(k_min, k_max + 1)]
    return robust_polygonal_union(copies)

def stitch_periodic_shape(geom, period=2 * np.pi):
    """Union a geometry with its seam-neighbour copies for periodic containment tests."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()
    return robust_polygonal_union([
        shift_polygon(geom, -period),
        geom,
        shift_polygon(geom, period),
    ])


def build_hierarchy(polygons, period=2 * np.pi, periodic_indices=None):
    from shapely import STRtree
    periodic_indices = set(periodic_indices or ())
    hierarchy = [[] for _ in polygons]
    if len(polygons) < 2:
        return hierarchy
    tree = STRtree(polygons)
    for i, outer_poly in enumerate(polygons):
        if not isinstance(outer_poly, Polygon) or outer_poly.is_empty:
            continue
        seen = set()
        if i in periodic_indices:
            stitched_outer = stitch_periodic_shape(outer_poly, period=period)
            if stitched_outer.is_empty:
                continue
            for j in tree.query(stitched_outer):
                j = int(j)
                if j == i or j in seen:
                    continue
                child = polygons[j]
                if not isinstance(child, Polygon) or child.is_empty:
                    continue
                try:
                    contained = stitched_outer.contains(child)
                except shapely.errors.GEOSException:
                    contained = False
                if not contained:
                    for child_shift in (-period, period):
                        shifted_child = shift_polygon(child, child_shift)
                        try:
                            contained = stitched_outer.contains(shifted_child)
                        except shapely.errors.GEOSException:
                            contained = False
                        if contained:
                            break
                if contained:
                    seen.add(j)
                    hierarchy[i].append(j)
            continue
        for shift in [0, -period, period]:
            shifted = shift_polygon(outer_poly, shift) if shift != 0 else outer_poly
            for j in tree.query(shifted, predicate='contains'):
                j = int(j)
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


def canonicalize_full_circle_loop(polygon, period=2 * np.pi, seam_angle=0.0):
    """Cut a full-circle loop at a real period seam without changing phase.

    The only allowed theta translation here is an integer-period wrap.  Rotating
    the vertex order at a seam crossing is fine; subtracting the first theta or
    theta_min is not, because that physically shifts the perimeter relative to
    the STL slice linework.
    """
    coords = np.asarray(polygon, dtype=float)
    if len(coords) < 2:
        return coords.copy()

    coords = coords.copy()
    if len(coords) >= 2 and np.linalg.norm(coords[0] - coords[-1]) <= 1e-9:
        coords = coords[:-1]
    if len(coords) < 2:
        return coords.copy()

    for _ in range(2):
        shifted_theta = coords[:, 0] - seam_angle
        wrapped_diffs = np.diff(shifted_theta, append=shifted_theta[0])
        wrapped_diffs = (wrapped_diffs + np.pi) % period - np.pi
        if wrapped_diffs.sum() >= 0:
            break
        coords = coords[::-1].copy()

    shifted_theta = coords[:, 0] - seam_angle
    wrapped_diffs = np.diff(shifted_theta, append=shifted_theta[0])
    wrapped_diffs = (wrapped_diffs + np.pi) % period - np.pi

    theta = [float(coords[0, 0])]
    for step in wrapped_diffs[:-1]:
        theta.append(theta[-1] + float(step))
    theta = np.asarray(theta, dtype=float)
    total = float(np.sum(wrapped_diffs))
    if abs(total) < period * 0.5:
        return coords.copy()

    closed_theta = np.concatenate([theta, [theta[0] + total]])
    closed_z = np.concatenate([coords[:, 1], [coords[0, 1]]])
    if closed_theta[-1] < closed_theta[0]:
        closed_theta = closed_theta[::-1]
        closed_z = closed_z[::-1]

    seam = seam_angle + math.ceil((closed_theta[0] - seam_angle - 1e-9) / period) * period
    if seam > closed_theta[-1] + 1e-9:
        seam -= period

    seam_i = None
    seam_z = None
    for i in range(len(closed_theta) - 1):
        a = closed_theta[i]
        b = closed_theta[i + 1]
        if min(a, b) - 1e-9 <= seam <= max(a, b) + 1e-9:
            if abs(b - a) <= 1e-12:
                seam_z = float(closed_z[i])
            else:
                t = (seam - a) / (b - a)
                seam_z = float(closed_z[i] + t * (closed_z[i + 1] - closed_z[i]))
            seam_i = i
            break

    if seam_i is None or seam_z is None:
        # Fall back to an integer-period wrap only. This preserves phase, even if
        # the caller later rejects the trace as not spanning a usable seam window.
        period_shift = math.floor((float(np.min(theta)) + 0.25) / period) * period
        out = coords.copy()
        out[:, 0] = theta - period_shift
        return out

    out = [[seam, seam_z]]
    for j in range(seam_i + 1, len(closed_theta)):
        out.append([closed_theta[j], closed_z[j]])
    for j in range(1, seam_i + 1):
        out.append([closed_theta[j] + period, closed_z[j]])
    out.append([seam + period, seam_z])

    out = np.asarray(out, dtype=float)
    out[:, 0] -= seam
    step = np.abs(np.diff(out, axis=0))
    keep = np.concatenate([[True], np.any(step > 1e-9, axis=1)])
    return out[keep]


def is_simple_full_circle_loop(polygon, tol=1e-6, period=2 * np.pi):
    """Return True for full-turn loops that stay monotone in theta.

    Simple wrap bands (tube walls, annuli) can be paired by concatenating one
    full-turn loop with the reversed opposite loop. Once features create local
    theta reversals within a loop, that shortcut becomes topologically wrong and
    can stitch unrelated z-bands together.
    """
    coords = np.asarray(polygon, dtype=float)
    if len(coords) < 4 or not is_full_circle(coords):
        return False
    theta = coords[:, 0]
    diffs = np.diff(theta, append=theta[0])
    diffs = (diffs + np.pi) % period - np.pi
    signs = np.sign(diffs[np.abs(diffs) > tol])
    if len(signs) == 0:
        return True
    return np.all(signs > 0) or np.all(signs < 0)


def reconstruct_complex_full_circle_regions(
    loops,
    mesh,
    radius,
    nozzle_size,
    debug_layer=None,
    debug_plot_radius=None,
    period=2 * np.pi,
    z_connect_tol=0.05,
):
    """Pair full-span STL boundary traces into wrap-spanning polygons.

    Complex full-circle slices are still made from boundary lines that each span
    one complete turn.  Normalize those lines into a shared [0, 2π] theta window,
    pair adjacent traces by Z, and build one polygon per pair.  Mesh occupancy
    sampling is intentionally not used here; the STL slice linework is the source
    of truth.
    """
    if not loops:
        return []
    _ = mesh

    @dataclass
    class FullCircleTrace:
        """Normalized one-turn boundary trace used for z-order pairing."""

        coords: np.ndarray
        z_samples: np.ndarray
        z_median: float
        source_index: int
        run_index: int

    trace_span_tol = 0.05
    trace_endpoint_tol = 0.25
    trace_internal_slack = 0.35
    coord_tol = 1e-8
    sample_thetas = np.linspace(0.0, period, 49)[1:-1]

    z_values = [np.asarray(poly, dtype=float)[:, 1] for poly in loops if len(poly) >= 2]
    if not z_values:
        return []
    z_loop_min = min(float(np.min(zs)) for zs in z_values)
    z_loop_max = max(float(np.max(zs)) for zs in z_values)

    def _dedupe_trace_coords(coords):
        """Remove only duplicate points that are adjacent in trace order.

        Trace extraction repeatedly clips and merges linework, which can leave
        consecutive identical vertices. Those zero-length steps confuse span
        checks and LineString construction but do not carry geometry. Non-adjacent
        repeats, including a final point equal to the first, can represent real
        trace connectivity and must be preserved.
        """
        coords = np.asarray(coords, dtype=float)
        if len(coords) < 2:
            return coords
        step = np.abs(np.diff(coords, axis=0))
        keep = np.concatenate([[True], np.any(step > coord_tol, axis=1)])
        return coords[keep]

    def _z_values_from_intersection(geom):
        """Collect z ordinates from a Shapely point/line/collection result.

        Vertical theta probes may intersect a trace at a point, overlap an edge,
        or return a nested collection. Return all available y coordinates so the
        caller can choose a robust representative z value.
        """
        zs = []
        if geom is None or geom.is_empty:
            return zs
        if geom.geom_type == "Point":
            zs.append(float(geom.y))
        elif geom.geom_type == "LineString":
            zs.extend(float(p[1]) for p in geom.coords)
        elif hasattr(geom, "geoms"):
            for part in geom.geoms:
                zs.extend(_z_values_from_intersection(part))
        return zs

    def _sample_trace_z(coords):
        """Sample a trace's z value at fixed interior theta positions.

        The samples provide a coarse shape profile for deciding whether two
        traces cross before they are paired. Missing/intersection-error samples
        are stored as NaN so callers can compare only reliable positions.
        """
        line = LineString(coords)
        z_pad = max(1.0, abs(z_loop_max - z_loop_min) + max(0.25, nozzle_size))
        sample_line_min = z_loop_min - z_pad
        sample_line_max = z_loop_max + z_pad
        samples = []
        for theta in sample_thetas:
            try:
                inter = line.intersection(LineString([(theta, sample_line_min), (theta, sample_line_max)]))
                zs = _z_values_from_intersection(inter)
            except shapely.errors.GEOSException:
                zs = []
            samples.append(float(np.median(zs)) if zs else np.nan)
        return np.asarray(samples, dtype=float)

    def _make_trace(coords, source_index, run_index):
        """Validate and normalize a candidate one-turn boundary trace.

        A valid trace must span exactly one angular period by endpoint or theta
        range. The trace is oriented toward increasing theta, shifted only by
        whole periods to preserve the STL slice phase, seam endpoints are snapped
        to 0/period, and implausible internal seam overshoots are rejected.
        """
        coords = _dedupe_trace_coords(coords)
        if len(coords) < 2:
            return None
        if coords[-1, 0] < coords[0, 0]:
            coords = coords[::-1].copy()

        theta_min = float(np.min(coords[:, 0]))
        theta_max = float(np.max(coords[:, 0]))
        theta_span = theta_max - theta_min
        endpoint_span = float(coords[-1, 0] - coords[0, 0])
        endpoint_is_full = abs(endpoint_span - period) <= trace_endpoint_tol
        range_is_full = abs(theta_span - period) <= trace_span_tol
        if not endpoint_is_full and not range_is_full:
            return None

        coords = coords.copy()
        period_shift = math.floor((float(np.min(coords[:, 0])) + trace_internal_slack) / period) * period
        if period_shift:
            coords[:, 0] -= period_shift
        # Preserve the slice phase.  Only whole-period wrapping is allowed here.
        # Do not rebase by coords[0] or theta_min; that creates a real angular
        # shift relative to the STL slice linework.
        coords[:, 0] = np.where(np.abs(coords[:, 0]) <= trace_span_tol, 0.0, coords[:, 0])
        coords[:, 0] = np.where(np.abs(coords[:, 0] - period) <= trace_span_tol, period, coords[:, 0])

        # Coarse, unsubdivided STL triangles can create small local theta
        # overshoots around the seam even when the boundary walk itself spans
        # exactly one turn. Keep those traces; reject only implausibly large
        # excursions that would make the pair polygon ambiguous.
        if (
            float(np.min(coords[:, 0])) < -trace_internal_slack
            or float(np.max(coords[:, 0])) > period + trace_internal_slack
        ):
            return None

        z_samples = _sample_trace_z(coords)
        finite = z_samples[np.isfinite(z_samples)]
        z_median = float(np.median(finite)) if len(finite) else float(np.median(coords[:, 1]))
        return FullCircleTrace(coords, z_samples, z_median, source_index, run_index)

    def _trace_signature(trace):
        """Build an orientation-independent rounded coordinate key for a trace."""
        coords = np.asarray(trace.coords, dtype=float)
        if len(coords) < 2:
            return ()
        rounded = tuple(map(tuple, np.round(coords, 7)))
        reversed_rounded = tuple(map(tuple, np.round(coords[::-1], 7)))
        return min(rounded, reversed_rounded)

    def _dedupe_traces(traces):
        """Keep the first instance of each geometric trace signature."""
        unique = []
        seen = set()
        for trace in traces:
            sig = _trace_signature(trace)
            if not sig or sig in seen:
                continue
            seen.add(sig)
            unique.append(trace)
        return unique

    def _windowed_full_span_traces(coords, source_index):
        """Extract full-period traces by clipping periodic loop copies.

        Complex loops can cross the theta seam or contain local theta reversals,
        so the raw coordinate order is not always a valid one-turn trace. Shift
        enough whole-period copies to cover the canonical [0, period] window,
        clip them to that window, merge compatible fragments, then validate each
        fragment directly or after splitting it into monotone theta runs.
        """
        coords = np.asarray(coords, dtype=float)
        if len(coords) < 2:
            return []
        if not np.allclose(coords[0], coords[-1], atol=coord_tol):
            coords = np.vstack([coords, coords[0]])

        theta_min = float(np.min(coords[:, 0]))
        theta_max = float(np.max(coords[:, 0]))
        z_pad = max(0.25, nozzle_size)
        clip = box(0.0, z_loop_min - z_pad, period, z_loop_max + z_pad)
        k_min = int(math.floor((0.0 - theta_max) / period)) - 1
        k_max = int(math.ceil((period - theta_min) / period)) + 1
        fragments = []
        for k in range(k_min, k_max + 1):
            shifted = coords.copy()
            shifted[:, 0] += k * period
            try:
                clipped = LineString(shifted).intersection(clip)
            except (ValueError, shapely.errors.GEOSException):
                continue
            _collect_lines(clipped, fragments, min_length=1e-9)

        traces = []
        candidate_fragments = []
        if fragments:
            fragments = _dedupe_lines(fragments, min_length=1e-9)
            try:
                merged = linemerge(MultiLineString(fragments))
                _collect_lines(merged, candidate_fragments, min_length=1e-9)
            except (ValueError, shapely.errors.GEOSException):
                candidate_fragments = []
        if not candidate_fragments:
            candidate_fragments = fragments

        for frag_i, fragment in enumerate(candidate_fragments):
            frag_coords = np.asarray(fragment.coords, dtype=float)
            if len(frag_coords) < 2:
                continue
            # Clipped full-wrap boundaries may come back as one long fragment, or
            # as monotone pieces when the raw loop locally reverses in theta.
            trace = _make_trace(frag_coords, source_index, frag_i)
            if trace is not None:
                traces.append(trace)
                continue
            traces.extend(_full_span_monotone_runs(frag_coords, source_index))
        return _dedupe_traces(traces)

    def _full_span_monotone_runs(coords, source_index):
        """Split a line into theta-monotone runs and keep full-span runs only."""
        runs = []
        try:
            split_runs = _split_theta_monotone_runs(LineString(coords), tol=coord_tol)
        except (ValueError, shapely.errors.GEOSException):
            split_runs = []
        for run_i, run in enumerate(split_runs):
            run_coords = np.asarray(run.coords, dtype=float)
            if len(run_coords) < 2:
                continue
            trace = _make_trace(run_coords, source_index, run_i)
            if trace is not None:
                runs.append(trace)
        return runs

    def _extract_traces_from_loop(loop, source_index):
        """Return all reconstructable full-period boundary traces from one loop.

        Extraction proceeds from least to most invasive: try the loop as-is,
        then clipped periodic copies, then a canonicalized loop and its monotone
        runs. This preserves straightforward traces while still recovering the
        separable full-circle boundaries inside more tangled STL slice loops.
        """
        coords = np.asarray(loop, dtype=float)
        if len(coords) < 4 or not is_full_circle(coords):
            return []

        whole_trace = _make_trace(coords, source_index, 0)
        if whole_trace is not None:
            return [whole_trace]

        windowed_traces = _windowed_full_span_traces(coords, source_index)
        if windowed_traces:
            return windowed_traces

        try:
            canonical = canonicalize_full_circle_loop(coords, period=period, seam_angle=0.0)
        except (ValueError, FloatingPointError):
            return []
        canonical = _dedupe_trace_coords(canonical)
        if len(canonical) < 2:
            return []

        split_traces = _full_span_monotone_runs(canonical, source_index)
        if split_traces:
            return split_traces

        canonical_trace = _make_trace(canonical, source_index, 0)
        return [canonical_trace] if canonical_trace is not None else []

    traces = []
    for source_index, loop in enumerate(loops):
        traces.extend(_extract_traces_from_loop(loop, source_index))
    traces = _dedupe_traces(traces)
    if len(traces) < 2:
        save_paired_full_circle_traces_layer_plot(
            debug_layer,
            debug_plot_radius if debug_plot_radius is not None else radius,
            [],
            traces,
        )
        if debug_layer in debug_layers:
            print(f"\nOnly found {len(traces)} full-circle trace(s) on layer {debug_layer}; skipping reconstruction.")
        return []

    def _trace_pair_crosses(lower, upper):
        """Return True when two candidate pair traces swap z order.

        Adjacent traces are paired by median z, but locally crossing profiles
        would create a self-intersecting or topologically wrong polygon. Compare
        only theta samples present on both traces and reject pairs whose lower
        decile violates the expected z ordering beyond the connection tolerance.
        """
        finite = np.isfinite(lower.z_samples) & np.isfinite(upper.z_samples)
        if np.count_nonzero(finite) < 3:
            return False
        dz = upper.z_samples[finite] - lower.z_samples[finite]
        if float(np.nanmedian(dz)) < 0.0:
            dz = -dz
        return float(np.nanpercentile(dz, 10)) < -max(z_connect_tol, 1e-4)

    sorted_traces = sorted(traces, key=lambda trace: trace.z_median)
    trace_pairs = []
    leftovers = []
    for pair_i in range(0, len(sorted_traces) - 1, 2):
        lower = sorted_traces[pair_i]
        upper = sorted_traces[pair_i + 1]
        if upper.z_median < lower.z_median:
            lower, upper = upper, lower
        if _trace_pair_crosses(lower, upper):
            leftovers.extend([lower, upper])
            continue
        trace_pairs.append((lower, upper))
    if len(sorted_traces) % 2:
        leftovers.append(sorted_traces[-1])

    save_paired_full_circle_traces_layer_plot(
        debug_layer,
        debug_plot_radius if debug_plot_radius is not None else radius,
        trace_pairs,
        leftovers,
    )
    if leftovers and debug_layer in debug_layers:
        print(f"\nSkipped {len(leftovers)} unpaired/crossing full-circle trace(s) on layer {debug_layer}.")
    if not trace_pairs:
        return []

    reconstructed = []
    for lower, upper in trace_pairs:

        ring = np.vstack([lower.coords+np.array([-2*np.pi,0]),lower.coords,lower.coords+np.array([2*np.pi,0]),(upper.coords+np.array([2*np.pi,0]))[::-1],upper.coords[::-1],(upper.coords+np.array([-2*np.pi,0]))[::-1]])
        try:
            geom = polygonal_only(shapely.make_valid(Polygon(ring)))
            geom = polygonal_only(
                fold_to_theta_window(
                    geom,
                    0.0,
                    period,
                    z_loop_min,
                    z_loop_max,
                    period=period,
                )
            )
        except (ValueError, shapely.errors.GEOSException):
            continue

        parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
        for part in parts:
            if isinstance(part, Polygon) and not part.is_empty and part.area > 1e-6:
                theta_width = float(part.bounds[2] - part.bounds[0])
                if theta_width < period * 0.75:
                    continue
                reconstructed.append(part)
    return reconstructed


def promote_periodic_reconstruction_holes(
    polygons,
    period=2 * np.pi,
    min_area=1e-6,
    layer_i=None,
    current_r=None,
):
    """Split seam-discovered holes out of reconstructed full-wrap polygons.

    A hole whose side lands exactly on the theta cut may only become a closed
    interior ring when neighbouring period copies are stitched together. Promote
    those rings into normal polygon candidates so hierarchy/hole subtraction sees
    them like any other STL loop.
    """
    if not polygons:
        return [], []

    source_parts = []
    shell_parts = []
    for poly in polygons:
        if not isinstance(poly, Polygon) or poly.is_empty:
            continue
        repaired = repair_polygonal_geometry(poly)
        parts = list(repaired.geoms) if hasattr(repaired, "geoms") else [repaired]
        for part in parts:
            if not isinstance(part, Polygon) or part.is_empty or part.area <= min_area:
                continue
            source_parts.append(part)
            try:
                shell = repair_polygonal_geometry(Polygon(part.exterior.coords))
            except (ValueError, shapely.errors.GEOSException):
                continue
            shell_geoms = list(shell.geoms) if hasattr(shell, "geoms") else [shell]
            for shell_part in shell_geoms:
                if isinstance(shell_part, Polygon) and not shell_part.is_empty and shell_part.area > min_area:
                    shell_parts.append(shell_part)

    if not source_parts:
        return [], []

    stitched = robust_polygonal_union(
        [
            shift_polygon(part, k * period) if k else part
            for part in source_parts
            for k in (0, 1)
        ]
    )

    def _polygon_parts(geom):
        geom = polygonal_only(geom)
        if geom.is_empty:
            return []
        return list(geom.geoms) if hasattr(geom, "geoms") else [geom]

    def _canonical_centroid_copy(poly):
        centroid_x = poly.centroid.x
        k = math.floor(-centroid_x / period)
        while centroid_x + k * period < 0.0:
            k += 1
        while centroid_x + k * period >= period:
            k -= 1
        return shift_polygon(poly, k * period) if k else poly

    def _periodic_signature(poly):
        coords = np.asarray(poly.exterior.coords[:-1], dtype=float)
        if len(coords) == 0:
            return None
        theta = np.mod(coords[:, 0], period)
        theta = np.where(np.isclose(theta, period, atol=1e-8), 0.0, theta)
        points = sorted((round(float(x), 7), round(float(z), 7)) for x, z in zip(theta, coords[:, 1]))
        return round(float(poly.area), 7), tuple(points)

    holes = []
    seen = set()
    for stitched_part in _polygon_parts(stitched):
        if not isinstance(stitched_part, Polygon) or stitched_part.is_empty:
            continue
        for interior in stitched_part.interiors:
            try:
                hole_geom = repair_polygonal_geometry(Polygon(interior.coords))
            except (ValueError, shapely.errors.GEOSException):
                continue
            for hole_part in _polygon_parts(hole_geom):
                if not isinstance(hole_part, Polygon) or hole_part.is_empty or hole_part.area <= min_area:
                    continue
                hole_part = _canonical_centroid_copy(hole_part)
                sig = _periodic_signature(hole_part)
                if sig is None or sig in seen:
                    continue
                seen.add(sig)
                holes.append(hole_part)

    promoted_shells = shell_parts or source_parts
    if holes and False:
        plot_radius = current_r if current_r is not None else 1.0
        fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
        title_layer = f" layer {layer_i}" if layer_i is not None else ""
        fig.suptitle(
            f"Periodic reconstruction hole promotion{title_layer}: "
            f"{len(holes)} promoted"
        )

        def _draw_geom(ax, geom, facecolor, edgecolor, label, alpha=0.25, copies=(0,)):
            parts = _polygon_parts(geom)
            for part in parts:
                if not isinstance(part, Polygon) or part.is_empty:
                    continue
                for copy_i, k in enumerate(copies):
                    draw_part = shift_polygon(part, k * period) if k else part
                    draw_label = label if copy_i == 0 else None
                    _plot_theta_z_polygon(
                        ax,
                        draw_part,
                        plot_radius,
                        facecolor,
                        edgecolor,
                        label=draw_label,
                        alpha=alpha,
                    )

        for part in source_parts:
            _draw_geom(axs[0], part, "#a8dadc", "#1d3557", "input", alpha=0.25, copies=(-1, 0, 1))
        axs[0].set_title("input + period copies")

        _draw_geom(axs[1], stitched, "#ffd166", "#bc6c25", "stitched union", alpha=0.30)
        for hole in holes:
            _draw_geom(axs[1], hole, "#ef476f", "#d00000", "promoted hole", alpha=0.60)
        axs[1].set_title("stitched union holes")

        for shell in promoted_shells:
            _draw_geom(axs[2], shell, "#90be6d", "#2d6a4f", "output shell", alpha=0.25)
        for hole in holes:
            _draw_geom(axs[2], hole, "#ef476f", "#d00000", "added hole candidate", alpha=0.65)
        axs[2].set_title("returned to normal hole pipeline")

        for ax in axs:
            ax.axvline(0.0, color="black", linewidth=0.8, alpha=0.6)
            ax.axvline(period * plot_radius, color="black", linewidth=0.8, alpha=0.6)
            ax.set_xlabel("arc length theta*r (mm)")
            ax.set_aspect("equal", adjustable="box")
            handles, labels = ax.get_legend_handles_labels()
            if labels:
                by_label = dict(zip(labels, handles))
                ax.legend(by_label.values(), by_label.keys(), loc="best")
        axs[0].set_ylabel("z (mm)")
        fig.tight_layout()
        for hole in holes:
            fig, ax = plt.subplots(1, 1, figsize=(18, 5), sharey=True)
            _plot_theta_z_polygon(ax,poly,current_r,"#00ff00","#ff00ff")
        plt.show()

    return promoted_shells, holes


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


def _dedupe_lines(lines, min_length=0.0):
    """Remove duplicate linework without polygon-unioning seam-sensitive regions."""
    lines = [line for line in lines if isinstance(line, LineString) and not line.is_empty]
    if len(lines) <= 1:
        return lines
    try:
        unioned = shapely.unary_union(MultiLineString(lines))
    except (ValueError, shapely.errors.GEOSException):
        return lines
    result = []
    _collect_lines(unioned, result, min_length=min_length)
    return result or lines


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
                part = robust_polygonal_buffer(part, r_target - r_current)
            inflated.append(part)
        if not inflated:
            return shapely.GeometryCollection()
        metric = shapely.unary_union(inflated)

    if merge_gap_mm > 0:
        metric = robust_polygonal_buffer(metric, merge_gap_mm / 2)
        metric = robust_polygonal_buffer(metric, -(merge_gap_mm / 2))

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
    folded = robust_polygonal_union([
        geom,
        shift_polygon(geom, period),
        shift_polygon(geom, -period),
    ])
    return robust_polygonal_intersection(folded, clip_box)


def fold_to_theta_window(geom, lower, upper, z_min, z_max, period=2 * np.pi):
    """Fold periodic (θ,z) geometry into an arbitrary one-turn theta window."""
    if geom is None or geom.is_empty:
        return geom
    geom = shapely.make_valid(geom)
    geom = shapely.set_precision(geom, 1e-8)
    minx, _, maxx, _ = geom.bounds
    k_min = int(math.floor((lower - maxx) / period)) - 1
    k_max = int(math.ceil((upper - minx) / period)) + 1
    copies = [shift_polygon(geom, k * period) if k else geom for k in range(k_min, k_max + 1)]
    folded = robust_polygonal_union(copies)
    return robust_polygonal_intersection(
        shapely.set_precision(folded, 1e-8),
        box(lower, z_min, upper, z_max),
    )


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
    theta_min, z_min, theta_max, z_max = bounds
    multi_turn = (theta_max - theta_min) > (2 * np.pi + 1e-6)

    # Anchor θ positions to r_base so lines stay at the same physical angles on every layer,
    # matching the strategy used by do_sparse_infill.
    r_base = slicing_params["r_base"]

    def _theta_grid(spacing):
        step = spacing / r_base
        origin = step / 2
        if multi_turn:
            k_start = int(np.ceil((theta_min - origin) / step))
            k_end = int(np.floor((theta_max - origin) / step))
            return origin + np.arange(k_start, k_end + 1) * step
        n = int(np.floor(2 * np.pi / step))
        return np.arange(n) * step + origin

    # Support strips are ~layer_height thick in z, so their infill intersections are
    # naturally short.  Only filter true zero-length degenerate geometry.
    min_length = 1e-6

    if direction in ("radial", "body_radial"):
        spacing = slicing_params["support_interface_spacing"] if direction == "radial" else slicing_params["support_body_spacing"]
        # Radial lines only (constant z, varying θ) — stitch seam-split arcs
        z_lo = np.floor(z_min / spacing) * spacing
        for z in np.arange(z_lo, z_max + spacing, spacing):
            ls = shapely.LineString([[theta_min, z], [theta_max, z]])
            inter = shapely.intersection(ls, support_geom)
            if hasattr(inter, 'geoms'):
                stitched = list(inter.geoms) if multi_turn else stitch_radial_seam(inter.geoms)
                for seg in stitched:
                    seg = seg if multi_turn else _snap_radial_to_full_circle(seg)
                    _collect_lines(seg, result, min_length)
            else:
                inter = inter if multi_turn else _snap_radial_to_full_circle(inter)
                _collect_lines(inter, result, min_length)

    elif direction in ("axial", "body_axial"):
        spacing = slicing_params["support_interface_spacing"] if direction == "axial" else slicing_params["support_body_spacing"]
        thetas = _theta_grid(spacing)
        # Axial lines only (constant θ, varying z)
        for theta in thetas:
            ls = shapely.LineString([[theta, z_min - spacing], [theta, z_max + spacing]])
            _collect_lines(shapely.intersection(ls, support_geom), result, min_length)

    else:  # "body"
        spacing = slicing_params["support_body_spacing"]
        thetas = _theta_grid(spacing)

        mode = slicing_params.get("support_infill_mode", "grid")

        if mode == "sine":
            amplitude_rad = slicing_params["support_sine_amplitude_mm"] / current_r
            period_mm = slicing_params["support_sine_period_mm"]
            # Extend z slightly so waves fully cover the geometry bounds
            z_lo = z_min - period_mm
            z_hi = z_max + period_mm
            for theta in thetas:
                if multi_turn:
                    shifts = [0.0]
                else:
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
                ls = shapely.LineString([[theta_min, z], [theta_max, z]])
                inter = shapely.intersection(ls, support_geom)
                if hasattr(inter, 'geoms'):
                    stitched = list(inter.geoms) if multi_turn else stitch_radial_seam(inter.geoms)
                    for seg in stitched:
                        seg = seg if multi_turn else _snap_radial_to_full_circle(seg)
                        _collect_lines(seg, result, min_length)
                else:
                    inter = inter if multi_turn else _snap_radial_to_full_circle(inter)
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
    theta_min, z_min, theta_max, z_max = bounds
    multi_turn = (theta_max - theta_min) > (2 * np.pi + 1e-6)

    def _collect(inter):
        if inter is None or inter.is_empty:
            return
        if isinstance(inter, LineString):
            if do_transform(inter).length <= 1e-6:
                return
            if plot_infill:
                plot_line(do_transform(inter), ax=ax, color="greenyellow", add_points=False)
            infill.append(inter)
        elif isinstance(inter, MultiLineString) or type(inter) == shapely.geometry.collection.GeometryCollection:
            for seg in inter.geoms:
                _collect(seg)
        elif type(inter) == shapely.geometry.multipoint.MultiPoint:
            return

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
from shapely.ops import linemerge, polygonize

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
    raw_coords = np.asarray(ring.coords, dtype=float)
    if len(raw_coords) < 2:
        return []

    period = upper - lower
    seam_slack = math.radians(5.0)
    raw_dtheta = np.diff(raw_coords[:, 0])
    raw_dz = np.diff(raw_coords[:, 1])
    raw_theta_span = raw_coords[:, 0].max() - raw_coords[:, 0].min()
    orthogonal_ring = (
        ring.is_ring
        and np.all((np.abs(raw_dtheta) <= tol) | (np.abs(raw_dz) <= tol))
    )

    if orthogonal_ring and raw_theta_span >= period - tol:
        z_levels = []
        for z in raw_coords[:, 1]:
            if not any(abs(z - existing) <= 1e-6 for existing in z_levels):
                z_levels.append(float(z))
        z_levels.sort()
        if 1 <= len(z_levels) <= 2:
            return [
                LineString([[lower, z], [upper, z]])
                for z in z_levels
            ]

    # Closed orthogonal rings (tube / band slices) can contain genuine full-width
    # constant-z edges. `_restore_theta_walk` would reinterpret those ±2π/±4π edges
    # as seam wraps and collapse them to near-zero span, dropping one side wall.
    # Preserve their raw walk; helical/diagonal rings still need theta restoration.
    if orthogonal_ring and raw_theta_span > period + tol:
        coords = raw_coords
    else:
        coords = _restore_theta_walk(raw_coords)
    theta_span = coords[:, 0].max() - coords[:, 0].min()
    if theta_span < 1e-3:
        return []

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

    unioned = shapely.unary_union(MultiLineString(pieces))
    if isinstance(unioned, LineString):
        merged = unioned
    else:
        merged_parts = []
        _collect_lines(unioned, merged_parts, min_length=0.0)
        if not merged_parts:
            return []
        merged = merged_parts[0] if len(merged_parts) == 1 else linemerge(MultiLineString(merged_parts))
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


def _canonical_full_wrap_strip_polygon(poly, lower=0.0, upper=2 * math.pi, seam_tol=0.1, fill_tol=0.995, ortho_tol=1e-6):
    """Snap a near-rectangular full-wrap band onto the exact canonical turn."""
    if not isinstance(poly, Polygon) or poly.is_empty or poly.interiors:
        return None

    minx, miny, maxx, maxy = poly.bounds
    period = upper - lower
    theta_span = maxx - minx
    z_span = maxy - miny
    if z_span <= 0:
        return None
    if minx > lower + seam_tol or maxx < upper - seam_tol:
        return None

    coords = np.asarray(poly.exterior.coords, dtype=float)
    if len(coords) < 4 or not poly.exterior.is_ring:
        return None
    dtheta = np.diff(coords[:, 0])
    dz = np.diff(coords[:, 1])
    if not np.all((np.abs(dtheta) <= ortho_tol) | (np.abs(dz) <= ortho_tol)):
        return None
    z_levels = []
    for z in coords[:, 1]:
        if not any(abs(z - existing) <= 1e-6 for existing in z_levels):
            z_levels.append(float(z))
    if len(z_levels) != 2:
        return None

    bbox_area = theta_span * z_span
    if bbox_area <= 0 or poly.area / bbox_area < fill_tol:
        return None

    return box(lower, miny, upper, maxy)


def _canonicalize_full_wrap_strips(geom, lower=0.0, upper=2 * math.pi):
    """Recursively snap full-wrap orthogonal strips to exact canonical bounds."""
    if geom is None or geom.is_empty:
        return geom
    if isinstance(geom, Polygon):
        snapped = _canonical_full_wrap_strip_polygon(geom, lower=lower, upper=upper)
        return snapped if snapped is not None else geom
    if isinstance(geom, shapely.MultiPolygon):
        parts = [_canonicalize_full_wrap_strips(g, lower=lower, upper=upper) for g in geom.geoms]
        parts = [g for g in parts if g is not None and not g.is_empty]
        if not parts:
            return shapely.GeometryCollection()
        return shapely.unary_union(parts)
    if isinstance(geom, shapely.GeometryCollection):
        parts = [_canonicalize_full_wrap_strips(g, lower=lower, upper=upper) for g in geom.geoms]
        parts = [g for g in parts if g is not None and not g.is_empty]
        if not parts:
            return shapely.GeometryCollection()
        return shapely.unary_union(parts)
    return geom


def _is_effectively_full_wrap_polygon(poly, lower=0.0, upper=2 * math.pi, seam_tol=0.15):
    """Return True when a polygon is one physical wrap with only a small seam gap.

    Full-wrap regions often come out of slice reconstruction or buffering with a
    narrow gap at the canonical seam instead of exact bounds at 0 and 2π. Those
    polygons still need periodic neighbour copies during inward offsets so we
    do not erode against synthetic seam endcaps.
    """
    if not isinstance(poly, Polygon) or poly.is_empty:
        return False

    period = upper - lower
    minx, _, maxx, _ = poly.bounds
    theta_span = maxx - minx
    if theta_span < period - seam_tol or theta_span > period + seam_tol:
        return False

    window = max(seam_tol, 1e-6)
    left = poly.intersection(box(lower, -1e15, lower + window, 1e15))
    right = poly.intersection(box(upper - window, -1e15, upper, 1e15))
    return (
        not left.is_empty
        and left.area > 1e-8
        and not right.is_empty
        and right.area > 1e-8
    )


def has_full_wrap_component(geom, lower=0.0, upper=2 * math.pi, seam_tol=0.15):
    """Return True if any polygon component physically spans the periodic seam."""
    geom = polygonal_only(geom)
    if geom.is_empty:
        return False
    parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
    period = upper - lower
    for part in parts:
        if not isinstance(part, Polygon) or part.is_empty:
            continue
        minx, _, maxx, _ = part.bounds
        if maxx - minx >= period - seam_tol and _is_effectively_full_wrap_polygon(part, lower=lower, upper=upper, seam_tol=seam_tol):
            return True
    return False


def _full_wrap_strip_lines(poly, lower=0.0, upper=2 * math.pi, seam_tol=0.1, fill_tol=0.995):
    """Return the two physical side-wall lines for a near-rectangular full-wrap band."""
    poly = _canonical_full_wrap_strip_polygon(poly, lower=lower, upper=upper, seam_tol=seam_tol, fill_tol=fill_tol)
    if poly is None:
        return None

    return [
        LineString([[lower, poly.bounds[1]], [upper, poly.bounds[1]]]),
        LineString([[lower, poly.bounds[3]], [upper, poly.bounds[3]]]),
    ]


def _periodic_boundary_lines(poly, lower=0.0, upper=2 * math.pi, tol=1e-9, seam_tol=0.15, cut_lower=None):
    """Extract physical boundary lines for a polygon spanning the periodic seam.

    A full-wrap polygon folded into [0, 2π] can have a planar exterior dominated
    by synthetic seam cuts. Build the ±one-turn periodic union first so seam cuts
    become interior edges, then clip the physical boundary back to one turn.
    """
    if not isinstance(poly, Polygon) or poly.is_empty:
        return None

    period = upper - lower
    minx, _, maxx, _ = poly.bounds
    theta_span = maxx - minx
    if (
        theta_span < period - seam_tol
        or theta_span > period + seam_tol
        or minx > lower + seam_tol
        or maxx < upper - seam_tol
    ):
        return None
    left_gap = max(0.0, minx - lower)
    right_gap = max(0.0, upper - maxx)
    snap_tol = max(left_gap, right_gap) + 1e-6
    if snap_tol > 0.025:
        return None

    def _snap_theta_to_seam(coords):
        coords = np.asarray(coords, dtype=float).copy()
        coords[:, 0] = np.where(coords[:, 0] <= lower + snap_tol, lower, coords[:, 0])
        coords[:, 0] = np.where(coords[:, 0] >= upper - snap_tol, upper, coords[:, 0])
        return coords

    poly = polygonal_only(shapely.make_valid(shapely.transform(poly, _snap_theta_to_seam)))
    if not isinstance(poly, Polygon) or poly.is_empty:
        return None

    try:
        periodic_geom = robust_polygonal_union([
            shift_polygon(poly, -period),
            poly,
            shift_polygon(poly, period),
        ])
        # Cut the printable boundary away from local features when the caller is
        # using the canonical window. For already-shifted print windows, the
        # caller can request the window start directly.
        if cut_lower is None:
            cut_lower = lower + period * 0.5
        cut_upper = cut_lower + period
        boundary = periodic_geom.boundary.intersection(box(cut_lower, -1e15, cut_upper, 1e15))
        raw_lines = []
        _collect_lines(boundary, raw_lines, min_length=1e-9)
        if not raw_lines:
            return None
        unioned = shapely.unary_union(MultiLineString(raw_lines))
    except shapely.errors.GEOSException:
        return None

    merged_lines = []
    _collect_lines(unioned, merged_lines, min_length=1e-9)
    runs = []
    for line in merged_lines:
        coords = np.asarray(line.coords, dtype=float)
        if len(coords) < 2:
            continue
        theta_span = coords[:, 0].max() - coords[:, 0].min()
        if line.is_ring and theta_span < period - tol:
            runs.append(line)
        else:
            runs.extend(_split_theta_monotone_runs(line, tol=tol))
    return [line for line in runs if not line.is_empty and line.length >= 1e-9]


def polygon_to_wrapped_lines(poly, lower=0.0, upper=2*math.pi, cut_lower=None):
    """
    Convert a polygon (or collection) to LineStrings for gcode/plotting.

    Polygonal sections may include ±2π seam-neighbour copies, or a single
    canonical polygon whose centroid falls slightly outside [0, 2π) because it
    straddles the seam.  Normalize polygonal input onto the canonical turn with
    overlap-based shifting before extracting boundary lines. Rings that fit
    inside one canonical turn are returned unchanged; rings that cross the seam
    or extend past one turn are clipped as boundary lines, shifted back into
    unwrapped theta space, and split only where their theta direction actually
    changes.
    """
    if poly is None or (hasattr(poly, 'is_empty') and poly.is_empty):
        return []

    # Recurse into collections / lists
    if isinstance(poly, list):
        ret = []
        for each in poly:
            ret += polygon_to_wrapped_lines(each, lower, upper, cut_lower=cut_lower)
        return _dedupe_lines(ret)

    if isinstance(poly, (LineString, shapely.MultiLineString)):
        # Already a line (e.g. infill passed through here) — flatten and return as-is
        ret = []
        _collect_lines(poly, ret, min_length=0.0)
        return ret

    if isinstance(poly, (shapely.MultiPolygon, shapely.GeometryCollection)):
        ret = []
        for geom in poly.geoms:
            ret += polygon_to_wrapped_lines(geom, lower, upper, cut_lower=cut_lower)
        return _dedupe_lines(ret)

    if isinstance(poly, Polygon):
        poly = keep_canonical_polygons(poly, lower=lower, upper=upper)
        if poly is None or poly.is_empty:
            return []
        if isinstance(poly, (shapely.MultiPolygon, shapely.GeometryCollection)):
            ret = []
            for geom in poly.geoms:
                ret += polygon_to_wrapped_lines(geom, lower, upper, cut_lower=cut_lower)
            return _dedupe_lines(ret)

    if not isinstance(poly, Polygon):
        return []

    periodic_lines = _periodic_boundary_lines(poly, lower=lower, upper=upper, cut_lower=cut_lower)
    if periodic_lines:
        return periodic_lines

    strip_lines = _full_wrap_strip_lines(poly, lower=lower, upper=upper)
    if strip_lines is not None:
        return strip_lines

    result = []
    for ring in [poly.exterior] + list(poly.interiors):
        result.extend(_ring_to_wrapped_lines(ring, lower=lower, upper=upper))
    return result


def _remove_periodic_seam_caps(
    lines,
    period=2 * math.pi,
    seam_origin=0.0,
    seam_tol=0.15,
    max_cap_theta=0.12,
    max_leg_z=6.0,
    coord_tol=1e-6,
):
    """Replace small rectangular caps centered on the periodic seam with a diagonal.

    Full-wrap polygons folded through a synthetic seam can leave tiny U-shaped
    caps after inward offsets. They show up as perfect 90-degree protrusions near
    theta=0 even though the physical boundary should pass straight through the
    seam. Only short caps near an integer-turn seam are touched.
    """
    if not lines:
        return lines

    original_lines = [line for line in lines if line is not None and not line.is_empty]
    edge_lines = []
    for line in original_lines:
        coords = np.asarray(line.coords, dtype=float)
        if len(coords) < 2:
            continue
        step = np.abs(np.diff(coords, axis=0))
        keep = np.concatenate([[True], np.any(step > coord_tol, axis=1)])
        coords = coords[keep]
        for p0, p1 in zip(coords[:-1], coords[1:]):
            if np.linalg.norm(p1 - p0) > coord_tol:
                edge_lines.append(LineString([p0, p1]))
    if not edge_lines:
        return lines
    lines = edge_lines

    def _coords(line):
        coords = np.asarray(line.coords, dtype=float)
        step = np.abs(np.diff(coords, axis=0))
        keep = np.concatenate([[True], np.any(step > coord_tol, axis=1)])
        return coords[keep]

    def _same_point(a, b):
        return abs(a[0] - b[0]) <= coord_tol and abs(a[1] - b[1]) <= coord_tol

    def _point_key(p):
        return (round(float(p[0]) / coord_tol), round(float(p[1]) / coord_tol))

    def _seam_distance(theta):
        return abs((theta - seam_origin + period / 2) % period - period / 2)

    entries = []
    for idx, line in enumerate(lines):
        if line is None or line.is_empty:
            continue
        coords = _coords(line)
        if len(coords) != 2:
            entries.append({"idx": idx, "line": line, "kind": "other"})
            continue
        dx = coords[1, 0] - coords[0, 0]
        dz = coords[1, 1] - coords[0, 1]
        if abs(dz) <= coord_tol and abs(dx) > coord_tol:
            kind = "theta_span"
        elif abs(dx) <= coord_tol and abs(dz) > coord_tol:
            kind = "z_span"
        else:
            kind = "other"
        entries.append({"idx": idx, "line": line, "coords": coords, "kind": kind})

    remove = set()
    replacements = []
    z_spans = [entry for entry in entries if entry.get("kind") == "z_span"]
    z_span_endpoints = defaultdict(list)
    for leg in z_spans:
        l0, l1 = leg["coords"]
        z_span_endpoints[_point_key(l0)].append((leg, l1, l0))
        z_span_endpoints[_point_key(l1)].append((leg, l0, l1))

    for cap in [entry for entry in entries if entry.get("kind") == "theta_span"]:
        cap_idx = cap["idx"]
        if cap_idx in remove:
            continue
        c0, c1 = cap["coords"]
        theta_width = abs(c1[0] - c0[0])
        theta_mid = (c0[0] + c1[0]) * 0.5
        if theta_width > max_cap_theta:
            continue
        if seam_tol is not None and _seam_distance(theta_mid) > seam_tol:
            continue

        legs = []
        for cap_endpoint in (c0, c1):
            matches = []
            for leg, free_endpoint, matched_endpoint in z_span_endpoints.get(_point_key(cap_endpoint), []):
                leg_idx = leg["idx"]
                if leg_idx == cap_idx or leg_idx in remove:
                    continue
                if _same_point(matched_endpoint, cap_endpoint):
                    matches.append((leg, free_endpoint))
            if len(matches) != 1:
                legs = []
                break
            legs.append(matches[0])
        if len(legs) != 2 or legs[0][0]["idx"] == legs[1][0]["idx"]:
            continue

        free0 = legs[0][1]
        free1 = legs[1][1]
        cap_z = c0[1]
        if not (cap_z < min(free0[1], free1[1]) - coord_tol or cap_z > max(free0[1], free1[1]) + coord_tol):
            continue
        if abs(free0[1] - cap_z) > max_leg_z or abs(free1[1] - cap_z) > max_leg_z:
            continue

        remove.update({cap_idx, legs[0][0]["idx"], legs[1][0]["idx"]})
        replacements.append(LineString([free0, free1]))

    if not remove:
        return original_lines

    cleaned = [line for idx, line in enumerate(lines) if idx not in remove]
    cleaned.extend(replacements)
    try:
        merged = linemerge(MultiLineString(cleaned))
    except (ValueError, shapely.errors.GEOSException):
        return cleaned
    merged_lines = []
    _collect_lines(merged, merged_lines, min_length=0.0)
    return merged_lines or cleaned


def _remove_local_orthogonal_caps(lines, period=2 * math.pi):
    """Remove small U-shaped perimeter shelves away from the periodic cut.

    At the first layers where a full-wrap ring appears, the full-wrap cell
    reconstruction can leave a tiny orthogonal shelf where the ring's axial edge
    meets a sloped cube face. It has the same three-edge topology as a seam cap,
    but it is local to the part geometry rather than the artificial theta cut.
    """
    return _remove_periodic_seam_caps(
        lines,
        period=period,
        seam_tol=None,
        max_cap_theta=0.08,
        max_leg_z=8.0,
    )


def _remove_periodic_seam_aligned_edges(
    lines,
    seam_origins=(0.0,),
    period=2 * math.pi,
    seam_tol=0.15,
    max_theta_span=0.25,
    min_z_span=1e-4,
):
    """Drop artificial cut edges that lie on a periodic theta seam.

    Full-wrap polygons need a flat cut edge so Shapely can store them as planar
    polygons, but that edge is not a physical perimeter on the cylinder.  When
    the true boundary winds through the seam, periodic unioning does not always
    cancel the cut edge because the two seam intervals can sit at different z
    ranges.  Remove only short-in-theta, z-changing segments centered on known
    seam origins; real boundary traces leave the seam with a larger theta span
    and are preserved.
    """
    if not lines:
        return []

    seam_origins = tuple(seam_origins or (0.0,))

    def _seam_distance(theta, origin):
        return abs((theta - origin + period / 2) % period - period / 2)

    def _is_cut_edge(p0, p1):
        dz = abs(float(p1[1] - p0[1]))
        if dz <= min_z_span:
            return False
        raw_span = abs(float(p1[0] - p0[0]))
        wrapped_span = abs((float(p1[0] - p0[0]) + period / 2) % period - period / 2)
        theta_span = min(raw_span, wrapped_span)
        if theta_span > max_theta_span:
            return False
        for origin in seam_origins:
            if max(_seam_distance(p0[0], origin), _seam_distance(p1[0], origin)) <= seam_tol:
                return True
        return False

    result = []
    for line in lines:
        if not isinstance(line, LineString) or line.is_empty:
            continue
        coords = np.asarray(line.coords, dtype=float)
        if len(coords) < 2:
            continue

        current = [coords[0]]
        for p1 in coords[1:]:
            p0 = current[-1]
            if _is_cut_edge(p0, p1):
                if len(current) >= 2:
                    result.append(LineString(current))
                current = [p1]
                continue
            current.append(p1)

        if len(current) >= 2:
            result.append(LineString(current))

    return _dedupe_lines(result, min_length=1e-9)


def _split_local_seam_closure_edges(
    lines,
    period=2 * math.pi,
    min_short_theta=1e-4,
    max_short_theta=0.45,
    z_tol=1e-3,
    min_length=1e-9,
):
    """Split local seam-crossing loop closures out of perimeter lines.

    Full-wrap perimeter extraction intentionally works in a shifted one-turn
    print window. A small lobe that crosses that window seam can come back as an
    open line whose closing chord is represented as almost a full positive turn
    (for example 6.39 -> 12.45 rad). That is topologically the short chord
    through the seam, but forward-only G-code would print the long way around.
    Keep true two-point full-circle boundaries, and split only multi-vertex
    local closures whose complement is a short angular span.
    """
    if not lines:
        return lines

    result = []
    for line in lines:
        if line is None or line.is_empty:
            continue
        coords = np.asarray(line.coords, dtype=float)
        if len(coords) < 2:
            continue
        current = [coords[0]]
        for p0, p1 in zip(coords[:-1], coords[1:]):
            dtheta = p1[0] - p0[0]
            dz = p1[1] - p0[1]
            split_edge = False

            if abs(dz) <= z_tol and dtheta > period / 2:
                short = period - dtheta
                if min_short_theta < short < max_short_theta:
                    if len(current) >= 2:
                        result.append(LineString(current))
                    result.append(LineString([p1 - [period, 0.0], p0]))
                    current = [p1]
                    split_edge = True
            elif abs(dz) <= z_tol and dtheta < -period / 2:
                short = period + dtheta
                if min_short_theta < short < max_short_theta:
                    if len(current) >= 2:
                        result.append(LineString(current))
                    result.append(LineString([p0 - [period, 0.0], p1]))
                    current = [p1]
                    split_edge = True

            if not split_edge:
                current.append(p1)

        if len(current) >= 2:
            result.append(LineString(current))

    return [line for line in result if line.length >= min_length]


def _split_path_local_seam_closure_edges(
    path,
    period=2 * math.pi,
    min_short_theta=1e-4,
    max_short_theta=0.45,
    x_tol=1e-3,
):
    """Split near-full-turn local closure edges after path rotation."""
    if path is None or not getattr(path, "points", None):
        return [path]

    pts = np.asarray(path.points, dtype=float)
    if len(pts) < 2:
        return [path]

    def _restore_open_thetas(theta_values):
        diffs = np.diff(theta_values)
        restored = []
        for diff in diffs:
            wrapped = (diff + math.pi) % period - math.pi
            k = round((diff - wrapped) / period)
            step = wrapped + k * period
            if step < -1e-9:
                step += math.ceil((-step - 1e-9) / period) * period
            restored.append(step)
        return np.concatenate([[theta_values[0]], theta_values[0] + np.cumsum(restored)])

    if getattr(path, "is_accumulated", False):
        emitted_theta = pts[:, 0].copy()
    elif path.closed:
        emitted_theta = np.unwrap(pts[:, 0])
    else:
        emitted_theta = _restore_open_thetas(pts[:, 0])

    parts = []
    current = [pts[0]]
    split_any = False

    for idx, (p0, p1) in enumerate(zip(pts[:-1], pts[1:])):
        dtheta = emitted_theta[idx + 1] - emitted_theta[idx]
        dx = p1[1] - p0[1]
        split_edge = False

        if abs(dx) <= x_tol and dtheta > period / 2:
            short = period - (dtheta % period)
            if min_short_theta < short < max_short_theta:
                if len(current) >= 2:
                    parts.append(current)
                parts.append([p1, p0])
                current = [p1]
                split_edge = True
                split_any = True
        elif abs(dx) <= x_tol and dtheta < -period / 2:
            short = period - ((-dtheta) % period)
            if min_short_theta < short < max_short_theta:
                if len(current) >= 2:
                    parts.append(current)
                parts.append([p0, p1])
                current = [p1]
                split_edge = True
                split_any = True

        if not split_edge:
            current.append(p1)

    if len(current) >= 2:
        parts.append(current)

    if not split_any:
        return [path]

    split_paths = []
    for coords in parts:
        arr = np.asarray(coords, dtype=float)
        if len(arr) < 2:
            continue
        arr[:, 0] = np.mod(arr[:, 0], period)
        new_path = path.__class__(arr.tolist(), False, path.kind, path.radius)
        split_paths.append(new_path)
    return split_paths or [path]


def keep_canonical_polygons(geom, lower=0.0, upper=2 * math.pi):
    """Shift polygon components onto the canonical turn and deduplicate copies.

    `external_peris` intentionally carries seam-neighbor copies at ±2π so support
    and boundary extraction can see across the wrap. Centroid-based filtering is
    not robust once a valid polygon itself spans across the seam, because its
    centroid may fall just outside `[0, 2π)` even though most of the geometry is
    the canonical physical copy. Instead, shift each polygon by integer turns to
    maximize its overlap with the canonical window, then union the results.
    """
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()

    period = upper - lower

    def _canonical_shift(poly):
        minx, _, maxx, _ = poly.bounds
        k_min = int(math.floor((lower - maxx) / period)) - 1
        k_max = int(math.ceil((upper - minx) / period)) + 1
        best_k = 0
        best_overlap = -1.0
        best_center_dist = float("inf")
        for k in range(k_min, k_max + 1):
            shifted_min = minx + k * period
            shifted_max = maxx + k * period
            overlap = max(0.0, min(shifted_max, upper) - max(shifted_min, lower))
            center = poly.centroid.x + k * period
            center_dist = abs(center - (lower + upper) * 0.5)
            if overlap > best_overlap + 1e-9 or (
                abs(overlap - best_overlap) <= 1e-9 and center_dist < best_center_dist
            ):
                best_k = k
                best_overlap = overlap
                best_center_dist = center_dist
        return best_k, best_overlap

    if isinstance(geom, Polygon):
        k, overlap = _canonical_shift(geom)
        if overlap <= 1e-9:
            return shapely.GeometryCollection()
        return shift_polygon(geom, k * period) if k else geom

    if isinstance(geom, shapely.MultiPolygon):
        polys = [keep_canonical_polygons(g, lower=lower, upper=upper) for g in geom.geoms if isinstance(g, Polygon)]
        polys = [g for g in polys if g is not None and not g.is_empty]
        if not polys:
            return shapely.GeometryCollection()
        if len(polys) == 1:
            return polys[0]
        return shapely.unary_union(polys)

    if isinstance(geom, shapely.GeometryCollection):
        kept = [keep_canonical_polygons(g, lower=lower, upper=upper) for g in geom.geoms]
        kept = [g for g in kept if g is not None and not g.is_empty]
        if not kept:
            return shapely.GeometryCollection()
        return shapely.unary_union(kept)

    return geom


def polygonal_only(geom):
    """Drop non-area fragments from a geometry tree."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()
    if isinstance(geom, (Polygon, shapely.MultiPolygon)):
        return geom
    if isinstance(geom, shapely.GeometryCollection):
        polys = [polygonal_only(g) for g in geom.geoms]
        polys = [g for g in polys if g is not None and not g.is_empty]
        if not polys:
            return shapely.GeometryCollection()
        return shapely.unary_union(polys)
    if hasattr(geom, "geoms"):
        polys = [polygonal_only(g) for g in geom.geoms]
        polys = [g for g in polys if g is not None and not g.is_empty]
        if not polys:
            return shapely.GeometryCollection()
        return shapely.unary_union(polys)
    return shapely.GeometryCollection()


def repair_polygonal_geometry(geom, grid_size=1e-8):
    """Normalize invalid overlay/buffer results into polygonal geometry."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()

    repaired = geom
    for repair in (
        lambda g: shapely.make_valid(g),
        lambda g: shapely.set_precision(g, grid_size),
        lambda g: g.buffer(0),
    ):
        try:
            repaired = repair(repaired)
        except (ValueError, shapely.errors.GEOSException):
            continue

    return polygonal_only(repaired)


def clean_polygonal_result(geom, grid_size=1e-8):
    """Keep valid polygonal output cheap; repair only invalid results."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()
    try:
        if geom.is_valid:
            return polygonal_only(geom)
    except (ValueError, shapely.errors.GEOSException):
        pass
    return repair_polygonal_geometry(geom, grid_size)


def robust_polygonal_buffer(geom, distance, grid_size=1e-8, **kwargs):
    """Buffer polygon geometry, repairing precision/topology failures."""
    if geom is None or geom.is_empty:
        return shapely.GeometryCollection()

    last_error = None
    try:
        return clean_polygonal_result(geom.buffer(distance, **kwargs), grid_size)
    except (ValueError, shapely.errors.GEOSException) as exc:
        last_error = exc

    candidate = repair_polygonal_geometry(geom, grid_size)
    if candidate is not None and not candidate.is_empty:
        try:
            return clean_polygonal_result(candidate.buffer(distance, **kwargs), grid_size)
        except (ValueError, shapely.errors.GEOSException) as exc:
            last_error = exc

    if last_error is not None:
        raise last_error
    return shapely.GeometryCollection()


def robust_polygonal_intersection(a, b, grid_size=1e-8):
    """Intersect polygon geometry with repair/retry for GEOS topology errors."""
    if a is None or b is None or a.is_empty or b.is_empty:
        return shapely.GeometryCollection()

    last_error = None
    candidate_pairs = [(a, b, None)]
    for left, right, precision in candidate_pairs:
        if left is None or right is None or left.is_empty or right.is_empty:
            continue
        try:
            if precision is None:
                result = shapely.intersection(left, right)
            else:
                result = shapely.intersection(left, right, grid_size=precision)
            return clean_polygonal_result(result, grid_size)
        except (ValueError, shapely.errors.GEOSException) as exc:
            last_error = exc

    for precision in (None, grid_size, grid_size * 10):
        repair_grid = precision or grid_size
        left = repair_polygonal_geometry(a, repair_grid)
        right = repair_polygonal_geometry(b, repair_grid)
        if left is None or right is None or left.is_empty or right.is_empty:
            continue
        try:
            if precision is None:
                result = shapely.intersection(left, right)
            else:
                result = shapely.intersection(left, right, grid_size=precision)
            return clean_polygonal_result(result, grid_size)
        except (ValueError, shapely.errors.GEOSException) as exc:
            last_error = exc

    if last_error is not None:
        raise last_error
    return shapely.GeometryCollection()


def robust_polygonal_difference(a, b, grid_size=1e-8):
    """Subtract polygon geometry with repair/retry for GEOS topology errors."""
    if a is None or a.is_empty:
        return shapely.GeometryCollection()
    if b is None or b.is_empty:
        return repair_polygonal_geometry(a, grid_size)

    last_error = None
    candidate_pairs = [(a, b, None)]
    for left, right, precision in candidate_pairs:
        if left is None or left.is_empty:
            return shapely.GeometryCollection()
        if right is None or right.is_empty:
            return clean_polygonal_result(left, grid_size)
        try:
            if precision is None:
                result = shapely.difference(left, right)
            else:
                result = shapely.difference(left, right, grid_size=precision)
            return clean_polygonal_result(result, grid_size)
        except (ValueError, shapely.errors.GEOSException) as exc:
            last_error = exc

    for precision in (None, grid_size, grid_size * 10):
        repair_grid = precision or grid_size
        left = repair_polygonal_geometry(a, repair_grid)
        right = repair_polygonal_geometry(b, repair_grid)
        if left is None or left.is_empty:
            return shapely.GeometryCollection()
        if right is None or right.is_empty:
            return clean_polygonal_result(left, grid_size)
        try:
            if precision is None:
                result = shapely.difference(left, right)
            else:
                result = shapely.difference(left, right, grid_size=precision)
            return clean_polygonal_result(result, grid_size)
        except (ValueError, shapely.errors.GEOSException) as exc:
            last_error = exc

    if last_error is not None:
        raise last_error
    return shapely.GeometryCollection()


def robust_polygonal_union(geoms, grid_size=1e-8):
    """Union polygon geometry with repair/retry for GEOS topology errors."""
    parts = [g for g in geoms if g is not None and not g.is_empty]
    parts = [g for g in parts if g is not None and not g.is_empty]
    if not parts:
        return shapely.GeometryCollection()
    if len(parts) == 1:
        return parts[0]

    last_error = None
    for precision in (None, grid_size, grid_size * 10):
        try:
            if precision is None:
                result = shapely.unary_union(parts)
            else:
                result = shapely.union_all(parts, grid_size=precision)
            return clean_polygonal_result(result, grid_size)
        except (ValueError, shapely.errors.GEOSException) as exc:
            last_error = exc
            parts = [repair_polygonal_geometry(g, (precision or grid_size) * 10) for g in parts]
            parts = [g for g in parts if g is not None and not g.is_empty]

    if last_error is not None:
        raise last_error
    return shapely.GeometryCollection()


def remove_tiny_regions(geom_theta_z, radius, min_area_mm2):
    """Drop polygon fragments too small to print after converting theta to mm."""
    geom_theta_z = polygonal_only(geom_theta_z)
    if geom_theta_z.is_empty or min_area_mm2 <= 0:
        return geom_theta_z

    parts = list(geom_theta_z.geoms) if hasattr(geom_theta_z, "geoms") else [geom_theta_z]
    kept = []
    for part in parts:
        if part.is_empty:
            continue
        metric = shapely.affinity.scale(part, xfact=radius, yfact=1.0, origin=(0, 0))
        if metric.area >= min_area_mm2:
            kept.append(part)
    if not kept:
        return shapely.GeometryCollection()
    return shapely.unary_union(kept)

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
    "surface_min_area": 0.05,        # mm² — discard numerical top/bottom slivers below this printable area
    "surface_min_solid_lines": 8,    # discard top/bottom skin islands too small to form a printable patch
    "slice_sample_offset": 0.25,     # fraction of layer height used to sample slice geometry within the layer interior
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
max_model_layer = args.max_layer
mesh = trimesh.load_mesh(args.mesh)
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

# Subdivide large faces so cylindrical cross-sections have enough sample points
# along long edges (avoids jitter from straight-line approximation of curved arcs).
#
# This used to run up to 8 unconditional loop-subdivision passes. Since every pass
# splits every face into 4, one very long edge could force the whole mesh to explode
# to tens of millions of faces before slicing started. Keep the same manifold-safe
# subdivision, but cap it with a face budget and make the knobs explicit.
_max_edge = float(os.environ.get(
    "SLICER_MAX_EDGE_MM",
    slicing_params["layer_height"] * 20,
))
_max_subdivide_faces = int(os.environ.get("SLICER_MAX_SUBDIVIDE_FACES", "2000000"))
_max_subdivide_passes = int(os.environ.get("SLICER_MAX_SUBDIVIDE_PASSES", "8"))
_disable_subdivide = os.environ.get("SLICER_DISABLE_SUBDIVIDE", "").lower() in {"1", "true", "yes"}
_orig_faces = len(mesh.faces)
_sub_verts, _sub_faces = mesh.vertices, mesh.faces
_subdivide_passes = 0
_cur_max = trimesh.Trimesh(_sub_verts, _sub_faces, process=False).edges_unique_length.max()
for _ in range(_max_subdivide_passes):
    if _disable_subdivide:
        break
    _cur_max = trimesh.Trimesh(_sub_verts, _sub_faces, process=False).edges_unique_length.max()
    if _cur_max <= _max_edge:
        break
    _next_face_count = len(_sub_faces) * 4
    if _next_face_count > _max_subdivide_faces:
        print(
            f"Subdivision capped at {len(_sub_faces)} faces "
            f"(next pass would be {_next_face_count}, cap {_max_subdivide_faces})."
        )
        break
    _sub_verts, _sub_faces = trimesh.remesh.subdivide(_sub_verts, _sub_faces)
    _subdivide_passes += 1
    gc.collect()
mesh = trimesh.Trimesh(_sub_verts, _sub_faces, process=False)
_final_max_edge = mesh.edges_unique_length.max()
print(
    f"Subdivided mesh: {_orig_faces} -> {len(mesh.faces)} faces "
    f"({_subdivide_passes} pass(es), target max edge {_max_edge:.2f} mm, "
    f"actual {_final_max_edge:.2f} mm)"
)
if _final_max_edge > _max_edge:
    print(
        "  warning: subdivision target not fully reached; raise "
        "SLICER_MAX_SUBDIVIDE_FACES/PASSES for smoother cylindrical slices."
    )
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
_layer_height = slicing_params["layer_height"]
_uncapped_slice_radii = np.arange(r_start, r_end, _layer_height)
_uncapped_layer_count = len(_uncapped_slice_radii)
if max_model_layer is not None:
    _requested_layer_count = max_model_layer + 1
    _slice_layer_count = min(_uncapped_layer_count, _requested_layer_count)
else:
    _slice_layer_count = _uncapped_layer_count
_slice_radii = _uncapped_slice_radii[:_slice_layer_count]
if len(_slice_radii):
    _cap_trimmed_layers = _slice_layer_count < _uncapped_layer_count
    r_effective_end = float(_slice_radii[-1] + _layer_height) if _cap_trimmed_layers else r_end
else:
    _cap_trimmed_layers = False
    r_effective_end = r_start
print(f"Mesh radial range: {np.min(pts[:,0]):.2f} – {np.max(pts[:,0]):.2f} mm")
print(f"Slicing range:     {r_start:.2f} – {r_effective_end:.2f} mm  ({_slice_layer_count} layers)")
if max_model_layer is not None:
    if _cap_trimmed_layers:
        print(f"Layer cap:         keeping model layers 0–{max_model_layer}; layers above {max_model_layer} ignored")
    else:
        print(f"Layer cap:         --max-layer {max_model_layer} is outside the model range; no layers trimmed")
if r_start >= r_end:
    raise ValueError(f"bed_radius ({printer_params['bed_radius']} mm) is at or beyond the mesh max radius ({np.max(pts[:,0]):.2f} mm) — nothing to slice.")
slicing_params["r_base"] = r_start
_total_layers = len(_slice_radii)
print(f"\nPass 1 — slicing geometry ({_total_layers} layers, r={r_start:.2f}–{r_effective_end:.2f} mm)")
for layer_i,current_r in enumerate(_slice_radii):
    times = [(time.perf_counter(), "")]
    print(f"  [{sliced_layer_count+1}/{_total_layers}] r={current_r:.3f} mm", end="\r", flush=True)
    # Sample within the radial interior of the deposited layer, not exactly on its
    # boundary. Boundary sampling is prone to degenerate cross-sections when the
    # slice coincides with mesh vertices or tangent transitions, which can fragment
    # one physical wrap band into many tiny loops on a single layer.
    _slice_r = current_r + slicing_params["layer_height"] * slicing_params["slice_sample_offset"]
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
    save_raw_stl_layer_plot(layer_i, current_r, _slice_r, poly_pts_raw)
    poly_pts = [p for p in poly_pts_raw if len(p) >= 3]
    # for poly in poly_pts:
    #     plt.plot(poly[:,0],poly[:,1],"-")
    # plt.show()
    full_circles = [is_full_circle(poly) for poly in poly_pts]
    fc_original_indices = [i for i, fc in enumerate(full_circles) if fc]
    complex_fc_loops = []
    simple_fc_indices = []
    for i in fc_original_indices:
        if is_simple_full_circle_loop(poly_pts[i]):
            simple_fc_indices.append(i)
        else:
            complex_fc_loops.append(poly_pts[i])
    save_classified_loops_layer_plot(
        layer_i,
        current_r,
        poly_pts,
        fc_original_indices,
        simple_fc_indices,
        complex_fc_loops,
    )

    # Once any full-turn boundary is complex, simple full-turn loops on the
    # same layer may still be part of the same periodic occupied region.  Pairing
    # or dropping them independently changes phase at the seam, so reconstruct
    # the whole full-circle set together from mesh occupancy.
    if complex_fc_loops:
        full_circle_pairs = []
        fc_loops_for_reconstruction = [poly_pts[i] for i in fc_original_indices]
    else:
        fc_loops_for_reconstruction = []
    if not complex_fc_loops and len(simple_fc_indices) >= 2:
        fc_z_mean = [np.mean(np.array(poly_pts[i])[:, 1]) for i in simple_fc_indices]
        sorted_order = np.argsort(fc_z_mean)
        sorted_fc = [simple_fc_indices[j] for j in sorted_order]
        full_circle_pairs = [np.array([sorted_fc[k], sorted_fc[k+1]]) for k in range(0, len(sorted_fc) - 1, 2)]
    elif not complex_fc_loops:
        full_circle_pairs = []
    full_circle_pairs = np.array(full_circle_pairs) if full_circle_pairs else np.empty((0, 2), dtype=int)

    complex_full_circle_parts = reconstruct_complex_full_circle_regions(
        fc_loops_for_reconstruction,
        mesh,
        _slice_r,
        slicing_params["nozzle_size"],
        debug_layer=layer_i,
        debug_plot_radius=current_r,
    )
    complex_full_circle_parts, promoted_periodic_holes = promote_periodic_reconstruction_holes(
        complex_full_circle_parts,
        layer_i=layer_i,
        current_r=current_r,
    )
    if promoted_periodic_holes:
        print(
            f"\nPromoted {len(promoted_periodic_holes)} seam-spanning periodic "
            f"hole candidate(s) on layer {layer_i}."
        )

    base_polys = [Polygon(p) for i,p in enumerate(poly_pts) if i not in fc_original_indices]
    base_polys.extend(promoted_periodic_holes)
    polys = list(base_polys)
    periodic_poly_indices = set()
    _complex_start = len(polys)
    polys.extend(complex_full_circle_parts)
    periodic_poly_indices.update(range(_complex_start, len(polys)))
    mapping = set()
    simple_pair_parts = []
    for pair in full_circle_pairs:
        p1 = canonicalize_full_circle_loop(poly_pts[pair[0]])
        p2 = canonicalize_full_circle_loop(poly_pts[pair[1]])
        final_poly = np.concatenate([p1, p2[::-1, :]])
        final_poly = shapely.make_valid(Polygon(final_poly))
        # make_valid can return a MultiPolygon when the strip self-intersects.
        # Flatten into individual
        # Polygons so the containment hierarchy loop (which only handles
        # Polygon instances) processes each component.
        fp_parts = list(final_poly.geoms) if hasattr(final_poly, "geoms") else [final_poly]
        for part in fp_parts:
            if not isinstance(part, Polygon) or part.is_empty:
                continue
            mapping.add(len(polys))
            periodic_poly_indices.add(len(polys))
            polys.append(part)
            simple_pair_parts.append(part)
    save_reconstructed_candidates_layer_plot(
        layer_i,
        current_r,
        base_polys,
        complex_full_circle_parts,
        simple_pair_parts,
    )
    # for poly in poly_pts:
    #     if is_full_circle(poly):
    #         plt.plot(poly[:,0],poly[:,1])
    # plt.show()
    hier = build_hierarchy(polys, periodic_indices=periodic_poly_indices)
    times.append((time.perf_counter(),"build_hierarchy"))

    flat_hier = [each2 for each in hier for each2 in each]
    levels = [0]*len(polys)
    _levels_i,_levels_count = np.unique(flat_hier,return_counts=True)
    for i,count in zip(_levels_i,_levels_count):
        levels[i] = count
    save_hierarchy_layer_plot(layer_i, current_r, polys, hier, levels)
    external_perimeters = []
    external_perimeters_print = []
    stitched_hole_debug_polys = []
    times.append((time.perf_counter(),"count levels"))


    plot_perimeters = False
    def do_transform(x):
        return shapely.transform(x, lambda p: p * [current_r, 1])
    def do_transform_inverse(x):
        return shapely.transform(x, lambda p: p * [1/current_r, 1])
    def do_transform_at_radius(x, radius):
        return shapely.transform(x, lambda p: p * [radius, 1])
    def do_transform_inverse_at_radius(x, radius):
        return shapely.transform(x, lambda p: p * [1/radius, 1])

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

    def _clean_polygon_shell(poly):
        """Return a polygonal shell with slicer-scale collinear cleanup applied."""
        if not isinstance(poly, Polygon) or poly.is_empty:
            return shapely.GeometryCollection()
        coords = remove_collinear_points_angle(
            np.asarray(poly.exterior.coords, dtype=float) * np.array([current_r, 1]),
            0.01,
        ) * np.array([1 / current_r, 1])
        if len(coords) < 4:
            return shapely.GeometryCollection()
        try:
            return polygonal_only(shapely.make_valid(Polygon(coords)))
        except (ValueError, shapely.errors.GEOSException):
            return shapely.GeometryCollection()

    def _is_contained_hole_copy(parent_geom, hole_geom):
        if hole_geom is None or hole_geom.is_empty:
            return False
        try:
            return parent_geom.contains(hole_geom) or parent_geom.covers(hole_geom)
        except shapely.errors.GEOSException:
            return False

    def _contained_hole_copies(parent_geom, child_poly, period=2 * np.pi):
        child_shell = _clean_polygon_shell(child_poly)
        if child_shell.is_empty:
            return []
        child_parts = list(child_shell.geoms) if hasattr(child_shell, "geoms") else [child_shell]
        hole_parts = []
        for child_part in child_parts:
            if not isinstance(child_part, Polygon) or child_part.is_empty:
                continue
            for shift in (-period, 0.0, period):
                candidate = shift_polygon(child_part, shift) if shift else child_part
                if _is_contained_hole_copy(parent_geom, candidate):
                    hole_parts.append(candidate)
        return hole_parts

    for poly_i,sub_polys in enumerate(hier):
        if levels[poly_i] % 2 == 1:
            continue
        if not isinstance(polys[poly_i], Polygon) or polys[poly_i].is_empty:
            continue
        _is_periodic_poly = poly_i in periodic_poly_indices
        parent_shell = _clean_polygon_shell(polys[poly_i])
        if parent_shell.is_empty:
            continue
        if _is_periodic_poly:
            parent_for_holes = stitch_periodic_shape(parent_shell)
        else:
            parent_for_holes = parent_shell
        if parent_for_holes.is_empty:
            continue

        hole_parts = []
        for sub_poly in sub_polys:
            if (
                levels[sub_poly] != levels[poly_i] + 1
                or not isinstance(polys[sub_poly], Polygon)
                or polys[sub_poly].is_empty
            ):
                continue
            hole_parts.extend(_contained_hole_copies(parent_for_holes, polys[sub_poly]))

        if hole_parts:
            poly_with_holes = robust_polygonal_difference(
                parent_for_holes,
                robust_polygonal_union(hole_parts),
            )
        else:
            poly_with_holes = parent_for_holes
        poly_with_holes = repair_polygonal_geometry(poly_with_holes)
        if poly_with_holes.is_empty:
            continue
        if layer_i in debug_layers:
            stitched_hole_debug_polys.append(poly_with_holes)

        if plot_perimeters:
            plot_polygon(do_transform(poly_with_holes),add_points=False,facecolor=(0, 0, 0, 0))

        _theta_min, _, _theta_max, _ = poly_with_holes.bounds
        _theta_span = _theta_max - _theta_min
        _buffer_source = poly_with_holes
        _needs_canonical_buffer_fold = False
        if _is_periodic_poly:
            _needs_canonical_buffer_fold = True
        elif (
            _theta_span <= 2 * np.pi + 1e-6
            and (
                (_theta_min <= 1e-6 and _theta_max >= 2 * np.pi - 1e-6)
                or _is_effectively_full_wrap_polygon(poly_with_holes)
            )
        ):
            # Full-wrap one-turn bands are periodic across the 0/2π seam.  Buffer
            # them with seam-neighbour copies so the inward wall offset does not
            # erode along synthetic seam endcaps.
            _buffer_source = stitch_periodic_shape(poly_with_holes)
            _needs_canonical_buffer_fold = True

        _buffer_source_metric = do_transform(_buffer_source)
        _buffer_source_metric = repair_polygonal_geometry(_buffer_source_metric)
        _buffered_metric = robust_polygonal_buffer(
            _buffer_source_metric,
            -slicing_params["nozzle_size"] / 2,
            join_style="mitre",
        )
        _buffered_metric = robust_polygonal_intersection(_buffered_metric, _buffer_source_metric)
        _poly_buffered_raw = do_transform_inverse(_buffered_metric)
        if _is_periodic_poly:
            # The periodic parent was already stitched from -2π/0/+2π copies for
            # hole assignment and buffering. Folding that stitched geometry by
            # unioning shifted copies can fill holes when any neighbour copy covers
            # the same canonical coordinate. Clip the stitched result directly to
            # each one-turn window so voids survive into skin and infill clipping.
            poly_with_holes = robust_polygonal_intersection(
                _poly_buffered_raw,
                theta_clip_box,
            )
            poly_with_holes_print = robust_polygonal_intersection(
                _poly_buffered_raw,
                box(math.pi, z_clip_min, math.pi + 2 * math.pi, z_clip_max),
            )
        elif _needs_canonical_buffer_fold:
            poly_with_holes = fold_to_canonical(_poly_buffered_raw, theta_clip_box)
            poly_with_holes_print = fold_to_theta_window(
                _poly_buffered_raw,
                math.pi,
                math.pi + 2 * math.pi,
                z_clip_min,
                z_clip_max,
            )
        else:
            poly_with_holes = _poly_buffered_raw
            poly_with_holes_print = _poly_buffered_raw
        poly_with_holes = robust_polygonal_intersection(poly_with_holes, axial_clip_box)
        poly_with_holes_print = robust_polygonal_intersection(poly_with_holes_print, axial_clip_box)
        if poly_with_holes.is_empty:
            continue

        # Normalize each physical perimeter into a canonical turn before adding
        # seam-neighbor copies. This keeps one stored copy selectable by the
        # centroid-based canonical filter even when the raw geometry spans many
        # full turns.
        _k_canon = _shift_periods_to_canonical_centroid(poly_with_holes)
        if _k_canon:
            poly_with_holes = shift_polygon(poly_with_holes, _k_canon * 2 * np.pi)
            if not _needs_canonical_buffer_fold:
                poly_with_holes_print = shift_polygon(poly_with_holes_print, _k_canon * 2 * np.pi)

        external_perimeters.append(poly_with_holes)
        if not poly_with_holes_print.is_empty:
            external_perimeters_print.append(poly_with_holes_print)
        _stored_theta_min, _, _stored_theta_max, _ = poly_with_holes.bounds
        _stored_theta_span = _stored_theta_max - _stored_theta_min
        # Once a physical perimeter already spans a full canonical turn, adding
        # ±2π seam-neighbour copies makes those copies overlap. `unary_union`
        # can then stitch the exterior out of mixed turns, which rotates local
        # features when separate islands merge into one wrap-touching polygon.
        # Only add seam-neighbour copies for genuinely sub-turn islands.
        if poly_i not in periodic_poly_indices and _stored_theta_span < 2 * np.pi - 1e-6:
            external_perimeters.append(shift_polygon(poly_with_holes,2*np.pi))
            external_perimeters.append(shift_polygon(poly_with_holes,-2*np.pi))
        if plot_perimeters:
            plot_polygon(do_transform(external_perimeters[-1]),add_points=False,color="green",facecolor=(0, 0, 0, 0))
    if plot_perimeters:
        plt.show()
    save_stitched_holes_layer_plot(layer_i, current_r, stitched_hole_debug_polys)
    times.append((time.perf_counter(),"form poly with holes"))
    
    final = robust_polygonal_union(external_perimeters)
    layer_data[layer_i]["external_peris_has_periodic"] = any(
        poly_i in periodic_poly_indices
        and levels[poly_i] % 2 == 0
        and isinstance(polys[poly_i], Polygon)
        and not polys[poly_i].is_empty
        for poly_i in range(len(polys))
    )
    layer_data[layer_i]["external_peris"] = final
    layer_data[layer_i]["external_peris_print"] = (
        robust_polygonal_union(external_perimeters_print)
        if external_perimeters_print
        else final
    )
    times.append((time.perf_counter(),"full union"))

    if layer_i == 0:
        layer_data[layer_i]["bottom_surfaces"] = layer_data[layer_i]["external_peris"]
        continue

    e = do_transform(final)

    prev_periodic = periodic_copies(layer_data[layer_i-1]["external_peris"], turns=2)
    e = robust_polygonal_difference(
        e,
        robust_polygonal_buffer(do_transform(prev_periodic), slicing_params["nozzle_size"] / 2),
    )

    if not e.is_empty:
        ns = do_transform_inverse(e)
        ns = clean_support_region(
            ns, current_r,
            min_area_mm2=slicing_params["support_min_area"],
            merge_gap_mm=slicing_params["support_merge_gap"],
        )
        # Fold shifted fragments into the canonical [0, 2π] range and remove
        # areas already occupied by the model at this layer.
        ns = keep_canonical_polygons(ns)
        if ns.area > 1e-6:
            layer_data[layer_i]["needs_support"] = ns
    times.append((time.perf_counter(),"support check"))

    _surface_cover_tol = slicing_params["nozzle_size"] * 3.0

    prev_r = layer_data[layer_i-1]["radius"]
    e = do_transform_at_radius(layer_data[layer_i-1]["external_peris"], prev_r)
    # Shrink the previous layer before testing for top surfaces so wall centerline
    # offsets are not classified as exposed top skin. Compare coverage in metric
    # space, with a nozzle-scale tolerance, so tiny phase differences across the
    # periodic seam do not create solid infill inside enclosed cube/ring volume.
    e = robust_polygonal_buffer(e, -slicing_params["nozzle_size"], join_style="mitre")
    current_cover = robust_polygonal_buffer(
        do_transform_at_radius(periodic_copies(final, turns=2), prev_r),
        _surface_cover_tol,
        join_style="mitre",
    )
    e = do_transform_inverse_at_radius(robust_polygonal_difference(e, current_cover), prev_r)
    e = remove_tiny_regions(
        e,
        layer_data[layer_i-1]["radius"],
        slicing_params.get("surface_min_area", 0.0),
    )
    times.append((time.perf_counter(),"top surface check"))

    layer_data[layer_i-1]["top_surfaces"] = e

    # Bottom surfaces of current layer: areas not covered by the layer below.
    # Use the same metric-space tolerant comparison as the top-surface test.
    e = do_transform(final)
    e = robust_polygonal_buffer(e, -slicing_params["nozzle_size"], join_style="mitre")
    prev_cover = robust_polygonal_buffer(
        do_transform(prev_periodic),
        _surface_cover_tol,
        join_style="mitre",
    )
    e = do_transform_inverse(robust_polygonal_difference(e, prev_cover))
    e = remove_tiny_regions(
        e,
        current_r,
        slicing_params.get("surface_min_area", 0.0),
    )
    layer_data[layer_i]["bottom_surfaces"] = e
    times.append((time.perf_counter(), "bottom surface check"))
    # for i in range(1,len(times)):
    #     print(f"{times[i][1]} took {times[i][0]-times[i-1][0]:.2f} seconds")



if not layer_data:
    limit_text = f" up to --max-layer {max_model_layer}" if max_model_layer is not None else ""
    raise ValueError(f"No non-empty model layers were sliced{limit_text}.")

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
        _prop = robust_polygonal_intersection(_ts, layer_data[_ik]["external_peris"])
        _prop = remove_tiny_regions(
            _prop,
            layer_data[_ik]["radius"],
            slicing_params.get("surface_min_area", 0.0),
        )
        if _prop.is_empty:
            break
        _existing = layer_data[_ik].get("top_surfaces", _empty_geom)
        layer_data[_ik]["top_surfaces"] = robust_polygonal_union([_existing, _prop])

for _pos, _key, _bs in _bot_sources:
    for _depth in range(1, slicing_params["bottom_layers"]):
        _outer_pos = _pos + _depth
        if _outer_pos >= len(_sorted_keys):
            break
        _ok = _sorted_keys[_outer_pos]
        _prop = robust_polygonal_intersection(_bs, layer_data[_ok]["external_peris"])
        _prop = remove_tiny_regions(
            _prop,
            layer_data[_ok]["radius"],
            slicing_params.get("surface_min_area", 0.0),
        )
        if _prop.is_empty:
            break
        _existing = layer_data[_ok].get("bottom_surfaces", _empty_geom)
        layer_data[_ok]["bottom_surfaces"] = robust_polygonal_union([_existing, _prop])

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
        metric_exp     = robust_polygonal_buffer(metric_ext, xy_dist_mm)
        _zb            = metric_ext.bounds          # (xmin, z_min, xmax, z_max)
        metric_exp     = robust_polygonal_intersection(metric_exp, shapely.geometry.box(-1e9, _zb[1], 1e9, _zb[3]))
        model_expanded = shapely.affinity.scale(metric_exp, xfact=1.0 / current_r, yfact=1.0, origin=(0, 0))
        support_needed = robust_polygonal_difference(support_needed, model_expanded)
    else:
        support_needed = robust_polygonal_difference(support_needed, current_ext)

    # Record what this layer must print as support.
    # Store with ±2π shifted copies (like external_peris) so seam-crossing
    # regions remain as one continuous polygon rather than two split pieces.
    if not support_needed.is_empty:
        support_canon = keep_canonical_polygons(support_needed)
        layer_data[i]["support"] = robust_polygonal_union([
            support_canon,
            shift_polygon(support_canon, 2*np.pi),
            shift_polygon(support_canon, -2*np.pi),
        ])

    # Add this layer's overhangs so inner layers know to support them
    if "needs_support" in layer_data[i]:
        support_needed = robust_polygonal_union([support_needed, layer_data[i]["needs_support"]])
        support_needed = keep_canonical_polygons(support_needed)

    # Organic taper: erode the support footprint so columns narrow toward the bed
    if not support_needed.is_empty :
        metric = shapely.affinity.scale(support_needed, xfact=current_r, yfact=1.0, origin=(0, 0))
        metric = robust_polygonal_buffer(metric, -taper_mm)
        if metric.is_empty:
            support_needed = shapely.GeometryCollection()
        else:
            support_needed = shapely.affinity.scale(metric, xfact=1.0 / current_r, yfact=1.0, origin=(0, 0))
            support_needed = keep_canonical_polygons(support_needed)

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
    gap_footprint = keep_canonical_polygons(
        robust_polygonal_union([
            support_needed,
            keep_canonical_polygons(layer_data[sorted_layer_keys[0]]["external_peris"]),
        ])
    )
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
                m = robust_polygonal_buffer(m, taper_mm * layers_from_top)
                tapered = shapely.affinity.scale(m, xfact=1.0 / r, yfact=1.0, origin=(0, 0))
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
keep_layer_temp = os.environ.get("SLICER_KEEP_LAYER_TEMP", "").lower() in {"1", "true", "yes"}

def on_key(event):
    """Ctrl-C or Q in any plot window exits the program cleanly."""
    if event.key in ('ctrl+c', 'q', 'escape'):
        plt.close('all')
        sys.exit(0)


def plot_line_endpoints_2d(line_list, ax, current_r, marker_size=8):
    """Plot start/end markers for generated line segments in metric space."""
    if not line_list:
        return
    for line in line_list:
        if line is None or line.is_empty:
            continue
        coords = np.array(line.coords)
        if len(coords) < 2:
            continue
        coords_metric = coords * np.array([current_r, 1])
        ax.plot(coords_metric[0, 0], coords_metric[0, 1], 'go',
                markersize=marker_size, alpha=0.7)
        ax.plot(coords_metric[-1, 0], coords_metric[-1, 1], 'ro',
                markersize=marker_size, alpha=0.7)


def plot_polygonal_2d(geom, ax, **kwargs):
    """Plot polygonal debug geometry, ignoring empty collection fragments."""
    geom = polygonal_only(geom)
    if geom is None or geom.is_empty:
        return
    parts = list(geom.geoms) if hasattr(geom, "geoms") else [geom]
    for part in parts:
        if isinstance(part, Polygon) and not part.is_empty:
            plot_polygon(part, ax=ax, add_points=False, **kwargs)


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
first_model_ext = keep_canonical_polygons(
    layer_data[sorted_layer_keys[0]]["external_peris"]
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
    gap_layer["support_infill"] = support_infill
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
    if not keep_layer_temp:
        gap_layer.pop("support_infill", None)

_sorted_model_keys = sorted(layer_data.keys())
for _model_idx, i in enumerate(_sorted_model_keys):
    real_layer_count += 1

    total = len(layer_data)
    print(f">>> Layer {real_layer_count}/{total} (layer_data key: {i})", end="\r",flush=True)
    export_debug_layer = i in debug_layers
    draw_layer_plot = plot_per_layer or export_debug_layer
    include_layer_3d = plot_per_layer or args.debug_plot_3d

    # Per-layer figure: 2D plus optional 3D for interactive/debug output.
    ax_2d_layer = None
    ax_3d_layer = None
    if draw_layer_plot:
        fig_layer = plt.figure(figsize=(14, 6) if include_layer_3d else (8, 6))
        ax_2d_layer = fig_layer.add_subplot(1, 2, 1) if include_layer_3d else fig_layer.add_subplot(1, 1, 1)
        if include_layer_3d:
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
    if draw_layer_plot:
        plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(current_ext_print))),ax=ax_2d_layer,add_points=False,color="blue")
    final_perimeter_raw = current_ext_print
    final_perimeter = do_transform(final_perimeter_raw)


    layer_data[i]["internal_peris"] = []
    for wall_number in range(1,slicing_params["wall_count"]):
        if has_full_wrap_component(current_ext_print):
            internal_source = periodic_copies(current_ext_print, turns=1)
            _internal_source_metric = repair_polygonal_geometry(do_transform(internal_source))
            internal_perimeter = robust_polygonal_buffer(
                _internal_source_metric,
                -slicing_params["nozzle_size"] * wall_number,
                join_style="mitre",
            )
            internal_perimeter = robust_polygonal_intersection(internal_perimeter, _internal_source_metric)
            internal_perimeter_raw = do_transform_inverse(internal_perimeter)
            internal_perimeter_raw = fold_to_canonical(internal_perimeter_raw, theta_clip_box)
            internal_perimeter_raw = _canonicalize_full_wrap_strips(polygonal_only(internal_perimeter_raw))
        else:
            _internal_source_metric = repair_polygonal_geometry(do_transform(current_ext_print))
            internal_perimeter = robust_polygonal_buffer(
                _internal_source_metric,
                -slicing_params["nozzle_size"] * wall_number,
                join_style="mitre",
            )
            internal_perimeter = robust_polygonal_intersection(internal_perimeter, _internal_source_metric)
            internal_perimeter_raw = _canonicalize_full_wrap_strips(do_transform_inverse(internal_perimeter))
        if draw_layer_plot:
            plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(internal_perimeter_raw))),ax=ax_2d_layer,add_points=False,color="green")
        layer_data[i]["internal_peris"].append(internal_perimeter_raw)
        final_perimeter_raw = internal_perimeter_raw
        final_perimeter = do_transform(final_perimeter_raw)



    final_top_surface = robust_polygonal_intersection(
        robust_polygonal_buffer(
            do_transform(periodic_copies_for_bounds(current_top_print, final_perimeter_raw)),
            slicing_params["nozzle_size"] * slicing_params["wall_count"],
            join_style="mitre",
        ),
        final_perimeter,
    )
    if draw_layer_plot:
        plot_polygonal_2d(final_top_surface, ax=ax_2d_layer, color="indigo")
    _solid_direction = "axial"
    _min_solid_lines = slicing_params.get("surface_min_solid_lines", 0)

    def _surface_solid_infill(surface_raw):
        surface_raw = polygonal_only(surface_raw)
        if surface_raw.is_empty:
            return []
        parts = list(surface_raw.geoms) if hasattr(surface_raw, "geoms") else [surface_raw]
        result = []
        for part in parts:
            part_solid = do_solid_infill(part,slicing_params,current_r,draw_layer_plot,ax=ax_2d_layer,direction=_solid_direction) or []
            if 0 < len(part_solid) < _min_solid_lines:
                continue
            result.extend(part_solid)
        return result

    final_top_surface_raw = _canonicalize_full_wrap_strips(do_transform_inverse(final_top_surface))
    top_solid_infill = _surface_solid_infill(final_top_surface_raw)

    final_bottom_surface = robust_polygonal_intersection(
        robust_polygonal_buffer(
            do_transform(periodic_copies_for_bounds(current_bottom_print, final_perimeter_raw)),
            slicing_params["nozzle_size"] * slicing_params["wall_count"],
            join_style="mitre",
        ),
        final_perimeter,
    )
    if draw_layer_plot:
        plot_polygonal_2d(final_bottom_surface, ax=ax_2d_layer, color="cyan")
    final_bottom_surface_raw = _canonicalize_full_wrap_strips(do_transform_inverse(final_bottom_surface))
    bottom_solid_infill = _surface_solid_infill(final_bottom_surface_raw)
    layer_data[i]["solid_infill"] = top_solid_infill + bottom_solid_infill

    # Support infill: interface layers (near model) alternate radial/axial direction across
    # support_interface_layers shells; body layers use the body pattern.
    support_geom = layer_data[i].get("support", shapely.GeometryCollection())
    if not support_geom.is_empty:
        interface_n = slicing_params.get("support_interface_layers", 2)
        gap_n = slicing_params.get("support_gap_layers", 0)
        support_print = keep_canonical_polygons(support_geom)

        # Check distance to upper model boundary (i+d) and lower model boundary (i-d) separately.
        # For bottom interface, depth 1 is the innermost layer (printed first) → radial.
        # For top interface, depth 1 is the outermost layer (printed last) → parity is flipped
        # so the first-printed layer is still radial.
        def _iface_depth(delta_sign):
            for d in range(1, interface_n + 1):
                ext = layer_data.get(i + delta_sign * d, {}).get("external_peris", shapely.GeometryCollection())
                if not ext.is_empty:
                    candidate = robust_polygonal_intersection(support_print, keep_canonical_polygons(ext))
                    if not candidate.is_empty:
                        return d
            return None

        top_depth = _iface_depth(+1)   # distance to model above (outward)
        bot_depth = _iface_depth(-1)   # distance to model below (inward)

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

            layer_data[i]["support_infill"] = do_support_infill(support_print, slicing_params, current_r, direction=direction)

        if draw_layer_plot:
            plot_polygonal_2d(do_transform(support_geom), ax=ax_2d_layer, color="magenta")
    else:
        layer_data[i]["support_infill"] = []

    if "needs_support" in layer_data[i] and draw_layer_plot:
        plot_polygonal_2d(do_transform(layer_data[i]["needs_support"]), ax=ax_2d_layer, color="red")

    # plt.savefig(f"frames/{i}.png",dpi=600)

    # Inset by half a nozzle so sparse infill lines land just inside the innermost
    # perimeter without overlapping it (its center line is at the boundary of final_perimeter).
    _half_noz = slicing_params["nozzle_size"] * 0.5
    if has_full_wrap_component(final_perimeter_raw):
        # Keep the artificial cut away from theta=0 while preparing sparse infill.
        # The cube in cubeandring sits on theta=0, so canonical buffering/opening
        # creates a visible seam gap through the cube interior.
        _periodic_final = periodic_copies(final_perimeter_raw, turns=1)
        _infill_boundary = robust_polygonal_buffer(
            do_transform(_periodic_final),
            -_half_noz,
            join_style="mitre",
        )
        _diff = robust_polygonal_difference(
            robust_polygonal_difference(
                _infill_boundary,
                do_transform(periodic_copies(final_top_surface_raw, turns=1)),
            ),
            do_transform(periodic_copies(final_bottom_surface_raw, turns=1)),
        )
        if not _diff.is_empty:
            _diff = robust_polygonal_buffer(_diff, -_half_noz, join_style="mitre")
            _diff = robust_polygonal_buffer(_diff, _half_noz, join_style="mitre")
        _diff_raw = do_transform_inverse(_diff)
        shape_to_infill = _canonicalize_full_wrap_strips(
            polygonal_only(
                fold_to_theta_window(
                    _diff_raw,
                    math.pi,
                    math.pi + 2 * math.pi,
                    z_clip_min,
                    z_clip_max,
                )
            )
        )
    else:
        _infill_boundary = robust_polygonal_buffer(final_perimeter, -_half_noz, join_style="mitre")
        _diff = robust_polygonal_difference(
            robust_polygonal_difference(_infill_boundary, final_top_surface),
            final_bottom_surface,
        )
        if not _diff.is_empty:
            # Morphological opening: removes any remaining floating-point slivers
            # narrower than half a nozzle width.
            _diff = robust_polygonal_buffer(_diff, -_half_noz, join_style="mitre")
            _diff = robust_polygonal_buffer(_diff, _half_noz, join_style="mitre")
        shape_to_infill = _canonicalize_full_wrap_strips(do_transform_inverse(_diff))

    fully_radial_infill, radial_infill, axial_infill = do_sparse_infill(shape_to_infill,slicing_params,current_r,draw_layer_plot,ax=ax_2d_layer)
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

    # --- Per-layer 3D (interactive or selected debug exports) ---
    if draw_layer_plot and include_layer_3d:
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

    # --- Finalize and optionally save/show per-layer figure ---
    if draw_layer_plot:
        plot_line_endpoints_2d(layer_data[i].get("solid_infill", []), ax_2d_layer, current_r, marker_size=6)
        plot_line_endpoints_2d(layer_data[i].get("radial_infill", []), ax_2d_layer, current_r, marker_size=6)
        plot_line_endpoints_2d(layer_data[i].get("axial_infill", []), ax_2d_layer, current_r, marker_size=6)
        plot_line_endpoints_2d(layer_data[i].get("fully_radial_infill", []), ax_2d_layer, current_r, marker_size=6)
        plot_line_endpoints_2d(layer_data[i].get("support_infill", []), ax_2d_layer, current_r, marker_size=5)
        ax_2d_layer.set_title(f'Layer {i}')
        ax_2d_layer.set_aspect('equal')
        fig_layer.tight_layout()
        fig_layer.canvas.draw()
        if export_debug_layer:
            debug_path = os.path.join(args.debug_dir, f"layer_{i:04d}.png")
            fig_layer.savefig(debug_path, dpi=200)
            print(f"\nSaved debug plot: {debug_path}")
        if export_debug_layer and defer_single_debug_layer_plots and not args.no_plot:
            fig_layer.canvas.draw_idle()
        elif (plot_per_layer or args.show_debug_layers) and not args.no_plot:
            print(f"Close layer {i} window to continue...")
            while plt.fignum_exists(fig_layer.number):
                plt.pause(0.1)
        else:
            plt.close(fig_layer)

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

    debug_emitted_segments = []
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
        _layer_has_periodic_external = bool(layer_data[i].get("external_peris_has_periodic", False))
        def _print_window_for_paths(geom):
            geom = keep_canonical_polygons(geom)
            if _layer_has_periodic_external or has_full_wrap_component(geom):
                lower = math.pi
                upper = lower + 2 * math.pi
                return (
                    fold_to_theta_window(periodic_copies(geom, turns=1), lower, upper, z_clip_min, z_clip_max),
                    lower,
                    upper,
                    lower,
                )
            return geom, 0.0, 2 * math.pi, None

        def _prep_section(subsection):
            if subsection == "external_peris":
                if _layer_has_periodic_external or has_full_wrap_component(current_ext_print):
                    geom = layer_data[i].get("external_peris_print", current_ext_print)
                    lower = math.pi
                    upper = lower + 2 * math.pi
                    cut = lower
                    is_full_wrap_section = True
                else:
                    geom, lower, upper, cut = _print_window_for_paths(current_ext_print)
                    is_full_wrap_section = cut is not None or has_full_wrap_component(geom, lower=lower, upper=upper)
                lines = polygon_to_wrapped_lines(geom, lower=lower, upper=upper, cut_lower=cut)
                lines = _remove_periodic_seam_caps(lines, seam_origin=cut if cut is not None else lower)
                if is_full_wrap_section:
                    # In the shifted full-wrap print window, 0 mod 2pi lies in
                    # the middle of the path window. Only clean the actual cut.
                    if cut is not None:
                        lines = _remove_periodic_seam_aligned_edges(
                            lines,
                            seam_origins=(cut,),
                            seam_tol=0.02,
                            max_theta_span=1e-6,
                        )
                    lines = _remove_local_orthogonal_caps(lines)
                lines = _split_local_seam_closure_edges(lines)
                raw = shapely_to_paths(lines, current_r)
                return [p for p in raw if p is not None]
            elif subsection == "internal_peris":
                raw = []
                if _layer_has_periodic_external or has_full_wrap_component(current_ext_print):
                    lower = math.pi
                    upper = lower + 2 * math.pi
                    source = layer_data[i].get("external_peris_print", current_ext_print)
                    source_metric = repair_polygonal_geometry(do_transform(source))
                    for wall_number in range(1, slicing_params["wall_count"]):
                        internal_metric = robust_polygonal_buffer(
                            source_metric,
                            -slicing_params["nozzle_size"] * wall_number,
                            join_style="mitre",
                        )
                        internal_metric = robust_polygonal_intersection(internal_metric, source_metric)
                        internal_raw = do_transform_inverse(internal_metric)
                        lines = polygon_to_wrapped_lines(
                            internal_raw,
                            lower=lower,
                            upper=upper,
                            cut_lower=lower,
                        )
                        lines = _remove_periodic_seam_caps(lines, seam_origin=lower)
                        # The canonical 0/2pi line is valid printable geometry
                        # here; only the shifted window cut is synthetic.
                        lines = _remove_periodic_seam_aligned_edges(
                            lines,
                            seam_origins=(lower,),
                            seam_tol=0.02,
                            max_theta_span=1e-6,
                        )
                        lines = _remove_local_orthogonal_caps(lines)
                        lines = _split_local_seam_closure_edges(lines)
                        raw.extend(
                            shapely_to_paths(lines, current_r)
                        )
                else:
                    for each in layer_data[i].get(subsection, []):
                        geom, lower, upper, cut = _print_window_for_paths(each)
                        lines = polygon_to_wrapped_lines(geom, lower=lower, upper=upper, cut_lower=cut)
                        lines = _remove_periodic_seam_caps(lines, seam_origin=cut if cut is not None else lower)
                        lines = _split_local_seam_closure_edges(lines)
                        raw.extend(
                            shapely_to_paths(lines, current_r)
                        )
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
                # For fully_radial_infill (full circles), rotate to start at current theta
                # to avoid unnecessary travel back to theta=0
                if subsection == "fully_radial_infill":
                    path = rotate_full_circle_path(path, last_pos)
                # For perimeters (external/internal), rotate to start at nearest point
                # Check by subsection name since paths from polygon_to_wrapped_lines are open LineStrings
                elif subsection in ("external_peris", "internal_peris"):
                    path = rotate_closed_path_to_nearest(path, last_pos)
                emit_paths = (
                    _split_path_local_seam_closure_edges(path)
                    if subsection in ("external_peris", "internal_peris")
                    else [path]
                )
                for emit_path in emit_paths:
                    f.write(f";TYPE:{section_name}\n")
                    _path_start_e = last_e
                    g,last_pos,last_e = emit_path.to_gcode(_path_params,last_pos,last_e)
                    if export_debug_layer:
                        debug_emitted_segments.extend(
                            extract_extrusion_segments_from_gcode(g, _path_start_e, section_name)
                        )
                    f.write(g)

    if export_debug_layer:
        save_generated_lines_layer_plot(i, current_r, debug_emitted_segments)

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

    if not keep_layer_temp:
        # These geometries are only needed while emitting this layer. Keeping them
        # for the full run can dominate memory on dense models.
        for _temp_key in (
            "internal_peris",
            "solid_infill",
            "support_infill",
            "fully_radial_infill",
            "radial_infill",
            "axial_infill",
        ):
            layer_data[i].pop(_temp_key, None)

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
