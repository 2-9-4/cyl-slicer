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
from togcode import shapely_to_paths
from plotutil import set_axes_equal

from vispy import scene
from vispy.scene import Line
from vispy import app, gloo
import solid2
def find_polygons(edges, num_points):
    # Build undirected edge map
    adjacency = defaultdict(list)
    degree = [0] * num_points

    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
        degree[a] += 1
        degree[b] += 1

    # Check for dead ends: all vertices must have even degree
    for i, d in enumerate(degree):
        if d % 2 != 0:
            raise ValueError(f"Deadend or unmatched edge at vertex {i}: degree={d}")

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
            prev = start
            visited_edges.add(ek)

            while curr != start:
                poly.append(curr)
                for next_node in adjacency[curr]:
                    if next_node == prev:
                        continue
                    ek2 = edge_key(curr, next_node)
                    if ek2 not in visited_edges:
                        visited_edges.add(ek2)
                        prev, curr = curr, next_node
                        break
                else:
                    raise ValueError(f"Dead-end or cycle break at vertex {curr}")

            polygons.append(poly)

    return polygons

def shift_polygon(poly, dx):
    return shapely.transform(poly,lambda x: x + [dx,0])

def build_hierarchy(polygons, period=2 * np.pi):

    hierarchy = [[] for _ in polygons]
    for i, outer_poly in enumerate(polygons):
        outer_shifts = [
            shift_polygon(outer_poly, 0),
            shift_polygon(outer_poly, -period),
            shift_polygon(outer_poly, +period)
        ]
        for j, inner_poly in enumerate(polygons):
            if i == j:
                continue
            # If the inner poly is contained in any of the shifted variants of the outer
            if any(outer.contains(inner_poly) for outer in outer_shifts):
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

def do_solid_infill(peri,slicing_params,current_r,plot_infill=False):
    if peri.is_empty:
        return
    def do_transform(x):
        return shapely.transform(x,lambda p: p * [current_r,1])
    def do_transform_inverse(x):
        return shapely.transform(x,lambda p: p * [1/current_r,1])
    infill = []
    bounds = shapely.bounds(peri)

    for current_line_h in np.arange(0,2*np.pi,slicing_params["nozzle_size"]/current_r):
        lp_x,lp_y = [current_line_h,current_line_h],[bounds[1],bounds[3]]
        ls = shapely.LineString(np.array([lp_x,lp_y]).swapaxes(0,1))
        inter = shapely.intersection(ls,peri)
        geoms = None
        if type(inter) == shapely.MultiLineString:
            geoms = inter.geoms
        elif type(inter) == shapely.geometry.collection.GeometryCollection:
            geoms = inter.geoms
        if geoms is not None:
            for line_segment in geoms:
                if plot_infill:
                    plot_line(do_transform(line_segment),color="greenyellow",add_points=False)
                infill.append(line_segment)
        else:
            if(type(inter) == shapely.geometry.multipoint.MultiPoint):
                print("got MultiPoint in solid infill")
                shapely.plotting.plot_points(do_transform(inter),color="greenyellow")
            else:
                if plot_infill:
                    plot_line(do_transform(inter),color="greenyellow",add_points=False)
                infill.append(inter)

    return infill

def do_sparse_infill(innermost_peri,slicing_params,current_r,plot_infill=False):
    if innermost_peri.is_empty:
        return [],[],[]
    def do_transform(x):
        return shapely.transform(x,lambda p: p * [current_r,1])
    def do_transform_inverse(x):
        return shapely.transform(x,lambda p: p * [1/current_r,1])
    fully_radial_infill = []
    radial_infill = []
    axial_infill = []
    lowest = shapely.bounds(innermost_peri)[1]
    lowest /= slicing_params["infill_spacing"]
    lowest = np.floor(lowest) * slicing_params["infill_spacing"]
    highest = shapely.bounds(innermost_peri)[3]
    highest /= slicing_params["infill_spacing"]
    highest = np.ceil(highest) * slicing_params["infill_spacing"]

    for current_line_v in np.arange(lowest,highest+1,slicing_params["infill_spacing"]):
        lp_x,lp_y = [0,2*np.pi],[current_line_v,current_line_v]
        ls = shapely.LineString(np.array([lp_x,lp_y]).swapaxes(0,1))
        inter = shapely.intersection(ls,innermost_peri)
        geoms = None
        if type(inter) == shapely.MultiLineString:
            geoms = inter.geoms
        elif type(inter) == shapely.geometry.collection.GeometryCollection:
            geoms = inter.geoms
        if geoms is not None:
            for line_segment in geoms:
                if plot_infill:
                    plot_line(do_transform(line_segment),color="orange",add_points=False)
                radial_infill.append(line_segment)
        else:
            if(type(inter) == shapely.geometry.multipoint.MultiPoint):
                print("got MultiPoint")
                shapely.plotting.plot_points(do_transform(inter),color="red")
            else:
                if plot_infill:
                    plot_line(do_transform(inter),color="red",add_points=False)
                fully_radial_infill.append(inter)
    for current_line_h in np.arange(0,2*np.pi,slicing_params["infill_spacing"]/current_r):
        lp_x,lp_y = [current_line_h,current_line_h],[lowest,highest]
        ls = shapely.LineString(np.array([lp_x,lp_y]).swapaxes(0,1))
        inter = shapely.intersection(ls,innermost_peri)
        geoms = None
        if type(inter) == shapely.MultiLineString:
            geoms = inter.geoms
        elif type(inter) == shapely.geometry.collection.GeometryCollection:
            geoms = inter.geoms
        
        if geoms is not None:
            for line_segment in geoms:
                if plot_infill:
                    plot_line(do_transform(line_segment),color="darkviolet",add_points=False)
                axial_infill.append(line_segment)
        else:
            if plot_infill:
                plot_line(do_transform(inter),color="darkviolet",add_points=False)
            axial_infill.append(inter)
    return fully_radial_infill, radial_infill, axial_infill
from shapely.geometry import Polygon, LineString, box, MultiLineString
from shapely.ops import linemerge

def polygon_to_wrapped_lines(poly: Polygon, lower=0.0, upper=2*math.pi):
    """
    Cut polygon at x=lower and x=upper, return open LineStrings for exterior
    and interiors. If an interior gets clipped, the pieces are wrapped
    around the seam and stitched together.
    """
    if isinstance(poly,list):
        ret = []
        for each in poly:
            ret += polygon_to_wrapped_lines(each,lower,upper)
        return ret
    elif isinstance(poly,shapely.MultiPolygon):
        ret = []
        for each in poly.geoms:
            ret += polygon_to_wrapped_lines(each,lower,upper)
        return ret
    
    if not isinstance(poly, Polygon):
        print(poly)
        print(type(poly))
        raise ValueError("Expected a shapely Polygon")

    clip_rect = box(lower, -1e12, upper, 1e12)
    result = []

    # --- exterior (no wrapping) ---
    boundary_line = LineString(poly.exterior.coords)
    clipped = boundary_line.intersection(clip_rect)


    if isinstance(clipped, LineString):
        coords = np.array(clipped.coords)
        if(np.any(coords[:,0] == 0)):
            coords[:,0] += 2*np.pi
        result.append(LineString(coords))
    elif not clipped.is_empty:
        for g in clipped.geoms:
            coords = np.array(g.coords)
            if(np.any(coords[:,0] == 0)):
                coords[:,0] += 2*np.pi
            result.append(LineString(coords))
        # result.extend([g for g in clipped.geoms if isinstance(g, LineString)])

    # --- interiors (with wrapping) ---
    for ring in poly.interiors:
        hole_line = LineString(ring.coords)
        clipped = hole_line.intersection(clip_rect)

        if clipped.is_empty:
            continue

        # convert to list of LineStrings
        if isinstance(clipped, LineString):
            lines = [clipped]
        else:
            lines = [g for g in clipped.geoms if isinstance(g, LineString)]

        wrapped_lines = []

        for line in lines:
            coords = np.array(line.coords)
            if(np.any(coords[:,0] == 0)):
                coords[:,0] += 2*np.pi

            wrapped_lines.append(LineString(coords))

        # stitch together the wrapped pieces
        merged = linemerge(MultiLineString(wrapped_lines))
        if isinstance(merged, LineString):
            result.append(merged)
        else:  # sometimes linemerge returns MultiLineString
            result.extend([l for l in merged.geoms if isinstance(l, LineString)])

    return result
    
def sum_poly(p):
    if isinstance(p, Polygon):
        return p 
    f = p[0]
    for e in p[1:]:
        f = f.union(e)
    return f


def remove_collinear_points_angle(coords, angle_tolerance=1e-6):
    """
    Remove collinear points using angle change method.
    More intuitive tolerance setting.
    
    Args:
        coords: List/array of (x, y) coordinates
        angle_tolerance: Tolerance in radians for angle change
    
    Returns:
        Simplified list of coordinates
    """
    # plt.plot(coords[:,0],coords[:,1],"o-")
    coords = np.asarray(coords)
    if len(coords) <= 2:
        return coords.tolist()
    
    keep = [True]  # Always keep first point
    last_kept = coords[0]
    for i in range(1, len(coords) - 1):
        p1 = last_kept
        p2 = coords[i]
        p3 = coords[i+1]
        
        # Calculate vectors
        v1 = p2 - p1
        v2 = p3 - p2
        
        # Calculate angle between vectors
        # Handle zero-length vectors
        norm1 = np.linalg.norm(v1)
        norm2 = np.linalg.norm(v2)
        
        if norm1 < 1e-10 or norm2 < 1e-10:
            keep.append(False)  # Remove duplicate/near-duplicate points
            continue
            
        # Normalize vectors
        v1_norm = v1 / norm1
        v2_norm = v2 / norm2

        # Calculate angle change (0 = straight line, π/2 = right angle, etc.)
        dot_product = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
        angle_change = np.arccos(abs(dot_product))
        # Keep point if angle change is significant
        keep.append(angle_change > angle_tolerance)
        if keep[-1]:
            last_kept = coords[i]
    
    keep.append(True)  # Always keep last point
    # plt.plot(coords[keep,0],coords[keep,1],"o-")
    # plt.gca().set_aspect('equal') 
    # plt.title(f"{angle_tolerance*180/np.pi:.2f} degrees")
    # plt.show()
    return coords[keep].tolist()


slicing_params = {
    "layer_height": 0.2,
    "infill_spacing" : 2,
    "nozzle_size": 0.4,
    "wall_count": 3,
    "feed_mm_per_min" : 6000,
    "travel_mm_per_min" : 6000,
    "extrusion_per_mm" : 0.05
}
max_real_layers = 8  # Limit for testing
mesh = trimesh.load_mesh("test.stl")
pts = cartesian_to_cylindrical(mesh.vertices)
max_z = np.max(pts[:,2])
min_z = np.min(pts[:,2])
part_height = (max_z-min_z)


# for current_r in np.arange(max_r-layer_height*35,max_r,layer_height):
layer_data = {}
sliced_layer_count = 0
for layer_i,current_r in enumerate( np.arange(np.min(pts[:,0])+slicing_params["layer_height"],np.max(pts[:,0])-1e-3,slicing_params["layer_height"])):
    # Limit slicing to 8 layers for testing
    if sliced_layer_count >= max_real_layers:
        print(f"\nSliced {sliced_layer_count} layers, stopping.")
        break
    
    times = []
    sections = 1000
    cyl = trimesh.creation.cylinder(
        radius=current_r,
        height=part_height*1.2,
        sections=sections
    )
    cyl.apply_translation([0,0,part_height*0.5+min_z])
    assert(cyl.is_volume)
    times.append((time.perf_counter(),""))
    good = False
    print(f"{current_r/(np.max(pts[:,0])-1e-3)*100:.2f}%")
    while not good:
        cut_mesh = cyl.intersection(mesh,engine="manifold")
        cut_mesh = fix_mesh(cut_mesh)
        good = cut_mesh.is_volume
        if not good:
            #this bs makes intersections behave better
            current_r+=np.random.uniform(-1e-3,1e-3)
            cyl = trimesh.creation.cylinder(
                radius=current_r,
                height=part_height*1.2,
                sections=sections
            )
            cyl.apply_translation([0,0,part_height*0.5+min_z])
            cyl = cyl.apply_transform(trimesh.transformations.rotation_matrix(np.random.uniform(0,2*np.pi/sections), [0,0,1]))

            print(f"retrying intersection on {layer_i}")

    layer_data[layer_i] = {"radius": current_r}
    sliced_layer_count += 1
    

    times.append((time.perf_counter(),"intersection"))
    assert( len(cut_mesh.vertices))
    layer_pts = cartesian_to_cylindrical(cut_mesh.vertices)
    mask = np.argwhere(np.isclose(layer_pts[:,0],current_r,1e-5)).squeeze()
    if(len(mask) < 3):
        print(f"empty layer for {layer_i}")
        continue
    times.append((time.perf_counter(),"mask generation"))


    triangles = fast_extract_faces(cut_mesh,mask)
    triangles = [tri for tri in triangles if tri[0] != tri[1] and tri[0] != tri[2] and tri[1] != tri[2]]
    times.append((time.perf_counter(),"tri list"))

    outer_edges = find_boundary_edges(triangles)
    # try:
    poly_idxs = find_polygons(outer_edges,len(layer_pts))
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
    full_circle_pairs = find_connected_shapes([poly for i,poly in enumerate(poly_idxs) if full_circles[i]],triangles)
    full_circle_pairs = np.array([np.argwhere(full_circles)[each].squeeze() for each in full_circle_pairs])
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
    
    for poly_i,sub_polys in enumerate(hier):
        if levels[poly_i] % 2 == 1:
            continue
        # if not len(sub_polys):
        #     continue
        # plt.plot(np.array(polys[poly_i].exterior.coords)[:,0],np.array(polys[poly_i].exterior.coords)[:,1])
        # plt.show()
        def do_transform(x):
            return shapely.transform(x,lambda p: p * [current_r,1])
        def do_transform_inverse(x):
            return shapely.transform(x,lambda p: p * [1/current_r,1])
        sub_poly_list = [np.array(polys[sub_poly].exterior.coords).astype(float) for sub_poly in sub_polys if levels[sub_poly] == levels[poly_i]+1]
        sub_poly_list = [[each - [2*np.pi,0],each,each + [2*np.pi,0]] for each in sub_poly_list]
        sub_poly_list = [remove_collinear_points_angle(each2*np.array([current_r,1]),0.01)*np.array([1/current_r,1]) for each in sub_poly_list for each2 in each]
        poly_exterior = remove_collinear_points_angle(polys[poly_i].exterior.coords*np.array([current_r,1]),0.01)*np.array([1/current_r,1])
        if not len(sub_poly_list):
            poly_with_holes = Polygon(poly_exterior)
        else:
            poly_with_holes = Polygon(poly_exterior,holes = sub_poly_list)

        if plot_perimeters:
            plot_polygon(do_transform(poly_with_holes),add_points=False,facecolor=(0, 0, 0, 0))

        poly_with_holes = do_transform_inverse(do_transform(poly_with_holes).buffer(-slicing_params["nozzle_size"]/2,join_style="mitre"))
        external_perimeters.append(poly_with_holes)
        if poly_i not in mapping:
            external_perimeters.append(shift_polygon(poly_with_holes,2*np.pi))
            external_perimeters.append(shift_polygon(poly_with_holes,-2*np.pi))
        if plot_perimeters:
            plot_polygon(do_transform(external_perimeters[-1]),add_points=False,color="green",facecolor=(0, 0, 0, 0))
    times.append((time.perf_counter(),"form poly with holes"))
    
    print(mapping)
    final = external_perimeters[0]

    for pi,p in enumerate(external_perimeters[1:]):
        final = final.union(p)

    layer_data[layer_i]["external_peris"] = final
    times.append((time.perf_counter(),"full union"))

    if layer_i == 0:
        continue

    e = do_transform(final)

    e -= do_transform(layer_data[layer_i-1]["external_peris"]).buffer(slicing_params["nozzle_size"]/2)
    e -= do_transform(shift_polygon(layer_data[layer_i-1]["external_peris"],-np.pi*2)).buffer(slicing_params["nozzle_size"]/2)

    
    layer_data[layer_i]["needs_support"] = do_transform_inverse(e)
    times.append((time.perf_counter(),"support check"))

    e = do_transform(layer_data[layer_i-1]["external_peris"])
        
    e = e.buffer(-slicing_params["nozzle_size"],join_style="mitre")# shrink the area of the top layer so we don't flag slicing_params["wall_count"] as needing top layers, TODO tweak this metric slicing_params["nozzle_size"]*slicing_params["wall_count"] is probably too little overlap between slicing_params["wall_count"]
    e = do_transform_inverse(e)

    e -= final
    e -= shift_polygon(final,-np.pi*2)
    times.append((time.perf_counter(),"top surface check"))

    layer_data[layer_i-1]["top_surfaces"] = e
    for i in range(1,len(times)):
        print(f"{times[i][1]} took {times[i][0]-times[i-1][0]:.2f} seconds")




layer_data[list(layer_data)[-1]]["top_surfaces"] = layer_data[list(layer_data)[-1]]["external_peris"]

plot = True  # Disabled - use gcode_viewer_gl.py for visualization
plot_3d = False 
ax = None
if plot_3d:
    # canvas = scene.SceneCanvas(keys='interactive', bgcolor='white', show=True)
    # view = canvas.central_widget.add_view()
    # view.camera = 'turntable'  # 3D interactive camera
    ax = plt.figure().add_subplot(projection='3d')

# output_scad = solid2.cube()
gcodes = []
last_pos = [0,0]
last_e = 0
real_layer_count = 0
print(f"\n=== LIMITING TO {max_real_layers} REAL LAYERS FOR TESTING ===\n")
with open("output.gcode","w") as f:
    f.write(";start of slice\n")
for i in sorted(layer_data.keys()):
    # Exit after 8 real layers for testing
    real_layer_count += 1
    if real_layer_count > max_real_layers:
        print(f"\nReached {max_real_layers} real layers, exiting.")
        break
    
    print(f"\n>>> Processing real layer {real_layer_count}/{max_real_layers} (layer_data key: {i})")
    
    current_r = layer_data[i]["radius"]
    def do_transform(x):
        return shapely.transform(x,lambda p: p * [current_r,1])
    def do_transform_inverse(x):
        return shapely.transform(x,lambda p: p * [1/current_r,1])
    if plot:
        plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(layer_data[i]["external_peris"]))),add_points=False,color="blue")
    final_perimeter = do_transform(layer_data[i]["external_peris"])
    layer_data[i]["internal_peris"] = []
    for wall_number in range(1,slicing_params["wall_count"]):
        internal_perimeter = do_transform(layer_data[i]["external_peris"]).buffer(-slicing_params["nozzle_size"]*(wall_number),join_style="mitre")
        if plot:
            plot_line(do_transform(sum_poly(polygon_to_wrapped_lines(do_transform_inverse(internal_perimeter)))),add_points=False,color="green")
        layer_data[i]["internal_peris"].append(do_transform_inverse(internal_perimeter))
        final_perimeter = internal_perimeter


    
    final_top_surface = do_transform(layer_data[i]["top_surfaces"]).buffer(slicing_params["nozzle_size"]*slicing_params["wall_count"],join_style="mitre").intersection(final_perimeter)
    if plot:
        plot_polygon(final_top_surface,add_points=False,color="indigo")
    layer_data[i]["solid_infill"] = do_solid_infill(do_transform_inverse(final_top_surface),slicing_params,current_r,plot)

    if "needs_support" in layer_data[i] and plot:
        plot_polygon(do_transform(layer_data[i]["needs_support"]),add_points=False,color="red")

    # plt.savefig(f"frames/{i}.png",dpi=600)

    shape_to_infill = do_transform_inverse(final_perimeter-final_top_surface)

    fully_radial_infill, radial_infill, axial_infill = do_sparse_infill(shape_to_infill,slicing_params,current_r,plot)
    layer_data[i]["fully_radial_infill"] = fully_radial_infill
    layer_data[i]["radial_infill"] = radial_infill
    layer_data[i]["axial_infill"] = axial_infill
    print(f"LAYER {i}")
    if plot:
        plt.show()

    traces = []
    # app.use_app('pyqt6')
    # Section names and colors for G-code comments and visualization
    sections = [
        ("external_peris", "EXTERNAL_PERIMETER"),
        ("internal_peris", "INTERNAL_PERIMETER"),
        ("solid_infill", "SOLID_INFILL"),
        ("fully_radial_infill", "INFILL_RADIAL_FULL"),
        ("radial_infill", "INFILL_RADIAL"),
        ("axial_infill", "INFILL_AXIAL")
    ]
    with open("output.gcode","a") as f:
        f.write(f";LAYER:{i}\n")
        f.write(f";Z:{current_r}\n")
        f.write(f"G0 Z{current_r}\n")

        for subsection, section_name in sections:
            if subsection == "external_peris" or subsection == "internal_peris":
                gcode_paths = shapely_to_paths(polygon_to_wrapped_lines(layer_data[i][subsection]), current_r)
            else:
                gcode_paths = shapely_to_paths(layer_data[i][subsection], current_r)
            for path in gcode_paths:
                if path is None:
                    continue
                f.write(f";TYPE:{section_name}\n")
                g,last_pos,last_e = path.to_gcode(slicing_params,last_pos,last_e)
                f.write(g)
    # for each
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
# plt.show()
    # plt.clf()
# if plot_3d:
#     canvas.app.run()

# mesh.show()
# output_scad.save_as_scad()
