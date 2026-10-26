"""Conservative variable-width centerline fill for perimeter-free narrow regions.

Only uses straight-ish polygons thinner than the nominal line width. This is
not a full medial-axis/Arachne contour solver; rejecting ambiguous geometry is
safer than filling across holes or silently overextruding sparse interiors.
All geometry is measured in the physical theta*radius / axial millimeter plane.
"""
from __future__ import annotations

import math

from shapely.geometry import LineString, Polygon, MultiPolygon, GeometryCollection
from shapely.ops import unary_union
from shapely import affinity

from .models import Toolpath


def _polygons(geom):
    if isinstance(geom, Polygon):
        yield geom
    elif isinstance(geom, (MultiPolygon, GeometryCollection)):
        for part in geom.geoms:
            yield from _polygons(part)


def _lines(geom):
    if isinstance(geom, LineString):
        yield geom
    elif hasattr(geom, "geoms"):
        for part in geom.geoms:
            yield from _lines(part)


def _compact(points, widths, radius):
    """Merge only equal-width collinear segments, preserving path length and E."""
    if len(points) < 3:
        return points, widths
    out_points = [points[0], points[1]]
    out_widths = [widths[0]]
    for point, width in zip(points[2:], widths[1:]):
        p0, p1 = out_points[-2:]
        ab = ((p1[0]-p0[0])*radius, p1[1]-p0[1])
        bc = ((point[0]-p1[0])*radius, point[1]-p1[1])
        cross = abs(ab[0]*bc[1]-ab[1]*bc[0])
        if (abs(width-out_widths[-1]) <= 1e-5
                and cross < 1e-8*max(1,math.hypot(*ab),math.hypot(*bc))
                and ab[0]*bc[0]+ab[1]*bc[1] > 0):
            out_points[-1] = point
        else:
            out_points.append(point)
            out_widths.append(width)
    return out_points, out_widths


def _narrow_polygon_paths(poly, radius, nominal_width, minimum_width, sample_step,
                          minimum_length):
    """Sample local perpendicular chord widths along a strip's long axis."""
    if not isinstance(poly, Polygon) or poly.is_empty or poly.interiors:
        return []
    rectangle = poly.minimum_rotated_rectangle
    if not isinstance(rectangle, Polygon):
        return []
    corners = list(rectangle.exterior.coords)[:4]
    lengths = [math.dist(corners[i],corners[(i+1)%4]) for i in range(4)]
    major_index = max(range(4),key=lambda i:lengths[i])
    major, minor = lengths[major_index], min(lengths)
    if (major < max(minimum_length, nominal_width*2.5)
            or minor > nominal_width*1.02
            or minor < minimum_width*0.45):
        return []
    edge_start,edge_end = corners[major_index],corners[(major_index+1)%4]
    ux,uy = (edge_end[0]-edge_start[0])/major, (edge_end[1]-edge_start[1])/major
    vx,vy = -uy,ux
    coords = list(poly.exterior.coords)
    low=min(x*ux+y*uy for x,y in coords)
    high=max(x*ux+y*uy for x,y in coords)
    count=max(2,math.ceil((high-low)/sample_step))
    # Bound pathological extremely long strips: controlled fallback, never
    # create tens of thousands of motion segments from one text/serif edge.
    if count > 2048:
        return []
    samples = []
    for j in range(count):
        t=low+(j+.5)*(high-low)/count
        extent=major+minor+2
        cross=LineString([
            (ux*t-vx*extent,uy*t-vy*extent),
            (ux*t+vx*extent,uy*t+vy*extent),
        ])
        segments=[seg for seg in _lines(poly.intersection(cross)) if seg.length>1e-8]
        # Cross-sections with multiple nontrivial intervals are branched/
        # ambiguous and must not become a connecting extrusion over a void.
        segments.sort(key=lambda segment:segment.length,reverse=True)
        if len(segments)>1 and segments[1].length>minimum_width*.3:
            samples.append(None)
            continue
        if not segments:
            samples.append(None)
            continue
        line=segments[0]
        width=line.length
        if width < minimum_width or width > nominal_width*1.02:
            samples.append(None)
            continue
        midpoint=line.interpolate(.5,normalized=True)
        samples.append(((midpoint.x,midpoint.y),min(nominal_width,width)))

    paths=[]
    current=[]
    widths=[]
    def emit():
        if len(current)<2:
            return
        # Local chords and centerline motion must remain in the original
        # polygon; no outward stroke or gap-spanning diagonal is permitted.
        pts=[(x/radius,y) for x,y in current]
        length=sum(math.dist(a,b) for a,b in zip(current,current[1:]))
        if length < minimum_length:
            return
        pts,segment_widths=_compact(pts,widths,radius)
        if len(pts)>=2 and len(segment_widths)==len(pts)-1:
            paths.append(Toolpath("GAP_FILL",radius,pts,
                                  segment_widths=segment_widths))

    for item in samples:
        if item is None:
            emit()
            current=[]
            widths=[]
            continue
        point,width=item
        if not current:
            current=[point]
            continue
        previous=current[-1]
        proposed=0.5*(width+max(minimum_width, samples[samples.index(item)-1][1])) if False else width
        # The narrower of the neighboring local measurements bounds the
        # rectangular swept track conservatively; segment E uses this width.
        # Explicitly evaluate the preceding sample rather than assuming
        # a constant nozzle-width bead.
        # A width of zero never prints. Reject an entire suspect segment.
        segment=LineString([previous,point])
        if segment.length < 1e-8:
            continue
        effective=min(width, previous_width)
        sweep=segment.buffer(effective*.5, cap_style="flat",join_style="mitre")
        outside=sweep.difference(poly)
        if outside.area > max(1e-8, segment.length*.006):
            emit()
            current=[point]
            widths=[]
            previous_width=width
            continue
        current.append(point)
        widths.append(effective)
        previous_width=width
    emit()
    return paths


def narrow_gap_fill_paths(model, perimeter_paths, radius, config):
    """Recover thin edge features missing when nominal-width walls disappear.

    Model-space candidate strips are restricted to the first half-bead edge
    band and subtract already printed perimeter footprints. This prevents
    accidentally turning intentionally empty sparse-infill cells into solid
    plastic, and preserves actual polygon holes.
    """
    if model is None or model.is_empty or radius<=0:
        return []
    physical=affinity.scale(model,xfact=radius,yfact=1,origin=(0,0))
    width=config.nozzle_width
    # Only the edge band is eligible: the interior may be deliberately sparse.
    core=physical.buffer(-width*.5)
    edge=physical.difference(core) if not core.is_empty else physical
    footprints=[]
    for path in perimeter_paths:
        if len(path.points)<2:
            continue
        line=LineString([(theta*radius,axial) for theta,axial in path.points])
        if line.length>1e-8:
            footprints.append(line.buffer(width*.5,cap_style="flat",join_style="mitre"))
    if footprints:
        edge=edge.difference(unary_union(footprints))
    if edge.is_empty:
        return []
    result=[]
    for poly in _polygons(edge):
        result.extend(_narrow_polygon_paths(
            poly,radius,width,config.gap_fill_min_width_mm,
            config.gap_fill_sample_step_mm,
            max(config.min_path_length,config.gap_fill_min_length_mm),
        ))
    return result
