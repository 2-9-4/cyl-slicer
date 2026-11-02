"""Conservative single-junction medial recovery for thin branched edge gaps.

Raster skeletons only nominate paths. Every extruded variable-width bead
sweep is checked against the exact polygon, and arms must not overlap.
"""
from __future__ import annotations
import itertools
import math
import numpy as np
import shapely
from scipy.ndimage import distance_transform_edt
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union


def branched_gap_paths(poly, radius, nominal_width, minimum_width,
                       sample_step, minimum_length):
    """Print a one-junction three- or four-arm tree with trimmed side arms."""
    if not isinstance(poly, Polygon) or poly.is_empty or poly.interiors:
        return []
    ring=list(poly.exterior.coords)[:-1]
    signed_area=sum(a[0]*b[1]-b[0]*a[1]
                    for a,b in zip(ring,ring[1:]+ring[:1]))
    if abs(signed_area)<1e-12:
        return []
    concave=0
    for a,b,c in zip(ring[-1:]+ring[:-1],ring,ring[1:]+ring[:1]):
        u=(b[0]-a[0],b[1]-a[1])
        v=(c[0]-b[0],c[1]-b[1])
        if (u[0]*v[1]-u[1]*v[0])*signed_area < -1e-9:
            concave+=1
    if concave<2:
        return []  # Ordinary strips avoid raster work.
    rectangle=poly.minimum_rotated_rectangle
    if not isinstance(rectangle,Polygon):
        return []
    corners=list(rectangle.exterior.coords)[:4]
    major=max(math.dist(corners[i],corners[(i+1)%4]) for i in range(4))
    if (major<max(minimum_length*2,nominal_width*3)
            or poly.area>nominal_width*major*2.2):
        return []
    pitch=min(minimum_width/4,sample_step/4,.07)
    if pitch<=0 or not math.isfinite(pitch):
        return []
    xmin,ymin,xmax,ymax=poly.bounds
    nx=int(math.ceil((xmax-xmin)/pitch))
    ny=int(math.ceil((ymax-ymin)/pitch))
    if nx*ny>350000 or nx<3 or ny<3:
        return []
    px=xmin+(np.arange(nx)+.5)*pitch
    py=ymin+(np.arange(ny)+.5)*pitch
    xx,yy=np.meshgrid(px,py)
    interior=shapely.contains_xy(poly,xx,yy)
    if not np.any(interior):
        return []
    dist=distance_transform_edt(interior)*pitch
    skeleton=np.pad(interior,1)
    # Two-phase Zhang-Suen thinning. Never use a partial skeleton.
    for _ in range(96):
        removed=0
        for phase in (0,1):
            core=skeleton[1:-1,1:-1]
            p2=skeleton[:-2,1:-1];p3=skeleton[:-2,2:]
            p4=skeleton[1:-1,2:];p5=skeleton[2:,2:]
            p6=skeleton[2:,1:-1];p7=skeleton[2:,:-2]
            p8=skeleton[1:-1,:-2];p9=skeleton[:-2,:-2]
            n=(p2,p3,p4,p5,p6,p7,p8,p9)
            count=sum(x.astype(np.int16) for x in n)
            changes=sum((~n[k]&n[(k+1)%8]).astype(np.int16)
                        for k in range(8))
            forbidden=((p2&p4&p6)|(p4&p6&p8) if phase==0
                       else (p2&p4&p8)|(p2&p6&p8))
            discard=core&(count>=2)&(count<=6)&(changes==1)&(~forbidden)
            removed+=int(np.count_nonzero(discard))
            core[discard]=False
        if not removed:
            break
    else:
        return []
    nodes=set(map(tuple,np.argwhere(skeleton)))
    if len(nodes)>12000:
        return []
    graph={}
    for r,c in nodes:
        near=set()
        for dy in (-1,0,1):
            for dx in (-1,0,1):
                candidate=(r+dy,c+dx)
                if (dx or dy) and candidate in nodes:
                    if dx and dy and ((r+dy,c) in nodes or (r,c+dx) in nodes):
                        continue  # Suppress diagonal shortcut at corner.
                    near.add(candidate)
        graph[(r,c)]=near
    junctions=[n for n,near in graph.items() if len(near)>2]
    ends=sorted(n for n,near in graph.items() if len(near)==1)
    if (len(junctions)!=1 or len(ends) not in (3,4)
            or len(graph[junctions[0]])!=len(ends)
            or any(len(near) not in (1,2,len(ends)) for near in graph.values())):
        return []  # Cyclic, multi-junction and ambiguous trees rejected.
    junction=junctions[0]
    arms=[]
    for endpoint in ends:
        path=[endpoint];prior=None;node=endpoint
        while node!=junction and len(path)<=len(nodes):
            nexts=[n for n in graph[node] if n!=prior]
            if len(nexts)!=1:
                return []
            prior,node=node,nexts[0]
            path.append(node)
        if node!=junction:
            return []
        arms.append(list(reversed(path)))
    if sum(len(path)-1 for path in arms)!=len(nodes)-1:
        return []  # Cyclic, disconnected or stray graph vertex.
    def xy(node):
        r,c=node
        return float(px[c-1]),float(py[r-1])
    pair=max(itertools.combinations(range(len(arms)),2),
             key=lambda ij:math.dist(xy(arms[ij[0]][-1]),xy(arms[ij[1]][-1])))
    trunk=list(reversed(arms[pair[0]]))+arms[pair[1]][1:]
    side_arms=[arms[i] for i in range(len(arms)) if i not in pair]
    def candidate(nodes):
        points=[xy(n) for n in nodes]
        widths=[]
        for (r0,c0),(r1,c1) in zip(nodes,nodes[1:]):
            # Conservative erosion of distance-transform widths; round down
            # to avoid extra material from pixel-edge quantization.
            raw=min(nominal_width,
                    2*min(dist[r0-1,c0-1],dist[r1-1,c1-1])-.8*pitch)
            w=math.floor((raw+1e-10)*100)/100
            if w<minimum_width-1e-9:
                return None
            widths.append(w)
        if not widths:
            return None
        sweeps=[LineString([a,b]).buffer(w*.5,cap_style="flat")
                for a,b,w in zip(points,points[1:],widths)]
        footprint=unary_union(sweeps)
        length=LineString(points).length
        if (length<minimum_length
                or footprint.difference(poly).area>max(1e-8,length*.001)):
            return None
        return points,widths,footprint
    main=candidate(trunk)
    if main is None:
        return []
    accepted=[main]
    for side in side_arms:
        choice=None
        # Width-aware footprint validation against every already printed
        # arm, not merely no-intersection of the path centerlines.
        for trim in range(1,len(side)-1):
            branch=candidate(side[trim:])
            if branch is None:
                continue
            if any(branch[2].intersection(other[2]).area>.001
                   for other in accepted):
                continue
            if (LineString(branch[0]).distance(LineString(main[0]))
                    >nominal_width*.8):
                continue
            choice=branch
            break
        if choice is None:
            return []
        accepted.append(choice)
    from .gap_fill import _compact
    from .models import Toolpath
    output=[]
    for points,widths,_ in accepted:
        coords=[(x/radius,y) for x,y in points]
        coords,widths=_compact(coords,widths,radius)
        output.append(Toolpath("GAP_FILL",radius,coords,segment_widths=widths))
    return output
