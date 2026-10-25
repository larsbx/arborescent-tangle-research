"""A5.2e2b4c2d2c: actual tetrahedral vertex stars at 2-valent K corners.

The previous edge-star kernels only covered proper *interiors* of K
segments. This constructs one genuine exact-Q³ triangulated convex
octahedral 3-ball whose INTERIOR VERTEX is an existing polygonal K
corner, with both adjacent source K directions as literal tetrahedral
mesh edges. Noncollinear bends and collinear straight continuations
have distinct guarded exact rational constructions.

The independent audit checks every oriented tetrahedron and PL vertex
link, TWO incident K edge links, geometric clearance from EVERY other
K segment, and strict clearance from all closed Conway cube solids.

This is a LOCAL knot-vertex 3-ball, not a glued edge/vertex tube atlas,
a whole-K conforming ambient, source Fig12.10 plumbing or R+/-.
"""
from __future__ import annotations

from collections import Counter,defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .exterior_conway_collars import (
    ExteriorCollaredBandCandidate,build_exterior_collared_band,
)
from .local_plumbing_slit import _volume6
from .rational_k_segment_stars import (
    _xyz,_minus,_plus,_scale,_cross,_dot,_frame,
)
from .spatial_twist_blocks import Point3

Vector=tuple[Fraction,Fraction,Fraction]
Tet=tuple[int,int,int,int]


def _generators(before:Point3,corner:Point3,after:Point3,
                radius:Fraction) -> tuple[Vector,Vector,Vector,int,int]:
    """Exact nondegenerate affine octahedral frame and marked K rays.

    For a corner bend, +u and +v are embedded in the two K edges.
    For a straight vertex, +u and -u are the two K rays, with a
    guarded pair of exact rational transverse basis vectors.
    """
    if (not isinstance(before,Point3) or
            not isinstance(corner,Point3) or
            not isinstance(after,Point3) or
            type(radius) is not Fraction or not 0<radius<=Fraction(1,4)):
        raise ValueError("typed K corner and exact positive radius required")
    d1=_minus(_xyz(before),_xyz(corner))
    d2=_minus(_xyz(after),_xyz(corner))
    if d1==(0,0,0) or d2==(0,0,0):
        raise ValueError("K vertex has collapsed incident segment")
    cross=_cross(d1,d2)
    if cross==(0,0,0):
        if _dot(d1,d2)>=0:
            raise ValueError("K arc doubles back on itself")
        _,v0,w0=_frame(corner,before)
        u,v,w=_scale(d1,radius),_scale(v0,radius),_scale(w0,radius)
        second=2
    else:
        u,v=_scale(d1,radius),_scale(d2,radius)
        normal=_scale(cross,1/max(abs(k) for k in cross))
        w=_scale(normal,radius)
        second=3
    if _dot(u,_cross(v,w))<=0:
        raise ValueError("nonpositive rational octahedral orientation")
    return u,v,w,1,second


def _vertices(corner:Point3,u:Vector,v:Vector,w:Vector
              ) -> tuple[Point3,...]:
    origin=_xyz(corner)
    return tuple(Point3(*point) for point in (
        origin,_plus(origin,u),_minus(origin,u),
        _plus(origin,v),_minus(origin,v),
        _plus(origin,w),_minus(origin,w)
    ))


def _tetrahedra(vertices:tuple[Point3,...]) -> tuple[Tet,...]:
    if len(vertices)!=7:
        raise ValueError("one center and six outer octahedron vertices required")
    result=[]
    for a,b,c in product((1,2),(3,4),(5,6)):
        tet=(0,a,b,c)
        det=_volume6(*(vertices[i] for i in tet))
        if det==0:
            raise ValueError("degenerate PL 3D corner tetrahedron")
        if det<0:
            tet=(0,a,c,b)
        result.append(tet)
    return tuple(result)


def _edge_link(tetrahedra:tuple[Tet,...],target:int) -> tuple[int,int]:
    """Independent link-circle check for center-to-marked-K-vertex edge."""
    if target not in range(1,7):
        raise ValueError("K edge target must be outer mesh vertex")
    pairs=[]
    for tet in tetrahedra:
        if 0 in tet and target in tet:
            pairs.append(frozenset(v for v in tet if v not in (0,target)))
    if len(pairs)!=4 or len(set(pairs))!=4:
        raise ValueError("K edge link is not four distinct triangles")
    graph:dict[int,set[int]]=defaultdict(set)
    for p in pairs:
        if len(p)!=2:
            raise ValueError("collapsed link edge")
        a,b=tuple(p)
        graph[a].add(b)
        graph[b].add(a)
    if len(graph)!=4 or any(len(nbr)!=2 for nbr in graph.values()):
        raise ValueError("K edge has no four-cycle PL link")
    remaining=set(graph)
    stack=[remaining.pop()]
    while stack:
        for nxt in graph[stack.pop()]:
            if nxt in remaining:
                remaining.remove(nxt)
                stack.append(nxt)
    if remaining:
        raise ValueError("K edge link disconnected")
    return 4,4


def _convex_clip_segment(a:Point3,b:Point3,corner:Point3,
                         generators:tuple[Vector,Vector,Vector]
                         ) -> tuple[Fraction,Fraction]|None:
    """Exactly clip a segment against the CLOSED affine octahedron.

    Transform by Cramer's rule into basis (u,v,w), then intersect the
    8 inequalities ±x±y±z≤1. No approximations or implicit tolerances.
    """
    u,v,w=generators
    det=_dot(u,_cross(v,w))
    if det<=0 or a==b:
        raise ValueError("nondegenerate exact oriented corner prism required")
    def transform(p:Point3)->Vector:
        x=_minus(_xyz(p),_xyz(corner))
        return (_dot(x,_cross(v,w))/det,
                _dot(x,_cross(w,u))/det,
                _dot(x,_cross(u,v))/det)
    p,q=transform(a),transform(b)
    lo,hi=Fraction(0),Fraction(1)
    for signs in product((-1,1),repeat=3):
        f0=sum((signs[j]*p[j] for j in range(3)),Fraction(0))
        velocity=sum((signs[j]*(q[j]-p[j]) for j in range(3)),
                     Fraction(0))
        if velocity==0:
            if f0>1:
                return None
            continue
        boundary=(Fraction(1)-f0)/velocity
        if velocity>0:
            hi=min(hi,boundary)
        else:
            lo=max(lo,boundary)
        if lo>hi:
            return None
    return lo,hi


def _separated_from_cube(vertices:tuple[Point3,...],cube)->bool:
    """Conservative exact sufficient certificate of closed-solid separation.

    Axis-aligned extrema of *all seven octahedral vertices* bound its
    convex hull. Strict separation on any axis proves no shared point.
    If this certificate is inconclusive, shrink radius or fail closed.
    """
    return any(
        max(Fraction(getattr(p,axis)) for p in vertices)<minimum or
        min(Fraction(getattr(p,axis)) for p in vertices)>maximum
        for axis,minimum,maximum in (
            ("x",Fraction(cube.west),Fraction(cube.east)),
            ("y",Fraction(-3),Fraction(3)),
            ("z",Fraction(-3),Fraction(3))
        )
    )


@dataclass(frozen=True)
class RationalKCornerStar:
    """One certified K vertex+two incident subedges inside a real 3-ball."""
    candidate:ExteriorCollaredBandCandidate
    arc:int
    corner:int
    radius:Fraction
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tet,...]
    status:Literal["LOCAL_CORNER_STAR_NOT_WHOLE_K_AMBIENT"]=(
        "LOCAL_CORNER_STAR_NOT_WHOLE_K_AMBIENT"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.candidate,ExteriorCollaredBandCandidate)
                or type(self.arc) is not int or type(self.corner) is not int
                or type(self.radius) is not Fraction
                or not 0<self.radius<=Fraction(1,4)
                or self.status!="LOCAL_CORNER_STAR_NOT_WHOLE_K_AMBIENT"):
            raise ValueError("supported actual K interior corner and rational radius needed")
        arcs=self.candidate.patches.strands.embedded_arcs()
        if not 0<=self.arc<len(arcs) or not 1<=self.corner<len(arcs[self.arc])-1:
            raise ValueError("proper interior two-valent K vertex required")
        before,p,after=arcs[self.arc][self.corner-1:self.corner+2]
        u,v,w,first,second=_generators(before,p,after,self.radius)
        if (self.vertices!=_vertices(p,u,v,w)
                or self.tetrahedra!=_tetrahedra(self.vertices)):
            raise ValueError("noncanonical physical K corner octahedral 3-ball")
        if self.vertices[0]!=p:
            raise ValueError("actual knot vertex not the tetrahedral center")
        _edge_link(self.tetrahedra,first)
        _edge_link(self.tetrahedra,second)
        summary=verify_connected_tetrahedral_ball(self.vertices,self.tetrahedra)
        expected_volume=Fraction(4,3)*_dot(u,_cross(v,w))
        if (summary.vertices!=7 or summary.tetrahedra!=8
                or summary.boundary_components!=1
                or summary.boundary_euler!=2
                or summary.volume!=expected_volume):
            raise ValueError("not a valid geometric PL 3-ball centered on K corner")
        for cube in self.candidate.patches.strands.cubes:
            if not _separated_from_cube(self.vertices,cube):
                raise ValueError("corner 3-ball lacks strict Conway cube separation")
        allowed=set()
        for ai,arc in enumerate(arcs):
            for si,(a,b) in enumerate(zip(arc,arc[1:])):
                hit=_convex_clip_segment(a,b,p,(u,v,w))
                if hit is None:
                    continue
                if (ai,si) not in (
                    (self.arc,self.corner-1),(self.arc,self.corner)):
                    raise ValueError("extraneous source K segment intersects corner 3-ball")
                if hit[0]>=hit[1]:
                    raise ValueError("tangential corner-star contact")
                # The original K segment's common endpoint P must lie
                # in the closed hit interval; the other endpoint must
                # remain strictly outside the chosen corner 3-ball.
                expected_t=Fraction(1) if si==self.corner-1 else Fraction(0)
                if (not hit[0]<=expected_t<=hit[1] or
                        hit[0]==0 and hit[1]==1):
                    raise ValueError("K corner fan does not contain a proper incident subedge")
                if expected_t==0 and not 0<hit[1]<1:
                    raise ValueError("outgoing corner star contains full K edge")
                if expected_t==1 and not 0<hit[0]<1:
                    raise ValueError("incoming corner star contains full K edge")
                allowed.add((ai,si))
        if allowed!={(self.arc,self.corner-1),(self.arc,self.corner)}:
            raise ValueError("both actual K incident edges must enter local 3-ball")
        return summary

    def require_global_source_ambient(self)->None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: a local corner octahedron does not "
            "construct glued K edge/vertex tubes or source Conway slit")


def build_rational_k_corner_star(
    atom:AtomicInstruction,arc:int,corner:int,
    max_radius:Fraction=Fraction(1,4)
)->RationalKCornerStar:
    """Deterministically shrink exact corner neighborhoods, fail closed."""
    if (not isinstance(atom,AtomicInstruction)
            or type(arc) is not int or type(corner) is not int
            or type(max_radius) is not Fraction
            or not 0<max_radius<=Fraction(1,4)):
        raise ValueError("typed source K interior corner and rational radius required")
    candidate=build_exterior_collared_band(atom)
    arcs=candidate.patches.strands.embedded_arcs()
    if not 0<=arc<len(arcs) or not 1<=corner<len(arcs[arc])-1:
        raise ValueError("existing proper interior K corner required")
    before,p,after=arcs[arc][corner-1:corner+2]
    for k in range(24):
        r=max_radius/Fraction(2**k)
        try:
            u,v,w,_,_=_generators(before,p,after,r)
            vertices=_vertices(p,u,v,w)
            return RationalKCornerStar(
                candidate,arc,corner,r,vertices,_tetrahedra(vertices))
        except ValueError:
            continue
    raise ValueError(
        "UNVERIFIED_GEOMETRY: bounded exact rational corner-star search "
        "failed; no nonembedding or source-isotopy conclusion follows")
