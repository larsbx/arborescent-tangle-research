"""A5.2e2b4c2d2a: actual K-edge-CONFORMING tetrahedral 3-ball tubes.

Two long straight RETURN-rail segments are already physically part
of the chosen r=0 K subset of the globally embedded ribbon boundary.
Around an exact rational proper SUBSEGMENT of each, construct a
tetrahedral 3-ball whose *interior mesh EDGE is literally that K
subsegment*. Verify all geometric source K segments, cube solids,
signed tetrahedral volumes, PL 3D vertex links, knot-edge links and
the exact occupied tube volume independently. No tolerance or snapping.

These are two DISJOINT local knot-conforming 3-balls. They do NOT
mesh the complete 2n arcs or the band/Conway spheres as conforming
subcomplexes, are not a connected knot-containing source ambient
(M,K), and DO NOT supply Figure 12.10 slitting or R+/- proofs.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .exterior_conway_collars import (
    ExteriorCollaredBandCandidate, build_exterior_collared_band,
)
from .local_plumbing_slit import _volume6
from .spatial_twist_blocks import Point3

Tet = tuple[int,int,int,int]
Edge = frozenset[int]


def _q(value: object) -> Fraction:
    if type(value) not in (int,Fraction):
        raise ValueError("exact integer or rational coordinate required")
    return Fraction(value)


def _clip_closed_box(
    a: Point3,b: Point3,
    limits: tuple[tuple[Fraction,Fraction],...]
) -> tuple[Fraction,Fraction] | None:
    """Exact convex clipping of ONE spatial segment to a closed cuboid."""
    if a==b or len(limits)!=3:
        raise ValueError("nondegenerate rational segment and 3D box required")
    lo,hi=Fraction(0),Fraction(1)
    for k,(minimum,maximum) in enumerate(limits):
        if not minimum<maximum:
            raise ValueError("nonempty exact box interval required")
        p=_q(getattr(a,("x","y","z")[k]))
        slope=_q(getattr(b,("x","y","z")[k]))-p
        if slope==0:
            if p<minimum or p>maximum:
                return None
            continue
        i,j=(minimum-p)/slope,(maximum-p)/slope
        lo=max(lo,min(i,j))
        hi=min(hi,max(i,j))
        if lo>hi:
            return None
    return lo,hi


def _long_rail_edge(arcs: tuple[tuple[Point3,...],...],
                    number: int) -> tuple[Point3,Point3]:
    if (number not in (0,1) or number>=len(arcs)
            or not isinstance(arcs[number],tuple)):
        raise ValueError("two original return rail arcs required")
    expected_y=-6 if number==0 else 6
    possibilities=[
        (a,b) for a,b in zip(arcs[number],arcs[number][1:])
        if (a.y==b.y==expected_y and a.z==b.z==0
            and a.x!=b.x)
    ]
    if len(possibilities)!=1:
        raise ValueError("original geometric K has no unique long return edge")
    return possibilities[0]


def _tube_vertices(start: Fraction,end: Fraction,
                   center_y: int,r: Fraction) -> tuple[Point3,...]:
    return tuple(Point3(x,y,z) for x in (start,end)
                 for y,z in (
                     (Fraction(center_y),Fraction(0)),
                     (Fraction(center_y)-r,-r),
                     (Fraction(center_y)+r,-r),
                     (Fraction(center_y)+r,r),
                     (Fraction(center_y)-r,r)
                 ))


def _fan_tetrahedra(vertices: tuple[Point3,...]) -> tuple[Tet,...]:
    """Global vertex-order staircase over four triangles of one 2-disk."""
    tets=[]
    for i in range(4):
        b,c=sorted((i+1,(i+1)%4+1))
        a,A,B,C=0,5,5+b,5+c
        for candidate in ((a,b,c,C),(a,b,B,C),(a,A,B,C)):
            det=_volume6(*(vertices[j] for j in candidate))
            if det==0:
                raise ValueError("degenerate tetrahedron from K fan")
            tet=candidate if det>0 else (
                candidate[0],candidate[2],candidate[1],candidate[3])
            tets.append(tet)
    return tuple(tets)


def _core_edge_link(tets: tuple[Tet,...]) -> tuple[int,int]:
    """Independently check the internal K core has one circular edge link."""
    cross=Edge((0,5))
    opposite=[]
    for tet in tets:
        if cross.issubset(tet):
            opposite.append(tuple(v for v in tet if v not in cross))
    links=Counter(frozenset(pair) for pair in opposite)
    if not opposite or any(m!=1 for m in links.values()):
        raise ValueError("K mesh edge star has duplicate/collapsed link edges")
    graph: dict[int,set[int]]=defaultdict(set)
    for pair in links:
        a,b=tuple(pair)
        graph[a].add(b)
        graph[b].add(a)
    if (len(graph)!=4 or len(opposite)!=4
            or any(len(nbr)!=2 for nbr in graph.values())):
        raise ValueError("K core edge has no four-cycle PL link")
    unseen=set(graph)
    frontier=[unseen.pop()]
    while frontier:
        for v in graph[frontier.pop()]:
            if v in unseen:
                unseen.remove(v)
                frontier.append(v)
    if unseen:
        raise ValueError("K core edge link disconnected")
    return len(graph),len(opposite)


@dataclass(frozen=True)
class KEdgeConformingTube:
    """Genuine 3D tetrahedral ball with an internal K 1-simplex edge."""
    candidate: ExteriorCollaredBandCandidate
    source_arc: int
    radius: Fraction
    vertices: tuple[Point3,...]
    tetrahedra: tuple[Tet,...]
    status: Literal["LOCAL_K_EDGE_CONFORMING_NOT_GLOBAL_AMBIENT"] = (
        "LOCAL_K_EDGE_CONFORMING_NOT_GLOBAL_AMBIENT"
    )

    def __post_init__(self) -> None:
        self.audit()

    def core(self) -> tuple[Point3,Point3]:
        return self.vertices[0],self.vertices[5]

    def bounds(self) -> tuple[tuple[Fraction,Fraction],...]:
        if not self.vertices:
            raise ValueError("no physical K tube")
        return tuple((
            min(_q(getattr(p,k)) for p in self.vertices),
            max(_q(getattr(p,k)) for p in self.vertices))
            for k in ("x","y","z"))

    def audit(self):
        if (not isinstance(self.candidate,ExteriorCollaredBandCandidate)
                or type(self.source_arc) is not int
                or self.source_arc not in (0,1)
                or type(self.radius) is not Fraction
                or not 0<self.radius<=Fraction(1,4)
                or self.status!="LOCAL_K_EDGE_CONFORMING_NOT_GLOBAL_AMBIENT"):
            raise ValueError("physical source K segment and rational radius required")
        source_arcs=self.candidate.patches.strands.embedded_arcs()
        original_a,original_b=_long_rail_edge(source_arcs,self.source_arc)
        x1,x2=sorted((_q(original_a.x),_q(original_b.x)))
        # Exactly one unit is trimmed at BOTH ends to avoid the turn
        # and to put core endpoints in the INTERIOR of the source rail.
        start,end=x1+1,x2-1
        cy=-6 if self.source_arc==0 else 6
        expected=_tube_vertices(start,end,cy,self.radius)
        if (not start<end
                or self.vertices!=expected
                or self.tetrahedra!=_fan_tetrahedra(expected)):
            raise ValueError("K-conforming 3D tube is not canonical")
        # The K subsegment is an ACTUAL internal 1-simplex in the
        # tetrahedral complex, not just a curve passing through cells.
        if not any({0,5}.issubset(tet) for tet in self.tetrahedra):
            raise ValueError("physical K core is not a tetrahedral mesh edge")
        _core_edge_link(self.tetrahedra)
        summary=verify_connected_tetrahedral_ball(self.vertices,self.tetrahedra)
        if (summary.vertices!=10 or summary.tetrahedra!=12
                or summary.boundary_euler!=2
                or summary.boundary_components!=1
                or summary.volume!=(end-start)*(2*self.radius)**2):
            raise ValueError("K-conforming tetrahedral 3-ball fails exact audit")

        limits=self.bounds()
        if (not (-7<limits[1][0]<limits[1][1]<-5)
                if self.source_arc==0 else
                not (5<limits[1][0]<limits[1][1]<7)):
            raise ValueError("K tube escaped physically separated outer return rail")
        # Independent exact segment-box clipping of EVERY physical K
        # edge. Only the SOURCE long rail edge can enter the closed
        # 3D tube, and it must contain the whole mesh core.
        allowed=0
        for ai,arc in enumerate(source_arcs):
            for a,b in zip(arc,arc[1:]):
                clipped=_clip_closed_box(a,b,limits)
                if clipped is None:
                    continue
                if ai!=self.source_arc or (a,b)!=(original_a,original_b):
                    raise ValueError("some other K segment crosses local 3D tube")
                allowed+=1
                if not (min(a.x,b.x)<start<end<max(a.x,b.x)):
                    raise ValueError("mesh core is not interior to its source K rail")
        if allowed!=1:
            raise ValueError("source K core edge not uniquely contained in 3D tube")
        # All existing Conway closed solid balls have |y|<=3.
        if any(not (limits[1][1]<-3 or limits[1][0]>3)
               for _ in self.candidate.patches.strands.cubes):
            raise ValueError("K tube intersects removed Conway ball")
        return summary

    def require_source_ambient(self):
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: local tetrahedral core subedges "
            "do not tetrahedralize ALL K, the supporting band, "
            "Conway spheres, or the source global ambient (M,K)")


@dataclass(frozen=True)
class KCoreTubePair:
    """Two physically disjoint knot edge stars, NOT one ambient manifold."""
    first: KEdgeConformingTube
    second: KEdgeConformingTube
    status: Literal["TWO_DISJOINT_K_TUBES_NOT_SOURCE_AMBIENT"] = (
        "TWO_DISJOINT_K_TUBES_NOT_SOURCE_AMBIENT"
    )

    def __post_init__(self) -> None:
        if (not isinstance(self.first,KEdgeConformingTube)
                or not isinstance(self.second,KEdgeConformingTube)
                or self.first.source_arc!=0 or self.second.source_arc!=1
                or self.first.candidate!=self.second.candidate
                or self.status!="TWO_DISJOINT_K_TUBES_NOT_SOURCE_AMBIENT"):
            raise ValueError("two distinct K-return-rail physical tubes required")
        self.first.audit()
        self.second.audit()
        x,y=self.first.bounds(),self.second.bounds()
        if not (x[1][1]<y[1][0]):
            raise ValueError("chosen 3D K tubes not spatially disjoint")

    def require_source_global_cut(self):
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: two K-edge local 3-balls are "
            "NOT a connected source 3-manifold or a Conway slit")


def build_k_core_tubes(
    atom: AtomicInstruction, radius: Fraction=Fraction(1,4)
) -> KCoreTubePair:
    if (not isinstance(atom,AtomicInstruction)
            or type(radius) is not Fraction
            or not 0<radius<=Fraction(1,4)):
        raise ValueError("supported source candidate and exact positive radius required")
    candidate=build_exterior_collared_band(atom)
    # Candidate constructor independently audits all 3D surface
    # intersections and closed Conway ball clearance.
    rails=candidate.patches.strands.embedded_arcs()
    result=[]
    for number,center_y in ((0,-6),(1,6)):
        p,q=_long_rail_edge(rails,number)
        lo,hi=sorted((_q(p.x),_q(q.x)))
        verts=_tube_vertices(lo+1,hi-1,center_y,radius)
        result.append(KEdgeConformingTube(
            candidate,number,radius,verts,_fan_tetrahedra(verts)))
    return KCoreTubePair(result[0],result[1])
