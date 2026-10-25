"""A5.2e2b4c2d2d3a: incoming positive-area cap on a REAL K-corner 3-ball.

The previous outgoing-only truncation cannot conform to the far end
of a preceding edge tube: an uncut incoming octahedral apex has ZERO
transverse area. This module truncates the INCOMING K ray at half of
the original octahedral generator. The five-vertex, four-triangle
actual exposed cap is a positive-area disk with a real source-K
1-simplex from its center toward the original corner P.

This is a source-unverified LOCAL 3D geometric cap. It is NOT yet
welded to an edge tube. Global (M,K), source Conway-slit, and R±
certification remain explicitly unavailable.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Literal

from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .local_plumbing_slit import _volume6
from .rational_k_corner_stars import RationalKCornerStar, _generators, build_rational_k_corner_star
from .truncated_k_corner_caps import _average, _oriented_tet, _cap_disk, _link_cycle, _tip
from .spatial_twist_blocks import Point3

Tet=tuple[int,int,int,int]
Tri=tuple[int,int,int]


def _incoming_geometry(
    source:RationalKCornerStar,
)->tuple[tuple[Point3,...],tuple[Tet,...],tuple[Tri,...],int,int,int]:
    """Clip the actual incoming +u octahedron apex at positive area."""
    if not isinstance(source,RationalKCornerStar):
        raise ValueError("an independently validated actual K corner is required")
    source.audit()
    tips=source.vertices
    outgoing,incoming,_=_tip(source)
    if incoming!=1:
        raise ValueError("unsupported original physical K incoming axis")
    # d1 = preceding K vertex - corner is always the +u direction.
    # Cross-sectional cyclic order viewed down +u uses +v,+w,-v,-w.
    neighbors=(3,5,4,6)
    kept=tuple(i for i in range(7) if i!=incoming)
    places={label:i for i,label in enumerate(kept)}
    points=[tips[i] for i in kept]
    center=len(points)
    points.append(_average(tips[0],tips[incoming]))
    mid={}
    for vertex in neighbors:
        mid[vertex]=len(points)
        points.append(_average(tips[vertex],tips[incoming]))
    vertices=tuple(points)
    cells=[]
    for a,b,c in product((1,2),(3,4),(5,6)):
        facet=(a,b,c)
        if incoming not in facet:
            cells.append(_oriented_tet(
                vertices,(places[0],places[a],places[b],places[c])))
        else:
            x,y=(j for j in facet if j!=incoming)
            cells.extend((
                _oriented_tet(vertices,(places[0],places[x],
                                       places[y],mid[y])),
                _oriented_tet(vertices,(places[0],places[x],
                                       mid[y],mid[x]))
            ))
    cap=[]
    for i in range(4):
        first=mid[neighbors[i]]
        second=mid[neighbors[(i+1)%4]]
        triangle=(center,first,second)
        cap.append(triangle)
        cells.append(_oriented_tet(vertices,(places[0],*triangle)))
    return vertices,tuple(cells),tuple(cap),places[outgoing],center,places[0]


@dataclass(frozen=True)
class IncomingTruncatedKCorner:
    """A verified local 3-ball with an actual positive-area incoming K disk."""
    original:RationalKCornerStar
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tet,...]
    cap_faces:tuple[Tri,...]
    outgoing_vertex:int
    cap_center:int
    k_center:int
    status:Literal["INCOMING_K_CAP_NOT_TUBE_GLUED"]=(
        "INCOMING_K_CAP_NOT_TUBE_GLUED"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.original,RationalKCornerStar)
                or self.status!="INCOMING_K_CAP_NOT_TUBE_GLUED"):
            raise ValueError("source-unverified incoming K corner required")
        expected=_incoming_geometry(self.original)
        actual=(self.vertices,self.tetrahedra,self.cap_faces,
                self.outgoing_vertex,self.cap_center,self.k_center)
        if actual!=expected:
            raise ValueError("tampered physical incoming corner cap")
        topology=verify_connected_tetrahedral_ball(self.vertices,self.tetrahedra)
        if (topology.vertices!=11 or topology.tetrahedra!=16
                or topology.boundary_euler!=2
                or topology.boundary_components!=1):
            raise ValueError("incoming truncated K cap is not a PL 3-ball")
        if _cap_disk(self.vertices,self.cap_faces,self.cap_center)!=(5,8,4):
            raise ValueError("not a physical 5-vertex / 4-face cap disk")
        incidences=Counter()
        for a,b,c,d in self.tetrahedra:
            for face in ((b,c,d),(a,d,c),(a,b,d),(a,c,b)):
                incidences[frozenset(face)]+=1
        if any(incidences[frozenset(face)]!=1 for face in self.cap_faces):
            raise ValueError("incoming cap face is not exposed on the 3-ball")
        # The incoming K-edge ends at the newly created boundary cap
        # center. The untouched outgoing ray is a genuine mesh edge.
        incoming_link=_link_cycle(
            self.tetrahedra,self.k_center,self.cap_center)
        outgoing_link=_link_cycle(
            self.tetrahedra,self.k_center,self.outgoing_vertex)
        outgoing_tip,_,_=_tip(self.original)
        if incoming_link!=(4,4) or outgoing_link!=(
            (4,4) if outgoing_tip==2 else (6,6)
        ):
            raise ValueError("K rays have wrong tetrahedral edge-link cycles")
        k_paths=self.original.candidate.patches.strands.embedded_arcs()
        before,p,after=k_paths[self.original.arc][
            self.original.corner-1:self.original.corner+2]
        u,v,w,_,_=_generators(before,p,after,self.original.radius)
        det=u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0])
        if topology.volume!=Fraction(5,4)*det:
            raise ValueError("incoming geometric truncation volume mismatch")
        if self.vertices[self.k_center]!=p:
            raise ValueError("original K corner was moved")
        expected_k_cap=Point3(*(
            Fraction(getattr(p,k))+
            self.original.radius*(
                Fraction(getattr(before,k))-Fraction(getattr(p,k)))/2
            for k in ("x","y","z")))
        if self.vertices[self.cap_center]!=expected_k_cap:
            raise ValueError("incoming cap point not on original physical K edge")
        return topology

    def require_full_source_plumbing(self)->None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: positive-area incoming cap is not "
            "a glued K tube, entire K ambient, 4-punctured Conway "
            "source slit, framed isotopy, or ring R+/- proof")


def truncate_incoming_k_corner(source:RationalKCornerStar)->IncomingTruncatedKCorner:
    return IncomingTruncatedKCorner(source,*_incoming_geometry(source))


def build_incoming_k_corner_cap(atom,arc:int,corner:int,
                                max_radius:Fraction=Fraction(1,4)
                                )->IncomingTruncatedKCorner:
    return truncate_incoming_k_corner(
        build_rational_k_corner_star(atom,arc,corner,max_radius))
