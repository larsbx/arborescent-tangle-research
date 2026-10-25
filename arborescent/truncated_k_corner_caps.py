"""A5.2e2b4c2d2d1: positive-area, tetrahedrally conforming K corner cap.

A previous convex octahedral K-corner star ends its marked outgoing
K ray at a SINGLE octahedron apex, so it cannot be glued to the
positive-area end of an edge tube without a nonmanifold point link.

This module TRUNCATES that actual corner octahedron at one half of
the outgoing axial coordinate. The newly exposed section is a
genuine five-vertex/four-triangle disk: its CENTER lies on the
outgoing K ray, and its four corners form an exact rational diamond.
A canonical 16-tetrahedron subdivision of the surviving convex
octahedron is independently checked as a genuine oriented PL 3-ball;
both K ray edges are real conforming 1-simplices. The exact removed
tip volume is 1/12 of the frame determinant.

The cap is a POSITIVE-AREA ATTACHMENT INTERFACE, not a glued tube.
No interface identification, source Fig.12.10 manifold, new Conway
boundary or ring R+/- is certified.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .local_plumbing_slit import _volume6
from .rational_k_corner_stars import (
    RationalKCornerStar, _generators, build_rational_k_corner_star,
)
from .spatial_twist_blocks import Point3

Tetrahedron = tuple[int,int,int,int]
Triangle = tuple[int,int,int]


def _average(a:Point3,b:Point3)->Point3:
    return Point3(*((Fraction(getattr(a,k))+Fraction(getattr(b,k)))/2
                     for k in ("x","y","z")))


def _oriented_tet(vertices:tuple[Point3,...],indices:Tetrahedron)->Tetrahedron:
    a,b,c,d=(vertices[i] for i in indices)
    det=_volume6(a,b,c,d)
    if det==0:
        raise ValueError("collapsed 3D truncated-corner simplex")
    if det<0:
        return (indices[0],indices[2],indices[1],indices[3])
    return indices


def _tip(original:RationalKCornerStar)->tuple[int,int,tuple[int,int,int,int]]:
    """Return outgoing apex, incoming apex and oriented cap-ring neighbor IDs."""
    arcs=original.candidate.patches.strands.embedded_arcs()
    a,p,b=arcs[original.arc][original.corner-1:original.corner+2]
    _,_,_,incoming,outgoing=_generators(a,p,b,original.radius)
    if (incoming,outgoing)==(1,3):
        # Cross-section normal points along +v:
        # oriented cap-ring order is +u,-w,-u,+w.
        return outgoing,incoming,(1,6,2,5)
    if (incoming,outgoing)==(1,2):
        # Straight case: outgoing normal points along -u, whose
        # positive cap-ring order is +v,-w,-v,+w.
        return outgoing,incoming,(3,6,4,5)
    raise ValueError("unsupported physical K corner axis choice")


def _canonical_geometry(original:RationalKCornerStar
                        )->tuple[tuple[Point3,...],tuple[Tetrahedron,...],
                                 tuple[Triangle,...],int,int,int]:
    original.audit()
    initial=original.vertices
    tip,incoming,neighbors=_tip(original)
    opposite=tip+1 if tip%2 else tip-1
    # Keep the center, the incoming/opposite apex, and the four
    # remaining octahedron tips. The outgoing tip itself is removed.
    kept=tuple(i for i in range(7) if i!=tip)
    coordinates=[initial[i] for i in kept]
    position={i:index for index,i in enumerate(kept)}
    cap_center=len(coordinates)
    coordinates.append(_average(initial[0],initial[tip]))
    cuts={}
    for i in neighbors:
        cuts[i]=len(coordinates)
        coordinates.append(_average(initial[i],initial[tip]))
    vertices=tuple(coordinates)
    tetra=[]
    # Eight original octahedron facets, each the choice of one
    # signed apex from each of its three axis pairs. Exactly four
    # facets contain the outgoing tip. Their clipped boundaries
    # become quadrilaterals, each subdivided into two triangles.
    for a,b,c in product((1,2),(3,4),(5,6)):
        face=(a,b,c)
        if tip not in face:
            tetra.append(_oriented_tet(vertices,(
                position[0],position[a],position[b],position[c])))
            continue
        x,y=(i for i in face if i!=tip)
        tetra.extend((
            _oriented_tet(vertices,(position[0],position[x],
                                    position[y],cuts[y])),
            _oriented_tet(vertices,(position[0],position[x],
                                    cuts[y],cuts[x])),
        ))
    # The new physical cap is a triangulated diamond disk with
    # a real K center point, NOT an apex-only shared vertex.
    cap=[]
    for i in range(4):
        p,q=cuts[neighbors[i]],cuts[neighbors[(i+1)%4]]
        face=(cap_center,p,q)
        cap.append(face)
        tetra.append(_oriented_tet(vertices,(position[0],*face)))
    return vertices,tuple(tetra),tuple(cap),position[incoming],cap_center,position[0]


def _cap_disk(
    vertices:tuple[Point3,...],triangles:tuple[Triangle,...],
    cap_center:int
)->tuple[int,int,int]:
    """Independent PL disk topology, oriented area and cap boundary."""
    edges=Counter()
    direction=Counter()
    used=set()
    for tri in triangles:
        if (len(tri)!=3 or len(set(tri))!=3
                or cap_center not in tri):
            raise ValueError("cap mesh is not a fan about the real K point")
        used.update(tri)
        for i,j in zip(tri,tri[1:]+tri[:1]):
            key=frozenset((i,j))
            edges[key]+=1
            direction[(i,j)]+=1
        if _volume6(vertices[0],*(vertices[i] for i in tri))==0:
            raise ValueError("cap face collapsed against K center")
    boundary=[e for e,c in edges.items() if c==1]
    if (len(triangles)!=4 or len(used)!=5 or len(edges)!=8
            or len(boundary)!=4 or len(used)-len(edges)+len(triangles)!=1
            or any(c not in (1,2) for c in edges.values())
            or any(direction[(a,b)]!=direction[(b,a)]
                   for e,c in edges.items() if c==2
                   for a,b in (tuple(e),))):
        raise ValueError("truncated K cap is not an oriented PL 2-disk")
    counts=Counter(v for e in boundary for v in e)
    if len(counts)!=4 or any(c!=2 for c in counts.values()):
        raise ValueError("cap boundary is not a simple four-cycle")
    return len(used),len(edges),len(triangles)


def _link_cycle(tets:tuple[Tetrahedron,...],
                center:int,target:int)->tuple[int,int]:
    """Independently certify a whole tetrahedral edge's link is a circle.

    Triangular boundary refinements can increase the degree, so
    unlike the source octahedron no hardcoded 'four' is assumed.
    """
    opposite=[]
    for tet in tets:
        if center in tet and target in tet:
            opposite.append(frozenset(v for v in tet
                                       if v not in (center,target)))
    if len(opposite)<3 or len(set(opposite))!=len(opposite):
        raise ValueError("missing or repeated K-edge link faces")
    graph={}
    for edge in opposite:
        if len(edge)!=2:
            raise ValueError("collapsed K edge link")
        a,b=tuple(edge)
        graph.setdefault(a,set()).add(b)
        graph.setdefault(b,set()).add(a)
    if (len(graph)!=len(opposite) or
            any(len(nbr)!=2 for nbr in graph.values())):
        raise ValueError("K edge link is not a simple cycle")
    remaining=set(graph)
    stack=[remaining.pop()]
    while stack:
        for next_vertex in graph[stack.pop()]:
            if next_vertex in remaining:
                remaining.remove(next_vertex)
                stack.append(next_vertex)
    if remaining:
        raise ValueError("K edge link is disconnected")
    return len(graph),len(opposite)


@dataclass(frozen=True)
class TruncatedKCorner:
    """Real K 3D corner with a positive-area outgoing attachment disk."""
    original:RationalKCornerStar
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tetrahedron,...]
    cap_faces:tuple[Triangle,...]
    incoming_vertex:int
    cap_center:int
    k_center:int
    status:Literal["TRUNCATED_K_CORNER_CAP_NOT_GLUED"]=(
        "TRUNCATED_K_CORNER_CAP_NOT_GLUED"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.original,RationalKCornerStar)
                or self.status!="TRUNCATED_K_CORNER_CAP_NOT_GLUED"):
            raise ValueError("typed source-UNVERIFIED truncated corner required")
        expected=_canonical_geometry(self.original)
        if (self.vertices,self.tetrahedra,self.cap_faces,self.incoming_vertex,
                self.cap_center,self.k_center)!=expected:
            raise ValueError("truncated K corner geometry or mesh was forged")
        topology=verify_connected_tetrahedral_ball(
            self.vertices,self.tetrahedra)
        if (topology.vertices!=11 or topology.tetrahedra!=16
                or topology.boundary_components!=1
                or topology.boundary_euler!=2):
            raise ValueError("truncated corner is not a tetrahedral PL 3-ball")
        _cap_disk(self.vertices,self.cap_faces,self.cap_center)
        # Independently reconstruct the actual external boundary
        # triangles from the tetrahedra: the four cap faces must each
        # occur EXACTLY once, i.e. be genuinely exposed 2-simplices.
        incidence=Counter()
        for tet in self.tetrahedra:
            a,b,c,d=tet
            for tri in ((b,c,d),(a,d,c),(a,b,d),(a,c,b)):
                incidence[frozenset(tri)]+=1
        if any(incidence[frozenset(face)]!=1 for face in self.cap_faces):
            raise ValueError("cap face is buried rather than boundary-attached")
        incoming=_link_cycle(self.tetrahedra,
                             self.k_center,self.incoming_vertex)
        outgoing=_link_cycle(self.tetrahedra,
                             self.k_center,self.cap_center)
        tip,_,_=_tip(self.original)
        # In a bend the incoming and outgoing axes are adjacent,
        # splitting two original adjacent facets, yielding 6 faces
        # around the incoming K edge. For a straight continuation
        # their tips are opposite, so its link remains a 4-cycle.
        incoming_expected=(4,4) if tip==2 else (6,6)
        if incoming!=incoming_expected or outgoing!=(4,4):
            raise ValueError("truncated K corner has wrong edge-link refinement")
        # Removing the half-size cone over one octahedron tip costs
        # det(u,v,w)/12, leaving (5/4)det(u,v,w) of volume.
        before,p,after=self.original.candidate.patches.strands.embedded_arcs()[
            self.original.arc][self.original.corner-1:self.original.corner+2]
        u,v,w,_,_=_generators(before,p,after,self.original.radius)
        det=(Fraction(u[0])*(v[1]*w[2]-v[2]*w[1])
             -u[1]*(v[0]*w[2]-v[2]*w[0])
             +u[2]*(v[0]*w[1]-v[1]*w[0]))
        if topology.volume!=Fraction(5,4)*det:
            raise ValueError("truncated octahedron volume conservation failed")
        return topology

    def require_source_tube_attachment(self):
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: the corner has a verified positive-area "
            "cap but no conforming PR63 edge-tube transition or global "
            "source Conway slitting/isotopy")


def truncate_source_k_corner(original:RationalKCornerStar)->TruncatedKCorner:
    if not isinstance(original,RationalKCornerStar):
        raise ValueError("verified actual K two-valent corner required")
    return TruncatedKCorner(original,*_canonical_geometry(original))


def build_truncated_k_corner(
    atom:AtomicInstruction,arc:int,corner:int,
    max_radius:Fraction=Fraction(1,4)
)->TruncatedKCorner:
    return truncate_source_k_corner(
        build_rational_k_corner_star(atom,arc,corner,max_radius))
