"""A5.2e2b4c2d2d3b: two real K corner 3-balls joined by one edge tube.

Take a previously independently checked outgoing K-corner/tube
40-tetrahedron PL 3-ball, and the actual NEXT K corner's positive-area
incoming cap. Weld its five-vertex diamond disk to the tube's far
five-vertex square using four actual triangular prisms (12 tetrahedra).
The 68 tetrahedra must form ONE embedded oriented PL 3-ball; every
pair of tetrahedra undergoes exact rational physical collision audit.

Only one original K edge and its two adjacent true K vertices are
connected. The *complete* 2n source K arc set, source ribbon/Conway
ambient, source Fig12.10 and ring R+/- remain UNVERIFIED.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .exact_tetrahedron_contacts import audit_embedded_tetrahedra
from .incoming_k_corner_caps import (
    IncomingTruncatedKCorner,build_incoming_k_corner_cap,
)
from .spatial_twist_blocks import Point3
from .welded_k_corner_edge import (
    WeldedKCornerToEdge,build_welded_k_corner_to_edge,
    _oriented,_edges,_boundary_faces,_clip_segment_tet,
    _closed_box_separates,_interior_source_parameter,
)

Tet=tuple[int,int,int,int]
Tri=tuple[int,int,int]


def _assemble(
    near:WeldedKCornerToEdge,far:IncomingTruncatedKCorner,
    rotation:int,reverse:bool,
)->tuple[tuple[Point3,...],tuple[Tet,...],tuple[Tet,...],
         tuple[Tri,...],tuple[Tri,...]]:
    if (type(rotation) is not int or rotation not in range(4)
            or type(reverse) is not bool):
        raise ValueError("typed cyclic/dihedral positive-area interface mapping required")
    all_vertices=list(near.vertices)
    ids={p:i for i,p in enumerate(all_vertices)}
    def add(point:Point3)->int:
        if point not in ids:
            ids[point]=len(all_vertices)
            all_vertices.append(point)
        return ids[point]
    corner_ids=tuple(add(v) for v in far.vertices)
    tube=near.tube
    near_ring=tuple(ids[tube.vertices[i]] for i in (6,7,8,9))
    near_center=ids[tube.vertices[5]]
    far_ring_original=tuple(
        corner_ids[face[1]] for face in far.cap_faces)
    if (len(set(far_ring_original))!=4 or
            len(set(near_ring))!=4 or
            near_center in near_ring):
        raise ValueError("invalid physical closed cap disks")
    far_ring=tuple(
        far_ring_original[(rotation+(-i if reverse else i))%4]
        for i in range(4))
    far_center=corner_ids[far.cap_center]
    vertices=tuple(all_vertices)
    bridge=[]
    first=[]
    last=[]
    for i in range(4):
        j=(i+1)%4
        a,b,c=near_center,near_ring[i],near_ring[j]
        A,B,C=far_center,far_ring[i],far_ring[j]
        first.append((a,b,c))
        last.append((A,B,C))
        for tet in ((a,b,c,C),(a,b,B,C),(a,A,B,C)):
            bridge.append(_oriented(vertices,tet))
    incoming_tets=tuple(
        tuple(corner_ids[i] for i in tet)
        for tet in far.tetrahedra)
    return (vertices,near.tetrahedra+tuple(bridge)+incoming_tets,
            tuple(bridge),tuple(first),tuple(last))


@dataclass(frozen=True)
class TwoEndedKVertexEdge:
    """True connected rational 3-ball covering an actual source K edge."""
    near:WeldedKCornerToEdge
    far:IncomingTruncatedKCorner
    rotation:int
    reverse:bool
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tet,...]
    bridge_tetrahedra:tuple[Tet,...]
    near_faces:tuple[Tri,...]
    far_faces:tuple[Tri,...]
    status:Literal["TWO_CORNER_K_EDGE_NOT_WHOLE_SOURCE_AMBIENT"]=(
        "TWO_CORNER_K_EDGE_NOT_WHOLE_SOURCE_AMBIENT"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.near,WeldedKCornerToEdge)
                or not isinstance(self.far,IncomingTruncatedKCorner)
                or self.status!="TWO_CORNER_K_EDGE_NOT_WHOLE_SOURCE_AMBIENT"):
            raise ValueError("one genuine two-corner source-candidate K edge required")
        self.near.audit()
        self.far.audit()
        ai,si=self.near.tube.arc,self.near.tube.segment
        if (self.far.original.candidate!=self.near.tube.candidate
                or self.far.original.arc!=ai
                or self.far.original.corner!=si+1):
            raise ValueError("incoming cap is not NEXT actual physical K corner")
        expected=_assemble(self.near,self.far,self.rotation,self.reverse)
        if ((self.vertices,self.tetrahedra,self.bridge_tetrahedra,
             self.near_faces,self.far_faces)!=expected):
            raise ValueError("forged or nonconforming two-sided tetrahedral glue")
        if (len(self.tetrahedra)!=68 or len(self.vertices)!=32
                or len(self.bridge_tetrahedra)!=12):
            raise ValueError("missing 3D tetrahedron in two-sided K attachment")
        summary=verify_connected_tetrahedral_ball(
            self.vertices,self.tetrahedra)
        if summary.boundary_euler!=2 or summary.boundary_components!=1:
            raise ValueError("two-ended K edge is not one PL 3-ball")
        collisions=audit_embedded_tetrahedra(
            self.vertices,self.tetrahedra)
        if collisions.pair_count!=68*67//2:
            raise ValueError("not all nonlocal tetrahedra independently checked")
        boundary=_boundary_faces(self.tetrahedra)
        if any(frozenset(tri) in boundary
               for tri in (*self.near_faces,*self.far_faces)):
            raise ValueError("diamond/square attachment face was not welded internally")
        paths=self.near.corner.original.candidate.patches.strands.embedded_arcs()
        path=paths[ai]
        p,q=path[si:si+2]
        original_k=(
            self.near.corner.vertices[self.near.corner.k_center],
            self.near.corner.vertices[self.near.corner.cap_center],
            self.near.tube.vertices[0],self.near.tube.vertices[5],
            self.far.vertices[self.far.cap_center],
            self.far.vertices[self.far.k_center],
        )
        if original_k[0]!=p or original_k[-1]!=q:
            raise ValueError("true K endpoints moved during tetrahedral welding")
        param=tuple(_interior_source_parameter(p,q,x)
                    for x in original_k[1:-1])
        if not (Fraction(0)<param[0]<param[1]<param[2]<param[3]<Fraction(1)):
            raise ValueError("two-ended tetrahedral K chain reverses/skips source")
        ids={v:i for i,v in enumerate(self.vertices)}
        mesh_edges=_edges(self.tetrahedra)
        if any(frozenset((ids[a],ids[b])) not in mesh_edges
               for a,b in zip(original_k,original_k[1:])):
            raise ValueError("K source centerline is not a true PL 1-subcomplex")
        # A separate geometric closed-segment/convex-tetrahedron clip
        # rejects unrelated K in the NEW transition piece. Already
        # vetted near and far component cells remain independently audited.
        for source_arc,arc in enumerate(paths):
            for source_segment,(a,b) in enumerate(zip(arc,arc[1:])):
                if (source_arc,source_segment)==(ai,si):
                    continue
                for tet in self.bridge_tetrahedra:
                    physical=tuple(self.vertices[i] for i in tet)
                    if _clip_segment_tet(a,b,physical) is not None:
                        raise ValueError("extraneous original K segment hits new bridge")
        for tet in self.bridge_tetrahedra:
            physical=tuple(self.vertices[i] for i in tet)
            if any(not _closed_box_separates(physical,cube)
                   for cube in self.near.tube.candidate.patches.strands.cubes):
                raise ValueError("two-sided bridge intersects Conway solid cube")
        return summary

    def require_source_plumbing_certificate(self)->None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: only one actual K edge and two corners "
            "are mesh-conforming; the complete knot/ribbon/Conway "
            "source ambient, four-puncture slit and ring R+/- remain open")


def build_two_ended_k_edge(
    atom:AtomicInstruction,arc:int,first_corner:int,
    max_radius:Fraction=Fraction(1,4),
)->TwoEndedKVertexEdge:
    """Bounded exact trial of radius/dihedral face maps. Failure=UNVERIFIED."""
    if (not isinstance(atom,AtomicInstruction)
            or type(arc) is not int or type(first_corner) is not int
            or type(max_radius) is not Fraction
            or not Fraction(0)<max_radius<=Fraction(1,4)):
        raise ValueError("typed actual K corners and rational radius required")
    for k in range(8):
        radius=max_radius/Fraction(2**k)
        try:
            near=build_welded_k_corner_to_edge(
                atom,arc,first_corner,radius)
            far=build_incoming_k_corner_cap(
                atom,arc,first_corner+1,radius)
        except (ValueError,NotImplementedError):
            continue
        for reverse in (False,True):
            for rotation in range(4):
                try:
                    geometry=_assemble(near,far,rotation,reverse)
                    return TwoEndedKVertexEdge(
                        near,far,rotation,reverse,*geometry)
                except ValueError:
                    continue
    raise ValueError(
        "UNVERIFIED_GEOMETRY: bounded rational two-ended knot "
        "attachment search found no valid embedded PL 3-ball; "
        "this is NOT a nonisotopy theorem")
