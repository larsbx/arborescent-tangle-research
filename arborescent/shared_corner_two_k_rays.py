"""A5.2e2b4c2e2: two incident K-edge tubes welded to ONE shared 3D corner.

The dual-cap knot core is the ONLY vertex neighborhood at an actual
two-valent K corner. Its incoming and outgoing four-face disks are
welded to the respective square sections of two distinct real
source-candidate K-segment tubes by 12 tetrahedra on EACH side.

No overlapping incoming-only/outgoing-only octahedra are allowed.
All physical tetrahedra and K mesh edges are exact Q^3 geometry.
This is a LOCAL two-edge K-vertex neighborhood, not the whole
arborescent K, source (M,K), Figure 12.10 plumbing or R+/- proof.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .doubly_truncated_k_corner import (
    DoublyTruncatedKCorner,build_doubly_truncated_k_corner,
)
from .exact_tetrahedron_contacts import audit_embedded_tetrahedra
from .rational_k_segment_stars import (
    RationalKSegmentStar,build_rational_k_segment_star,
)
from .welded_k_corner_edge import (
    _oriented,_boundary_faces,_edges,_clip_segment_tet,
    _closed_box_separates,_interior_source_parameter,
)
from .spatial_twist_blocks import Point3

Tet=tuple[int,int,int,int]
Tri=tuple[int,int,int]


def _permute_ring(ring:tuple[int,...],rotation:int,reverse:bool)->tuple[int,...]:
    if (len(ring)!=4 or len(set(ring))!=4
            or type(rotation) is not int or rotation not in range(4)
            or type(reverse) is not bool):
        raise ValueError("physical oriented four-vertex cap ring required")
    return tuple(ring[(rotation+(-i if reverse else i))%4]
                 for i in range(4))


def _assemble(
    center:DoublyTruncatedKCorner,
    incoming:RationalKSegmentStar,
    outgoing:RationalKSegmentStar,
    incoming_rotation:int,outgoing_rotation:int,
    incoming_reverse:bool,outgoing_reverse:bool,
)->tuple[tuple[Point3,...],tuple[Tet,...],tuple[Tet,...],
         tuple[Tet,...],tuple[Tri,...],tuple[Tri,...],
         tuple[Tri,...],tuple[Tri,...]]:
    """Exact shared physical vertices; prism cells use uniform side diagonals."""
    for rot,rev in ((incoming_rotation,incoming_reverse),
                    (outgoing_rotation,outgoing_reverse)):
        _permute_ring((0,1,2,3),rot,rev)
    points=list(center.vertices)
    indexes={p:i for i,p in enumerate(points)}
    def add(p:Point3)->int:
        if p not in indexes:
            indexes[p]=len(points)
            points.append(p)
        return indexes[p]
    ids_in=tuple(add(p) for p in incoming.vertices)
    ids_out=tuple(add(p) for p in outgoing.vertices)
    vertices=tuple(points)

    def bridge(
        near_center:int,near_ring:tuple[int,...],
        far_center:int,far_ring:tuple[int,...],
    )->tuple[tuple[Tet,...],tuple[Tri,...],tuple[Tri,...]]:
        cells=[];faces_near=[];faces_far=[]
        for i in range(4):
            j=(i+1)%4
            a,b,c=near_center,near_ring[i],near_ring[j]
            A,B,C=far_center,far_ring[i],far_ring[j]
            faces_near.append((a,b,c))
            faces_far.append((A,B,C))
            for cell in ((a,b,c,C),(a,b,B,C),(a,A,B,C)):
                cells.append(_oriented(vertices,cell))
        return tuple(cells),tuple(faces_near),tuple(faces_far)

    # Original incoming K direction is from its tube FAR section
    # toward the existing shared source knot-corner cap.
    tube_far=ids_in[5]
    tube_far_ring=tuple(ids_in[i] for i in (6,7,8,9))
    incoming_cap_ring=tuple(face[1] for face in center.incoming_faces)
    incoming_cells,in_square,in_diamond=bridge(
        tube_far,tube_far_ring,
        center.incoming_center,
        _permute_ring(incoming_cap_ring,
                      incoming_rotation,incoming_reverse),
    )

    outgoing_cap_ring=tuple(face[1] for face in center.outgoing_faces)
    tube_near=ids_out[0]
    tube_near_ring=tuple(ids_out[i] for i in (1,2,3,4))
    outgoing_cells,out_diamond,out_square=bridge(
        center.outgoing_center,outgoing_cap_ring,
        tube_near,
        _permute_ring(tube_near_ring,
                      outgoing_rotation,outgoing_reverse),
    )

    cells=(
        tuple(tuple(ids_in[i] for i in tet) for tet in incoming.tetrahedra)
        +incoming_cells+center.tetrahedra+outgoing_cells+
        tuple(tuple(ids_out[i] for i in tet) for tet in outgoing.tetrahedra)
    )
    return (vertices,cells,incoming_cells,outgoing_cells,
            in_square,in_diamond,out_diamond,out_square)


@dataclass(frozen=True)
class SharedCornerTwoKDirections:
    """ONE verified PL 3-ball around K edge-in → shared vertex → edge-out."""
    corner:DoublyTruncatedKCorner
    incoming:RationalKSegmentStar
    outgoing:RationalKSegmentStar
    incoming_rotation:int
    outgoing_rotation:int
    incoming_reverse:bool
    outgoing_reverse:bool
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tet,...]
    incoming_bridge:tuple[Tet,...]
    outgoing_bridge:tuple[Tet,...]
    incoming_square:tuple[Tri,...]
    incoming_diamond:tuple[Tri,...]
    outgoing_diamond:tuple[Tri,...]
    outgoing_square:tuple[Tri,...]
    status:Literal["TWO_K_RAYS_LOCAL_NOT_SOURCE_AMBIENT"]=(
        "TWO_K_RAYS_LOCAL_NOT_SOURCE_AMBIENT"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.corner,DoublyTruncatedKCorner)
                or not isinstance(self.incoming,RationalKSegmentStar)
                or not isinstance(self.outgoing,RationalKSegmentStar)
                or self.status!="TWO_K_RAYS_LOCAL_NOT_SOURCE_AMBIENT"):
            raise ValueError("shared K corner/tube types and source gate required")
        c=self.corner.original
        if (self.incoming.candidate!=c.candidate
                or self.outgoing.candidate!=c.candidate
                or self.incoming.arc!=c.arc or self.outgoing.arc!=c.arc
                or self.incoming.segment!=c.corner-1
                or self.outgoing.segment!=c.corner):
            raise ValueError("K segments must be ACTUAL incident edges of shared corner")
        self.corner.audit()
        self.incoming.audit()
        self.outgoing.audit()
        canonical=_assemble(
            self.corner,self.incoming,self.outgoing,
            self.incoming_rotation,self.outgoing_rotation,
            self.incoming_reverse,self.outgoing_reverse)
        actual=(self.vertices,self.tetrahedra,
                self.incoming_bridge,self.outgoing_bridge,
                self.incoming_square,self.incoming_diamond,
                self.outgoing_diamond,self.outgoing_square)
        if canonical!=actual:
            raise ValueError("noncanonical or forged rational 3D two-ray seam")
        if (len(self.vertices)!=35 or len(self.tetrahedra)!=72 or
                len(self.incoming_bridge)!=12 or
                len(self.outgoing_bridge)!=12):
            raise ValueError("incomplete doubly welded exact K two-edge mesh")
        topology=verify_connected_tetrahedral_ball(
            self.vertices,self.tetrahedra)
        if (topology.vertices!=35 or topology.tetrahedra!=72
                or topology.boundary_components!=1
                or topology.boundary_euler!=2):
            raise ValueError("two K rays through shared corner are not ONE PL 3-ball")
        geometry=audit_embedded_tetrahedra(
            self.vertices,self.tetrahedra)
        if geometry.pair_count!=72*71//2:
            raise ValueError("not all tetrahedron physical intersections certified")
        exposed=_boundary_faces(self.tetrahedra)
        for faces in (self.incoming_square,self.incoming_diamond,
                      self.outgoing_diamond,self.outgoing_square):
            if len(faces)!=4 or any(
                    frozenset(face) in exposed for face in faces):
                raise ValueError("joined K attachment fan remains on exterior boundary")
        # The knot's real consecutive 3D source segments must be
        # literally ordered tetrahedral 1-subchains through the SAME P.
        physical_arcs=c.candidate.patches.strands.embedded_arcs()
        path=physical_arcs[c.arc]
        p,q,r=path[c.corner-1:c.corner+2]
        if self.corner.vertices[0]!=q:
            raise ValueError("middle original physical K corner moved")
        incoming_points=(
            self.incoming.vertices[0],self.incoming.vertices[5],
            self.corner.vertices[self.corner.incoming_center],q
        )
        outgoing_points=(
            q,self.corner.vertices[self.corner.outgoing_center],
            self.outgoing.vertices[0],self.outgoing.vertices[5]
        )
        incoming_params=tuple(
            _interior_source_parameter(p,q,x)
            for x in incoming_points[:-1])
        outgoing_params=tuple(
            _interior_source_parameter(q,r,x)
            for x in outgoing_points[1:])
        if (not incoming_params[0]<incoming_params[1]<incoming_params[2]<1
                or not 0<outgoing_params[0]<outgoing_params[1]
                                  <outgoing_params[2]):
            raise ValueError("source K ray geometric parameter ordering invalid")
        positions={v:i for i,v in enumerate(self.vertices)}
        edges=_edges(self.tetrahedra)
        for polyline in (incoming_points,outgoing_points):
            if any(frozenset((positions[a],positions[b])) not in edges
                   for a,b in zip(polyline,polyline[1:])):
                raise ValueError("physical consecutive K rays not real mesh edges")
        # Every original source K segment OTHER than the appropriate
        # bridge's own incident segment must be disjoint from its
        # CLOSED 3D transition simplices (including near tangencies).
        for bridge,target in (
            (self.incoming_bridge,c.corner-1),
            (self.outgoing_bridge,c.corner)
        ):
            for arc_index,arc in enumerate(physical_arcs):
                for seg_index,(a,b) in enumerate(zip(arc,arc[1:])):
                    if (arc_index,seg_index)==(c.arc,target):
                        continue
                    for cell in bridge:
                        tet=tuple(self.vertices[i] for i in cell)
                        if _clip_segment_tet(a,b,tet) is not None:
                            raise ValueError("unrelated physical K strand intersects bridge")
            for cell in bridge:
                tet=tuple(self.vertices[i] for i in cell)
                if any(not _closed_box_separates(
                        tet,cube) for cube in c.candidate.patches.strands.cubes):
                    raise ValueError("corner tube bridge intersects closed Conway ball")
        return topology

    def require_source_full_k(self)->None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: one K vertex with two incident proper "
            "subedges is not the whole K, ribbon/Conway ambient, "
            "source slitting, framed isotopy or ring R+/-")


def build_shared_corner_two_k_directions(
    atom:AtomicInstruction,arc:int,corner:int,
    max_radius:Fraction=Fraction(1,4),
)->SharedCornerTwoKDirections:
    """Bounded exact search: lack of positive witness is UNVERIFIED."""
    if (not isinstance(atom,AtomicInstruction)
            or type(arc) is not int or type(corner) is not int
            or type(max_radius) is not Fraction
            or not 0<max_radius<=Fraction(1,4)):
        raise ValueError("typed actual K corner / guarded exact radius required")
    if atom.ring_count!=0 or atom.valence<2:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: ringed or valence-one source "
            "K neighborhoods are not physically realized")
    for exponent in range(8):
        radius=max_radius/Fraction(2**exponent)
        try:
            center=build_doubly_truncated_k_corner(
                atom,arc,corner,radius)
            first=build_rational_k_segment_star(
                atom,arc,corner-1,radius)
            second=build_rational_k_segment_star(
                atom,arc,corner,radius)
        except (ValueError,NotImplementedError):
            continue
        # Each independently instantiated mesh may shrink further.
        # The construction must be explicit about the chosen source
        # geometry and all three exact radii.
        if not (center.original.radius==first.radius==second.radius==radius):
            continue
        for ir in (False,True):
            for ix in range(4):
                for orev in (False,True):
                    for ox in range(4):
                        try:
                            args=_assemble(center,first,second,
                                           ix,ox,ir,orev)
                            return SharedCornerTwoKDirections(
                                center,first,second,ix,ox,ir,orev,*args)
                        except ValueError:
                            continue
    raise ValueError(
        "UNVERIFIED_GEOMETRY: finite rational search found no "
        "independently certified shared corner/two-edge PL 3-ball; "
        "this is not a nonembedding or nonisotopy theorem")
