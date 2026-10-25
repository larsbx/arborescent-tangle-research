"""A5.2e2b4c2e3: ONE actual K tube between TWO shared dual-cap vertices.

A5.2e2b4c2e2's individual double-capped K-corner models share a
positive-volume full tube if naively unioned across their common
source K edge. This module uses only ONE actual physical K segment
tube and ONE doubly truncated PL 3-ball at EACH of its two true
two-valent K endpoints. Two four-sector transition prisms weld both
ends via complete positive-area physical mesh triangle interfaces.

Exact Q^3, independent full PL vertex and edge links, closed convex
tetrahedron/tetrahedron collisions, and physical original K subchain
checked on 84 tetrahedra. Not the source complete ambient (M,K).
"""
from __future__ import annotations

from collections import Counter
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
from .truncated_k_corner_caps import _link_cycle
from .shared_corner_two_k_rays import _permute_ring
from .welded_k_corner_edge import (
    _oriented,_boundary_faces,_edges,_clip_segment_tet,
    _closed_box_separates,_interior_source_parameter,
)
from .spatial_twist_blocks import Point3

Tet=tuple[int,int,int,int]
Tri=tuple[int,int,int]


def _assemble(
    left:DoublyTruncatedKCorner,
    tube:RationalKSegmentStar,
    right:DoublyTruncatedKCorner,
    left_rotation:int,right_rotation:int,
    left_reverse:bool,right_reverse:bool,
)->tuple[tuple[Point3,...],tuple[Tet,...],tuple[Tet,...],
         tuple[Tet,...],tuple[Tri,...],tuple[Tri,...],
         tuple[Tri,...],tuple[Tri,...]]:
    """Build one physical global vertex table, TWO welded four-prism seams."""
    for rotation,reverse in ((left_rotation,left_reverse),
                             (right_rotation,right_reverse)):
        _permute_ring((0,1,2,3),rotation,reverse)
    vertices=list(left.vertices)
    ids={p:i for i,p in enumerate(vertices)}
    def add(p:Point3)->int:
        if p not in ids:
            ids[p]=len(vertices)
            vertices.append(p)
        return ids[p]
    tube_ids=tuple(add(p) for p in tube.vertices)
    right_ids=tuple(add(p) for p in right.vertices)
    points=tuple(vertices)

    def prism(center1:int,ring1:tuple[int,...],
              center2:int,ring2:tuple[int,...]
             )->tuple[tuple[Tet,...],tuple[Tri,...],tuple[Tri,...]]:
        if (len(ring1)!=4 or len(ring2)!=4
                or len(set(ring1))!=4 or len(set(ring2))!=4):
            raise ValueError("nondegenerate conforming 2D section rings needed")
        cells=[];first=[];last=[]
        for i in range(4):
            j=(i+1)%4
            a,b,c=center1,ring1[i],ring1[j]
            A,B,C=center2,ring2[i],ring2[j]
            first.append((a,b,c))
            last.append((A,B,C))
            for tet in ((a,b,c,C),(a,b,B,C),(a,A,B,C)):
                cells.append(_oriented(points,tet))
        return tuple(cells),tuple(first),tuple(last)

    left_ring=tuple(face[1] for face in left.outgoing_faces)
    tube_near=tuple(tube_ids[i] for i in (1,2,3,4))
    near,near_diamond,near_square=prism(
        left.outgoing_center,left_ring,tube_ids[0],
        _permute_ring(tube_near,left_rotation,left_reverse)
    )
    tube_far=tuple(tube_ids[i] for i in (6,7,8,9))
    right_ring=tuple(right_ids[face[1]] for face in right.incoming_faces)
    far,far_square,far_diamond=prism(
        tube_ids[5],tube_far,right_ids[right.incoming_center],
        _permute_ring(right_ring,right_rotation,right_reverse)
    )
    all_tets=(
        left.tetrahedra+near+
        tuple(tuple(tube_ids[i] for i in tet) for tet in tube.tetrahedra)
        +far+
        tuple(tuple(right_ids[i] for i in tet) for tet in right.tetrahedra)
    )
    return (points,all_tets,near,far,
            near_diamond,near_square,far_square,far_diamond)


@dataclass(frozen=True)
class SharedTubeBetweenKVertices:
    """A genuine embedded 3-ball: dual-cap vertex–one K tube–dual-cap vertex."""
    left:DoublyTruncatedKCorner
    tube:RationalKSegmentStar
    right:DoublyTruncatedKCorner
    left_rotation:int
    right_rotation:int
    left_reverse:bool
    right_reverse:bool
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tet,...]
    near_bridge:tuple[Tet,...]
    far_bridge:tuple[Tet,...]
    near_diamond:tuple[Tri,...]
    near_square:tuple[Tri,...]
    far_square:tuple[Tri,...]
    far_diamond:tuple[Tri,...]
    status:Literal["ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT"]=(
        "ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.left,DoublyTruncatedKCorner)
                or not isinstance(self.right,DoublyTruncatedKCorner)
                or not isinstance(self.tube,RationalKSegmentStar)
                or self.status!="ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT"):
            raise ValueError("actual source-candidate K mesh types required")
        L,R=self.left.original,self.right.original
        if (L.candidate!=R.candidate or L.candidate!=self.tube.candidate
                or L.arc!=R.arc or L.arc!=self.tube.arc
                or R.corner!=L.corner+1 or self.tube.segment!=L.corner):
            raise ValueError("adjacent genuine K source vertices and same tube required")
        # Cheap canonical-geometry comparison comes BEFORE expensive
        # independently replayed component/collision validation. This
        # prevents forged fields from forcing repeated O(N^2) audits;
        # every accepted candidate still receives ALL audits below.
        expected=_assemble(
            self.left,self.tube,self.right,
            self.left_rotation,self.right_rotation,
            self.left_reverse,self.right_reverse)
        actual=(self.vertices,self.tetrahedra,self.near_bridge,self.far_bridge,
                self.near_diamond,self.near_square,self.far_square,self.far_diamond)
        if actual!=expected:
            raise ValueError("forged physical K mesh or incorrect cap mapping")
        self.left.audit();self.tube.audit();self.right.audit()
        if (len(self.vertices)!=40 or len(self.tetrahedra)!=84 or
                len(self.near_bridge)!=12 or len(self.far_bridge)!=12):
            raise ValueError("missing shared tube / corners / 3D transition cells")
        stats=verify_connected_tetrahedral_ball(self.vertices,self.tetrahedra)
        if (stats.boundary_components,stats.boundary_euler)!=(1,2):
            raise ValueError("joined dual-cap K vertices are not one PL 3-ball")
        collision=audit_embedded_tetrahedra(self.vertices,self.tetrahedra)
        if collision.pair_count!=3486:
            raise ValueError("not all 3486 tetrahedron physical contacts checked")
        physical_boundary=_boundary_faces(self.tetrahedra)
        for faces in (self.near_diamond,self.near_square,
                      self.far_square,self.far_diamond):
            if len(faces)!=4 or any(frozenset(face) in physical_boundary
                                    for face in faces):
                raise ValueError("true positive-area K face failed to weld internally")
        # The *unused* left incoming cap and right outgoing cap
        # must stay exposed for assembling a larger K arc later.
        indexes={p:i for i,p in enumerate(self.vertices)}
        for core,cap in ((self.left,self.left.incoming_faces),
                         (self.right,self.right.outgoing_faces)):
            for face in cap:
                physical=frozenset(indexes[core.vertices[i]] for i in face)
                if physical not in physical_boundary:
                    raise ValueError("unused exterior corner cap was buried")
        arc=L.candidate.patches.strands.embedded_arcs()[L.arc]
        p,q=arc[L.corner:L.corner+2]
        centerline=(
            self.left.vertices[0],
            self.left.vertices[self.left.outgoing_center],
            self.tube.vertices[0],
            self.tube.vertices[5],
            self.right.vertices[self.right.incoming_center],
            self.right.vertices[0]
        )
        if centerline[0]!=p or centerline[-1]!=q:
            raise ValueError("original source K endpoints moved")
        times=tuple(_interior_source_parameter(p,q,v)
                    for v in centerline[1:-1])
        if not 0<times[0]<times[1]<times[2]<times[3]<1:
            raise ValueError("nonmonotone physical K chain")
        actual_edges=_edges(self.tetrahedra)
        for a,b in zip(centerline,centerline[1:]):
            ia,ib=indexes[a],indexes[b]
            if frozenset((ia,ib)) not in actual_edges:
                raise ValueError("physical source K chain not a mesh 1-skeleton")
            if _link_cycle(self.tetrahedra,ia,ib)[0]<4:
                raise ValueError("K interior edge has noncircular tetrahedral link")
        # Both corners' OTHER source K rays are retained as mesh
        # edges to their UNUSED exposed attachment caps.
        for core,end in ((self.left,self.left.incoming_center),
                         (self.right,self.right.outgoing_center)):
            i,j=indexes[core.vertices[0]],indexes[core.vertices[end]]
            if frozenset((i,j)) not in actual_edges:
                raise ValueError("unused original source K corner ray was lost")
        paths=L.candidate.patches.strands.embedded_arcs()
        for transition in (self.near_bridge,self.far_bridge):
            for ai,original in enumerate(paths):
                for si,(a,b) in enumerate(zip(original,original[1:])):
                    if (ai,si)==(L.arc,L.corner):
                        continue
                    for tet in transition:
                        if _clip_segment_tet(a,b,tuple(self.vertices[i] for i in tet)
                                             ) is not None:
                            raise ValueError("other physical source K strand hits tube weld")
            for tet in transition:
                physical=tuple(self.vertices[i] for i in tet)
                if any(not _closed_box_separates(physical,cube)
                       for cube in L.candidate.patches.strands.cubes):
                    raise ValueError("welded transition touches closed Conway solid")
        return stats

    def require_source_ambient(self)->None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: only one physical K edge between TWO "
            "source-candidate shared corners is conforming; whole K, "
            "the source ribbon/Conway ambient, four-puncture slit, "
            "framed isotopy and ring R+/- remain open")


def build_one_shared_k_tube(
    atom:AtomicInstruction,arc:int,first_corner:int,
    max_radius:Fraction=Fraction(1,4),
)->SharedTubeBetweenKVertices:
    if (not isinstance(atom,AtomicInstruction)
            or type(arc) is not int or type(first_corner) is not int
            or type(max_radius) is not Fraction
            or not Fraction(0)<max_radius<=Fraction(1,4)):
        raise ValueError("actual source K endpoints and rational radius required")
    if atom.ring_count!=0 or atom.valence<2:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: ringed / valence-one source K geometry unavailable")
    for step in range(8):
        radius=max_radius/Fraction(2**step)
        try:
            left=build_doubly_truncated_k_corner(
                atom,arc,first_corner,radius)
            right=build_doubly_truncated_k_corner(
                atom,arc,first_corner+1,radius)
            tube=build_rational_k_segment_star(
                atom,arc,first_corner,radius)
        except (ValueError,NotImplementedError):
            continue
        if not (left.original.radius==tube.radius==right.original.radius==radius):
            continue
        for rev_a in (False,True):
            for rot_a in range(4):
                for rev_b in (False,True):
                    for rot_b in range(4):
                        try:
                            mesh=_assemble(left,tube,right,rot_a,rot_b,rev_a,rev_b)
                            # Reject invalid 3D interface arrangements
                            # BEFORE rerunning component audits + the
                            # 3,486-pair nonlocal collision proof. A
                            # successful final object STILL checks all.
                            verify_connected_tetrahedral_ball(mesh[0],mesh[1])
                            exposed=_boundary_faces(mesh[1])
                            if any(frozenset(face) in exposed
                                   for faces in mesh[4:8] for face in faces):
                                continue
                            return SharedTubeBetweenKVertices(
                                left,tube,right,rot_a,rot_b,rev_a,rev_b,*mesh)
                        except ValueError:
                            continue
    raise ValueError(
        "UNVERIFIED_GEOMETRY: finite exact rational shared-tube search "
        "found no independently verified PL 3-ball, not a proof "
        "that the source K has no embedding/isotopy")
