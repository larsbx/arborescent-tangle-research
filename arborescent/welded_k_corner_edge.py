"""A5.2e2b4c2d2d2: physically welded exact-Q3 K corner-to-edge 3-ball.

A PR66 truncated corner has a real 4-triangle diamond disk. A PR63
edge tube has a real 4-triangle square disk. This module triangulates
FOUR corresponding 3D triangular prisms, then welds the corner,
transition and tube at IDENTICAL actual physical mesh triangles.
A generic exact tetrahedron-contact kernel separately rejects any
nonconforming interior collision or tangency.

This positive result concerns *one* existing source-candidate K
corner and its outgoing proper subedge, not the whole-K global
ambient, full band/conway face conformance or source Fig12.10.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .exact_tetrahedron_contacts import audit_embedded_tetrahedra, _planes
from .local_plumbing_slit import _volume6
from .rational_k_segment_stars import (
    RationalKSegmentStar,build_rational_k_segment_star,
    _xyz,_minus,_cross,_dot,
)
from .truncated_k_corner_caps import (
    TruncatedKCorner,build_truncated_k_corner,
)
from .spatial_twist_blocks import Point3

Tet=tuple[int,int,int,int]
Triangle=tuple[int,int,int]


def _oriented(vertices:tuple[Point3,...],tet:Tet)->Tet:
    det=_volume6(*(vertices[i] for i in tet))
    if det==0:
        raise ValueError("transition tetrahedron has zero spatial volume")
    return tet if det>0 else (tet[0],tet[2],tet[1],tet[3])


def _clip_segment_tet(
    a:Point3,b:Point3,tet:tuple[Point3,Point3,Point3,Point3]
)->tuple[Fraction,Fraction]|None:
    """Exact closed segment clipped against all tetrahedron faces."""
    if a==b:
        raise ValueError("proper physical K segment required")
    lo,hi=Fraction(0),Fraction(1)
    pa,pb=_xyz(a),_xyz(b)
    delta=_minus(pb,pa)
    for normal,bound in _planes(tet):
        f0=_dot(normal,pa)-bound
        slope=_dot(normal,delta)
        if slope==0:
            if f0>0:
                return None
            continue
        cutoff=-f0/slope
        if slope>0:
            hi=min(hi,cutoff)
        else:
            lo=max(lo,cutoff)
        if lo>hi:
            return None
    return lo,hi


def _as_meshes(
    corner:TruncatedKCorner,
    tube:RationalKSegmentStar,
    rotation:int,
)->tuple[tuple[Point3,...],tuple[Tet,...],tuple[Tet,...],
         tuple[Triangle,...],tuple[Triangle,...]]:
    """Exact vertex welding and 4-prism staircase triangulation."""
    if type(rotation) is not int or rotation not in range(4):
        raise ValueError("four cyclic attachment alignments supported")
    vertices=list(corner.vertices)
    indexes={p:i for i,p in enumerate(vertices)}
    def add(p:Point3)->int:
        if p not in indexes:
            indexes[p]=len(vertices)
            vertices.append(p)
        return indexes[p]
    tube_ids=tuple(add(p) for p in tube.vertices)
    near_center=corner.cap_center
    far_center=tube_ids[0]
    near_ring=tuple(face[1] for face in corner.cap_faces)
    if len(set(near_ring))!=4:
        raise ValueError("physical cap must be four distinct ring corners")
    far_ring=tuple(tube_ids[1+(i+rotation)%4] for i in range(4))
    if len(set(far_ring))!=4:
        raise ValueError("tube cross-section cannot collapse")
    all_vertices=tuple(vertices)
    tets=list(corner.tetrahedra)
    bridge=[]
    near=[]
    far=[]
    for i in range(4):
        j=(i+1)%4
        a,b,c=near_center,near_ring[i],near_ring[j]
        A,B,C=far_center,far_ring[i],far_ring[j]
        near.append((a,b,c))
        far.append((A,B,C))
        # Full 3D 3-simplex staircase in one triangular prism.
        for cell in ((a,b,c,C),(a,b,B,C),(a,A,B,C)):
            bridge.append(_oriented(all_vertices,cell))
    tets.extend(bridge)
    tets.extend(tuple(tube_ids[i] for i in cell)
                for cell in tube.tetrahedra)
    return all_vertices,tuple(tets),tuple(bridge),tuple(near),tuple(far)


def _closed_box_separates(vertices:tuple[Point3,...],cube)->bool:
    """Sufficient strict convex hull/cube separation on one axis."""
    return any(
        max(Fraction(getattr(p,k)) for p in vertices)<lo
        or min(Fraction(getattr(p,k)) for p in vertices)>hi
        for k,lo,hi in (
            ("x",Fraction(cube.west),Fraction(cube.east)),
            ("y",Fraction(-3),Fraction(3)),
            ("z",Fraction(-3),Fraction(3)),
        )
    )


def _edges(tets:tuple[Tet,...])->set[frozenset[int]]:
    result=set()
    for a,b,c,d in tets:
        for i,j in ((a,b),(a,c),(a,d),(b,c),(b,d),(c,d)):
            result.add(frozenset((i,j)))
    return result


def _boundary_faces(tets:tuple[Tet,...])->set[frozenset[int]]:
    c=Counter()
    for a,b,c0,d in tets:
        for face in ((b,c0,d),(a,d,c0),(a,b,d),(a,c0,b)):
            c[frozenset(face)]+=1
    if any(n not in (1,2) for n in c.values()):
        raise ValueError("nonmanifold welded tetrahedral face incidence")
    return {face for face,n in c.items() if n==1}


def _interior_source_parameter(a:Point3,b:Point3,x:Point3)->Fraction:
    u=_minus(_xyz(b),_xyz(a))
    v=_minus(_xyz(x),_xyz(a))
    if _cross(u,v)!=(0,0,0):
        raise ValueError("fabricated knot vertex does not lie on source edge")
    i=next(i for i in range(3) if u[i])
    t=v[i]/u[i]
    if not 0<t<1:
        raise ValueError("knot mesh edge midpoint lies outside original K segment")
    return t


@dataclass(frozen=True)
class WeldedKCornerToEdge:
    """A REAL connected 3-ball with a continuous conforming K 1-chain."""
    corner:TruncatedKCorner
    tube:RationalKSegmentStar
    rotation:int
    vertices:tuple[Point3,...]
    tetrahedra:tuple[Tet,...]
    bridge_tetrahedra:tuple[Tet,...]
    near_faces:tuple[Triangle,...]
    far_faces:tuple[Triangle,...]
    status:Literal["ONE_K_CORNER_EDGE_GLUE_NOT_SOURCE_AMBIENT"]=(
        "ONE_K_CORNER_EDGE_GLUE_NOT_SOURCE_AMBIENT"
    )

    def __post_init__(self)->None:
        self.audit()

    def audit(self):
        if (not isinstance(self.corner,TruncatedKCorner)
                or not isinstance(self.tube,RationalKSegmentStar)
                or self.status!="ONE_K_CORNER_EDGE_GLUE_NOT_SOURCE_AMBIENT"
                or self.corner.original.candidate!=self.tube.candidate
                or self.corner.original.arc!=self.tube.arc
                or self.corner.original.corner!=self.tube.segment):
            raise ValueError("physical outgoing source K segment mismatch")
        self.corner.audit()
        self.tube.audit()
        expected=_as_meshes(self.corner,self.tube,self.rotation)
        if (self.vertices,self.tetrahedra,self.bridge_tetrahedra,
                self.near_faces,self.far_faces)!=expected:
            raise ValueError("changed welded 3D tetrahedral mesh")
        if len(self.bridge_tetrahedra)!=12:
            raise ValueError("incomplete conforming four-sector prism")
        topology=verify_connected_tetrahedral_ball(
            self.vertices,self.tetrahedra)
        if (topology.tetrahedra!=40 or topology.vertices!=21
                or topology.boundary_components!=1
                or topology.boundary_euler!=2):
            raise ValueError("connected corner-to-edge 3-ball topology fails")
        # Deliberately *independent* spatial contact proof.
        geometry=audit_embedded_tetrahedra(
            self.vertices,self.tetrahedra)
        if geometry.pair_count!=len(self.tetrahedra)*(
                len(self.tetrahedra)-1)//2:
            raise ValueError("not all 3D geometric collision pairs were checked")
        exposed=_boundary_faces(self.tetrahedra)
        if any(frozenset(f) in exposed for f in
               (*self.near_faces,*self.far_faces)):
            raise ValueError("cap failed to weld as an INTERNAL tetrahedral face")
        # The ORIGINAL K centerline, not a made-up diagram arc, must
        # be a continuous 1-skeleton path through all three 3-balls.
        arcs=self.corner.original.candidate.patches.strands.embedded_arcs()
        a=arcs[self.tube.arc][self.tube.segment]
        b=arcs[self.tube.arc][self.tube.segment+1]
        cap=self.corner.vertices[self.corner.cap_center]
        tube_start=self.tube.vertices[0]
        tube_end=self.tube.vertices[5]
        center=self.corner.vertices[self.corner.k_center]
        if center!=a:
            raise ValueError("physical source K corner differs")
        times=tuple(_interior_source_parameter(
            a,b,x) for x in (cap,tube_start,tube_end))
        if not times[0]<times[1]<times[2]:
            raise ValueError("K centerline reverses or skips at physical weld")
        node={p:i for i,p in enumerate(self.vertices)}
        required=(center,cap,tube_start,tube_end)
        mesh_edges=_edges(self.tetrahedra)
        if any(frozenset((node[p],node[q])) not in mesh_edges
               for p,q in zip(required,required[1:])):
            raise ValueError("physical K core chain not contained in 3D 1-skeleton")
        # Independent contact checker against every original source
        # K segment, including the incoming ray to this K corner.
        for ai,arc in enumerate(arcs):
            for si,(p,q) in enumerate(zip(arc,arc[1:])):
                if (ai,si)==(self.tube.arc,self.tube.segment):
                    continue
                # Incoming K edge belongs to the CORNER, not the
                # transition or outgoing segment tube.
                for tet in (*self.bridge_tetrahedra,
                            *(tuple(self.tube.vertices[i] for i in t)
                              for t in self.tube.tetrahedra)):
                    if _clip_segment_tet(p,q,tuple(
                        self.vertices[i] for i in tet)
                        if all(type(i) is int for i in tet) else tet
                    ) is not None:
                        raise ValueError("unrelated K intersects glued bridge/tube")
        for tet in self.bridge_tetrahedra:
            physical=tuple(self.vertices[i] for i in tet)
            if any(not _closed_box_separates(physical,cube)
                   for cube in self.corner.original.candidate.patches.strands.cubes):
                raise ValueError("bridge lacks strict Conway solid-ball clearance")
        return topology

    def require_source_ambient(self)->None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: only one K corner and outgoing edge "
            "are conforming; full source ambient, 4-punctured cutting "
            "and ring R+/- remain open")


def build_welded_k_corner_to_edge(
    atom:AtomicInstruction,arc:int,corner:int,
    max_radius:Fraction=Fraction(1,4)
)->WeldedKCornerToEdge:
    """Try exact rational radii/4 cyclic mappings, fail closed if none."""
    if (not isinstance(atom,AtomicInstruction)
            or type(arc) is not int or type(corner) is not int
            or type(max_radius) is not Fraction
            or not 0<max_radius<=Fraction(1,4)):
        raise ValueError("source candidate K corner/rational radius required")
    for shift in range(10):
        r=max_radius/Fraction(2**shift)
        try:
            corner_obj=build_truncated_k_corner(atom,arc,corner,r)
            tube_obj=build_rational_k_segment_star(atom,arc,corner,r)
        except ValueError:
            continue
        for rotate in range(4):
            try:
                mesh=_as_meshes(corner_obj,tube_obj,rotate)
                return WeldedKCornerToEdge(
                    corner_obj,tube_obj,rotate,*mesh)
            except ValueError:
                continue
    raise ValueError(
        "UNVERIFIED_GEOMETRY: bounded corner-edge 3D welding search "
        "found no positive PL embedding; NOT proof of nonisotopy")
