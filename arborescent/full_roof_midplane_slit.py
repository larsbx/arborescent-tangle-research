"""A5.2e2b4c2d1: complete two-sided PL split of ONE connected 3-ball.

A5.2e2b4c2b's actual tetrahedral roof slab has a full internal
midplane disk z=4. This module cuts the *ENTIRE* separating disk,
including connector intervals, duplicates every cut edge/vertex
by 3D side, and checks the two resulting oriented 3-balls and exact
inverse 3D reglue. Unlike earlier 'face barriers', no cut vertex
or link remains shared between the two halves.

THEOREM SCOPE: ONLY this slab. Its z-range is >= 15/4 and its
existing proper K has z<=1, so the full cut has ZERO K punctures.
It is not a source-identified plumbing patch, not a new 4-punctured
Conway sphere, and not the global ambient (M,K) or an R+/- isotopy.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import (
    ConnectedRoofAmbientSlab, _face_incidence, _connected,
    verify_connected_tetrahedral_ball, build_connected_roof_ambient_slab,
)
from .local_plumbing_slit import _tet_faces, _cross, _vec
from .spatial_twist_blocks import Point3

Side = Literal["LOWER","UPPER","EXTERIOR"]
CutTri = tuple["SlitVertex","SlitVertex","SlitVertex"]
CutTet = tuple["SlitVertex","SlitVertex","SlitVertex","SlitVertex"]


@dataclass(frozen=True)
class SlitVertex:
    original: int
    side: Side

    def __post_init__(self) -> None:
        if (type(self.original) is not int or self.original<0
                or self.side not in ("LOWER","UPPER","EXTERIOR")):
            raise ValueError("invalid physical cut vertex/side")


def _tet_layer(slab: ConnectedRoofAmbientSlab, tet: tuple[int,int,int,int]
               ) -> Literal["LOWER","UPPER"]:
    mid=Fraction(4)
    heights=[Fraction(slab.vertices[i].z) for i in tet]
    if min(heights)>=mid and max(heights)>mid:
        return "UPPER"
    if max(heights)<=mid and min(heights)<mid:
        return "LOWER"
    raise ValueError("tetrahedron crosses cut disk interior")


def _tag(slab: ConnectedRoofAmbientSlab, i: int, side: str) -> SlitVertex:
    z=Fraction(slab.vertices[i].z)
    if z==4:
        return SlitVertex(i,side)
    if (side=="LOWER" and z<4) or (side=="UPPER" and z>4):
        return SlitVertex(i,"EXTERIOR")
    raise ValueError("inconsistent cut-side vertex position")


def _canonical_cells(slab: ConnectedRoofAmbientSlab
                     ) -> tuple[tuple[CutTet,...],tuple[CutTet,...]]:
    result: dict[str,list[CutTet]]={"LOWER":[],"UPPER":[]}
    for tet in slab.tetrahedra:
        side=_tet_layer(slab,tet)
        result[side].append(tuple(_tag(slab,i,side) for i in tet))
    return tuple(result["LOWER"]),tuple(result["UPPER"])


def _outward_faces(cells: tuple[CutTet,...]) -> tuple[CutTri,...]:
    incidences: dict[frozenset[SlitVertex],list[CutTri]]=defaultdict(list)
    for tet in cells:
        for face in _tet_faces(tet):
            incidences[frozenset(face)].append(face)
    if any(len(occurrences)>2 for occurrences in incidences.values()):
        raise ValueError("cut 3D cell face has excessive incidence")
    return tuple(occurrences[0] for occurrences in incidences.values()
                 if len(occurrences)==1)


def _midplane_boundary(slab: ConnectedRoofAmbientSlab,
                       cells: tuple[CutTet,...]) -> tuple[CutTri,...]:
    faces=[face for face in _outward_faces(cells)
           if all(slab.vertices[v.original].z==4 for v in face)]
    return tuple(sorted(faces,key=lambda f:tuple(sorted(v.original for v in f))))


def _side_topology(slab: ConnectedRoofAmbientSlab,
                   cells: tuple[CutTet,...]):
    """Derive real compact 3D triangulation FROM supplied cut cells."""
    nodes={v for tet in cells for v in tet}
    if not nodes:
        raise ValueError("empty cut 3D component")
    # Build own per-half vertex table, not original 3D global indices:
    # midplane vertices must not remain identified between halves.
    ordered=tuple(sorted(nodes,key=lambda v:(v.original,v.side)))
    lookup={v:i for i,v in enumerate(ordered)}
    verts=tuple(slab.vertices[v.original] for v in ordered)
    tets=tuple(tuple(lookup[v] for v in tet) for tet in cells)
    return verify_connected_tetrahedral_ball(verts,tets)


def _disk_topology(cut: tuple[CutTri,...]) -> tuple[int,int,int]:
    """Independent simplicial disk check on the actual cut face table.

    Returns V,E,F. Its boundary must be a single simple closed
    polygonal circuit and the interior a connected triangulated disk.
    """
    if not cut:
        raise ValueError("no cut disk faces")
    edges=Counter()
    triangle_links: dict[frozenset[SlitVertex],list[int]]=defaultdict(list)
    vertices=set()
    for fi,face in enumerate(cut):
        if len(set(face))!=3:
            raise ValueError("degenerate cut triangle")
        vertices.update(face)
        for a,b in zip(face,face[1:]+face[:1]):
            key=frozenset((a,b))
            edges[key]+=1
            triangle_links[key].append(fi)
    if any(m not in (1,2) for m in edges.values()):
        raise ValueError("nonmanifold cut-disk edge")
    dual: dict[int,list[int]]=defaultdict(list)
    boundary_graph: dict[SlitVertex,list[SlitVertex]]=defaultdict(list)
    for edge,faces in triangle_links.items():
        if len(faces)==2:
            a,b=faces
            dual[a].append(b)
            dual[b].append(a)
        else:
            a,b=tuple(edge)
            boundary_graph[a].append(b)
            boundary_graph[b].append(a)
    if (_connected(dual,set(range(len(cut))))!=1
            or len(vertices)-len(edges)+len(cut)!=1
            or not boundary_graph
            or any(len(nbr)!=2 for nbr in boundary_graph.values())
            or _connected(boundary_graph,set(boundary_graph))!=1):
        raise ValueError("cut disk not a connected PL 2-disk")
    return len(vertices),len(edges),len(cut)


@dataclass(frozen=True)
class FullRoofSlitStats:
    lower_tetrahedra: int
    upper_tetrahedra: int
    duplicated_midplane_vertices: int
    cut_triangles_per_side: int
    each_ball_euler: int
    cut_disk_euler: int
    knot_piercings: int
    source_cut_four_punctured_spheres: int
    volume_lower: Fraction
    volume_upper: Fraction


@dataclass(frozen=True)
class FullRoofSlit:
    """Actual edge- AND vertex-split oriented PL cut of a roof 3-ball."""
    slab: ConnectedRoofAmbientSlab
    lower: tuple[CutTet,...]
    upper: tuple[CutTet,...]
    lower_faces: tuple[CutTri,...]
    upper_faces: tuple[CutTri,...]
    status: Literal["FULL_ROOF_PL_SLIT_NOT_SOURCE_CONWAY"] = (
        "FULL_ROOF_PL_SLIT_NOT_SOURCE_CONWAY"
    )

    def __post_init__(self) -> None:
        self.audit()

    def audit(self) -> FullRoofSlitStats:
        if (not isinstance(self.slab,ConnectedRoofAmbientSlab)
                or self.status!="FULL_ROOF_PL_SLIT_NOT_SOURCE_CONWAY"
                or not isinstance(self.lower,tuple)
                or not isinstance(self.upper,tuple)
                or not isinstance(self.lower_faces,tuple)
                or not isinstance(self.upper_faces,tuple)):
            raise ValueError("invalid local/global PL slit types or status")
        original=self.slab.audit()
        expected_lower,expected_upper=_canonical_cells(self.slab)
        if (self.lower,self.upper)!=(expected_lower,expected_upper):
            raise ValueError("side-tagged tetrahedral cells disagree with geometry")
        m=len(self.slab.x_slices)-1
        if len(self.lower)!=6*m or len(self.upper)!=6*m:
            raise ValueError("missing 3D half-slab cells")
        if self.lower_faces!=_midplane_boundary(self.slab,self.lower):
            raise ValueError("lower cut triangles are not actual exterior faces")
        if self.upper_faces!=_midplane_boundary(self.slab,self.upper):
            raise ValueError("upper cut triangles are not actual exterior faces")
        face_groups=(
            {frozenset(v.original for v in f) for f in self.lower_faces},
            {frozenset(v.original for v in f) for f in self.upper_faces}
        )
        source_mid={
            frozenset(face) for face,occ in
            _face_incidence(self.slab.tetrahedra).items()
            if len(occ)==2 and all(self.slab.vertices[i].z==4 for i in face)
        }
        if (face_groups[0]!=face_groups[1]
                or face_groups[0]!=source_mid
                or len(source_mid)!=2*m):
            raise ValueError("wrong midplane face pairing")
        for lower,upper in zip(self.lower_faces,self.upper_faces):
            rawlo=tuple(v.original for v in lower)
            rawhi=tuple(v.original for v in upper)
            if (frozenset(rawlo)!=frozenset(rawhi)
                    or not all(v.side=="LOWER" for v in lower)
                    or not all(v.side=="UPPER" for v in upper)):
                raise ValueError("cut disks share vertex copies or wrongly paired")
            normal_lo=_cross(_vec(self.slab.vertices[rawlo[0]],
                                  self.slab.vertices[rawlo[1]]),
                             _vec(self.slab.vertices[rawlo[0]],
                                  self.slab.vertices[rawlo[2]]))
            normal_hi=_cross(_vec(self.slab.vertices[rawhi[0]],
                                  self.slab.vertices[rawhi[1]]),
                             _vec(self.slab.vertices[rawhi[0]],
                                  self.slab.vertices[rawhi[2]]))
            if not normal_lo[2]>0 or not normal_hi[2]<0:
                raise ValueError("cut faces have incorrect outward normals")
        low_disk=_disk_topology(self.lower_faces)
        high_disk=_disk_topology(self.upper_faces)
        if (low_disk!=high_disk or low_disk!=(2*(m+1),4*m+1,2*m)):
            raise ValueError("cut plane failed simplicial disk topology")
        low=_side_topology(self.slab,self.lower)
        high=_side_topology(self.slab,self.upper)
        if (low.boundary_euler!=2 or high.boundary_euler!=2
                or low.boundary_components!=1 or high.boundary_components!=1
                or low.volume+high.volume!=original.volume
                or low.volume!=high.volume):
            raise ValueError("cut parts are not complementary oriented PL 3-balls")
        mid_nodes={v.original for f in self.lower_faces for v in f}
        lower_nodes={v for tet in self.lower for v in tet}
        upper_nodes={v for tet in self.upper for v in tet}
        if (len(mid_nodes)!=2*(m+1)
                or any(SlitVertex(i,"LOWER") not in lower_nodes
                       or SlitVertex(i,"UPPER") not in upper_nodes for i in mid_nodes)
                or lower_nodes & upper_nodes):
            raise ValueError("slit did not duplicate all physical cut vertices")
        # Exact inverse cell-level reglue. All coordinates are original
        # and the split side tags are eliminated ONLY after checking.
        rejoined=tuple(tuple(v.original for v in tet)
                       for tet in (*self.lower,*self.upper))
        if (set(map(frozenset,rejoined))!=
                set(map(frozenset,self.slab.tetrahedra))
                or len(rejoined)!=len(self.slab.tetrahedra)):
            raise ValueError("reglue does not reproduce exact 3D cells")
        # This entire slab was independently proven to avoid K
        # (z>=15/4 vs source-candidate K z<=1). The new cut plane
        # cannot create source four-punctured Conway boundaries.
        return FullRoofSlitStats(
            len(self.lower),len(self.upper),len(mid_nodes),
            len(self.lower_faces),2,1,0,0,low.volume,high.volume)

    def require_source_conway_certificate(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: an edge/vertex-complete PL slit of "
            "the roof 3-ball has ZERO K punctures and is NOT the "
            "Bonahon-Siebenmann source slit in the ambient (M,K)"
        )


def slit_complete_roof_midplane(
    atom: AtomicInstruction, thickness: Fraction=Fraction(1,16)
) -> FullRoofSlit:
    """Construct both true PL halves of an audited chosen-roof 3-ball."""
    slab=build_connected_roof_ambient_slab(atom,thickness)
    lower,upper=_canonical_cells(slab)
    return FullRoofSlit(
        slab,lower,upper,
        _midplane_boundary(slab,lower),
        _midplane_boundary(slab,upper)
    )
