"""A5.2e2b4c2b: a genuinely CONNECTED triangulated ambient 3-ball slab.

Given the physically embedded r=0 candidate band, construct one
exact rational PL 3-ball enclosing all n roof plumbing-chart candidates
and the unoccupied intervals between them. The 3D complex has actual
shared tetrahedron faces across the former n disjoint product charts.

A separate 3D cell verifier checks oriented positive tetrahedra,
internal faces, boundary sphere, individual vertex-link 2-disks, dual
connectedness and exact occupied prism volume. The n selected oriented
midplane patch faces are audited as two-sided internal tetra faces.

This is NOT the global S^3 / (M,K) ambient triangulation in the source:
K lies strictly below the slab, so each patch has ZERO K piercings.
Removing dual-face adjacencies is not a full PL slit at the patch
edges (no split-vertex links are produced). It is therefore forbidden
to claim four-punctured Conway boundaries, source isotopy or R+/-.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction, Port
from .disjoint_local_slit_atlas import (
    DisjointLocalSlitAtlas, build_disjoint_local_slit_atlas,
)
from .local_plumbing_slit import _product_tets, _tet_faces, _volume6
from .spatial_twist_blocks import Point3

Tet = tuple[int, int, int, int]
Face = tuple[int, int, int]


def _connected(graph: dict, nodes: set) -> int:
    remaining = set(nodes)
    count = 0
    while remaining:
        count += 1
        stack = [remaining.pop()]
        while stack:
            for neighbor in graph.get(stack.pop(), ()):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
    return count


def _oriented_opposite(p: Face, q: Face) -> bool:
    return q in ((p[1], p[0], p[2]), (p[0], p[2], p[1]),
                 (p[2], p[1], p[0]))


def _edge(a: int, b: int) -> frozenset[int]:
    return frozenset((a, b))


def _tet_map(slab: int) -> tuple[int, ...]:
    """Local 12 vertices -> one GLOBAL grid with actually shared faces."""
    local = []
    for layer in range(3):
        local.extend((
            6*slab + 2*layer,
            6*(slab+1) + 2*layer,
            6*(slab+1) + 2*layer + 1,
            6*slab + 2*layer + 1,
        ))
    return tuple(local)


def _face_incidence(tets: tuple[Tet, ...]
                   ) -> dict[frozenset[int], list[tuple[int, Face]]]:
    result: dict[frozenset[int],list[tuple[int,Face]]] = defaultdict(list)
    for i, tet in enumerate(tets):
        for face in _tet_faces(tet):
            result[frozenset(face)].append((i, face))
    return result


@dataclass(frozen=True)
class Connected3BallAudit:
    vertices: int
    tetrahedra: int
    boundary_triangles: int
    boundary_euler: int
    boundary_components: int
    interior_face_pairs: int
    volume: Fraction


def verify_connected_tetrahedral_ball(
    vertices: tuple[Point3, ...], tetrahedra: tuple[Tet, ...]
) -> Connected3BallAudit:
    """Independent, coordinate-only PL 3-ball checker for this grid.

    Validates individual vertex links as 2-disks or 2-spheres,
    not merely face incidence/Euler characteristic.
    This does NOT detect arbitrary nonlocal self-overlap of tetrahedra;
    the canonical disjoint x-slab construction is checked separately.
    """
    if (not isinstance(vertices,tuple) or not vertices
            or any(not isinstance(p,Point3) for p in vertices)
            or len(vertices)!=len(set(vertices))
            or not isinstance(tetrahedra,tuple) or not tetrahedra):
        raise ValueError("exact unique 3D vertices and nonempty tetrahedra required")
    incidence = _face_incidence(tetrahedra)
    tet_dual: dict[int,list[int]] = defaultdict(list)
    used = set()
    volume = Fraction(0)
    for i,tet in enumerate(tetrahedra):
        if (not isinstance(tet,tuple) or len(tet)!=4
                or any(type(v) is not int or not 0<=v<len(vertices) for v in tet)
                or len(set(tet))!=4 or frozenset(tet) in used):
            raise ValueError("bad or repeated tetrahedron")
        used.add(frozenset(tet))
        det = _volume6(*(vertices[v] for v in tet))
        if det<=0:
            raise ValueError("non-positive oriented tetrahedral volume")
        volume += det/6
    for pair in incidence.values():
        if len(pair) not in (1,2):
            raise ValueError("nonmanifold tetrahedron face incidence")
        if len(pair)==2:
            (a,f),(b,g)=pair
            if a==b or not _oriented_opposite(f,g):
                raise ValueError("internal 3D faces lack opposite orientations")
            tet_dual[a].append(b)
            tet_dual[b].append(a)
    if _connected(tet_dual,set(range(len(tetrahedra))))!=1:
        raise ValueError("ambient 3D cell complex is disconnected")

    boundary = [pair[0][1] for pair in incidence.values() if len(pair)==1]
    boundary_edges = Counter()
    bdgraph: dict[int,set[int]] = defaultdict(set)
    for f in boundary:
        for u,v in zip(f,f[1:]+f[:1]):
            boundary_edges[_edge(u,v)]+=1
            bdgraph[u].add(v)
            bdgraph[v].add(u)
    bdv={v for f in boundary for v in f}
    if (any(n!=2 for n in boundary_edges.values())
            or len(bdv)-len(boundary_edges)+len(boundary)!=2
            or _connected(bdgraph,bdv)!=1):
        raise ValueError("ambient boundary is not a connected Euler-2 2-sphere")

    # Independent vertex links: for each tetrahedral vertex v,
    # the link has one triangle from each incident tetrahedron.
    # The link must be a connected disk (boundary vertex) or
    # sphere (interior vertex), with no edge used >2 times.
    for v in range(len(vertices)):
        links=[tuple(u for u in tet if u!=v)
               for tet in tetrahedra if v in tet]
        if not links:
            raise ValueError("isolated ambient vertex")
        edges=Counter()
        faces_dual: dict[int,list[int]] = defaultdict(list)
        shared_edges: dict[frozenset[int],list[int]] = defaultdict(list)
        nodes=set()
        for fid,tri in enumerate(links):
            if len(set(tri))!=3:
                raise ValueError("collapsed vertex link")
            nodes.update(tri)
            for a,b in zip(tri,tri[1:]+tri[:1]):
                key=_edge(a,b)
                edges[key]+=1
                shared_edges[key].append(fid)
        if any(count not in (1,2) for count in edges.values()):
            raise ValueError("nonmanifold ambient vertex link edge")
        for faces in shared_edges.values():
            if len(faces)==2:
                a,b=faces
                faces_dual[a].append(b)
                faces_dual[b].append(a)
        if _connected(faces_dual,set(range(len(links))))!=1:
            raise ValueError("pinched or disconnected 3D vertex star")
        chi=len(nodes)-len(edges)+len(links)
        boundary_link_edges=[key for key,count in edges.items() if count==1]
        if boundary_link_edges:
            graph: dict[int,list[int]]=defaultdict(list)
            for edge in boundary_link_edges:
                a,b=tuple(edge)
                graph[a].append(b)
                graph[b].append(a)
            if (chi!=1 or any(len(nbr)!=2 for nbr in graph.values())
                    or _connected(graph,set(graph))!=1):
                raise ValueError("boundary vertex link is not one PL disk")
        elif chi!=2:
            raise ValueError("interior vertex link is not a PL sphere")
    return Connected3BallAudit(
        len(vertices),len(tetrahedra),len(boundary),2,1,
        sum(len(pair)==2 for pair in incidence.values()),volume)


@dataclass(frozen=True)
class PairedRoofCutFace:
    """Both sides of one actual INTERNAL 3D triangular patch face."""
    site: int
    port: Port
    lower_tetra: int
    upper_tetra: int
    lower_face: Face
    upper_face: Face


@dataclass(frozen=True)
class RoofFaceBarriers:
    """A *dual-adjacency* detachment, NOT full edge-aware PL cutting."""
    paired_faces: tuple[PairedRoofCutFace, ...]
    dual_components_after_detachment: int
    actual_knot_piercings: int
    status: Literal["DUAL_FACE_BARRIERS_NOT_PL_MANIFOLD_SLIT"] = (
        "DUAL_FACE_BARRIERS_NOT_PL_MANIFOLD_SLIT"
    )


@dataclass(frozen=True)
class ConnectedRoofAmbientSlab:
    """Single PL 3-ball around the roofs, but NOT the source ambient S³."""
    atlas: DisjointLocalSlitAtlas
    x_slices: tuple[Fraction, ...]
    vertices: tuple[Point3, ...]
    tetrahedra: tuple[Tet, ...]
    status: Literal["CONNECTED_ROOF_SLAB_NOT_SOURCE_AMBIENT"] = (
        "CONNECTED_ROOF_SLAB_NOT_SOURCE_AMBIENT"
    )

    def __post_init__(self) -> None:
        self.audit()

    def topology(self) -> Connected3BallAudit:
        return verify_connected_tetrahedral_ball(self.vertices,self.tetrahedra)

    def face_barriers(self) -> RoofFaceBarriers:
        """Recompute all two-sided marked faces FROM stored physical tets."""
        incidences=_face_incidence(self.tetrahedra)
        expected=[]
        detached=set()
        for site_index,site in enumerate(self.atlas.sites):
            patch=site.solid.patch
            for tri in patch.triangles():
                marks=tuple(self.vertices.index(p) for p in tri)
                pair=incidences.get(frozenset(marks),())
                if len(pair)!=2:
                    raise ValueError("source candidate patch is not a 3D internal face")
                lower=[(tid,f) for tid,f in pair
                       if all(self.vertices[v].z<=patch.corners[0].z
                              for v in self.tetrahedra[tid])]
                upper=[(tid,f) for tid,f in pair
                       if all(self.vertices[v].z>=patch.corners[0].z
                              for v in self.tetrahedra[tid])]
                if (len(lower)!=1 or len(upper)!=1
                        or not _oriented_opposite(lower[0][1],upper[0][1])):
                    raise ValueError("oriented patch face does not divide top/bottom tets")
                lo,hi=lower[0],upper[0]
                expected.append(PairedRoofCutFace(
                    site_index,site.port,lo[0],hi[0],lo[1],hi[1]))
                detached.add(frozenset(marks))

        # Only the specific 2n internal faces are detached from the
        # tetrahedron DUAL graph. No edge/vertex splitting yet.
        dual: dict[int,list[int]]=defaultdict(list)
        for face,pair in incidences.items():
            if len(pair)==2 and face not in detached:
                a,b=pair[0][0],pair[1][0]
                dual[a].append(b)
                dual[b].append(a)
        comps=_connected(dual,set(range(len(self.tetrahedra))))
        return RoofFaceBarriers(tuple(expected),comps,0)

    def audit(self) -> Connected3BallAudit:
        if (not isinstance(self.atlas,DisjointLocalSlitAtlas)
                or self.status!="CONNECTED_ROOF_SLAB_NOT_SOURCE_AMBIENT"):
            raise ValueError("unsupported ambient source promotion")
        self.atlas.audit()
        h=self.atlas.thickness
        xs=[]
        for site in self.atlas.sites:
            lo,hi=site.spatial_bounds()[0]
            xs.extend((lo,hi))
        if len(xs)!=2*len(self.atlas.sites) or any(a>=b
                for a,b in zip(xs,xs[1:])):
            raise ValueError("connected source-site slab order invalid")
        if self.x_slices!=tuple(xs):
            raise ValueError("physical grid x-axis or source patch positions altered")
        expected_vertices=tuple(
            Point3(x,y,z)
            for x in self.x_slices
            for z in (Fraction(4)-h,Fraction(4),Fraction(4)+h)
            for y in (-2,0)
        )
        if self.vertices!=expected_vertices:
            raise ValueError("3D ambient grid does not match actual site coordinates")
        expected_tets=tuple(
            tuple(_tet_map(s)[i] for i in local)
            for s in range(len(xs)-1)
            for local in _product_tets()
        )
        if self.tetrahedra!=expected_tets:
            raise ValueError("connected 3D ambient tetrahedra changed")
        stats=self.topology()
        if stats.tetrahedra!=12*(2*len(self.atlas.sites)-1):
            raise ValueError("missing connected tetrahedral bridge")
        expected_volume=(xs[-1]-xs[0])*2*(2*h)
        if stats.volume!=expected_volume:
            raise ValueError("connected ambient slab has wrong exact volume")
        if (Fraction(4)-h<=3
                or max(Fraction(p.z) for arc in
                       self.atlas.candidate.patches.strands.embedded_arcs()
                       for p in arc) >= Fraction(4)-h):
            raise ValueError("ambient roof slab falsely intersects knot or Conway balls")
        barriers=self.face_barriers()
        if (len(barriers.paired_faces)!=2*len(self.atlas.sites)
                or barriers.dual_components_after_detachment!=1
                or barriers.actual_knot_piercings!=0):
            raise ValueError("partial dual roof cuts not independently supported")
        return stats

    def require_source_global_slit(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: connected ROOF 3-ball slab omits K and "
            "most of ambient S3; removing interior dual adjacencies is "
            "not an edge-split manifold cut or a four-punctured Conway sphere"
        )


def build_connected_roof_ambient_slab(
    atom: AtomicInstruction, thickness: Fraction=Fraction(1,16)
) -> ConnectedRoofAmbientSlab:
    atlas=build_disjoint_local_slit_atlas(atom,thickness)
    xs=tuple(x for site in atlas.sites for x in site.spatial_bounds()[0])
    vertices=tuple(Point3(x,y,z) for x in xs
                   for z in (Fraction(4)-thickness,
                             Fraction(4),Fraction(4)+thickness)
                   for y in (-2,0))
    tetrahedra=tuple(
        tuple(_tet_map(s)[i] for i in local)
        for s in range(len(xs)-1)
        for local in _product_tets()
    )
    return ConnectedRoofAmbientSlab(atlas,xs,vertices,tetrahedra)
