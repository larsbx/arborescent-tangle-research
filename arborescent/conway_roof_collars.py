"""A5.2e2b2: physically coincident PL cap collars on marked Conway cubes.

This converts the n separate A5.2e2a embedded ribbon disks into a
CONNECTED exact-rational 2D triangulated surface by adding U-shaped
collars *on the boundary of each removed Conway cube*. Each collar
connects the west cap of the incoming disk to the east cap of the
outgoing disk, using real Q^3 vertices, not abstract quotient maps.

The result is a verified topological/coordinate candidate **only in the
chosen cubical geometry**. A collar lying on a Conway sphere is not
the source-specific transverse plumbing patch or a source-canonical
boundary parametrization. No source Figure 12.10 isotopy, R+/- move,
or r>0 ring is certified. A5.1 remains UNVERIFIED_GEOMETRY.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .marked_conway_geometry import MarkedConwayCube
from .pl_ribbon_patches import (
    TriangulatedRibbon, ZeroRingRibbonPatchSystem,
    build_zero_ring_ribbon_patches,
)
from .spatial_twist_blocks import Point3

Triangle = tuple[Point3, Point3, Point3]


def _difference(a: Point3, b: Point3) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(Fraction(getattr(b, d))-Fraction(getattr(a, d))
                 for d in ("x","y","z"))


def _cross(a: tuple[Fraction,...], b: tuple[Fraction,...]
           ) -> tuple[Fraction, Fraction, Fraction]:
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _connected(graph: dict, nodes: set) -> int:
    count = 0
    pending = set(nodes)
    while pending:
        count += 1
        queue = [pending.pop()]
        while queue:
            for node in graph.get(queue.pop(), ()):
                if node in pending:
                    pending.remove(node)
                    queue.append(node)
    return count


@dataclass(frozen=True)
class ExactMeshTopology:
    components: int
    euler: int
    boundary_components: int
    orientable: bool
    vertices: int
    edges: int
    faces: int


def verify_spatial_surface(triangles: tuple[Triangle, ...]) -> ExactMeshTopology:
    """Independent point-welded finite 2-manifold checker in actual Q^3.

    The input is just triangle coordinates; no source atom, gluing
    metadata or expected topological invariants are trusted. This
    checks *combinatorial* PL 2-manifoldness; isolated geometric
    triangle interior intersections require separate clearance proofs.
    """
    if not isinstance(triangles, tuple) or not triangles:
        raise ValueError("explicit nonempty exact triangle tuple required")
    vertices = set()
    edges: dict[frozenset[Point3], list[tuple[int,Point3,Point3]]] = defaultdict(list)
    corners: dict[Point3, set[int]] = defaultdict(set)
    prior_faces = set()
    for fid, tri in enumerate(triangles):
        if (not isinstance(tri, tuple) or len(tri)!=3
                or not all(isinstance(p, Point3) for p in tri)
                or len(set(tri))!=3
                or _cross(_difference(tri[0],tri[1]),
                          _difference(tri[0],tri[2]))==(0,0,0)):
            raise ValueError("degenerate or non-PL spatial triangle")
        face = frozenset(tri)
        if face in prior_faces:
            raise ValueError("duplicate coincident mesh face")
        prior_faces.add(face)
        vertices.update(tri)
        for p in tri:
            corners[p].add(fid)
        for a,b in zip(tri,tri[1:]+tri[:1]):
            edges[frozenset((a,b))].append((fid,a,b))
    if any(len(inc) not in (1,2) for inc in edges.values()):
        raise ValueError("nonmanifold edge shared by more than two faces")

    face_graph: dict[int,list[int]] = defaultdict(list)
    signed: dict[int,list[tuple[int,int]]] = defaultdict(list)
    links: dict[Point3,dict[int,list[int]]] = defaultdict(lambda:defaultdict(list))
    boundary_graph: dict[Point3,list[Point3]] = defaultdict(list)
    for edge,inc in edges.items():
        a,b=tuple(edge)
        if len(inc)==1:
            boundary_graph[a].append(b)
            boundary_graph[b].append(a)
            continue
        (f,u,v),(g,w,x)=inc
        factor = -1 if (u,v)==(w,x) else 1
        if (u,v)!=(w,x) and (u,v)!=(x,w):
            raise ValueError("edge not consistently welded")
        face_graph[f].append(g)
        face_graph[g].append(f)
        signed[f].append((g,factor))
        signed[g].append((f,factor))
        for point in edge:
            links[point][f].append(g)
            links[point][g].append(f)

    # Each physical vertex link must be a single interval or circle.
    for p, incident_faces in corners.items():
        adjacency = links[p]
        if any(len(adjacency[f])>2 for f in incident_faces):
            raise ValueError("branching vertex link")
        if _connected(adjacency, incident_faces)!=1:
            raise ValueError("pinched vertex / disconnected star")
        degrees=[len(adjacency[f]) for f in incident_faces]
        if not (all(deg==2 for deg in degrees)
                or degrees==[0]
                or (degrees.count(1)==2 and degrees.count(0)==0
                    and all(deg in (1,2) for deg in degrees))):
            raise ValueError("vertex link is neither an interval nor a circle")
    if any(len(neighbors)!=2 for neighbors in boundary_graph.values()):
        raise ValueError("boundary is not a union of circles")

    orientations={}
    orientable=True
    for base in range(len(triangles)):
        if base in orientations:
            continue
        orientations[base]=1
        queue=[base]
        while queue:
            now=queue.pop()
            for other,ratio in signed[now]:
                want=orientations[now]*ratio
                if other in orientations:
                    if orientations[other]!=want:
                        orientable=False
                else:
                    orientations[other]=want
                    queue.append(other)
    return ExactMeshTopology(
        _connected(face_graph,set(range(len(triangles)))),
        len(vertices)-len(edges)+len(triangles),
        _connected(boundary_graph,set(boundary_graph)),
        orientable,len(vertices),len(edges),len(triangles)
    )


@dataclass(frozen=True)
class ConwayRoofCollar:
    """A real triangulated U-shaped strip contained in a Conway cube boundary.

    Start: cap at x=west, z=0. Follow up the west face, over
    z=3 roof, and down the east face to x=east, z=0.
    Width sections y=-2,0,2 use exact spatial coordinates.
    """
    cube: MarkedConwayCube
    vertices: tuple[Point3,...]
    faces: tuple[tuple[int,int,int],...]
    status: Literal["GEOMETRIC_BOUNDARY_COLLAR"] = "GEOMETRIC_BOUNDARY_COLLAR"

    def __post_init__(self) -> None:
        self.audit()

    def triangles(self) -> tuple[Triangle,...]:
        return tuple(tuple(self.vertices[i] for i in face)
                     for face in self.faces)  # type: ignore[return-value]

    def west_cap(self) -> tuple[Point3,Point3,Point3]:
        return self.vertices[:3]  # type: ignore[return-value]

    def east_cap(self) -> tuple[Point3,Point3,Point3]:
        return self.vertices[-3:]  # type: ignore[return-value]

    def audit(self) -> None:
        if not isinstance(self.cube,MarkedConwayCube) or self.status!="GEOMETRIC_BOUNDARY_COLLAR":
            raise ValueError("invalid source-boundary collar")
        expected=tuple(Point3(x,y,z) for x,z in (
            (self.cube.west,0),(self.cube.west,3),
            (self.cube.east,3),(self.cube.east,0))
            for y in (-2,0,2))
        if self.vertices!=expected:
            raise ValueError("collar not on exact predeclared Conway cube faces")
        face_indices=[]
        for k in range(3):
            for j in range(2):
                a=3*k+j
                b=3*(k+1)+j
                face_indices.extend(((a,b,b+1),(a,b+1,a+1)))
        if self.faces!=tuple(face_indices):
            raise ValueError("collar faces fail exact geometric 2D embedding contract")
        for tri in self.triangles():
            if _cross(_difference(tri[0],tri[1]),
                      _difference(tri[0],tri[2]))==(0,0,0):
                raise ValueError("degenerate collar triangle")
            # Convex triangle lies on one cube boundary face when its
            # three vertices all have x=west / x=east / z=3.
            if not any(
                all(getattr(p,axis)==fixed for p in tri)
                for axis,fixed in (("x",self.cube.west),
                                   ("x",self.cube.east),("z",3))
            ):
                raise ValueError("collar does not stay on cube boundary")


def roof_collar(cube: MarkedConwayCube) -> ConwayRoofCollar:
    if not isinstance(cube,MarkedConwayCube):
        raise ValueError("Conway cube required")
    vertices=tuple(Point3(x,y,z) for x,z in (
        (cube.west,0),(cube.west,3),
        (cube.east,3),(cube.east,0)) for y in (-2,0,2))
    faces=[]
    for k in range(3):
        for j in range(2):
            a=3*k+j
            b=3*(k+1)+j
            faces.extend(((a,b,b+1),(a,b+1,a+1)))
    return ConwayRoofCollar(cube,vertices,tuple(faces))


def _actual_cap(ribbon: TriangulatedRibbon, end: str) -> tuple[Point3,...]:
    cap=ribbon.boundary_paths()[2 if end=="start" else 3]
    if len(cap)!=3:
        raise ValueError("physical ribbon cap must have exactly three vertices")
    return tuple(sorted(cap,key=lambda p:p.y))


@dataclass(frozen=True)
class SpatialCollaredBandCandidate:
    patches: ZeroRingRibbonPatchSystem
    collars: tuple[ConwayRoofCollar,...]
    status: Literal["COLLARED_BAND_SOURCE_UNVERIFIED"] = "COLLARED_BAND_SOURCE_UNVERIFIED"

    def __post_init__(self) -> None:
        self.audit()

    def triangles(self) -> tuple[Triangle,...]:
        return (tuple(face for r in self.patches.ribbons for face in r.triangles())
                +tuple(face for c in self.collars for face in c.triangles()))

    def topology(self) -> ExactMeshTopology:
        return verify_spatial_surface(self.triangles())

    def audit(self) -> None:
        if (not isinstance(self.patches,ZeroRingRibbonPatchSystem)
                or not isinstance(self.collars,tuple)
                or self.status!="COLLARED_BAND_SOURCE_UNVERIFIED"):
            raise ValueError("invalid collared band candidate")
        self.patches.audit()
        cubes=self.patches.strands.cubes
        ribbons=self.patches.ribbons
        if len(self.collars)!=len(cubes):
            raise ValueError("must have one physical collar per Conway boundary")
        for i,(cube,collar) in enumerate(zip(cubes,self.collars)):
            if not isinstance(collar,ConwayRoofCollar) or collar!=roof_collar(cube):
                raise ValueError("collar not canonical on its actual cube")
            incoming=ribbons[0] if i==0 else ribbons[i]
            outgoing=ribbons[0] if i==len(cubes)-1 else ribbons[i+1]
            if (_actual_cap(incoming,"end")!=collar.west_cap()
                    or _actual_cap(outgoing,"start")!=collar.east_cap()):
                raise ValueError("physical ribbon caps do not match collar coordinates")
        # Exact constructive clearance: all source patches are independently
        # cube-interior-free and meet each closed cube only at cap chords
        # (A5.2e2a audit); collar triangles are on that cube boundary.
        # The cubes are in disjoint closed x-slabs, so distinct collars
        # are disjoint; other patches are already separately audited.
        topo=self.topology()
        if topo.components!=1 or topo.euler!=0:
            raise ValueError("spatial collar arrangement not a connected chi=0 band")
        if (topo.orientable and topo.boundary_components!=2
                or not topo.orientable and topo.boundary_components!=1):
            raise ValueError("spatial collar topology is not annulus/Mobius")
        # No source framing or PL isotopy certificate is inferred.

    def canonical_coordinates(self) -> tuple:
        return (self.patches.strands.canonical_coordinates(),
                tuple(c.cube.center_x for c in self.collars),
                tuple(tuple(p.serial() for p in c.vertices) for c in self.collars),
                self.topology())

    def require_source_atomic_realization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: collared Q^3 band in chosen cubical frame "
            "has no source Fig.12.10 boundary/framing isotopy certificate"
        )


def build_spatial_collared_band(atom: AtomicInstruction) -> SpatialCollaredBandCandidate:
    """Create an exact physically connected, source-unverified band candidate."""
    patches=build_zero_ring_ribbon_patches(atom)
    return SpatialCollaredBandCandidate(
        patches,tuple(roof_collar(cube) for cube in patches.strands.cubes)
    )
