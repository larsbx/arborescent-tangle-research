"""A5.2e2b3b: exact off-Conway-boundary PL collar candidates in Q^3.

Unlike the A5.2e2b2 roof surfaces lying on whole Conway sphere faces,
these collars meet closed Conway balls ONLY on the two designated
three-point end-cap chords. Their open surface lies in the complement
of ALL the removed Conway balls.

The new independent exact triangle/triangle intersection kernel from
A5.2e2b3a audits the entire physical joined surface. This proves a
valid embedded *chosen-coordinate* band candidate, NOT correspondence
to Bonahon--Siebenmann Fig.12.10's source plumbing maps, framing,
boundary marking, or relative-boundary isotopy. Source certification
remains UNVERIFIED_GEOMETRY.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .conway_roof_collars import ExactMeshTopology, verify_spatial_surface
from .marked_conway_geometry import MarkedConwayCube
from .pl_ribbon_patches import (
    ZeroRingRibbonPatchSystem, _clip_cube, _in_segment,
    build_zero_ring_ribbon_patches,
)
from .spatial_twist_blocks import Point3
from .triangle_intersections import MeshIntersectionAudit, audit_embedded_triangle_mesh

Triangle = tuple[Point3,Point3,Point3]


@dataclass(frozen=True)
class ExteriorConwayCollar:
    """An actual PL ribbon collar with interior outside the solid cube.

    Longitudinal x-z profile from west cap to east cap:
       (west,0) -> (west-1,4) -> (east+1,4) -> (east,0).
    Each cross-section has y=-2,0,2. The first/last chords are
    exactly the existing ribbon patch cap meshes on marked spheres.
    All other collar points are outside the CLOSED solid cube.
    """
    cube: MarkedConwayCube
    vertices: tuple[Point3,...]
    faces: tuple[tuple[int,int,int],...]
    status: Literal["OFF_BOUNDARY_CANDIDATE"] = "OFF_BOUNDARY_CANDIDATE"

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
        if (not isinstance(self.cube,MarkedConwayCube)
                or self.status!="OFF_BOUNDARY_CANDIDATE"):
            raise ValueError("unsupported Conway collar type/status")
        w,e=self.cube.west,self.cube.east
        expected=tuple(Point3(x,y,z)
                       for x,z in ((w,0),(w-1,4),(e+1,4),(e,0))
                       for y in (-2,0,2))
        if self.vertices!=expected:
            raise ValueError("off-boundary collar coordinates or framing mutated")
        faces=[]
        for station in range(3):
            for half in range(2):
                a=3*station+half
                b=3*(station+1)+half
                faces.extend(((a,b,b+1),(a,b+1,a+1)))
        if self.faces!=tuple(faces):
            raise ValueError("unsupported off-boundary collar triangulation")
        # Exact whole-triangle-vs-solid-ball clipping, not edge sampling.
        # Each attached cap is the only legal cube contact.
        for tri in self.triangles():
            clipped=_clip_cube(tri,self.cube)
            if not all(
                any(_in_segment(point,a,b)
                    for cap in (self.west_cap(),self.east_cap())
                    for a,b in zip(cap,cap[1:]))
                for point in clipped
            ):
                raise ValueError("interior of collar crosses a Conway ball")


def exterior_collar(cube: MarkedConwayCube) -> ExteriorConwayCollar:
    if not isinstance(cube,MarkedConwayCube):
        raise ValueError("Conway sphere required")
    w,e=cube.west,cube.east
    vertices=tuple(Point3(x,y,z)
                   for x,z in ((w,0),(w-1,4),(e+1,4),(e,0))
                   for y in (-2,0,2))
    faces=[]
    for i in range(3):
        for j in range(2):
            a=3*i+j
            b=3*(i+1)+j
            faces.extend(((a,b,b+1),(a,b+1,a+1)))
    return ExteriorConwayCollar(cube,vertices,tuple(faces))


def _cap(ribbon,side: str) -> tuple[Point3,...]:
    path=ribbon.boundary_paths()[2 if side=="start" else 3]
    return tuple(sorted(path,key=lambda p:p.y))


@dataclass(frozen=True)
class ExteriorCollaredBandCandidate:
    """Full physical ribbon *candidate*, with source orientation unknown."""
    patches: ZeroRingRibbonPatchSystem
    collars: tuple[ExteriorConwayCollar,...]
    status: Literal["OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED"] = "OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED"

    def __post_init__(self) -> None:
        self.audit()

    def triangles(self) -> tuple[Triangle,...]:
        return (tuple(tri for ribbon in self.patches.ribbons
                      for tri in ribbon.triangles())
                +tuple(tri for collar in self.collars
                       for tri in collar.triangles()))

    def topology(self) -> ExactMeshTopology:
        return verify_spatial_surface(self.triangles())

    def geometric_intersections(self) -> MeshIntersectionAudit:
        return audit_embedded_triangle_mesh(self.triangles())

    def audit(self) -> None:
        if (not isinstance(self.patches,ZeroRingRibbonPatchSystem)
                or not isinstance(self.collars,tuple)
                or self.status!="OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED"):
            raise ValueError("unverified off-boundary candidate required")
        self.patches.audit()
        cubes=self.patches.strands.cubes
        ribbons=self.patches.ribbons
        if len(self.collars)!=len(cubes):
            raise ValueError("one physical collar required for every Conway cube")
        for i,(cube,collar) in enumerate(zip(cubes,self.collars)):
            if not isinstance(collar,ExteriorConwayCollar) or collar!=exterior_collar(cube):
                raise ValueError("incorrect source-port collar / physical mesh")
            incoming=ribbons[0] if i==0 else ribbons[i]
            outgoing=ribbons[0] if i==len(cubes)-1 else ribbons[i+1]
            if (_cap(incoming,"end")!=collar.west_cap()
                    or _cap(outgoing,"start")!=collar.east_cap()):
                raise ValueError("physical cap coordinates do not weld")
            # Not only the collar's own ball: test against EVERY removed
            # Conway ball using exact triangle-solid-cube clipping.
            for tri in collar.triangles():
                for other in cubes:
                    clipped=_clip_cube(tri,other)
                    if not clipped:
                        continue
                    if other!=cube or not all(
                        any(_in_segment(p,a,b)
                            for cap in (collar.west_cap(),collar.east_cap())
                            for a,b in zip(cap,cap[1:]))
                        for p in clipped
                    ):
                        raise ValueError("collar illegally intersects another Conway ball")
        topology=self.topology()
        if topology.components!=1 or topology.euler!=0:
            raise ValueError("off-boundary candidate is not one chi-zero surface")
        if (topology.orientable and topology.boundary_components!=2
                or not topology.orientable and topology.boundary_components!=1):
            raise ValueError("bad orientable/nonorientable boundary topology")
        # CRITICAL independent check: a manifold can still self-intersect
        # away from its mesh vertices/edges.
        self.geometric_intersections()

    def canonical_coordinates(self) -> tuple:
        return (self.patches.strands.canonical_coordinates(),
                tuple(c.cube.port for c in self.collars),
                tuple(tuple(p.serial() for p in c.vertices) for c in self.collars),
                self.topology())

    def require_source_atomic_realization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: physically embedded off-boundary PL band "
            "candidate has no source Fig.12.10 framing/Conway parametrization "
            "or relative-boundary isotopy witness"
        )


def build_exterior_collared_band(atom: AtomicInstruction) -> ExteriorCollaredBandCandidate:
    """Build exact non-source-certified global PL band and audit all faces."""
    patches=build_zero_ring_ribbon_patches(atom)
    collars=tuple(exterior_collar(cube) for cube in patches.strands.cubes)
    return ExteriorCollaredBandCandidate(patches,collars)
