"""A5.2e2b4c2c: exact K-versus-Conway-TRIANGLE intersection audits.

Unlike relying on the marked arc endpoints supplied in A5.2e1,
this kernel derives all physical point contacts with the actual
oriented Conway cube surface triangles using exact Q arithmetic.
It independently distinguishes a transverse pierce at a source
mark from tangential in-plane contact or an unintended fifth pierce.

This validates the n ALREADY-CONSTRUCTED cubical candidate spheres.
It does NOT create the four-punctured Conway spheres by slitting a
source-identified global ambient manifold. Those are still absent.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_return_strands import Arc, independent_pl_arc_audit
from .marked_conway_geometry import MarkedConwayCube
from .atomic_plumbing_contract import Port
from .spatial_twist_blocks import Point3

Vec = tuple[Fraction,Fraction,Fraction]
Triangle = tuple[Point3,Point3,Point3]


def _coords(p: Point3) -> Vec:
    return Fraction(p.x),Fraction(p.y),Fraction(p.z)


def _minus(a: Vec,b: Vec) -> Vec:
    return a[0]-b[0],a[1]-b[1],a[2]-b[2]


def _cross(a: Vec,b: Vec) -> Vec:
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _dot(a: Vec,b: Vec) -> Fraction:
    return sum((x*y for x,y in zip(a,b)),Fraction(0))


@dataclass(frozen=True)
class SegmentTriangleHit:
    kind: Literal["DISJOINT","TRANSVERSE","COPLANAR_POINT","COPLANAR_SEGMENT"]
    point: Point3 | None = None


def segment_triangle_hit(a: Point3,b: Point3,tri: Triangle) -> SegmentTriangleHit:
    """Exact closed segment vs exact oriented triangle, no epsilon.

    A transverse plane hit may occur at the segment endpoint. It is
    a valid source boundary piercing ONLY if it is a declared K
    terminal at a face-interior Conway marking.
    """
    if (not isinstance(a,Point3) or not isinstance(b,Point3) or a==b
            or not isinstance(tri,tuple) or len(tri)!=3
            or any(not isinstance(p,Point3) for p in tri)):
        raise ValueError("nondegenerate typed spatial segment/triangle needed")
    p,q,r=map(_coords,tri)
    normal=_cross(_minus(q,p),_minus(r,p))
    if normal==(0,0,0):
        raise ValueError("degenerate Conway surface triangle")
    u,v=_coords(a),_coords(b)
    da,db=_dot(normal,_minus(u,p)),_dot(normal,_minus(v,p))
    direction=_minus(v,u)
    def within(x: Vec) -> bool:
        return all(
            _dot(normal,_cross(_minus(t,s),_minus(x,s)))>=0
            for s,t in ((p,q),(q,r),(r,p))
        )
    if da==0 and db==0:
        # Exact rational interval clipping of the coplanar segment
        # against each oriented triangle interior half-plane.
        lo,hi=Fraction(0),Fraction(1)
        for s,t in ((p,q),(q,r),(r,p)):
            edge=_minus(t,s)
            f0=_dot(normal,_cross(edge,_minus(u,s)))
            f1=_dot(normal,_cross(edge,_minus(v,s)))
            if f0==f1:
                if f0<0:
                    return SegmentTriangleHit("DISJOINT")
                continue
            threshold=-f0/(f1-f0)
            if f1>f0:
                lo=max(lo,threshold)
            else:
                hi=min(hi,threshold)
            if lo>hi:
                return SegmentTriangleHit("DISJOINT")
        if lo==hi:
            pt=Point3(*(u[i]+lo*direction[i] for i in range(3)))
            return SegmentTriangleHit("COPLANAR_POINT",pt)
        return SegmentTriangleHit("COPLANAR_SEGMENT")
    if (da>0 and db>0) or (da<0 and db<0):
        return SegmentTriangleHit("DISJOINT")
    t=da/(da-db)
    if not 0<=t<=1:
        return SegmentTriangleHit("DISJOINT")
    x=tuple(u[i]+t*direction[i] for i in range(3))
    if not within(x):
        return SegmentTriangleHit("DISJOINT")
    return SegmentTriangleHit("TRANSVERSE",Point3(*x))


@dataclass(frozen=True)
class PhysicalConwayPiercing:
    port: Port
    puncture: str
    point: Point3
    knot_arc: int
    at_start: bool


@dataclass(frozen=True)
class ExistingConwayPiercingAudit:
    """Contacts verified on EXISTING spheres, not source-slit-generated."""
    piercings: tuple[PhysicalConwayPiercing,...]
    existing_spheres_checked: int
    source_cut_generated_spheres: int = 0
    status: Literal["EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED"] = (
        "EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED"
    )

    def __post_init__(self) -> None:
        if (self.status!="EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED"
                or self.source_cut_generated_spheres!=0
                or type(self.existing_spheres_checked) is not int
                or self.existing_spheres_checked<1
                or len(self.piercings)!=4*self.existing_spheres_checked):
            raise ValueError("source-cut certificate cannot be invented")

    def require_source_conway_boundaries(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: these 4-punctured cubical sphere "
            "intersections are not boundaries DERIVED from an ambient "
            "Figure 12.5 plumbing slit")


def audit_existing_cube_piercings(
    cubes: tuple[MarkedConwayCube,...],
    knot_arcs: tuple[Arc,...],
) -> ExistingConwayPiercingAudit:
    """Compute K piercings FROM face/segment geometry, then compare marks.

    Source-labelled points are NOT used to locate physical hits; only
    after all spatial segment/triangle intersections have been
    enumerated do we require exact agreement with the named marks.
    """
    if (not isinstance(cubes,tuple) or not cubes
            or any(not isinstance(c,MarkedConwayCube) or c.port is None
                   for c in cubes)
            or not isinstance(knot_arcs,tuple)
            or len(knot_arcs)!=2*len(cubes)):
        raise ValueError("source-port spheres and 2n physical K arcs required")
    records=[]
    for cube in cubes:
        marks={p:name for name,p in cube.punctures()}
        hits: dict[Point3,set[tuple[int,int]]] = {}
        for ai,arc in enumerate(knot_arcs):
            if not isinstance(arc,tuple) or len(arc)<2:
                raise ValueError("polyline K arc expected")
            for si,(a,b) in enumerate(zip(arc,arc[1:])):
                for tri in cube.triangles():
                    hit=segment_triangle_hit(a,b,tri)
                    if hit.kind.startswith("COPLANAR"):
                        raise ValueError("K lies in / touches Conway boundary tangentially")
                    if hit.kind=="TRANSVERSE":
                        assert hit.point is not None
                        hits.setdefault(hit.point,set()).add((ai,si))
        if set(hits)!=set(marks):
            raise ValueError("computed physical K/sphere intersections != four marks")
        if len(hits)!=4:
            raise ValueError("a source Conway sphere requires four unique piercings")
        for p,name in marks.items():
            matches=hits[p]
            if len(matches)!=1:
                raise ValueError("Conway piercing meets multiple K segments")
            ai,si=next(iter(matches))
            arc=knot_arcs[ai]
            at_start=p==arc[0] and si==0
            at_end=p==arc[-1] and si==len(arc)-2
            if at_start==at_end:
                raise ValueError("piercing must occur at one proper K endpoint")
            # Independently certify nonzero outward face-normal component.
            near=arc[1] if at_start else arc[-2]
            if p.x==cube.west:
                if not near.x<cube.west:
                    raise ValueError("west-face K piercing is tangent or enters solid")
            elif p.x==cube.east:
                if not near.x>cube.east:
                    raise ValueError("east-face K piercing is tangent or enters solid")
            else:
                raise ValueError("mark is not in an intended transverse face")
            if not (-3<p.y<3 and -3<p.z<3):
                raise ValueError("Conway piercing occurred on a mesh face edge")
            records.append(PhysicalConwayPiercing(
                cube.port,name,p,ai,at_start))
    # Independent ambient-arc check ensures no hidden through-cube
    # intervals between contact sites, overlapping arcs or duplicate
    # source terminals. It does not replace the face-by-face contacts.
    independent_pl_arc_audit(cubes,knot_arcs)
    return ExistingConwayPiercingAudit(tuple(records),len(cubes))


def audit_existing_band_conway_piercings(band) -> ExistingConwayPiercingAudit:
    from .exterior_conway_collars import ExteriorCollaredBandCandidate
    if not isinstance(band,ExteriorCollaredBandCandidate):
        raise ValueError("typed independently checked chosen-band candidate required")
    band.audit()
    return audit_existing_cube_piercings(
        band.patches.strands.cubes,band.patches.strands.embedded_arcs())
