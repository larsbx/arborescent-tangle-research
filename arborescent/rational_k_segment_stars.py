"""A5.2e2b4c2d2b: exact rational 3D K-edge stars in ANY direction.

Positive theorem: a proper interior subsegment of one previously
independently verified source-candidate K edge can be made a literal
interior tetrahedral 1-simplex, including vertical and signed-braid
skew edges. Exact rational transversals, shrinking clearance,
independent segment/tube intersections and 3D manifold links.

This is a LOCAL edge star; segment endpoints/corners and the entire
K/band/Conway ambient conforming triangulation remain unavailable.
No source-equivalent isotopy, 4-punctured Conway slit or R+/- proof.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .connected_roof_ambient_slab import verify_connected_tetrahedral_ball
from .exterior_conway_collars import (
    ExteriorCollaredBandCandidate, build_exterior_collared_band,
)
from .k_conforming_tetrahedral_tubes import (
    _fan_tetrahedra, _core_edge_link,
)
from .spatial_twist_blocks import Point3

Vector = tuple[Fraction,Fraction,Fraction]
Axis = tuple[Vector,Vector,Vector]
UnitAxis = ((Fraction(1),Fraction(0),Fraction(0)),
            (Fraction(0),Fraction(1),Fraction(0)),
            (Fraction(0),Fraction(0),Fraction(1)))


def _xyz(p: Point3) -> Vector:
    return Fraction(p.x),Fraction(p.y),Fraction(p.z)


def _minus(a: Vector,b: Vector) -> Vector:
    return tuple(x-y for x,y in zip(a,b))  # type: ignore[return-value]


def _plus(a: Vector,b: Vector) -> Vector:
    return tuple(x+y for x,y in zip(a,b))  # type: ignore[return-value]


def _scale(a: Vector,s: Fraction) -> Vector:
    return tuple(s*x for x in a)  # type: ignore[return-value]


def _cross(a: Vector,b: Vector) -> Vector:
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _dot(a: Vector,b: Vector) -> Fraction:
    return sum((x*y for x,y in zip(a,b)),Fraction(0))


def _frame(a: Point3,b: Point3) -> Axis:
    """Rational nondegenerate right-handed (longitudinal,u,v) frame.

    Transverse vectors have max absolute coordinate exactly 1;
    no sqrt, normalization approximation, or angle is involved.
    """
    d=_minus(_xyz(b),_xyz(a))
    if d==(0,0,0):
        raise ValueError("zero source knot segment cannot carry a tube")
    u=next((_cross(d,e) for e in UnitAxis
            if _cross(d,e)!=(0,0,0)),None)
    assert u is not None
    u=_scale(u,1/max(abs(x) for x in u))
    v=_cross(d,u)
    v=_scale(v,1/max(abs(x) for x in v))
    if (_dot(d,u)!=0 or _dot(d,v)!=0 or _dot(u,v)!=0
            or _dot(d,_cross(u,v))<=0):
        raise ValueError("rational transverse frame degenerate")
    return d,u,v


def _axes_vertices(a: Point3,b: Point3,r: Fraction) -> tuple[Point3,...]:
    d,u,v=_frame(a,b)
    base=_xyz(a)
    vertices=[]
    for t in (Fraction(1,4),Fraction(3,4)):
        x=_plus(base,_scale(d,t))
        vertices.append(Point3(*x))
        for su,sv in ((-1,-1),(1,-1),(1,1),(-1,1)):
            vertices.append(Point3(*_plus(x,_plus(_scale(u,su*r),
                                                _scale(v,sv*r)))))
    return tuple(vertices)


def _tube_corners(vertices: tuple[Point3,...]) -> tuple[Vector,...]:
    if len(vertices)!=10:
        raise ValueError("K 3D tube needs ten rational vertices")
    return tuple(_xyz(vertices[i]) for i in (1,2,3,4,6,7,8,9))


def _separator(points_a: tuple[Vector,...],points_b: tuple[Vector,...],
               axes: tuple[Vector,...]) -> bool:
    """Exact separating axis certificate for two convex 3D polytopes.

    Strict inequalities prove the CLOSED solids are disjoint, not
    merely that their open interiors are nonoverlapping.
    """
    for axis in axes:
        if axis==(0,0,0):
            continue
        x=tuple(_dot(p,axis) for p in points_a)
        y=tuple(_dot(p,axis) for p in points_b)
        if max(x)<min(y) or max(y)<min(x):
            return True
    return False


def _box_disjoint_from_cube(vertices: tuple[Point3,...],cube) -> bool:
    """Exact parallelepiped vs axis-aligned closed Conway cube SAT."""
    corners=_tube_corners(vertices)
    a,b,c=(_minus(corners[j],corners[0]) for j in (4,1,3))
    directions=(a,b,c)
    normals=tuple(_cross(directions[i],directions[(i+1)%3])
                  for i in range(3))
    axes=normals+UnitAxis+tuple(_cross(d,e) for d,e in
                                         product(directions,UnitAxis))
    return _separator(corners,tuple(_xyz(p) for p in cube.vertices()),axes)


def _clip_k_segment(
    a: Point3,b: Point3,
    source_a: Point3,source_b: Point3,
    radius: Fraction,
) -> tuple[Fraction,Fraction] | None:
    """Exact closed segment vs oblique K tube via 3x3 Cramer coords."""
    d,u,v=_frame(source_a,source_b)
    det=_dot(d,_cross(u,v))
    start=_xyz(source_a)
    def transform(p: Point3) -> Vector:
        w=_minus(_xyz(p),start)
        return (_dot(w,_cross(u,v))/det,
                _dot(w,_cross(v,d))/det,
                _dot(w,_cross(d,u))/det)
    p,q=transform(a),transform(b)
    lo,hi=Fraction(0),Fraction(1)
    for j,(minimum,maximum) in enumerate((
        (Fraction(1,4),Fraction(3,4)),(-radius,radius),
        (-radius,radius)
    )):
        velocity=q[j]-p[j]
        if velocity==0:
            if p[j]<minimum or p[j]>maximum:
                return None
            continue
        t0,t1=(minimum-p[j])/velocity,(maximum-p[j])/velocity
        lo=max(lo,min(t0,t1))
        hi=min(hi,max(t0,t1))
        if lo>hi:
            return None
    return lo,hi


@dataclass(frozen=True)
class RationalKSegmentStar:
    """One PHYSICALLY verified source K subedge, not source global ambient."""
    candidate: ExteriorCollaredBandCandidate
    arc: int
    segment: int
    radius: Fraction
    vertices: tuple[Point3,...]
    tetrahedra: tuple[tuple[int,int,int,int],...]
    status: Literal["ARBITRARY_EDGE_STAR_NOT_GLOBAL_K_AMBIENT"] = (
        "ARBITRARY_EDGE_STAR_NOT_GLOBAL_K_AMBIENT"
    )

    def __post_init__(self) -> None:
        self.audit()

    def audit(self):
        if (not isinstance(self.candidate,ExteriorCollaredBandCandidate)
                or type(self.arc) is not int or type(self.segment) is not int
                or type(self.radius) is not Fraction or not 0<self.radius<=Fraction(1,4)
                or self.status!="ARBITRARY_EDGE_STAR_NOT_GLOBAL_K_AMBIENT"):
            raise ValueError("typed exact local K star and radius required")
        arcs=self.candidate.patches.strands.embedded_arcs()
        if not 0<=self.arc<len(arcs) or not 0<=self.segment<len(arcs[self.arc])-1:
            raise ValueError("chosen K arc segment does not exist")
        a,b=arcs[self.arc][self.segment:self.segment+2]
        expected=_axes_vertices(a,b,self.radius)
        if self.vertices!=expected or self.tetrahedra!=_fan_tetrahedra(expected):
            raise ValueError("noncanonical exact K 3D tube or triangulation")
        _core_edge_link(self.tetrahedra)
        checked=verify_connected_tetrahedral_ball(self.vertices,self.tetrahedra)
        direction,u,v=_frame(a,b)
        expected_volume=(Fraction(1,2)*
                         _dot(direction,_cross(u,v))*
                         (2*self.radius)**2)
        if (checked.vertices!=10 or checked.tetrahedra!=12 or
                checked.boundary_components!=1 or checked.boundary_euler!=2 or
                checked.volume!=expected_volume):
            raise ValueError("K source segment 3D PL ball/volume invalid")
        if not all(_box_disjoint_from_cube(self.vertices,cube)
                   for cube in self.candidate.patches.strands.cubes):
            raise ValueError("K tube intersects removed solid Conway ball")
        allowed=0
        for ai,arc in enumerate(arcs):
            for si,(p,q) in enumerate(zip(arc,arc[1:])):
                interval=_clip_k_segment(p,q,a,b,self.radius)
                if interval is None:
                    continue
                if (ai,si)!=(self.arc,self.segment):
                    raise ValueError("an extraneous original K edge intersects tube")
                if interval!=(Fraction(1,4),Fraction(3,4)):
                    raise ValueError("tube core is not interior proper subsegment")
                allowed+=1
        if allowed!=1:
            raise ValueError("K core original edge is not unique")
        return checked

    def require_source_ambient(self):
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: a proper-interior K edge tube "
            "has no K vertex/corner stars, whole-K triangulation "
            "or source 4-punctured slit/isotopy witness")


def build_rational_k_segment_star(
    atom: AtomicInstruction, arc: int, segment: int,
    max_radius: Fraction=Fraction(1,4),
) -> RationalKSegmentStar:
    """Exact deterministic search for a physically separated edge star.

    Candidate radii halve geometrically. A failed bounded search is
    UNVERIFIED (ValueError), not proof of knot obstruction.
    """
    if (not isinstance(atom,AtomicInstruction)
            or type(max_radius) is not Fraction
            or not 0<max_radius<=Fraction(1,4)
            or type(arc) is not int or type(segment) is not int):
        raise ValueError("typed source K segment and rational radius required")
    candidate=build_exterior_collared_band(atom)
    paths=candidate.patches.strands.embedded_arcs()
    if not 0<=arc<len(paths) or not 0<=segment<len(paths[arc])-1:
        raise ValueError("selected existing K edge out of range")
    a,b=paths[arc][segment:segment+2]
    for k in range(24):
        radius=max_radius/Fraction(2**k)
        vertices=_axes_vertices(a,b,radius)
        try:
            return RationalKSegmentStar(
                candidate,arc,segment,radius,vertices,
                _fan_tetrahedra(vertices))
        except ValueError:
            continue
    raise ValueError(
        "UNVERIFIED_GEOMETRY: 24 exact radius refinements did not "
        "produce a separated K-edge star; not a nonembedding proof")
