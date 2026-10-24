"""A5.2e2b3a: independent exact-rational triangle/triangle contact kernel.

Checks actual points of intersection in Q^3, including skew, shared
vertices, full shared edges, coplanar transverse contacts, positive-area
coplanar overlaps, and noncoplanar penetrating segments.

This module checks geometry, not merely mesh-edge incidence. It does
not derive source Figure 12.10 Conway framing or produce an isotopy
certificate. Its positive result concerns *only the supplied mesh*.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Literal

from .spatial_twist_blocks import Point3

Tri = tuple[Point3, Point3, Point3]
Vec = tuple[Fraction, Fraction, Fraction]
Pt2 = tuple[Fraction, Fraction]


def _key(p: Point3) -> tuple[Fraction, Fraction, Fraction]:
    return Fraction(p.x), Fraction(p.y), Fraction(p.z)


def _minus(a: Vec, b: Vec) -> Vec:
    return a[0]-b[0],a[1]-b[1],a[2]-b[2]


def _cross(a: Vec, b: Vec) -> Vec:
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _dot(a: Vec, b: Vec) -> Fraction:
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]


def _normal(tri: Tri) -> Vec:
    p,q,r=(_key(t) for t in tri)
    normal=_cross(_minus(q,p),_minus(r,p))
    if normal==(0,0,0):
        raise ValueError("zero-area triangle")
    return normal


def _validate(tri: Tri) -> None:
    if (not isinstance(tri,tuple) or len(tri)!=3
            or not all(isinstance(p,Point3) for p in tri)
            or len(set(tri))!=3):
        raise ValueError("triangle requires three distinct exact Q3 points")
    _normal(tri)


def _bbox_disjoint(a: Tri, b: Tri) -> bool:
    return any(
        max(Fraction(getattr(p,d)) for p in a) <
        min(Fraction(getattr(p,d)) for p in b)
        or max(Fraction(getattr(p,d)) for p in b) <
        min(Fraction(getattr(p,d)) for p in a)
        for d in ("x","y","z")
    )


def _distinct(points: list[Point3]) -> tuple[Point3,...]:
    seen=set()
    result=[]
    for point in points:
        if point not in seen:
            seen.add(point)
            result.append(point)
    return tuple(result)


def _slice_by_plane(tri: Tri, plane_normal: Vec, plane_base: Vec
                    ) -> tuple[Point3,...]:
    pts=tuple(_key(p) for p in tri)
    d=tuple(_dot(plane_normal,_minus(p,plane_base)) for p in pts)
    if min(d)>0 or max(d)<0:
        return ()
    found=[tri[i] for i in range(3) if d[i]==0]
    for i in range(3):
        j=(i+1)%3
        if (d[i]<0<d[j]) or (d[j]<0<d[i]):
            f=d[i]/(d[i]-d[j])
            p=Point3(*(pts[i][k]+f*(pts[j][k]-pts[i][k])
                        for k in range(3)))
            found.append(p)
    return _distinct(found)


def _cross2(a: Pt2,b: Pt2) -> Fraction:
    return a[0]*b[1]-a[1]*b[0]


def _diff2(a: Pt2,b: Pt2) -> Pt2:
    return a[0]-b[0],a[1]-b[1]


def _signed_area(poly: list[Pt2]) -> Fraction:
    return sum((a[0]*b[1]-b[0]*a[1]
                for a,b in zip(poly,poly[1:]+poly[:1])),Fraction(0))


def _coplanar_contact(a: Tri,b: Tri,normal: Vec
                      ) -> tuple[Point3,...]:
    """Clip exact convex 2D triangles and lift back to their Q3 plane."""
    drop=next(i for i in range(3) if normal[i]!=0)
    keep=tuple(i for i in range(3) if i!=drop)
    def project(p: Point3) -> Pt2:
        k=_key(p)
        return k[keep[0]],k[keep[1]]
    subj=[project(p) for p in a]
    clip=[project(p) for p in b]
    if _signed_area(clip)<0:
        clip.reverse()
    for side_start,side_end in zip(clip,clip[1:]+clip[:1]):
        if not subj:
            break
        axis=_diff2(side_end,side_start)
        def signed(p: Pt2) -> Fraction:
            return _cross2(axis,_diff2(p,side_start))
        output=[]
        for p,q in zip(subj[-1:]+subj[:-1],subj):
            vp,vq=signed(p),signed(q)
            if (vp<0 and vq>=0) or (vp>=0 and vq<0):
                f=vp/(vp-vq)
                output.append((p[0]+f*(q[0]-p[0]),
                               p[1]+f*(q[1]-p[1])))
            if vq>=0:
                output.append(q)
        subj=[]
        for p in output:
            if p not in subj:
                subj.append(p)
    base=_key(a[0])
    lifted=[]
    for xy in subj:
        p=list(base)
        p[keep[0]]=xy[0]
        p[keep[1]]=xy[1]
        p[drop]=base[drop]-sum(
            (normal[k]*(p[k]-base[k]) for k in keep),
            Fraction(0))/normal[drop]
        lifted.append(Point3(*p))
    return _distinct(lifted)


def _farthest(points: tuple[Point3,...]) -> tuple[Point3,Point3]:
    return max(
        ((p,q) for p,q in combinations(points,2)),
        key=lambda pq:sum(((_key(pq[0])[i]-_key(pq[1])[i])**2
                           for i in range(3)),Fraction(0)),
    )


@dataclass(frozen=True)
class TriangleContact:
    kind: Literal["DISJOINT","POINT","SEGMENT","AREA"]
    points: tuple[Point3,...]

    def __post_init__(self) -> None:
        if self.kind not in ("DISJOINT","POINT","SEGMENT","AREA"):
            raise ValueError("invalid intersection classification")
        if self.kind=="DISJOINT" and self.points:
            raise ValueError("disjoint triangles have no common points")
        if self.kind=="POINT" and len(self.points)!=1:
            raise ValueError("point contact requires one point")
        if self.kind=="SEGMENT" and (len(self.points)!=2
                                      or self.points[0]==self.points[1]):
            raise ValueError("segment contact requires distinct endpoints")
        if self.kind=="AREA" and len(self.points)<3:
            raise ValueError("positive-area overlap needs polygon vertices")


def _result(points: tuple[Point3,...],coplanar: bool,
            plane_normal: Vec|None = None) -> TriangleContact:
    if not points:
        return TriangleContact("DISJOINT",())
    points=tuple(sorted(points,key=_key))
    if len(points)==1:
        return TriangleContact("POINT",points)
    if coplanar:
        assert plane_normal is not None
        # Shoelace in the nondegenerate coordinate projection.
        drop=next(i for i in range(3) if plane_normal[i]!=0)
        keep=tuple(i for i in range(3) if i!=drop)
        xy=[(_key(p)[keep[0]],_key(p)[keep[1]]) for p in points]
        # Sorting loses polygon order; recompute hull-independent
        # positive area by finding any noncollinear triple.
        for p,q,r in combinations(xy,3):
            if _cross2(_diff2(q,p),_diff2(r,p))!=0:
                return TriangleContact("AREA",points)
    p,q=_farthest(points)
    return TriangleContact("SEGMENT",tuple(sorted((p,q),key=_key)))


def triangle_contact(a: Tri,b: Tri) -> TriangleContact:
    """Classify intersection of two closed triangles in exact Q3.

    Robust for coplanar overlap, skew crossings, coincident corners,
    shared edges and penetration along a segment; no float tolerance.
    """
    _validate(a)
    _validate(b)
    if _bbox_disjoint(a,b):
        return TriangleContact("DISJOINT",())
    an,bn=_normal(a),_normal(b)
    ap,bp=_key(a[0]),_key(b[0])
    da=[_dot(bn,_minus(_key(p),bp)) for p in a]
    db=[_dot(an,_minus(_key(p),ap)) for p in b]
    if (min(da)>0 or max(da)<0 or min(db)>0 or max(db)<0):
        return TriangleContact("DISJOINT",())
    line=_cross(an,bn)
    if line==(0,0,0):
        if any(d!=0 for d in da):
            return TriangleContact("DISJOINT",())
        return _result(_coplanar_contact(a,b,an),True,an)
    ac=_slice_by_plane(a,bn,bp)
    bc=_slice_by_plane(b,an,ap)
    if not ac or not bc:
        return TriangleContact("DISJOINT",())
    axis=next(i for i in range(3) if line[i]!=0)
    lower=max(min(_key(p)[axis] for p in ac),
              min(_key(p)[axis] for p in bc))
    upper=min(max(_key(p)[axis] for p in ac),
              max(_key(p)[axis] for p in bc))
    if lower>upper:
        return TriangleContact("DISJOINT",())
    def point_at(points: tuple[Point3,...],value: Fraction) -> Point3:
        if len(points)==1:
            return points[0]
        lo,hi=min(points,key=lambda p:_key(p)[axis]),max(points,key=lambda p:_key(p)[axis])
        if lo==hi:
            return lo
        factor=(value-_key(lo)[axis])/(_key(hi)[axis]-_key(lo)[axis])
        return Point3(*(_key(lo)[k]+factor*(_key(hi)[k]-_key(lo)[k])
                        for k in range(3)))
    if lower==upper:
        return TriangleContact("POINT",(point_at(ac,lower),))
    return TriangleContact("SEGMENT",tuple(sorted(
        (point_at(ac,lower),point_at(ac,upper)),key=_key)))


@dataclass(frozen=True)
class MeshIntersectionAudit:
    checked_pairs: int
    bounding_box_rejections: int
    allowed_shared_vertices: int
    allowed_shared_edges: int


def audit_embedded_triangle_mesh(triangles: tuple[Tri,...]) -> MeshIntersectionAudit:
    """Reject non-embedded contact between *any* two closed triangles.

    Legal pairs are disjoint, touch at exactly one COMMON input
    vertex, or meet precisely along an entire COMMON input edge.
    The topology checker separately verifies that shared edges have
    the correct incidence and that vertex stars are manifolds.

    In particular, identical endpoints and an apparently legal
    finite 2-complex do not let triangle interiors cross in Q3.
    """
    if not isinstance(triangles,tuple) or not triangles:
        raise ValueError("typed nonempty triangle table required")
    for tri in triangles:
        _validate(tri)
    rejected=vertex_contacts=edge_contacts=0
    checked=0
    for i,a in enumerate(triangles):
        for j in range(i+1,len(triangles)):
            b=triangles[j]
            checked+=1
            if _bbox_disjoint(a,b):
                rejected+=1
                continue
            contact=triangle_contact(a,b)
            common=set(a)&set(b)
            if contact.kind=="DISJOINT":
                if common:
                    raise ValueError("mesh shares physical vertices without geometric contact")
                continue
            if (len(common)==1 and contact.kind=="POINT"
                    and contact.points[0] in common):
                vertex_contacts+=1
                continue
            if (len(common)==2 and contact.kind=="SEGMENT"
                    and set(contact.points)==common):
                edge_contacts+=1
                continue
            raise ValueError(
                f"triangle interiors intersect illegally: faces {i},{j} "
                f"contact={contact.kind}, shared vertices={len(common)}"
            )
    return MeshIntersectionAudit(checked,rejected,vertex_contacts,edge_contacts)
