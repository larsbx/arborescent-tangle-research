"""Independent exact-Q3 geometric intersection kernel for tetrahedral PL meshes.

This checks *geometric* contact, NOT just combinatorial vertex links.
For two closed tetrahedra, intersect their eight oriented rational
halfspaces: every extreme point of the bounded intersection is an
intersection of three linearly independent supporting planes. A
tetrahedron pair may meet only along the convex hull of its ACTUAL
shared mesh vertices (empty/vertex/edge/face). Unexpected point
contacts, face penetration and full-dimensional overlap fail closed.

Finite exact rational arithmetic only; no float or tolerances.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations

from .spatial_twist_blocks import Point3

Vec=tuple[Fraction,Fraction,Fraction]
Tetra=tuple[Point3,Point3,Point3,Point3]
Plane=tuple[Vec,Fraction]


def _xyz(p:Point3)->Vec:
    return Fraction(p.x),Fraction(p.y),Fraction(p.z)


def _sub(a:Vec,b:Vec)->Vec:
    return a[0]-b[0],a[1]-b[1],a[2]-b[2]


def _cross(a:Vec,b:Vec)->Vec:
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _dot(a:Vec,b:Vec)->Fraction:
    return sum((a[i]*b[i] for i in range(3)),Fraction(0))


def _planes(tetra:Tetra)->tuple[Plane,...]:
    if (not isinstance(tetra,tuple) or len(tetra)!=4 or
            any(not isinstance(v,Point3) for v in tetra) or
            len(set(tetra))!=4):
        raise ValueError("four distinct exact 3D vertices required")
    points=tuple(_xyz(p) for p in tetra)
    faces=[]
    for omit in range(4):
        ids=[i for i in range(4) if i!=omit]
        p,q,r=(points[i] for i in ids)
        normal=_cross(_sub(q,p),_sub(r,p))
        if normal==(0,0,0):
            raise ValueError("zero area tetrahedral face")
        if _dot(normal,_sub(points[omit],p))>0:
            normal=tuple(-x for x in normal)  # type: ignore[assignment]
        if _dot(normal,_sub(points[omit],p))>=0:
            raise ValueError("zero-volume tetrahedron")
        faces.append((normal,_dot(normal,p)))
    return tuple(faces)


def _aabb_separate(a:Tetra,b:Tetra)->bool:
    return any(
        max(getattr(v,axis) for v in a)<min(getattr(v,axis) for v in b)
        or max(getattr(v,axis) for v in b)<min(getattr(v,axis) for v in a)
        for axis in ("x","y","z")
    )


def _triple_intersection(a:Plane,b:Plane,c:Plane)->Vec|None:
    na,da=a;nb,db=b;nc,dc=c
    x=_cross(nb,nc)
    y=_cross(nc,na)
    z=_cross(na,nb)
    det=_dot(na,x)
    if det==0:
        return None
    return tuple((da*x[k]+db*y[k]+dc*z[k])/det
                 for k in range(3))  # type: ignore[return-value]


def tetrahedron_intersection_vertices(a:Tetra,b:Tetra)->tuple[Point3,...]:
    """Return ALL extreme points of exact closed tetrahedron intersection.

    Complete even for lower-dimensional edge/point contacts.
    """
    halfspaces=_planes(a)+_planes(b)
    if _aabb_separate(a,b):
        return ()
    corners=set()
    for triple in combinations(halfspaces,3):
        p=_triple_intersection(*triple)
        if p is not None and all(_dot(n,p)<=d for n,d in halfspaces):
            corners.add(Point3(*p))
    return tuple(sorted(corners,key=_xyz))


def _on_shared_face(p:Point3,common:tuple[Point3,...])->bool:
    a,b,c=(_xyz(v) for v in common)
    normal=_cross(_sub(b,a),_sub(c,a))
    q=_xyz(p)
    if _dot(normal,_sub(q,a))!=0:
        return False
    return all(_dot(normal,_cross(_sub(v,u),_sub(q,u)))>=0
               for u,v in ((a,b),(b,c),(c,a)))


def _within_shared_simplex(point:Point3,common:tuple[Point3,...])->bool:
    if not common:
        return False
    if len(common)==1:
        return point==common[0]
    if len(common)==2:
        a,b=(_xyz(v) for v in common)
        d=_sub(b,a);e=_sub(_xyz(point),a)
        if _cross(d,e)!=(0,0,0):
            return False
        return 0<=_dot(e,d)<=_dot(d,d)
    if len(common)==3:
        return _on_shared_face(point,common)
    raise ValueError("coincident or invalid tetrahedral mesh simplex")


@dataclass(frozen=True)
class TetrahedralCollisionAudit:
    pair_count:int
    strict_aabb_rejections:int
    legitimate_vertex_pairs:int
    legitimate_edge_pairs:int
    legitimate_face_pairs:int

    def __post_init__(self)->None:
        if (any(type(x) is not int or x<0
                for x in (self.pair_count,self.strict_aabb_rejections,
                          self.legitimate_vertex_pairs,
                          self.legitimate_edge_pairs,
                          self.legitimate_face_pairs))
                or (self.strict_aabb_rejections+
                    self.legitimate_vertex_pairs+
                    self.legitimate_edge_pairs+
                    self.legitimate_face_pairs>self.pair_count)):
            raise ValueError("bad geometric tetrahedron collision metrics")


def audit_embedded_tetrahedra(
    vertices:tuple[Point3,...],
    tetrahedra:tuple[tuple[int,int,int,int],...],
)->TetrahedralCollisionAudit:
    """Reject nonconforming geometric intersections for every pair.

    This is independent of the tetrahedron-manifold topology checker:
    it directly computes every closed convex 3D simplex intersection.
    """
    if (not isinstance(vertices,tuple) or not isinstance(tetrahedra,tuple)
            or not tetrahedra or len(vertices)!=len(set(vertices))
            or any(not isinstance(v,Point3) for v in vertices)):
        raise ValueError("physical distinct Q3 vertices and cells required")
    cells=[]
    for tet in tetrahedra:
        if (not isinstance(tet,tuple) or len(tet)!=4 or
                len(set(tet))!=4 or
                any(type(i) is not int or not 0<=i<len(vertices) for i in tet)):
            raise ValueError("malformed tetrahedral cell")
        obj=tuple(vertices[i] for i in tet)
        _planes(obj)
        cells.append(obj)
    metrics=[0,0,0,0,0]
    for i,a in enumerate(cells):
        for b in cells[i+1:]:
            metrics[0]+=1
            if _aabb_separate(a,b):
                metrics[1]+=1
                continue
            common=tuple(sorted(set(a)&set(b),key=_xyz))
            if len(common)==4:
                raise ValueError("two duplicate tetrahedra occupy the same solid")
            points=tetrahedron_intersection_vertices(a,b)
            if points and (not common or any(
                    not _within_shared_simplex(p,common) for p in points)):
                raise ValueError("unmatched physical tetrahedra cross or overlap")
            if common:
                # The intersection includes every actual common vertex;
                # since both convex tets contain their convex hull,
                # bounds on all extreme points imply exact equality.
                if not points or not set(common).issubset(points):
                    raise ValueError("physical 3D mesh has an invalid shared simplex")
                metrics[len(common)+1]+=1
    return TetrahedralCollisionAudit(*metrics)
