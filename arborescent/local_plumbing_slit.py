"""A5.2e2b4c1: exact local ambient-3D PL product slitting.

Bonahon--Siebenmann §12.1 Figure 12.5 (printed p.140) marks TWO
transverse arrows on a surface plumbing rectangle. Here an actual
two-triangle rectangle on an already embedded candidate ribbon is
given an exact oriented (u,v) chart. A 3D product ball P x [-h,h] is
triangulated into 12 tetrahedra, then split ALONG P x {0} into two
side-labelled 3D pieces. A separate checker replays the exact cut and
inverse reglue, including volumes, face incidence and oriented cuts.

LOCAL ONLY: no ambient S3 cut, no global 4-punctured Conway sphere,
no source-specified attachment/framing, no ambient isotopy, no R+/-.
All public statuses prohibit source-topology certification.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .spatial_twist_blocks import Point3
from .exterior_conway_collars import ExteriorConwayCollar

Triangle = tuple[int,int,int]
Tet = tuple[int,int,int,int]
Vec = tuple[Fraction,Fraction,Fraction]


def _vec(a: Point3,b: Point3) -> Vec:
    return tuple(Fraction(getattr(b,k))-Fraction(getattr(a,k))
                 for k in ("x","y","z"))  # type: ignore[return-value]


def _cross(u: Vec,v: Vec) -> Vec:
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def _dot(u: Vec,v: Vec) -> Fraction:
    return sum((a*b for a,b in zip(u,v)),Fraction(0))


def _volume6(a: Point3,b: Point3,c: Point3,d: Point3) -> Fraction:
    return _dot(_cross(_vec(a,b),_vec(a,c)),_vec(a,d))


def _tet_faces(t: Tet) -> tuple[Triangle,...]:
    a,b,c,d=t
    return (b,c,d),(a,d,c),(a,b,d),(a,c,b)


def _product_tets() -> tuple[Tet,...]:
    """Canonical consistent staircase tetrahedra for two triangular prisms."""
    result=[]
    for low,top in ((0,4),(4,8)):
        for a,b,c in ((0,1,2),(0,2,3)):
            i,j,k=(low+a,low+b,low+c)
            I,J,K=(top+a,top+b,top+c)
            result.extend(((i,j,k,K),(i,J,j,K),(i,I,J,K)))
    return tuple(result)


@dataclass(frozen=True)
class LocalPlumbingPatch:
    """Chosen, source-UNVERIFIED oriented rectangle on a physical PL band."""
    corners: tuple[Point3,Point3,Point3,Point3]
    origin: Point3
    u_tip: Point3
    v_tip: Point3
    status: Literal["LOCAL_ARROW_PATCH_SOURCE_UNVERIFIED"] = (
        "LOCAL_ARROW_PATCH_SOURCE_UNVERIFIED"
    )

    def __post_init__(self) -> None:
        if (not isinstance(self.corners,tuple) or len(self.corners)!=4
                or any(not isinstance(p,Point3) for p in self.corners)
                or any(not isinstance(p,Point3)
                           for p in (self.origin,self.u_tip,self.v_tip))
                or self.status!="LOCAL_ARROW_PATCH_SOURCE_UNVERIFIED"):
            raise ValueError("exact oriented PL patch required")
        a,b,c,d=self.corners
        if not (a.x<b.x and a.y<d.y and a.y==b.y and c.y==d.y
                and a.x==d.x and b.x==c.x and a.z==b.z==c.z==d.z):
            raise ValueError("nondegenerate counterclockwise xy rectangle required")
        center=Point3((Fraction(a.x)+Fraction(c.x))/2,
                      (Fraction(a.y)+Fraction(c.y))/2,a.z)
        expected_u=Point3(center.x+(Fraction(b.x)-Fraction(a.x))/4,
                          center.y,center.z)
        expected_v=Point3(center.x,
                          center.y+(Fraction(d.y)-Fraction(a.y))/4,
                          center.z)
        if (self.origin,self.u_tip,self.v_tip)!=(center,expected_u,expected_v):
            raise ValueError("source-inspired transverse arrows invalid")
        u,v=_vec(center,self.u_tip),_vec(center,self.v_tip)
        if _dot(u,v)!=0 or _cross(u,v)[2]<=0:
            raise ValueError("arrow vectors must be orthogonal with normal +z")

    def triangles(self) -> tuple[tuple[Point3,Point3,Point3],...]:
        a,b,c,d=self.corners
        return ((a,b,c),(a,c,d))

    def require_source_plumbing_patch(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: roof patch is a local chart, "
            "not source-identified Figure 12.5 plumbing")


def collar_plumbing_chart(collar: ExteriorConwayCollar) -> LocalPlumbingPatch:
    """Read a real two-triangle rectangle from the existing band mesh."""
    if not isinstance(collar,ExteriorConwayCollar):
        raise ValueError("physical exterior Conway collar required")
    collar.audit()
    a,b,c,d=(collar.vertices[i] for i in (3,6,7,4))
    center=Point3((Fraction(a.x)+Fraction(c.x))/2,
                  (Fraction(a.y)+Fraction(c.y))/2,a.z)
    patch=LocalPlumbingPatch(
        (a,b,c,d),center,
        Point3(center.x+(Fraction(b.x)-Fraction(a.x))/4,center.y,center.z),
        Point3(center.x,center.y+(Fraction(d.y)-Fraction(a.y))/4,center.z))
    mesh={frozenset(f) for f in collar.triangles()}
    if any(frozenset(t) not in mesh for t in patch.triangles()):
        raise ValueError("rectangle not actually triangulated on physical collar")
    return patch


@dataclass(frozen=True)
class TetrahedralProductNeighborhood:
    patch: LocalPlumbingPatch
    half_thickness: Fraction
    vertices: tuple[Point3,...]
    tetrahedra: tuple[Tet,...]
    status: Literal["LOCAL_PRODUCT_3BALL_ONLY"] = "LOCAL_PRODUCT_3BALL_ONLY"

    def __post_init__(self) -> None:
        self.audit()

    def midplane_faces(self) -> tuple[Triangle,Triangle]:
        return (4,5,6),(4,6,7)

    def audit(self) -> None:
        if (not isinstance(self.patch,LocalPlumbingPatch)
                or type(self.half_thickness) is not Fraction
                or self.half_thickness<=0
                or self.status!="LOCAL_PRODUCT_3BALL_ONLY"):
            raise ValueError("positive exact local product neighborhood required")
        h=self.half_thickness
        expected=tuple(
            Point3(p.x,p.y,Fraction(p.z)+delta)
            for delta in (-h,Fraction(0),h)
            for p in self.patch.corners)
        if self.vertices!=expected or len(set(self.vertices))!=12:
            raise ValueError("bad Q3 local product vertices")
        if self.tetrahedra!=_product_tets():
            raise ValueError("canonical 3-ball tetrahedra changed")
        incidence=Counter()
        orientations=defaultdict(list)
        volume=Fraction(0)
        for tet in self.tetrahedra:
            if len(set(tet))!=4 or any(type(i) is not int or not 0<=i<12
                                       for i in tet):
                raise ValueError("invalid tetrahedron")
            det=_volume6(*(self.vertices[i] for i in tet))
            if det<=0:
                raise ValueError("non-positive exact tetrahedron volume")
            volume+=det/6
            for face in _tet_faces(tet):
                incidence[frozenset(face)]+=1
                orientations[frozenset(face)].append(face)
        a,b,_,d=self.patch.corners
        expected_volume=(Fraction(b.x)-Fraction(a.x))*(
            Fraction(d.y)-Fraction(a.y))*2*h
        if volume!=expected_volume:
            raise ValueError("tetrahedra do not fill the 3D product volume")
        if any(n not in (1,2) for n in incidence.values()):
            raise ValueError("bad tetrahedral face incidence")
        for faces in orientations.values():
            if len(faces)==2:
                p,q=faces
                if q not in ((p[1],p[0],p[2]),(p[0],p[2],p[1]),
                              (p[2],p[1],p[0])):
                    raise ValueError("interior tetrahedral faces not oppositely oriented")
        if any(incidence[frozenset(face)]!=2
               for face in self.midplane_faces()):
            raise ValueError("midplane is not an internal cut disk")
        exterior=[faces[0] for faces in orientations.values()
                  if len(faces)==1]
        edges=Counter()
        for i,j,k in exterior:
            for u,v in ((i,j),(j,k),(k,i)):
                edges[frozenset((u,v))]+=1
        boundary_vertices={i for f in exterior for i in f}
        if (any(n!=2 for n in edges.values())
                or len(boundary_vertices)-len(edges)+len(exterior)!=2):
            raise ValueError("local product boundary is not Euler-2 closed")


def product_neighborhood(
    patch: LocalPlumbingPatch, half_thickness: Fraction=Fraction(1,16)
) -> TetrahedralProductNeighborhood:
    if (not isinstance(patch,LocalPlumbingPatch)
            or type(half_thickness) is not Fraction or half_thickness<=0):
        raise ValueError("positive rational thickness and framed patch required")
    points=tuple(Point3(p.x,p.y,Fraction(p.z)+delta)
                 for delta in (-half_thickness,Fraction(0),half_thickness)
                 for p in patch.corners)
    return TetrahedralProductNeighborhood(
        patch,half_thickness,points,_product_tets())


@dataclass(frozen=True,order=True)
class CutVertex:
    original_index: int
    side: Literal["LOWER","UPPER","EXTERIOR"]

    def __post_init__(self) -> None:
        if (type(self.original_index) is not int or
                not 0<=self.original_index<12
                or self.side not in ("LOWER","UPPER","EXTERIOR")):
            raise ValueError("invalid cut vertex/side")


def _audit_cut_halves(
    tetrahedra: tuple[tuple[CutVertex,CutVertex,CutVertex,CutVertex],...],
    lower_faces: tuple[tuple[CutVertex,CutVertex,CutVertex],...],
    upper_faces: tuple[tuple[CutVertex,CutVertex,CutVertex],...]
) -> bool:
    """Independently certify TWO cut 3-ball boundary 2-spheres.

    Count actual cut-cell incidences instead of trusting source
    midplane labels, including every face and closed boundary edge.
    """
    for side,declared in (("LOWER",lower_faces),("UPPER",upper_faces)):
        tets=[t for t in tetrahedra
              if any(v.side==side for v in t)]
        if len(tets)!=6 or any(any(v.side not in (side,"EXTERIOR")
                                  for v in tet) for tet in tets):
            return False
        face_count=Counter()
        oriented=defaultdict(list)
        for tet in tets:
            for face in _tet_faces(tet):
                face_count[frozenset(face)]+=1
                oriented[frozenset(face)].append(face)
        if any(n not in (1,2) for n in face_count.values()):
            return False
        for pair in oriented.values():
            if len(pair)==2:
                p,q=pair
                if q not in ((p[1],p[0],p[2]),(p[0],p[2],p[1]),
                              (p[2],p[1],p[0])):
                    return False
        boundary=[v[0] for v in oriented.values() if len(v)==1]
        if len(boundary)!=12:
            return False
        actual_cut={frozenset(f) for f in boundary
                    if all(v.side==side for v in f)}
        if actual_cut != {frozenset(f) for f in declared}:
            return False
        vertices={v for f in boundary for v in f}
        edges=Counter()
        for i,j,k in boundary:
            for a,b in ((i,j),(j,k),(k,i)):
                edges[frozenset((a,b))]+=1
        if (len(vertices)!=8 or any(n!=2 for n in edges.values())
                or len(vertices)-len(edges)+len(boundary)!=2):
            return False
    return True


@dataclass(frozen=True)
class LocalSlitWitness:
    """Independent replay object for the *local* ambient PL 3-ball cut."""
    solid: TetrahedralProductNeighborhood
    cut_tetrahedra: tuple[tuple[CutVertex,CutVertex,CutVertex,CutVertex],...]
    lower_cut: tuple[tuple[CutVertex,CutVertex,CutVertex],...]
    upper_cut: tuple[tuple[CutVertex,CutVertex,CutVertex],...]
    status: Literal["LOCAL_SLIT_REGLUE_ONLY"] = "LOCAL_SLIT_REGLUE_ONLY"

    def __post_init__(self) -> None:
        if not self.verify():
            raise ValueError("invalid exact local slitting certificate")

    def verify(self) -> bool:
        if (not isinstance(self.solid,TetrahedralProductNeighborhood)
                or self.status!="LOCAL_SLIT_REGLUE_ONLY"):
            return False
        try:
            self.solid.audit()
            check=[]
            for tet in self.solid.tetrahedra:
                lower=all(i<8 for i in tet)
                upper=all(i>=4 for i in tet)
                if lower==upper:
                    return False
                side="LOWER" if lower else "UPPER"
                check.append(tuple(CutVertex(i,side if 4<=i<8 else "EXTERIOR")
                                   for i in tet))
            if self.cut_tetrahedra!=tuple(check):
                return False
            mid=self.solid.midplane_faces()
            low_faces=tuple(tuple(CutVertex(i,"LOWER") for i in f) for f in mid)
            upper_faces=tuple(tuple(CutVertex(i,"UPPER") for i in
                                      (f[0],f[2],f[1])) for f in mid)
            if self.lower_cut!=low_faces or self.upper_cut!=upper_faces:
                return False
            if not _audit_cut_halves(self.cut_tetrahedra,
                                     self.lower_cut,self.upper_cut):
                return False
            for lo,hi in zip(self.lower_cut,self.upper_cut):
                if frozenset(self.solid.vertices[v.original_index] for v in lo
                        )!=frozenset(self.solid.vertices[v.original_index] for v in hi):
                    return False
            if set(v for face in self.lower_cut for v in face) & set(
                    v for face in self.upper_cut for v in face):
                return False
            # Independently reglue the distinct side-labelled vertices
            # back to the original 3D tetrahedral incidence table.
            return tuple(tuple(v.original_index for v in tet)
                         for tet in self.cut_tetrahedra)==self.solid.tetrahedra
        except (ValueError,TypeError,IndexError,AttributeError):
            return False

    def require_source_conway_slitting(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: local product-ball split does not "
            "construct global four-punctured Conway boundary spheres "
            "or source-relative 3D tangle isotopy")


def slit_product(solid: TetrahedralProductNeighborhood) -> LocalSlitWitness:
    if not isinstance(solid,TetrahedralProductNeighborhood):
        raise ValueError("typed exact local product 3-ball required")
    cut=[]
    for tet in solid.tetrahedra:
        side="LOWER" if all(i<8 for i in tet) else "UPPER"
        cut.append(tuple(CutVertex(i,side if 4<=i<8 else "EXTERIOR")
                         for i in tet))
    mid=solid.midplane_faces()
    lower=tuple(tuple(CutVertex(i,"LOWER") for i in f) for f in mid)
    upper=tuple(tuple(CutVertex(i,"UPPER") for i in
                          (f[0],f[2],f[1])) for f in mid)
    return LocalSlitWitness(solid,tuple(cut),lower,upper)
