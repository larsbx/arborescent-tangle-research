"""A5.2e2a: genuine triangulated PL ribbon patches on r=0 returning arcs.

For each pair of the 2n previously independently embedded strands, this
adds a 2-dimensional strip with exact rational vertices and oriented
triangles. All n patches are disjoint from Conway ball *interiors* and
one another under the chosen source-UNVERIFIED coordinate convention.

This is a ribbon PATCH SYSTEM, not the source's unique global twisted
band, source parametrizations or an ambient-isotopy equivalence proof.
No ring construction, 3D Conway sewing or R+/- certification follows.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_return_strands import (
    Arc, ZeroRingReturnStrandCandidate, realize_zero_ring_return_strands,
)
from .atomic_plumbing_contract import AtomicInstruction
from .marked_conway_geometry import MarkedConwayCube
from .spatial_twist_blocks import Point3


def _coord(p: Point3, axis: str) -> Fraction:
    return Fraction(getattr(p, axis))


def _diff(p: Point3, q: Point3) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(_coord(q, a)-_coord(p,a) for a in ("x","y","z"))


def _cross(a: tuple[Fraction,...], b: tuple[Fraction,...]
          ) -> tuple[Fraction,Fraction,Fraction]:
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _mid(a: Point3, b: Point3, height: int | None = None) -> Point3:
    p=Point3(*((_coord(a,d)+_coord(b,d))/2 for d in ("x","y","z")))
    return p if height is None else Point3(p.x,p.y,height)


def _in_segment(p: Point3, a: Point3, b: Point3) -> bool:
    """Exact point membership in a nondegenerate closed 3D segment."""
    d=_diff(a,b)
    e=_diff(a,p)
    if d==(0,0,0) or _cross(d,e)!=(0,0,0):
        return False
    index=next(k for k in range(3) if d[k])
    t=e[index]/d[index]
    return 0<=t<=1


def _clip_cube(triangle: tuple[Point3,Point3,Point3],
               cube: MarkedConwayCube) -> tuple[Point3,...]:
    """Independently clip a 3D triangle against a closed convex cube.

    Each clipping step computes exact rational intersections against
    a bounding half-space; returns polygon / segment / point / empty.
    """
    polygon=list(triangle)
    for axis,bound,sgn in (
        ("x",cube.west,1),("x",cube.east,-1),
        ("y",-3,1),("y",3,-1),("z",-3,1),("z",3,-1)
    ):
        if not polygon:
            break
        def value(p: Point3) -> Fraction:
            return sgn*(_coord(p,axis)-Fraction(bound))
        new=[]
        for i,q in enumerate(polygon):
            p=polygon[i-1]
            fp,fq=value(p),value(q)
            if (fp<0 and fq>=0) or (fp>=0 and fq<0):
                t=fp/(fp-fq)
                cut=Point3(*(_coord(p,k)+t*(_coord(q,k)-_coord(p,k))
                             for k in ("x","y","z")))
                new.append(cut)
            if fq>=0:
                new.append(q)
        polygon=[]
        for p in new:
            if not polygon or p!=polygon[-1]:
                polygon.append(p)
        if len(polygon)>1 and polygon[0]==polygon[-1]:
            polygon.pop()
    return tuple(polygon)


@dataclass(frozen=True)
class TriangulatedRibbon:
    """An explicit PL disk, with boundary rails and end-cap chords.

    Coordinates, ridge, AND triangle-index table are stored, not
    inferred from metadata at verification time.
    """
    kind: Literal["EXTERIOR_RETURN","SECTOR_CORRIDOR"]
    first: Arc
    second: Arc
    ridge: Arc
    faces: tuple[tuple[int,int,int], ...]

    def __post_init__(self) -> None:
        self.audit()

    @property
    def vertices(self) -> tuple[Point3,...]:
        return tuple(p for triple in zip(self.first,self.ridge,self.second)
                     for p in triple)

    def triangles(self) -> tuple[tuple[Point3,Point3,Point3],...]:
        v=self.vertices
        return tuple((v[a],v[b],v[c]) for a,b,c in self.faces)

    def boundary_paths(self) -> tuple[Arc,Arc,Arc,Arc]:
        """Physical boundary reconstructed from typed edge coordinates."""
        return (self.first,self.second,
                (self.first[0],self.ridge[0],self.second[0]),
                (self.first[-1],self.ridge[-1],self.second[-1]))

    def audit(self) -> None:
        if (self.kind not in ("EXTERIOR_RETURN","SECTOR_CORRIDOR")
                or not isinstance(self.first,tuple)
                or not isinstance(self.second,tuple)
                or not isinstance(self.ridge,tuple)
                or not 2<=len(self.first)==len(self.second)==len(self.ridge)
                or not all(isinstance(p,Point3)
                           for p in (*self.first,*self.second,*self.ridge))
                or not isinstance(self.faces,tuple)):
            raise ValueError("typed equally sampled PL ribbon rails required")
        n=len(self.first)
        v=self.vertices
        if len(set(v))!=len(v):
            raise ValueError("PL strip has coincident cross-section vertices")
        # Validate the ACTUAL triangulated surface, not simply its edges.
        directed: Counter[tuple[int,int]]=Counter()
        undirected: Counter[frozenset[int]]=Counter()
        adjacency: dict[int,set[int]]=defaultdict(set)
        for face in self.faces:
            if (not isinstance(face,tuple) or len(face)!=3
                    or any(type(i) is not int or not 0<=i<len(v) for i in face)
                    or len(set(face))!=3):
                raise ValueError("invalid oriented triangle indices")
            a,b,c=(v[i] for i in face)
            if _cross(_diff(a,b),_diff(a,c))==(0,0,0):
                raise ValueError("zero-area PL ribbon triangle")
            for i,j in zip(face,face[1:]+face[:1]):
                directed[i,j]+=1
                undirected[frozenset((i,j))]+=1
                adjacency[i].add(j)
                adjacency[j].add(i)
        if (len(self.faces)!=4*(n-1)
                or any(k not in (1,2) for k in undirected.values())
                or any(directed[i,j]!=directed[j,i]
                       for e,count in undirected.items() if count==2
                       for i,j in (tuple(e),))
                or len(v)-len(undirected)+len(self.faces)!=1):
            raise ValueError("PL disk has invalid triangulation, orientation or Euler number")
        seen={0}
        queue=[0]
        while queue:
            for other in adjacency[queue.pop()]:
                if other not in seen:
                    seen.add(other)
                    queue.append(other)
        if len(seen)!=len(v):
            raise ValueError("disconnected ribbon mesh")
        boundary={e for e,k in undirected.items() if k==1}
        want=set()
        for path in self.boundary_paths():
            indices=[]
            for point in path:
                indices.append(v.index(point))
            want.update(frozenset((a,b)) for a,b in zip(indices,indices[1:]))
        if (boundary!=want or len(boundary)!=2*n+2
                or any(sum(i in e for e in boundary)!=2
                       for i in set().union(*boundary))):
            raise ValueError("PL strip boundary does not follow both physical rails")
        # The vertex link and Euler constraints alone do not prove
        # geometric embedding. These are checked by the model audit.

    def canonical_mesh(self) -> tuple:
        return (self.kind,tuple(p.serial() for p in self.vertices),self.faces)


def _ribbon(kind: Literal["EXTERIOR_RETURN","SECTOR_CORRIDOR"],
            first: Arc, second: Arc, ridge: Arc) -> TriangulatedRibbon:
    faces=[]
    for i in range(len(first)-1):
        a=3*i
        b=a+3
        faces.extend(((a,a+1,b+1),(a,b+1,b),
                      (a+1,a+2,b+2),(a+1,b+2,b+1)))
    return TriangulatedRibbon(kind,first,second,ridge,tuple(faces))


@dataclass(frozen=True)
class ZeroRingRibbonPatchSystem:
    """Exact PL band-PATCH candidate; source twisted-band equivalence open."""
    strands: ZeroRingReturnStrandCandidate
    ribbons: tuple[TriangulatedRibbon,...]
    status: Literal["RIBBON_PATCHES_SOURCE_UNVERIFIED"] = "RIBBON_PATCHES_SOURCE_UNVERIFIED"

    def __post_init__(self) -> None:
        self.audit()

    def audit(self) -> None:
        if (not isinstance(self.strands,ZeroRingReturnStrandCandidate)
                or self.status!="RIBBON_PATCHES_SOURCE_UNVERIFIED"
                or not isinstance(self.ribbons,tuple)
                or len(self.ribbons)!=self.strands.atom.valence):
            raise ValueError("invalid partial twisted-band patch system")
        self.strands.audit()
        physical=self.strands.embedded_arcs()
        all_cubes=self.strands.cubes
        for i,ribbon in enumerate(self.ribbons):
            if (not isinstance(ribbon,TriangulatedRibbon)
                    or ribbon.kind!=("EXTERIOR_RETURN" if i==0 else "SECTOR_CORRIDOR")
                    or (ribbon.first,ribbon.second)!=physical[2*i:2*i+2]):
                raise ValueError("ribbon boundary does not equal actual checked strands")
            ribbon.audit()
            expected=tuple(_mid(a,b,8 if i==0 and j in (2,3) else None)
                           for j,(a,b) in enumerate(zip(ribbon.first,ribbon.second)))
            if ribbon.ridge!=expected:
                raise ValueError("unproved noncanonical ribbon ridge or twist inserted")
            for triangle in ribbon.triangles():
                for cube in all_cubes:
                    clipped=_clip_cube(triangle,cube)
                    if not clipped:
                        continue
                    cap0=ribbon.boundary_paths()[2]
                    cap1=ribbon.boundary_paths()[3]
                    if not all(
                        any(_in_segment(p,a,b) for cap in (cap0,cap1)
                            for a,b in zip(cap,cap[1:]))
                        for p in clipped
                    ):
                        raise ValueError("ribbon surface enters a Conway ball")
        self._prove_patch_separation()

    def _prove_patch_separation(self) -> None:
        """Exact slab + height enclosure test for chosen fixed ribbons.

        Other ribbons occupy mutually disjoint x-gaps between cubes,
        |y|<=2, |z|<=1. Exterior-return triangles outside the arch
        lie strictly before the first cube or after the last cube.
        The long arch is the two triangles per side at rows 2->3:
        z=8-4|y|/3, so for |y|<=2 it has z>=16/3>1.
        This is a *construction-specific* algebraic disjointness proof.
        """
        outer=self.ribbons[0]
        if len(outer.first)<7:
            raise ValueError("return band lacks both external U-turns")
        arch_faces=set(range(4*2,4*3))
        for j,face in enumerate(outer.triangles()):
            if j in arch_faces:
                ys=[_coord(p,"y") for p in face]
                if not (all(y<=0 for y in ys) or all(y>=0 for y in ys)):
                    raise ValueError("arch triangle crosses its ridge side")
                if any(_coord(p,"z")!=8-Fraction(4,3)*abs(_coord(p,"y"))
                       for p in face):
                    raise ValueError("return arch clearance inequality unproved")
            else:
                lo=min(_coord(p,"x") for p in face)
                hi=max(_coord(p,"x") for p in face)
                for inner in self.ribbons[1:]:
                    ilo=min(_coord(p,"x") for p in inner.vertices)
                    ihi=max(_coord(p,"x") for p in inner.vertices)
                    if not (hi<ilo or lo>ihi):
                        raise ValueError("return strip overlaps an interior longitudinal slab")
        for i,inner in enumerate(self.ribbons[1:]):
            if any(abs(_coord(p,"y"))>2 or abs(_coord(p,"z"))>1
                   for p in inner.vertices):
                raise ValueError("inner ribbon exceeds verified braid corridor box")
            if any(not all(_coord(a,"x")<_coord(b,"x")
                           for a,b in zip(rail,rail[1:]))
                   for rail in (inner.first,inner.ridge,inner.second)):
                raise ValueError("inner ribbon lacks strictly monotone x slices")
            # At a fixed x, the width vector from first to second
            # never vanishes inside an interval of the PL construction.
            for j in range(len(inner.first)-1):
                d0=_diff(inner.first[j],inner.second[j])
                d1=_diff(inner.first[j+1],inner.second[j+1])
                if d0==(0,0,0) or d1==(0,0,0):
                    raise ValueError("collapsed ribbon cross section")
                # An affine vector interpolation can become zero only
                # when d0 and d1 point in opposite collinear directions.
                if _cross(d0,d1)==(0,0,0):
                    dot=sum(a*b for a,b in zip(d0,d1))
                    if dot<=0:
                        raise ValueError("ribbon cross-section collapses or reverses")
            if i:
                prev=self.ribbons[i]
                if max(_coord(p,"x") for p in prev.vertices)>=min(
                        _coord(p,"x") for p in inner.vertices):
                    raise ValueError("different corridor surfaces share an x slab")
        # Outer return's near ends are exactly its source candidate arcs.
        # The arch's |y|<=2 region stays above every corridor (z<=1).

    def canonical_meshes(self) -> tuple:
        return tuple(r.canonical_mesh() for r in self.ribbons)

    def require_source_atomic_realization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: triangulated ribbon PATCHES are not "
            "source-certified Fig.12.8 global band/framing or Conway gluing"
        )


def build_zero_ring_ribbon_patches(atom: AtomicInstruction
                                  ) -> ZeroRingRibbonPatchSystem:
    """Real 2D PL strip candidates with independent face/cube audits."""
    candidate=realize_zero_ring_return_strands(atom)
    arcs=candidate.embedded_arcs()
    patches=[]
    for i in range(candidate.atom.valence):
        a,b=arcs[2*i:2*i+2]
        ridge=tuple(_mid(p,q,8 if i==0 and j in (2,3) else None)
                    for j,(p,q) in enumerate(zip(a,b)))
        patches.append(_ribbon("EXTERIOR_RETURN" if i==0 else "SECTOR_CORRIDOR",
                               a,b,ridge))
    return ZeroRingRibbonPatchSystem(candidate,tuple(patches))
