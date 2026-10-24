"""A5.2e1: exact global return STRANDS of an r=0 atomic-tangle candidate.

Source: Bonahon--Siebenmann Fig. 12.8 (p.142) and Fig. 12.29
(p.158), which have n four-punctured Conway boundaries, n signed
half-twist sectors, and long outer return tracks.

Unlike A5.2d, there is NO auxiliary anchor sphere and NO unjoined
strand endpoint. Two full return strands are embedded in Q^3,
around exactly n actual marked Conway-sphere boundaries.

But a pair of return *strands* is not the source's whole twisted
*band*: the ribbon surface, framing, source-compatible Conway maps,
independent source-diagram equivalence, and any r>0 rings are absent.
This is a candidate, not a certified source realization or R± proof.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction, Port
from .conway_braid_diagrams import Braid2Diagram
from .marked_conway_geometry import MarkedConwayCube, TwoStrandCorridor
from .spatial_twist_blocks import Point3, SpatialTwoBraid, spatial_two_braid

Arc = tuple[Point3, ...]
Pairing = tuple[tuple[Port, str], tuple[Port, str]]


def _vector(p: Point3, q: Point3) -> tuple[Fraction, Fraction, Fraction]:
    return tuple(Fraction(getattr(q, k)) - Fraction(getattr(p, k))
                 for k in ("x", "y", "z"))  # type: ignore[return-value]


def _cross(a: tuple[Fraction, ...], b: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _dot(a: tuple[Fraction, ...], b: tuple[Fraction, ...]) -> Fraction:
    return sum((x*y for x,y in zip(a,b)), Fraction(0))


def _segment_contact(p: Point3, q: Point3, r: Point3, s: Point3
                     ) -> tuple[str, Point3 | None]:
    """Independent exact-Q³ closed-segment contact classification.

    Returns none / point / overlap. Handles skew, parallel, collinear,
    endpoint and rational-coordinate cases, without a floating epsilon.
    A positive-length shared subsegment is ALWAYS an invalid overlap.
    """
    if p == q or r == s:
        raise ValueError("zero-length spatial segment")
    for dim in ("x","y","z"):
        if max(getattr(p,dim),getattr(q,dim)) < min(getattr(r,dim),getattr(s,dim)):
            return "none", None
        if max(getattr(r,dim),getattr(s,dim)) < min(getattr(p,dim),getattr(q,dim)):
            return "none", None
    a = _vector(p,q)
    b = _vector(r,s)
    w = _vector(p,r)
    cross = _cross(a,b)
    if cross != (0,0,0):
        if _dot(w,cross) != 0:
            return "none", None
        denom = _dot(cross,cross)
        t = _dot(_cross(w,b), cross)/denom
        u = _dot(_cross(w,a), cross)/denom
        if not (0 <= t <= 1 and 0 <= u <= 1):
            return "none", None
        return "point", Point3(*(Fraction(getattr(p,k))+t*a[i]
                                 for i,k in enumerate(("x","y","z"))))
    if _cross(w,a) != (0,0,0):
        return "none", None
    axis = next(i for i in range(3) if a[i] != 0)
    t1 = w[axis]/a[axis]
    t2 = (Fraction(getattr(s,("x","y","z")[axis]))
          -Fraction(getattr(p,("x","y","z")[axis])))/a[axis]
    lo, hi = max(Fraction(0),min(t1,t2)), min(Fraction(1),max(t1,t2))
    if lo > hi:
        return "none", None
    if lo < hi:
        return "overlap", None
    return "point", Point3(*(Fraction(getattr(p,k))+lo*a[i]
                             for i,k in enumerate(("x","y","z"))))


def _segment_cube_contact(a: Point3, b: Point3, cube: MarkedConwayCube
                          ) -> tuple[str, Point3 | None]:
    """Exact segment/closed-solid-cube clipping.

    Any interval of positive length inside a solid cube is forbidden.
    Single face contact is permitted ONLY at a marked terminal point.
    """
    low = Fraction(0)
    high = Fraction(1)
    bounds = ((cube.west,cube.east),(-3,3),(-3,3))
    for axis,(lo,hi) in zip(("x","y","z"),bounds):
        p = Fraction(getattr(a,axis))
        d = Fraction(getattr(b,axis))-p
        if d == 0:
            if not lo <= p <= hi:
                return "none", None
            continue
        start,end = Fraction(lo-p,d),Fraction(hi-p,d)
        low = max(low,min(start,end))
        high = min(high,max(start,end))
        if low > high:
            return "none", None
    if low < high:
        return "overlap", None
    point = Point3(*(Fraction(getattr(a,k))
                     + low*(Fraction(getattr(b,k))-Fraction(getattr(a,k)))
                     for k in ("x","y","z")))
    return "point", point


def independent_pl_arc_audit(cubes: tuple[MarkedConwayCube, ...],
                             arcs: tuple[Arc, ...]) -> tuple[Pairing, ...]:
    """Independently check global PL embed, not constructor's layout metadata.

    Requires *exactly* 4n marked endpoints, each used once, all
    interiors disjoint from all closed Conway balls and each other,
    and no self-intersections. Returns all checked boundary pairings.
    """
    if (not isinstance(cubes, tuple) or not cubes
            or not isinstance(arcs, tuple)):
        raise ValueError("typed sphere and arc tuples required")
    coordinate_labels: dict[Point3, tuple[Port,str]] = {}
    for cube in cubes:
        if not isinstance(cube, MarkedConwayCube) or cube.port is None:
            raise ValueError("all cubes must have genuine A5.1 ports")
        cube.audit()
        for name,point in cube.punctures():
            if point in coordinate_labels:
                raise ValueError("two punctures at one spatial coordinate")
            coordinate_labels[point] = (cube.port,name)
    if len(arcs)*2 != len(coordinate_labels):
        raise ValueError("wrong total boundary strand count")

    segments: list[tuple[int,int,Point3,Point3]] = []
    endpoint_counts: Counter[Point3] = Counter()
    pairings: list[Pairing] = []
    for i,arc in enumerate(arcs):
        if (not isinstance(arc, tuple) or len(arc) < 2
                or any(not isinstance(p, Point3) for p in arc)):
            raise ValueError("invalid embedded polyline")
        if arc[0] not in coordinate_labels or arc[-1] not in coordinate_labels:
            raise ValueError("all global arcs must terminate on Conway marks")
        if arc[0] == arc[-1]:
            raise ValueError("arcs may not close at one puncture")
        endpoint_counts.update((arc[0],arc[-1]))
        pairings.append((coordinate_labels[arc[0]],coordinate_labels[arc[-1]]))
        for j,(p,q) in enumerate(zip(arc,arc[1:])):
            if p == q:
                raise ValueError("degenerate line segment")
            for cube in cubes:
                kind,point = _segment_cube_contact(p,q,cube)
                if kind == "overlap" or (kind == "point" and
                    (point not in (arc[0],arc[-1])
                     or point not in coordinate_labels)):
                    raise ValueError("arc penetrates or touches Conway sphere at nonendpoint")
            segments.append((i,j,p,q))
    if (len(endpoint_counts) != len(coordinate_labels)
            or set(endpoint_counts) != set(coordinate_labels)
            or any(c != 1 for c in endpoint_counts.values())):
        raise ValueError("Conway puncture used missing or multiple times")

    for idx,(ai,si,p,q) in enumerate(segments):
        for bi,sj,r,s in segments[idx+1:]:
            kind,point = _segment_contact(p,q,r,s)
            if kind == "none":
                continue
            # Adjacent segments of the same embedded polyline must
            # meet at exactly their common enumerated vertex.
            if (kind == "point" and ai == bi and sj == si+1
                    and point == q == r):
                continue
            raise ValueError("two embedded PL segments intersect")
    return tuple(pairings)


def _local_twist(weight: int) -> SpatialTwoBraid:
    return spatial_two_braid(Braid2Diagram(
        (1,)*weight if weight >= 0 else (-1,)*(-weight)))


def _entry_strands(twist: SpatialTwoBraid) -> tuple[Arc,Arc]:
    return tuple(tuple(Point3(p.y,p.x,p.z) for p in path)
                 for path in (twist.nw,twist.ne))  # type: ignore[return-value]


def _rail(cap_sign: int, right_x: int) -> Arc:
    if cap_sign not in (-1,1):
        raise ValueError("return rail must be upper or lower")
    y = 2*cap_sign
    far = 6*cap_sign
    return (Point3(0,y,0),Point3(-3,y,0),Point3(-3,far,0),
            Point3(right_x+3,far,0),Point3(right_x+3,y,0),
            Point3(right_x,y,0))


@dataclass(frozen=True)
class ZeroRingReturnStrandCandidate:
    """n actual Conway spheres, n signed braid sectors, 2n open arcs.

    Return rails complete the *strands* but no globally twisted
    supporting ribbon/framing has been constructed.
    """
    atom: AtomicInstruction
    cubes: tuple[MarkedConwayCube, ...]
    entry: SpatialTwoBraid
    internal: tuple[TwoStrandCorridor, ...]
    return_rails: tuple[Arc,Arc]
    status: Literal["RETURN_STRANDS_SOURCE_UNVERIFIED"] = "RETURN_STRANDS_SOURCE_UNVERIFIED"

    def __post_init__(self) -> None:
        self.audit()

    def embedded_arcs(self) -> tuple[Arc, ...]:
        entry = _entry_strands(self.entry)
        returns = tuple(tuple(reversed(self.return_rails[i])) + entry[i][1:]
                        for i in (0,1))
        inner = tuple(arc for corridor in self.internal
                      for arc in corridor.strands())
        return returns + inner

    def boundary_pairings(self) -> tuple[Pairing,...]:
        return independent_pl_arc_audit(self.cubes,self.embedded_arcs())

    def canonical_coordinates(self) -> tuple:
        return (self.atom.address,self.atom.sign,self.atom.sector_weights,
                tuple((cube.center_x,cube.port,
                       tuple((name,p.serial()) for name,p in cube.punctures()))
                      for cube in self.cubes),
                tuple(tuple(p.serial() for p in arc) for arc in self.embedded_arcs()))

    def audit(self) -> None:
        if (not isinstance(self.atom,AtomicInstruction)
                or self.atom.ring_count != 0 or self.atom.valence < 2
                or self.status != "RETURN_STRANDS_SOURCE_UNVERIFIED"
                or len(self.cubes) != self.atom.valence
                or len(self.internal) != self.atom.valence-1
                or not isinstance(self.entry,SpatialTwoBraid)
                or not isinstance(self.return_rails,tuple)
                or len(self.return_rails) != 2):
            raise ValueError("invalid r=0 return-strand candidate")
        c = 3+2*abs(self.atom.sector_weights[0])
        expected_centers = [c]
        for weight in self.atom.sector_weights[1:]:
            c += 5+2*abs(weight)
            expected_centers.append(c)
        if tuple(cube.center_x for cube in self.cubes) != tuple(expected_centers):
            raise ValueError("sphere configuration changed")
        if tuple(cube.port for cube in self.cubes) != tuple(
            Port(self.atom.address,i) for i in range(self.atom.valence)):
            raise ValueError("incorrect source A5.1 port sequence")
        self.entry.audit()
        if self.entry.word != _local_twist(self.atom.sector_weights[0]).word:
            raise ValueError("left signed twist does not match a0")
        for i,corridor in enumerate(self.internal,1):
            corridor.audit()
            if (corridor.left != self.cubes[i-1]
                    or corridor.right != self.cubes[i]
                    or corridor.sector_index != i
                    or corridor.signed_weight != self.atom.sector_weights[i]):
                raise ValueError("internal twist corridor/source mismatch")
        if (self.return_rails != (_rail(-1,self.cubes[-1].east),
                                 _rail(1,self.cubes[-1].east))):
            raise ValueError("return rails differ from certified layout")
        if (tuple(arc[-1] for arc in _entry_strands(self.entry))
                != (Point3(self.cubes[0].west,-2,0),
                    Point3(self.cubes[0].west,2,0))
                and set(arc[-1] for arc in _entry_strands(self.entry))
                != {Point3(self.cubes[0].west,-2,0),
                    Point3(self.cubes[0].west,2,0)}):
            raise ValueError("entry twist fails to hit first sphere")
        if (tuple(arc[0] for arc in _entry_strands(self.entry))
                != (Point3(0,-2,0),Point3(0,2,0))):
            raise ValueError("left cap endpoints altered")
        pairings=self.boundary_pairings()
        if len(pairings) != 2*self.atom.valence:
            raise ValueError("not all Conway endpoints used")

    def require_source_atomic_realization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: return STRANDS only, not source-verified "
            "twisted ribbon, Conway framing, 3D gluing or R± isotopy"
        )


def realize_zero_ring_return_strands(atom: AtomicInstruction
                                     ) -> ZeroRingReturnStrandCandidate:
    """Build source-shaped r=0 candidate arcs with no auxiliary sphere."""
    if not isinstance(atom,AtomicInstruction):
        raise ValueError("A5.1 atomic instruction required")
    if atom.ring_count:
        raise NotImplementedError("UNVERIFIED_GEOMETRY: ring embeddings unavailable")
    if atom.valence < 2:
        raise NotImplementedError("UNVERIFIED_GEOMETRY: valence-one model unavailable")
    centers = [3+2*abs(atom.sector_weights[0])]
    for weight in atom.sector_weights[1:]:
        centers.append(centers[-1]+5+2*abs(weight))
    cubes = tuple(MarkedConwayCube(center,Port(atom.address,i))
                  for i,center in enumerate(centers))
    entry = _local_twist(atom.sector_weights[0])
    internal = tuple(
        TwoStrandCorridor(cubes[i-1],cubes[i],i,weight,_local_twist(weight))
        for i,weight in enumerate(atom.sector_weights) if i>0)
    return ZeroRingReturnStrandCandidate(
        atom,cubes,entry,internal,
        (_rail(-1,cubes[-1].east),_rail(1,cubes[-1].east)))
