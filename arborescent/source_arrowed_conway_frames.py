"""A5.2e2b4a — arrowed Conway marking and signed mapping-class obligations.

Bonahon–Siebenmann §12.1 printed p.141, Figs 12.6–12.7:
four endpoints are NOT sufficient for a parameterized Conway sphere.
An ORIENTED equatorial arc between two adjacent punctures specifies
the zero-slope arrow; the two arrows for plumbing emanate from one
common point along distinct adjacent arcs. §12.2 pp.142–143 then
uses these parametrizations when plumbing the atomic twisted band.

This implements exact SOURCE-INSPIRED boundary *combinatorics*, a
CHOSEN cubical spatial representative of its arrow, and a full
signed GL(2,Z) algebraic torus-cover boundary map (including the
orientation-reversing determinant -1 cases unavailable in A5.2b).
No class here claims a homeomorphism was geometrically realized on
our cubes, nor source-relative isotopy of the global atomic band.

The state is always ARROW_FRAMING_SOURCE_UNVERIFIED.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Literal

from .atomic_plumbing_contract import Port
from .conway_mapping_class import PureConwayMappingClass, BRANCH_COORDINATES
from .marked_conway_geometry import MarkedConwayCube
from .spatial_twist_blocks import Point3

PUNCTURE_NAMES: tuple[str,...] = ("P00", "P10", "P11", "P01")
ZERO_SLOPE_ARROW: tuple[str,str] = ("P00","P10")


def _primitive(a: int, b: int) -> tuple[int,int]:
    if (type(a) is not int or type(b) is not int or
            (a==0 and b==0) or gcd(a,b)!=1):
        raise ValueError("primitive nonzero integer slope required")
    return (-a,-b) if a<0 or a==0 and b<0 else (a,b)


@dataclass(frozen=True)
class ArrowedConwayBoundary:
    """Exact oriented puncture-to-puncture arrow with cyclic order.

    A source *standard* mark may be represented by the four P-labels,
    but no arbitrary physical cube mapping is thereby source-certified.
    This arrow only selects an adjacent equatorial arc; its complete
    isotopy class is additional structure not fixed by its endpoints.
    """
    port: Port
    arrow_start: str = ZERO_SLOPE_ARROW[0]
    arrow_end: str = ZERO_SLOPE_ARROW[1]
    cyclic: tuple[str,str,str,str] = PUNCTURE_NAMES
    status: Literal["ARROW_FRAMING_SOURCE_UNVERIFIED"] = (
        "ARROW_FRAMING_SOURCE_UNVERIFIED"
    )

    def __post_init__(self) -> None:
        if (not isinstance(self.port, Port)
                or self.cyclic != PUNCTURE_NAMES
                or self.status != "ARROW_FRAMING_SOURCE_UNVERIFIED"
                or self.arrow_start not in self.cyclic
                or self.arrow_end not in self.cyclic):
            raise ValueError("invalid source-arrow Conway frame")
        a,b=self.cyclic.index(self.arrow_start),self.cyclic.index(self.arrow_end)
        if (a-b)%4 not in (1,3):
            raise ValueError("equatorial arrow must connect cyclic neighbors")

    @property
    def source_zero_slope_arrow(self) -> bool:
        return (self.arrow_start,self.arrow_end)==ZERO_SLOPE_ARROW

    def require_source_parametrization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: arrow endpoints and chosen cubic arc "
            "do not certify a source Fig 12.6 Conway parametrization"
        )


@dataclass(frozen=True)
class CubicalArrowRepresentative:
    """The *chosen* physical meridian-side arc on an A5.2d cube.

    P00=NW, P10=NE, P11=SE, P01=SW in the existing (chosen)
    cubical coordinates. The arrow from P00 to P10 lies on west,
    y=-3 and east faces; it never goes THROUGH the closed solid cube.
    """
    cube: MarkedConwayCube
    frame: ArrowedConwayBoundary
    arc: tuple[Point3,...]
    status: Literal["CUBICAL_CHOICE_SOURCE_UNVERIFIED"] = (
        "CUBICAL_CHOICE_SOURCE_UNVERIFIED"
    )

    def __post_init__(self) -> None:
        if (not isinstance(self.cube,MarkedConwayCube)
                or not isinstance(self.frame,ArrowedConwayBoundary)
                or self.cube.port is None or self.frame.port != self.cube.port
                or self.status != "CUBICAL_CHOICE_SOURCE_UNVERIFIED"):
            raise ValueError("cube and arrow frame must share one physical port")
        pts=self.cube.punctures()
        marks=dict(zip(PUNCTURE_NAMES,(p for _,p in pts)))
        expected=(marks["P00"],
                  Point3(self.cube.west,-3,0),
                  Point3(self.cube.east,-3,0),
                  marks["P10"])
        if not self.frame.source_zero_slope_arrow or self.arc!=expected:
            raise ValueError("unverified cube arrow, reversed arc or wrong framing")
        # Check EVERY open segment is on an exposed *face*, not
        # through the solid ball; endpoints are exactly source marks.
        for a,b in zip(self.arc,self.arc[1:]):
            if a==b or not any(
                getattr(a,axis)==getattr(b,axis)==value
                for axis,value in (("x",self.cube.west),
                                   ("x",self.cube.east),("y",-3))
            ):
                raise ValueError("arrow leaves cubical Conway boundary")
        if not (self.arc[0]==marks["P00"] and
                self.arc[-1]==marks["P10"]):
            raise ValueError("arrow puncture endpoints changed")


def chosen_cubical_arrow(cube: MarkedConwayCube) -> CubicalArrowRepresentative:
    if not isinstance(cube,MarkedConwayCube) or cube.port is None:
        raise ValueError("genuine A5.1 four-puncture cube port required")
    frame=ArrowedConwayBoundary(cube.port)
    labels=dict(zip(PUNCTURE_NAMES,(p for _,p in cube.punctures())))
    return CubicalArrowRepresentative(
        cube,frame,(labels["P00"],Point3(cube.west,-3,0),
                    Point3(cube.east,-3,0),labels["P10"])
    )


@dataclass(frozen=True)
class SignedConwayBoundaryClass:
    """Algebraic signed 4-puncture self-map on the branched torus cover.

    An exact GL2(Z) lift of determinant +1 or -1, together with a
    binary translation of the four 2-torsion punctures. The linear
    reduction mod2 may permute them (not only V4). This is a proposed
    boundary *mapping class*, NOT a realized spatial 3D seam map.
    """
    a: int
    b: int
    c: int
    d: int
    tx: int = 0
    ty: int = 0

    def __post_init__(self) -> None:
        if (any(type(x) is not int for x in
                (self.a,self.b,self.c,self.d,self.tx,self.ty))
                or self.determinant not in (-1,1)
                or self.tx not in (0,1) or self.ty not in (0,1)):
            raise ValueError("signed Conway map needs GL2(Z) and 2-torsion shift")

    @property
    def determinant(self) -> int:
        return self.a*self.d-self.b*self.c

    def puncture_permutation(self) -> tuple[int,int,int,int]:
        """Compute all FOUR actual affine 2-torsion images independently."""
        result=[]
        for x,y in BRANCH_COORDINATES:
            image=((self.a*x+self.b*y+self.tx)%2,
                   (self.c*x+self.d*y+self.ty)%2)
            result.append(BRANCH_COORDINATES.index(image))
        if len(set(result))!=4:
            raise AssertionError("invertible mod 2 map must permute all punctures")
        return tuple(result)

    def image_puncture(self, label: str) -> str:
        if label not in PUNCTURE_NAMES:
            raise ValueError("unrecognized source puncture")
        return PUNCTURE_NAMES[self.puncture_permutation()[
            PUNCTURE_NAMES.index(label)]]

    def image_slope(self, p: int, q: int) -> tuple[int,int]:
        p,q=_primitive(p,q)
        return _primitive(self.a*p+self.b*q,
                          self.c*p+self.d*q)

    def inverse(self) -> "SignedConwayBoundaryClass":
        det=self.determinant
        a,b,c,d=(det*self.d,-det*self.b,-det*self.c,det*self.a)
        return SignedConwayBoundaryClass(
            a,b,c,d,(-a*self.tx-b*self.ty)%2,
            (-c*self.tx-d*self.ty)%2
        )

    def compose(self, other: "SignedConwayBoundaryClass"
                ) -> "SignedConwayBoundaryClass":
        if not isinstance(other,SignedConwayBoundaryClass):
            raise ValueError("typed signed Conway map required")
        return SignedConwayBoundaryClass(
            self.a*other.a+self.b*other.c,
            self.a*other.b+self.b*other.d,
            self.c*other.a+self.d*other.c,
            self.c*other.b+self.d*other.d,
            (self.a*other.tx+self.b*other.ty+self.tx)%2,
            (self.c*other.tx+self.d*other.ty+self.ty)%2
        )

    def projective_key(self) -> tuple[tuple[int,int,int,int],int,int]:
        x=(self.a,self.b,self.c,self.d)
        if next(v for v in x if v) <0:
            x=tuple(-v for v in x)
        return x,self.tx,self.ty


@dataclass(frozen=True)
class ArrowedConwaySeamObligation:
    """Source Fig12.7 adjacent arrow gate, plus signed GL2 map data.

    The arrows on the two source-side frames must share their starting
    puncture under the mapping and end along DIFFERENT adjacent
    equatorial arcs. This is a *necessary* source condition, not a
    sufficient 3D gluing or isotopy witness.
    """
    incoming: ArrowedConwayBoundary
    outgoing: ArrowedConwayBoundary
    mapping: SignedConwayBoundaryClass
    status: Literal["ARROWED_SEAM_SOURCE_UNVERIFIED"] = (
        "ARROWED_SEAM_SOURCE_UNVERIFIED"
    )

    def __post_init__(self) -> None:
        if (not isinstance(self.incoming,ArrowedConwayBoundary)
                or not isinstance(self.outgoing,ArrowedConwayBoundary)
                or self.incoming.port==self.outgoing.port
                or not isinstance(self.mapping,SignedConwayBoundaryClass)
                or self.mapping.determinant != -1
                or self.status!="ARROWED_SEAM_SOURCE_UNVERIFIED"):
            raise ValueError("seam needs distinct frames and orientation-reversing map")
        # In source Figure 12.7 the TWO oriented arrows emanate
        # from the same puncture on adjacent, distinct arcs.
        mapped_start=self.mapping.image_puncture(self.incoming.arrow_start)
        mapped_end=self.mapping.image_puncture(self.incoming.arrow_end)
        if (mapped_start != self.outgoing.arrow_start
                or mapped_end == self.outgoing.arrow_end):
            raise ValueError("source arrows must emanate from one point along distinct arcs")
        idx=self.outgoing.cyclic.index(mapped_start)
        target_idx=self.outgoing.cyclic.index(mapped_end)
        if (idx-target_idx)%4 not in (1,3):
            raise ValueError("mapped arrow not adjacent on equator")

    def source_alignment_record(self) -> tuple:
        return (self.incoming.port,self.outgoing.port,
                self.mapping.puncture_permutation(),
                self.mapping.projective_key(),
                (self.mapping.image_slope(1,0),
                 self.mapping.image_slope(0,1)))

    def require_ambient_isotopy(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: source Fig12.7 arrows, GL2 map and "
            "four punctures do not embed an orientation-reversing "
            "Conway 3D seam or certify source-band isotopy"
        )


def proposed_arrowed_seam(
    source: ArrowedConwayBoundary,
    target: ArrowedConwayBoundary,
    mapping: SignedConwayBoundaryClass,
) -> ArrowedConwaySeamObligation:
    return ArrowedConwaySeamObligation(source,target,mapping)
