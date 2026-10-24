"""A5.2c: exact PL embeddings of the *local two-strand twist blocks*.

Bonahon--Siebenmann Fig. 12.8 (printed p. 142) and Fig. 12.29
(printed p. 158) include signed integral half-twist sectors a_i.
This module realizes those individual braid boxes in Q^3, not the
entire atomic ring tangle, its punctured 3-manifold, or its seams.

y increases north-to-south, x is the projected horizontal coordinate,
and z is height toward the observer. Four endpoints are held fixed.
Only the precisely audited PL normal form is supported.

IMPORTANT: sector blocks alone DO NOT realize Fig. 12.29. The long
return strands, Conway-sphere holes, twisted bands, ring embeddings
and orientation-reversing sewing maps remain missing. No source
R+/R- topology is certified by this local-geometry module.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction
from .conway_braid_diagrams import Braid2Diagram, Braid2IsotopyCertificate, RIIEdit


def _exact(value: object) -> bool:
    return type(value) is int or type(value) is Fraction


@dataclass(frozen=True)
class Point3:
    x: int | Fraction
    y: int | Fraction
    z: int | Fraction

    def __post_init__(self) -> None:
        if not all(_exact(v) for v in (self.x, self.y, self.z)):
            raise ValueError("PL coordinates must be exact rational numbers")

    def serial(self) -> tuple[str, str, str]:
        return tuple(str(Fraction(v)) for v in (self.x, self.y, self.z))


@dataclass(frozen=True)
class SpatialTwoBraid:
    """Two disjoint, explicitly embedded polygonal arcs in exact Q^3.

    The prescribed PL normal form has one projected transverse crossing
    per generator and no others. The audit checks actual strand points,
    including height order, rather than trusting a claimed diagram.
    """
    word: tuple[int, ...]
    nw: tuple[Point3, ...]
    ne: tuple[Point3, ...]
    geometry_status: Literal["VERIFIED_LOCAL_PL_BRAID"] = "VERIFIED_LOCAL_PL_BRAID"

    def __post_init__(self) -> None:
        if (not isinstance(self.word, tuple)
                or any(type(a) is not int or a not in (-1, 1) for a in self.word)
                or self.geometry_status != "VERIFIED_LOCAL_PL_BRAID"):
            raise ValueError("invalid local braid word or geometry status")
        self.audit()

    def audit(self) -> None:
        """Verify the complete PL normal form, crossing heights and endpoints.

        In each crossing cell, one strand goes x=-2 to x=+2
        and the other x=+2 to x=-2. Their only projected collision
        is at x=0, y=2i+1, where their z values are +1 and -1.
        Between these levels, x values are distinct by linearity,
        and terminal vertical segments are also disjoint.
        """
        if (not isinstance(self.nw, tuple) or not isinstance(self.ne, tuple)
                or not all(isinstance(p, Point3) for p in (*self.nw, *self.ne))):
            raise ValueError("strands must be exact polygonal point tuples")
        n = len(self.word)
        if len(self.nw) != 2*n + 2 or len(self.ne) != 2*n + 2:
            raise ValueError("wrong number of polygonal vertices")
        tracks = {"NW": self.nw, "NE": self.ne}
        if (self.nw[0] != Point3(-2, 0, 0)
                or self.ne[0] != Point3(2, 0, 0)):
            raise ValueError("northern marked boundary points changed")
        left, right = "NW", "NE"
        for i, sign in enumerate(self.word):
            a, b = tracks[left], tracks[right]
            if (a[2*i] != Point3(-2, 2*i, 0)
                    or b[2*i] != Point3(2, 2*i, 0)):
                raise ValueError("cell must enter with separated framed endpoints")
            if (a[2*i+1] != Point3(0, 2*i+1, sign)
                    or b[2*i+1] != Point3(0, 2*i+1, -sign)):
                raise ValueError("projected crossing or overpass is wrong")
            if (a[2*i+2] != Point3(2, 2*i+2, 0)
                    or b[2*i+2] != Point3(-2, 2*i+2, 0)):
                raise ValueError("strand endpoints do not swap across crossing")
            left, right = right, left
        if (tracks[left][-2:] != (Point3(-2, 2*n, 0), Point3(-2, 2*n+1, 0))
                or tracks[right][-2:] != (Point3(2, 2*n, 0), Point3(2, 2*n+1, 0))):
            raise ValueError("southern boundary/framing altered")

    @property
    def crossing_count(self) -> int:
        return len(self.word)

    @property
    def writhe(self) -> int:
        return sum(self.word)

    def endpoint_matching(self) -> tuple[tuple[str, str], tuple[str, str]]:
        def south(path: tuple[Point3, ...]) -> str:
            return "SW" if path[-1].x == -2 else "SE"
        return (("NW", south(self.nw)), ("NE", south(self.ne)))

    def crossing_records(self) -> tuple[tuple[int, str, str], ...]:
        """Extract height order from actual 3D coordinates."""
        records = []
        for i in range(len(self.word)):
            p, q = self.nw[2*i+1], self.ne[2*i+1]
            if p.z == q.z:
                raise ValueError("actual over/under separation missing")
            records.append((i, "NW" if p.z > q.z else "NE",
                            "NE" if p.z > q.z else "NW"))
        return tuple(records)

    def canonical_coordinates(self) -> tuple:
        return (tuple(self.word),
                tuple(p.serial() for p in self.nw),
                tuple(p.serial() for p in self.ne))

    def require_source_atomic_realization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: a local braid is not a Fig.12.29 "
            "ringed atomic tangle with Conway boundary spheres"
        )


def spatial_two_braid(diagram: Braid2Diagram) -> SpatialTwoBraid:
    """Construct two polygonal arcs in 3-space for a typed two-braid."""
    if not isinstance(diagram, Braid2Diagram):
        raise ValueError("typed 2-strand braid diagram required")
    points = {"NW": [Point3(-2, 0, 0)], "NE": [Point3(2, 0, 0)]}
    left, right = "NW", "NE"
    for i, sign in enumerate(diagram.generators):
        points[left].extend((Point3(0, 2*i+1, sign), Point3(2, 2*i+2, 0)))
        points[right].extend((Point3(0, 2*i+1, -sign), Point3(-2, 2*i+2, 0)))
        left, right = right, left
    n = len(diagram.generators)
    points[left].append(Point3(-2, 2*n+1, 0))
    points[right].append(Point3(2, 2*n+1, 0))
    return SpatialTwoBraid(diagram.generators,
                           tuple(points["NW"]), tuple(points["NE"]))


@dataclass(frozen=True)
class AtomicSectorBlock:
    """Only one sector's Fig.12.29 a_i half-twist box, NOT a whole tile."""
    atom_address: tuple[int, ...]
    sector_index: int
    signed_half_twists: int
    geometry: SpatialTwoBraid
    source_scope: Literal["LOCAL_SECTOR_ONLY"] = "LOCAL_SECTOR_ONLY"

    def __post_init__(self) -> None:
        if (not isinstance(self.atom_address, tuple)
                or type(self.sector_index) is not int or self.sector_index < 0
                or type(self.signed_half_twists) is not int
                or not isinstance(self.geometry, SpatialTwoBraid)
                or self.source_scope != "LOCAL_SECTOR_ONLY"):
            raise ValueError("invalid source-sector block")
        expected = ((1,) * self.signed_half_twists if self.signed_half_twists >= 0
                    else (-1,) * (-self.signed_half_twists))
        if self.geometry.word != expected:
            raise ValueError("geometric half-twists disagree with sector weight")


def realize_unringed_sector_blocks(atom: AtomicInstruction) -> tuple[AtomicSectorBlock, ...]:
    """Construct each signed a_i box; do not claim the whole atomic tangle.

    Reject ring-bearing data rather than silently eliding its geometry.
    The source's valence-one closed atoms also need a separate model.
    """
    if not isinstance(atom, AtomicInstruction):
        raise ValueError("A5.1 AtomicInstruction required")
    if atom.ring_count:
        raise NotImplementedError("UNVERIFIED_GEOMETRY: r>0 Fig.12.29 rings not embedded")
    if atom.valence < 2:
        raise NotImplementedError("UNVERIFIED_GEOMETRY: closed/low-valence atom unsupported")
    return tuple(
        AtomicSectorBlock(
            atom.address, i, w,
            spatial_two_braid(Braid2Diagram((1,)*w if w >= 0 else (-1,)*(-w))),
        )
        for i, w in enumerate(atom.sector_weights)
    )


def verify_spatial_rii_certificate(cert: Braid2IsotopyCertificate) -> bool:
    """Replay braid RII and check actual Q^3 embeddings at every step.

    Known RII is a relative-boundary two-braid isotopy. This checker
    does NOT establish any ambient isotopy of ringed Conway tangles.
    """
    if not isinstance(cert, Braid2IsotopyCertificate):
        return False
    try:
        if (not isinstance(cert.start, Braid2Diagram)
                or not isinstance(cert.target, Braid2Diagram)):
            return False
        state = cert.start
        prev = spatial_two_braid(state)
        for edit in cert.edits:
            if not isinstance(edit, RIIEdit):
                return False
            child = edit.apply(state)
            nxt = spatial_two_braid(child)
            if (nxt.writhe != prev.writhe
                    or nxt.endpoint_matching() != prev.endpoint_matching()
                    or abs(nxt.crossing_count - prev.crossing_count) != 2):
                return False
            for d in (prev, nxt):
                if tuple((e.level, e.over, e.under)
                         for e in Braid2Diagram(d.word).crossing_events()) != d.crossing_records():
                    return False
            state, prev = child, nxt
        return state == cert.target
    except (ValueError, TypeError, AttributeError):
        return False
