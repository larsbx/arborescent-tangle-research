"""A5.2a: a *real*, but very small, four-ended braid-diagram kernel.

Unlike A5.1's abstract port/seam graph this module constructs genuine
two-strand braid diagrams as planar crossing projections with an
over/under assignment at every crossing. The north-west/north-east
inputs run monotonically to the south-west/south-east outputs.

A cancellation of adjacent sigma and sigma^{-1} is an explicit
boundary-relative Reidemeister-II certificate for THIS diagram class.
Every two-strand braid word reduces to sigma^w or sigma^{-w}; the
number w is a complete invariant ONLY of two-strand BRAID ISOTOPY,
not of all tangles in a 3-ball. Do not infer a general nonisotopy
statement from different braid exponents without further theorems.

The V4 action on four ordered punctures is a chosen, typed
double-transposition convention. It must be checked independently
against the actual Conway 4-punctured-sphere parametrizations in
Bonahon--Siebenmann §§12.5/12.7; it is NOT yet a verified seam map.
This kernel does NOT realize §12.5 Fig.12.29 atomic tangles with rings
and CANNOT establish soundness of any modified R± move.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .atomic_plumbing_contract import AtomicPlumbingPlan, Port

# Clockwise around a rectangular 4-puncture boundary.
PUNCTURES: tuple[str, ...] = ("NW", "NE", "SE", "SW")
# XI flips left/right, ETA flips top/bottom, ZETA does both.
# These specific geometric identifications are DECLARED conventions,
# not a source-certified interpretation of the labelled tangle flips.
PERMUTATIONS: tuple[tuple[int, int, int, int], ...] = (
    (0, 1, 2, 3),  # identity
    (1, 0, 3, 2),  # XI
    (3, 2, 1, 0),  # ETA
    (2, 3, 0, 1),  # ZETA
)


@dataclass(frozen=True)
class FourPunctureFrame:
    """A single abstract four-punctured boundary with named punctures.

    The frame is *not* a 3D sphere or a Conway gluing map.
    """
    port: Port
    cyclic_punctures: tuple[str, str, str, str] = PUNCTURES

    def __post_init__(self) -> None:
        if not isinstance(self.port, Port):
            raise ValueError("frame requires an A5.1 atomic port")
        if self.cyclic_punctures != PUNCTURES:
            raise ValueError("only declared clockwise NW/NE/SE/SW frame supported")


@dataclass(frozen=True)
class PuncturePermutation:
    """One V4 double-transposition as an *explicit puncture bijection*."""
    v4: int

    def __post_init__(self) -> None:
        if type(self.v4) is not int or self.v4 not in range(4):
            raise ValueError("puncture action must be an exact V4 element")

    @property
    def image_indices(self) -> tuple[int, int, int, int]:
        return PERMUTATIONS[self.v4]

    def apply(self, puncture: str) -> str:
        if puncture not in PUNCTURES:
            raise ValueError("unknown puncture")
        return PUNCTURES[self.image_indices[PUNCTURES.index(puncture)]]

    def compose(self, other: "PuncturePermutation") -> "PuncturePermutation":
        if not isinstance(other, PuncturePermutation):
            raise ValueError("both maps must be typed puncture permutations")
        composed = tuple(
            self.image_indices[other.image_indices[i]] for i in range(4)
        )
        if composed not in PERMUTATIONS:
            raise AssertionError("V4 puncture maps not closed under composition")
        return PuncturePermutation(PERMUTATIONS.index(composed))

    def inverse(self) -> "PuncturePermutation":
        return self  # the selected V4 actions are involutions


@dataclass(frozen=True)
class AbstractPunctureMatching:
    """Four named boundary puncture matches, NOT a validated 3D seam glue."""
    parent: FourPunctureFrame
    child: FourPunctureFrame
    mapping: PuncturePermutation
    status: Literal["ABSTRACT_ONLY"] = "ABSTRACT_ONLY"

    def __post_init__(self) -> None:
        if (not isinstance(self.parent, FourPunctureFrame)
                or not isinstance(self.child, FourPunctureFrame)
                or not isinstance(self.mapping, PuncturePermutation)
                or self.status != "ABSTRACT_ONLY"
                or self.parent.port == self.child.port):
            raise ValueError("puncture matching must remain an abstract bijection")

    def pairs(self) -> tuple[tuple[str, str], ...]:
        """Exact 4-puncture bijection, no strand seaming asserted."""
        return tuple((s, self.mapping.apply(s)) for s in PUNCTURES)

    def certify_seam_isotopy(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: no Conway seam gluing, puncture "
            "framing or ambient isotopy verifier is implemented"
        )


def frames_for_plan(plan: AtomicPlumbingPlan) -> tuple[FourPunctureFrame, ...]:
    """Expose one declared four-puncture frame for every A5.1 seam/free port."""
    if not isinstance(plan, AtomicPlumbingPlan):
        raise ValueError("A5.1 typed plumbing plan required")
    return tuple(
        FourPunctureFrame(Port(atom.address, slot))
        for atom in plan.atoms for slot in range(atom.valence)
    )


@dataclass(frozen=True)
class Crossing:
    """Projection crossing with exact over/under input-strand identity.

    'left_at_top'/'right_at_top' are strand labels originating at
    NW and NE; sign +1 means the left strand passes *over* right.
    """
    level: int
    sign: int
    left_at_top: str
    right_at_top: str
    over: str
    under: str

    def __post_init__(self) -> None:
        if (type(self.level) is not int or self.level < 0
                or type(self.sign) is not int or self.sign not in (-1, 1)
                or {self.left_at_top, self.right_at_top} != {"NW", "NE"}
                or self.left_at_top == self.right_at_top
                or {self.over, self.under} != {"NW", "NE"}
                or self.over == self.under):
            raise ValueError("invalid oriented two-braid crossing")
        if self.over != (self.left_at_top if self.sign == 1
                         else self.right_at_top):
            raise ValueError("overpass disagrees with signed braid convention")


@dataclass(frozen=True)
class Braid2Diagram:
    """A planar projection word of explicit signed two-strand crossings."""
    generators: tuple[int, ...]
    boundary: tuple[str, str, str, str] = ("NW", "NE", "SW", "SE")

    def __post_init__(self) -> None:
        if (not isinstance(self.generators, tuple)
                or any(type(a) is not int or a not in (-1, 1)
                       for a in self.generators)
                or self.boundary != ("NW", "NE", "SW", "SE")):
            raise ValueError("braid generators must be an exact tuple of ±1")

    @property
    def exponent(self) -> int:
        """Sum of crossing signs; B2 braid-isotopy classification ONLY."""
        return sum(self.generators)

    @property
    def crossing_count(self) -> int:
        return len(self.generators)

    def crossing_events(self) -> tuple[Crossing, ...]:
        """Compute actual left/right and over/under for every crossing."""
        left, right = "NW", "NE"
        events = []
        for i, sign in enumerate(self.generators):
            over = left if sign == 1 else right
            under = right if sign == 1 else left
            events.append(Crossing(i, sign, left, right, over, under))
            left, right = right, left
        return tuple(events)

    def endpoint_matching(self) -> tuple[tuple[str, str], tuple[str, str]]:
        """Two disjoint physical strand routes between the four endpoints."""
        if self.crossing_count % 2:
            return (("NW", "SE"), ("NE", "SW"))
        return (("NW", "SW"), ("NE", "SE"))

    def projection_paths(self) -> tuple[tuple[tuple[float, float], ...], ...]:
        """Exact crossing-level polygonal paths (over/under in crossing_events).

        x is 0 for left, 2 for right, y grows from top to bottom.
        Crossings occur at (1, i+0.5) for level i. Both planar paths
        share that projected point; their Z-order is the signed
        over/under metadata, not encoded by fake planar separation.
        """
        paths = {"NW": [(0.0, 0.0)], "NE": [(2.0, 0.0)]}
        left, right = "NW", "NE"
        for i in range(self.crossing_count):
            paths[left].extend(((1.0, i + .5), (2.0, i + 1.0)))
            paths[right].extend(((1.0, i + .5), (0.0, i + 1.0)))
            left, right = right, left
        # Extend each physical strand to its actual SOUTH boundary
        # puncture. Without this segment a zero-crossing picture would
        # contain only the two northern points and no full boundary arcs.
        y_end = float(self.crossing_count + 1)
        paths[left].append((0.0, y_end))
        paths[right].append((2.0, y_end))
        return tuple(tuple(paths[name]) for name in ("NW", "NE"))

    def render_svg(self) -> str:
        """Small deterministic 2D projection. The overpass is drawn last.

        This is a bona fide projection of a two-strand braid, NOT the
        Bonahon--Siebenmann atomic-ring tangle of Fig.12.29.
        """
        height = 60*(self.crossing_count+1)
        width = 160
        paths = self.projection_paths()
        # Draw through-strands with a white crossing gap on underpasses.
        # For inspection we also emit signed crossing metadata in SVG.
        pieces = [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            'role="img" aria-label="two-strand braid diagram">',
        ]
        for idx, route in enumerate(paths):
            pts = " ".join(f"{40+40*x:.1f},{20+60*y:.1f}" for x,y in route)
            pieces.append(f'<polyline points="{pts}" fill="none" '
                          f'stroke="{"#1e40af" if idx == 0 else "#b91c1c"}" '
                          'stroke-width="3"/>')
        for event in self.crossing_events():
            yy = 20+60*(event.level+.5)
            # Four-pixel under-gap at the exact crossing, followed by
            # the straight overstrand segment in the same projection.
            pieces.append(f'<circle cx="80" cy="{yy:.1f}" r="4" fill="white"/>')
            over_left = event.over == event.left_at_top
            if over_left:
                a,b = (70, yy-15), (90, yy+15)
            else:
                a,b = (90, yy-15), (70, yy+15)
            pieces.append(
                f'<path d="M {a[0]} {a[1]:.1f} L {b[0]} {b[1]:.1f}" '
                'stroke="#111827" stroke-width="3.5" fill="none"/>'
            )
            pieces.append(
                f'<desc>level={event.level} sign={event.sign} '
                f'over={event.over} under={event.under}</desc>'
            )
        pieces.append('</svg>')
        return "".join(pieces)


@dataclass(frozen=True)
class RIIEdit:
    """A single checked adjacent inverse pair deletion/insertion."""
    kind: Literal["delete", "insert"]
    index: int
    first_sign: int

    def __post_init__(self) -> None:
        if (self.kind not in ("delete", "insert")
                or type(self.index) is not int or self.index < 0
                or type(self.first_sign) is not int
                or self.first_sign not in (-1, 1)):
            raise ValueError("invalid typed Reidemeister-II braid edit")

    def apply(self, diagram: Braid2Diagram) -> Braid2Diagram:
        if not isinstance(diagram, Braid2Diagram):
            raise ValueError("braid diagram required")
        a = diagram.generators
        if self.kind == "delete":
            if (self.index + 1 >= len(a)
                    or a[self.index:self.index+2]
                    != (self.first_sign, -self.first_sign)):
                raise ValueError("RII deletion requires adjacent inverse crossings")
            return Braid2Diagram(a[:self.index] + a[self.index+2:])
        if self.index > len(a):
            raise ValueError("RII insertion index exceeds braid length")
        return Braid2Diagram(
            a[:self.index]
            + (self.first_sign, -self.first_sign)
            + a[self.index:]
        )

    def inverse(self) -> "RIIEdit":
        return RIIEdit(
            "insert" if self.kind == "delete" else "delete",
            self.index, self.first_sign,
        )


@dataclass(frozen=True)
class Braid2IsotopyCertificate:
    """Constructive *relative-boundary two-braid* RII-equivalence witness."""
    start: Braid2Diagram
    target: Braid2Diagram
    edits: tuple[RIIEdit, ...]

    def verify(self) -> bool:
        if not isinstance(self.start, Braid2Diagram) or not isinstance(self.target, Braid2Diagram):
            return False
        state = self.start
        for edit in self.edits:
            if not isinstance(edit, RIIEdit):
                return False
            child = edit.apply(state)
            if (child.exponent != state.exponent
                    or child.endpoint_matching() != state.endpoint_matching()):
                return False
            state = child
        return state == self.target


def reduce_two_braid(diagram: Braid2Diagram) -> Braid2IsotopyCertificate:
    """Independent RII-only stack reduction to the exact B2 normal word."""
    if not isinstance(diagram, Braid2Diagram):
        raise ValueError("two-strand braid required")
    current = diagram
    edits: list[RIIEdit] = []
    i = 0
    while i < len(current.generators)-1:
        a = current.generators
        if a[i] == -a[i+1]:
            step = RIIEdit("delete", i, a[i])
            current = step.apply(current)
            edits.append(step)
            i = max(0, i-1)
        else:
            i += 1
    expected = ((1,)*diagram.exponent if diagram.exponent >= 0
                else (-1,)*(-diagram.exponent))
    if current.generators != expected:
        raise AssertionError("RII stack reduction did not reach braid normal form")
    witness = Braid2IsotopyCertificate(diagram, current, tuple(edits))
    if not witness.verify():
        raise AssertionError("braid reduction certificate did not replay")
    return witness


def prove_two_braid_equivalence(
    source: Braid2Diagram, target: Braid2Diagram
) -> Braid2IsotopyCertificate:
    """Build an *actual positive* RII diagram witness when exponent agrees.

    If exponent differs, reject the request as not proved by this
    braid-isotopy method. Do not claim general tangle nonisotopy.
    """
    if not isinstance(source, Braid2Diagram) or not isinstance(target, Braid2Diagram):
        raise ValueError("two braid diagrams required")
    if source.exponent != target.exponent:
        raise ValueError("no braid-isotopy proof: exponents differ")
    a, b = reduce_two_braid(source), reduce_two_braid(target)
    edits = a.edits + tuple(edit.inverse() for edit in reversed(b.edits))
    cert = Braid2IsotopyCertificate(source, target, edits)
    if not cert.verify():
        raise AssertionError("constructed relative-boundary RII witness is invalid")
    return cert
