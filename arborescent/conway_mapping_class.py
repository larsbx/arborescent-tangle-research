"""A5.2b: why four Conway puncture labels do NOT determine a sewing map.

A four-punctured Conway sphere S^2\{4 points} has an INFINITE pure
mapping-class group (despite only four marked punctures and finite V4
label permutations). The branched double cover by a torus supplies
a concrete algebraic family: matrices in SL(2,Z), congruent to I mod 2,
fix the four branch points but can twist primitive homology slopes.

For every n, T_n=[[1,2n],[0,1]] is a distinct pure mapping class
modulo the central +/-I ambiguity. All have the SAME identity action
on the four marked punctures, including the V4 label.

Therefore A5.2a's explicit four-puncture bijection is only an
ABSTRACT skeleton for a Conway seam: no V4 label or puncture
permutation can encode the required mapping-class / twisting data.

This is an algebraic obstruction, NOT an implementation of
Bonahon--Siebenmann Fig.12.29 or an ambient tangle isotopy verifier.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import gcd

from .atomic_plumbing_contract import AtomicPlumbingPlan, Port
from .conway_braid_diagrams import (
    AbstractPunctureMatching, proposed_seam_puncture_matching,
)


def _primitive_slope(p: int, q: int) -> tuple[int, int]:
    if (type(p) is not int or type(q) is not int
            or (p == 0 and q == 0)
            or gcd(p, q) != 1):
        raise ValueError("slope must be an exact primitive nonzero integer pair")
    if p < 0 or (p == 0 and q < 0):
        return -p, -q
    return p, q


@dataclass(frozen=True)
class PureConwayMappingClass:
    """An SL(2,Z) lift of a puncture-fixing 4-punctured-sphere class.

    The projective +/-I ambiguity is normalized by choosing the sign
    of the first nonzero matrix entry in row-major order. We retain
    ONLY the pure congruence subgroup Gamma(2): M == I (mod 2).

    This is an *orientation-preserving self-map class* of a marked
    4-punctured sphere (via its torus double cover). It does not
    specify the orientation-reversing gluing of two boundary spheres.
    """
    a: int
    b: int
    c: int
    d: int

    def __post_init__(self) -> None:
        if any(type(x) is not int for x in (self.a, self.b, self.c, self.d)):
            raise ValueError("mapping class entries must be exact integers")
        if self.a * self.d - self.b * self.c != 1:
            raise ValueError("mapping class lift must have determinant +1")
        if (self.a % 2, self.b % 2, self.c % 2, self.d % 2) != (1, 0, 0, 1):
            raise ValueError("pure Conway subgroup must act trivially mod 2")

    @property
    def matrix(self) -> tuple[tuple[int, int], tuple[int, int]]:
        return (self.a, self.b), (self.c, self.d)

    def projective_key(self) -> tuple[int, int, int, int]:
        """Identify M and -M, the same sphere class under deck involution."""
        x = (self.a, self.b, self.c, self.d)
        return x if next(k for k in x if k != 0) > 0 else tuple(-k for k in x)

    def inverse(self) -> "PureConwayMappingClass":
        return PureConwayMappingClass(self.d, -self.b, -self.c, self.a)

    def compose(self, other: "PureConwayMappingClass") -> "PureConwayMappingClass":
        """Group product M * N, acting on column vectors by N then M."""
        if not isinstance(other, PureConwayMappingClass):
            raise ValueError("typed pure Conway mapping class required")
        return PureConwayMappingClass(
            self.a*other.a + self.b*other.c,
            self.a*other.b + self.b*other.d,
            self.c*other.a + self.d*other.c,
            self.c*other.b + self.d*other.d,
        )

    def image_slope(self, p: int, q: int) -> tuple[int, int]:
        """Projective primitive slope image of (p,q) in the torus cover."""
        p, q = _primitive_slope(p, q)
        return _primitive_slope(
            self.a*p + self.b*q,
            self.c*p + self.d*q,
        )

    def branch_point_action(
        self,
    ) -> tuple[tuple[int, int], tuple[int, int],
               tuple[int, int], tuple[int, int]]:
        """Act explicitly on the four 2-torsion classes of the torus cover.

        These are the four double-cover branch points. Every Gamma(2)
        matrix MUST fix each point, independent of our chosen V4
        coordinate convention for the Conway boundary.
        """
        return tuple(
            ((self.a*x + self.b*y) % 2,
             (self.c*x + self.d*y) % 2)
            for x, y in ((0, 0), (1, 0), (0, 1), (1, 1))
        )

    def fixes_four_punctures(self) -> bool:
        """Independently compute action on the double-cover branch points."""
        return self.branch_point_action() == (
            (0, 0), (1, 0), (0, 1), (1, 1)
        )

    def require_verified_conway_gluing(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: a pure sphere mapping class is not "
            "an orientation-reversing Conway seam gluing or "
            "boundary-parametrized arborescent tangle isotope"
        )


def pure_twist(power: int) -> PureConwayMappingClass:
    """Pure twist family [[1,2n],[0,1]], with n any signed integer."""
    if type(power) is not int:
        raise ValueError("twist power must be an exact integer")
    return PureConwayMappingClass(1, 2*power, 0, 1)


def algebraic_slope_det(
    first: tuple[int, int], second: tuple[int, int]
) -> int:
    """Absolute SL(2,Z)-invariant determinant of primitive homology slopes.

    Not by itself an isotopy oracle for any 3D tangle.
    """
    if (not isinstance(first, tuple) or len(first) != 2
            or not isinstance(second, tuple) or len(second) != 2):
        raise ValueError("two typed primitive slope pairs required")
    p, q = _primitive_slope(*first)
    r, s = _primitive_slope(*second)
    return abs(p*s - q*r)


# Coordinates for A5.2a's clockwise puncture enumeration:
# NW=(0,0), NE=(1,0), SE=(1,1), SW=(0,1) in torus 2-torsion.
# An affine half-translation by (tx/2,ty/2) acts as V4 on these points.
BRANCH_COORDINATES = ((0, 0), (1, 0), (1, 1), (0, 1))


@dataclass(frozen=True)
class AffineConwayBoundaryClass:
    """A pure Gamma(2) linear part plus a 2-torsion affine translation.

    The map z |-> M z + t/2 commutes with z |-> -z on T² whenever
    t=(tx,ty) has binary entries, because 2(t/2) is integral.
    It therefore descends to the four-branch-point sphere.

    Gamma(2) reduces to identity mod 2, so the induced labelled
    puncture permutation depends ONLY on the half-translation
    and matches the explicitly chosen A5.2a V4 action.

    This is still an ORIENTATION-PRESERVING boundary self-map class,
    NOT an oriented Conway gluing between two sphere boundaries.
    """
    pure: PureConwayMappingClass
    v4: int

    def __post_init__(self) -> None:
        if not isinstance(self.pure, PureConwayMappingClass):
            raise ValueError("exact pure matrix part required")
        if type(self.v4) is not int or self.v4 not in range(4):
            raise ValueError("exact V4 affine translation required")

    @property
    def translation(self) -> tuple[int, int]:
        return self.v4 & 1, (self.v4 >> 1) & 1

    def puncture_permutation(self) -> tuple[int, int, int, int]:
        """Compute branch-point action from affine torus data, not a lookup."""
        tx, ty = self.translation
        images = []
        for x, y in BRANCH_COORDINATES:
            mapped = (
                (self.pure.a*x+self.pure.b*y+tx) % 2,
                (self.pure.c*x+self.pure.d*y+ty) % 2,
            )
            images.append(BRANCH_COORDINATES.index(mapped))
        return tuple(images)

    def projective_key(self) -> tuple[tuple[int, int, int, int], int]:
        return self.pure.projective_key(), self.v4

    def image_slope(self, p: int, q: int) -> tuple[int, int]:
        """Half-translation does not change torus homology slopes."""
        return self.pure.image_slope(p, q)

    def compose(self, other: "AffineConwayBoundaryClass"
                ) -> "AffineConwayBoundaryClass":
        if not isinstance(other, AffineConwayBoundaryClass):
            raise ValueError("typed affine Conway boundary mapping class required")
        # Since each pure part is identity mod 2, M*u + t = u+t mod 2.
        return AffineConwayBoundaryClass(
            self.pure.compose(other.pure), self.v4 ^ other.v4
        )

    def inverse(self) -> "AffineConwayBoundaryClass":
        return AffineConwayBoundaryClass(self.pure.inverse(), self.v4)


@dataclass(frozen=True)
class ConwaySeamAmbiguityWitness:
    """Two different PURE mapping classes with identical puncture bijection."""
    matching: AbstractPunctureMatching
    first: PureConwayMappingClass
    second: PureConwayMappingClass
    distinguishing_slope: tuple[int, int] = (0, 1)

    def verify(self) -> bool:
        if (not isinstance(self.matching, AbstractPunctureMatching)
                or not isinstance(self.first, PureConwayMappingClass)
                or not isinstance(self.second, PureConwayMappingClass)):
            return False
        try:
            slope = _primitive_slope(*self.distinguishing_slope)
        except (ValueError, TypeError):
            return False
        return (
            self.first.fixes_four_punctures()
            and self.second.fixes_four_punctures()
            and self.first.projective_key() != self.second.projective_key()
            and self.first.image_slope(*slope) != self.second.image_slope(*slope)
            and len(set(b for _, b in self.matching.pairs())) == 4
            and self.matching.status == "ABSTRACT_ONLY"
        )


def seam_ambiguity_witness(
    plan: AtomicPlumbingPlan, parent: Port, child: Port,
    v4: int, twist_power: int = 1
) -> ConwaySeamAmbiguityWitness:
    """Give exact counterexample to seam uniqueness from 4-point + V4 data.

    The chosen V4 puncture matching is held FIXED. Composing it with
    a pure map T_n would not alter any individually marked puncture,
    yet changes the slope of the nontrivial simple closed curve
    represented by (0,1), for n != 0. We don't actually glue spheres.
    """
    if type(twist_power) is not int or twist_power == 0:
        raise ValueError("nonzero exact integer twist power required")
    matching = proposed_seam_puncture_matching(plan, parent, child, v4)
    witness = ConwaySeamAmbiguityWitness(
        matching, pure_twist(0), pure_twist(twist_power)
    )
    if not witness.verify():
        raise AssertionError("pure Conway seam ambiguity failed verification")
    return witness
