"""A4.10: Gaussian-binomial distance distribution for A4.1 *formal* ring slides.

The strict formal ring calculus allows adjacent swaps R O <-> O R,
where R is a fully recognized pendant ring stem and O is an ordinary
bond in a base's rooted, nonanchor cyclic incident list. There are
NO R-R or O-O swaps. It is not a 3D isotopy calculus.

For o ordinary attachments and r identical ring caps, the exact
distance-to-stable-form polynomial is the Gaussian binomial:
    G[o,r](q) = G[o,r-1](q) + q**r * G[o-1,r](q)
with G[0,r]=G[o,0]=1.

For a max-valence-four path with k>=2 independently ring-bearing bases,
the full indexed formal-slide orbit rank polynomial is
    (1+q)**2 * (1+q+q**2)**(k-2).

For two rooted, same-index/same-bond presentations in a fixed strict
ring-slide orbit, minimal distance is the sum over indexed ring stems
of the absolute difference between the numbers of *ordinary* bonds
before each ring stem. The constructive geodesic below replays exactly
that number of typed A4.1 steps. Arbitrary graph isomorphisms, flip
moves, F3' and true R± are explicitly NOT treated as ring slides.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from .ring_normalizer import _supported, align_rooted_cyclic_slots
from .ring_relocation import RingSlide, RingTrace, apply_ring_slide
from .tangle_labels import V4LabeledTree


def _validate(o: int, r: int) -> None:
    if (type(o) is not int or o < 0
            or type(r) is not int or r < 0):
        raise ValueError("o and r must be nonnegative integers")


def gaussian_ring_polynomial(o: int, r: int) -> tuple[int, ...]:
    """Exact q-binomial coefficients, validating types BEFORE cache lookup.

    Python bool and int have equal hashes (True == 1, False == 0);
    a public lru_cache decorator could silently return valid integer
    results for explicitly forbidden bool parameters. The cache lives
    only on the private helper after strict public validation.
    """
    _validate(o, r)
    return _cached_gaussian_ring_polynomial(o, r)


@lru_cache(maxsize=None)
def _cached_gaussian_ring_polynomial(o: int, r: int) -> tuple[int, ...]:
    if o == 0 or r == 0:
        return (1,)
    first = _cached_gaussian_ring_polynomial(o, r-1)
    second = (0,) * r + _cached_gaussian_ring_polynomial(o-1, r)
    return tuple(
        (first[i] if i < len(first) else 0)
        + (second[i] if i < len(second) else 0)
        for i in range(max(len(first), len(second)))
    )


def _product(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    coefficients = [0] * (len(a)+len(b)-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            coefficients[i+j] += ai*bj
    return tuple(coefficients)


def path_rank_polynomial(bases: int) -> tuple[int, ...]:
    """Full A4.1 formal distance-to-normal distribution on k-base path."""
    if type(bases) is not int or bases < 2:
        raise ValueError("at least two base vertices required")
    out = _product((1, 1), (1, 1))
    for _ in range(bases-2):
        out = _product(out, (1, 1, 1))
    return out


def _ordinary_before_rings(slots: tuple[str, ...], ring_stems: set[str]
                           ) -> tuple[int, ...]:
    seen_ordinary = 0
    positions = []
    for bond in slots[1:]:
        if bond in ring_stems:
            positions.append(seen_ordinary)
        else:
            seen_ordinary += 1
    return tuple(positions)


def _same_indexed_slide_class(
    source: V4LabeledTree, target: V4LabeledTree, sign: int,
    root: str
) -> tuple[V4LabeledTree, V4LabeledTree, dict[int, dict[str, int]]]:
    """Fail closed unless all differences are eligible R/O cyclic positions."""
    if not isinstance(source, V4LabeledTree) or not isinstance(target, V4LabeledTree):
        raise ValueError("V4-labelled source and target required")
    left = align_rooted_cyclic_slots(source, root)
    right = align_rooted_cyclic_slots(target, root)
    first = _supported(left, sign)
    second = _supported(right, sign)
    if (len(left.tree.vertices) != len(right.tree.vertices)
            or left.labels != right.labels or first != second):
        raise ValueError("presentations do not share indexed signed cap data")
    for i, (u, v) in enumerate(zip(left.tree.vertices, right.tree.vertices)):
        if i not in first:
            if u != v:
                raise ValueError("nonbase weighted presentation data differ")
            continue
        stems = set(first[i])
        if (u.slots[0] != v.slots[0]
                or set(u.slots) != set(v.slots)
                or u.sectors != v.sectors
                or any(u.sectors)):
            raise ValueError("base anchor, bonds or sectors differ")
        if (tuple(s for s in u.slots[1:] if s not in stems)
                != tuple(s for s in v.slots[1:] if s not in stems)
                or tuple(s for s in u.slots[1:] if s in stems)
                != tuple(s for s in v.slots[1:] if s in stems)):
            raise ValueError("ordinary or ring-stem relative order differs")
    return left, right, first


def formal_slide_distance(
    source: V4LabeledTree, target: V4LabeledTree, sign: int,
    root: str = "f0"
) -> int:
    """Exact geodesic distance inside one strictly guarded indexed slide orbit."""
    left, right, groups = _same_indexed_slide_class(source, target, sign, root)
    total = 0
    for i, stems in groups.items():
        first = _ordinary_before_rings(left.tree.vertices[i].slots, set(stems))
        second = _ordinary_before_rings(right.tree.vertices[i].slots, set(stems))
        if len(first) != len(second):
            raise AssertionError("cap multiplicity changed")
        total += sum(abs(a-b) for a, b in zip(first, second))
    return total


@dataclass(frozen=True)
class FormalSlideGeodesic:
    """Rooted cyclic presentation isomorphisms + minimal formal slide trace."""
    source: V4LabeledTree
    target: V4LabeledTree
    aligned_source: V4LabeledTree
    aligned_target: V4LabeledTree
    sign: int
    root: str
    trace: RingTrace
    distance: int

    def verify(self) -> bool:
        if self.aligned_source != align_rooted_cyclic_slots(self.source, self.root):
            return False
        if self.aligned_target != align_rooted_cyclic_slots(self.target, self.root):
            return False
        if self.trace.start != self.aligned_source or self.trace.target != self.aligned_target:
            return False
        if len(self.trace.moves) != self.distance or not self.trace.verify():
            return False
        return formal_slide_distance(self.source, self.target,
                                     self.sign, self.root) == self.distance


def formal_slide_geodesic(
    source: V4LabeledTree, target: V4LabeledTree, sign: int,
    root: str = "f0"
) -> FormalSlideGeodesic:
    """Construct exact shortest sequence of A4.1 typed formal ring slides."""
    left, right, groups = _same_indexed_slide_class(source, target, sign, root)
    expected = formal_slide_distance(left, right, sign, root)
    current = left
    moves: list[RingSlide] = []
    for i in sorted(groups):
        stems = groups[i]
        desired = right.tree.vertices[i].slots
        for j in range(1, len(desired)):
            slots = current.tree.vertices[i].slots
            if slots[j] == desired[j]:
                continue
            # Stable R/R and O/O subsequence preservation ensures that
            # every crossing in this transport is between unlike types.
            position = slots.index(desired[j])
            if position <= j:
                raise AssertionError("inconsistent stable relative order")
            while position > j:
                slots = current.tree.vertices[i].slots
                before, moving = slots[position-1:position+1]
                if before in stems and moving not in stems:
                    move = RingSlide(sign, stems[before], moving)
                elif moving in stems and before not in stems:
                    move = RingSlide(sign, stems[moving], before)
                else:
                    raise AssertionError("requested move crosses same-type bonds")
                current = apply_ring_slide(current, move)
                moves.append(move)
                position -= 1
    if current != right or len(moves) != expected:
        raise AssertionError("constructed ring path is not a minimal geodesic")
    result = FormalSlideGeodesic(
        source, target, left, right, sign, root,
        RingTrace(left, right, tuple(moves)), expected,
    )
    if not result.verify():
        raise AssertionError("failed formal geodesic replay certificate")
    return result
