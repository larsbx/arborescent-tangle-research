"""A4.2: terminating, replayable normal form for *formal* ring slides.

Restricted to A4.1's signed, disjoint ±2--±1--±2 caps at all-zero-sector
bases, no cap anchored at slot zero and no remaining unrecognized
ring-signature candidates. This is a canonicalization of one *specified
presentation rewrite system*, NOT canonical tangle isotopy.

Every allowed local slide swaps R O -> O R or O R -> R O at one base,
where R denotes a recognized cap stem and O an ordinary bond. The
stable normalizer moves all R symbols right of all O symbols while
preserving the relative order of ordinary bonds and indexed cap stems.

The inversion count # {(R,O): R occurs before O} decreases by exactly
one per recorded step. Confluence on the restricted class follows
because all words of fixed ordered O/R subsequences have a unique
O...O R...R normal form. No completeness claim outside this class.
"""
from __future__ import annotations

from dataclasses import dataclass

from .ring_abbreviation import strict_ring_caps
from .ring_signatures import ring_signatures
from .ring_relocation import RingSlide, apply_ring_slide
from .tangle_labels import V4LabeledTree


def _supported(source: V4LabeledTree, sign: int):
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled presentation required")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("sign must be ±1")
    candidates = ring_signatures(source)
    if any(candidate.sign != sign for candidate in candidates):
        raise ValueError("mixed-sign ring candidates are unsupported")
    caps = strict_ring_caps(source, sign)
    if len(caps) != len(candidates):
        raise ValueError("incomplete strict recognition; refusing normalization")
    center_set = {cap.center for cap in caps}
    if any(cap.base in center_set for cap in caps):
        raise ValueError("adjacent ring centers need exceptional handling")
    groups: dict[int, dict[str, int]] = {}
    for cap in caps:
        base = source.tree.vertices[cap.base]
        if any(base.sectors):
            raise ValueError("nonzero base sector is unsupported")
        if base.slots[0] == cap.stem:
            raise ValueError("ring stem occupies the anchored base slot")
        groups.setdefault(cap.base, {})[cap.stem] = cap.center
    return groups


def ring_position_inversions(source: V4LabeledTree, sign: int) -> int:
    """Nonnegative integer termination measure for the strict ring system."""
    groups = _supported(source, sign)
    inversions = 0
    for base_idx, ring_map in groups.items():
        ordinary_to_right = 0
        for bond in reversed(source.tree.vertices[base_idx].slots[1:]):
            if bond in ring_map:
                inversions += ordinary_to_right
            else:
                ordinary_to_right += 1
    return inversions


@dataclass(frozen=True)
class RingNormalForm:
    """Shortest monotone bubble-sort certificate, indexed by source bonds."""
    start: V4LabeledTree
    normal: V4LabeledTree
    sign: int
    moves: tuple[RingSlide, ...]

    def verify(self) -> bool:
        current = self.start
        measure = ring_position_inversions(current, self.sign)
        for move in self.moves:
            child = apply_ring_slide(current, move)
            next_measure = ring_position_inversions(child, self.sign)
            if next_measure != measure - 1:
                return False
            current, measure = child, next_measure
        return current == self.normal and measure == 0


def normalize_ring_positions(
    source: V4LabeledTree, sign: int, max_moves: int | None = None
) -> RingNormalForm:
    """Move every strict recognized ring to the end of its anchored base.

    The deterministic presentation-level algorithm is a stable sort,
    and each recorded move is an actual typed A4.1 formal slide.
    Does not quotient flips, F3', general R±, or genuine isotopy.
    """
    if (max_moves is not None
            and (type(max_moves) is not int or max_moves < 0)):
        raise ValueError("max_moves must be nonnegative or None")
    group = _supported(source, sign)
    initial_inversions = ring_position_inversions(source, sign)
    if max_moves is not None and initial_inversions > max_moves:
        raise ValueError("move budget is insufficient; no false normal form")

    state = source
    moves: list[RingSlide] = []
    for base_idx in sorted(group):
        stems = group[base_idx]
        while True:
            slots = state.tree.vertices[base_idx].slots
            next_move = None
            # Move only R O -> O R, rightward; cap stems never cross.
            for j in range(len(slots) - 2, 0, -1):
                if slots[j] in stems and slots[j + 1] not in stems:
                    next_move = RingSlide(sign, stems[slots[j]], slots[j+1])
                    break
            if next_move is None:
                break
            state = apply_ring_slide(state, next_move)
            moves.append(next_move)
    certificate = RingNormalForm(source, state, sign, tuple(moves))
    if len(moves) != initial_inversions or not certificate.verify():
        raise AssertionError("normalization failed its inversion certificate")
    return certificate


def strict_ring_normal_key(source: V4LabeledTree, sign: int) -> tuple:
    """A marked rooted *presentation* code modulo supported slides.

    Matching keys are not certificates of topological isotopy or full R±.
    """
    return normalize_ring_positions(source, sign).normal.labelled_rooted_code()
