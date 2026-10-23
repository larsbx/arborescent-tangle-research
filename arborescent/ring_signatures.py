"""A2.5: necessary local signatures for modified Bonahon--Siebenmann R±.

This module does NOT implement a ring *move*. It finds a deliberately
permissive local necessary signature extracted from the two modified-ring
diagrams at §12.4, printed p.155 (PDF page index 166):

R+: a valence-three vertex of total weight +1 with two adjacent
     valence-one vertices of weight +2.
R-: corresponding -1, -2, -2 configuration.

Cyclic sector placement, boundary/marking constraints and the full
pictured local replacement must be checked before any actual R± move
can be authorized. Therefore signature presence never certifies that
a ring move is valid; signature ABSENCE rules out that local pattern.

For F1/F2/strict-F3' (the current partial code), the total integer
weight at each vertex is unchanged. This yields an exact obstruction:
a presentation with no vertices of total weight ±1 cannot acquire
a modified-ring signature via any existing partial moves.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .tangle_labels import V4LabeledTree
from .weighted_trees import WeightedPlaneTree


@dataclass(frozen=True)
class RingSignature:
    """Necessary motif, *not* a proof that R± is available."""
    sign: int
    center: int
    leaf_vertices: tuple[int, int]
    stem_bond: str

    def __post_init__(self) -> None:
        if self.sign not in (-1, 1) or type(self.sign) is not int:
            raise ValueError("sign must be ±1")
        if self.center < 0 or self.leaf_vertices[0] == self.leaf_vertices[1]:
            raise ValueError("invalid ring signature")


def _tree(obj: WeightedPlaneTree | V4LabeledTree) -> WeightedPlaneTree:
    if isinstance(obj, V4LabeledTree):
        return obj.tree
    if isinstance(obj, WeightedPlaneTree):
        return obj
    raise ValueError("expected weighted plane tree or V4-labelled presentation")


def vertex_weights(obj: WeightedPlaneTree | V4LabeledTree) -> tuple[int, ...]:
    """Vertex-indexed total weights, preserved by implemented flips."""
    return tuple(v.weight for v in _tree(obj).vertices)


def ring_signatures(obj: WeightedPlaneTree | V4LabeledTree) -> tuple[RingSignature, ...]:
    """Enumerate candidate *necessary* ±(1;2,2) local configurations.

    Finds a valence-three center of total weight sign*1, with at least
    two internal incident edges ending at valence-one tips of sign*2.
    The remaining third incident slot is an uninterpreted stem bond.

    This is intentionally a SUPERSET of actual R± applicability,
    because the source diagram's oriented angular weights, branch
    attachment, and boundary-label hypotheses are not yet encoded.
    """
    tree = _tree(obj)
    owner: dict[str, list[int]] = {}
    for i, v in enumerate(tree.vertices):
        for slot in v.slots:
            owner.setdefault(slot, []).append(i)
    matches: set[RingSignature] = set()
    for i, v in enumerate(tree.vertices):
        sign = v.weight
        if len(v.slots) != 3 or sign not in (-1, 1):
            continue
        eligible: list[tuple[str, int]] = []
        for bond in v.slots:
            if not bond.startswith("e"):
                continue
            j = next(k for k in owner[bond] if k != i)
            leaf = tree.vertices[j]
            if len(leaf.slots) == 1 and leaf.weight == 2 * sign:
                eligible.append((bond, j))
        for (bond_a, leaf_a), (bond_b, leaf_b) in combinations(eligible, 2):
            remain = [b for b in v.slots if b not in (bond_a, bond_b)]
            if len(remain) != 1:
                raise AssertionError("valence-three stem not unique")
            matches.add(RingSignature(sign, i, tuple(sorted((leaf_a, leaf_b))), remain[0]))
    return tuple(sorted(matches, key=lambda z: (z.center, z.sign, z.leaf_vertices, z.stem_bond)))


def ring_absence_by_unit_weight(obj: WeightedPlaneTree | V4LabeledTree) -> bool:
    """Certificate: if True, *no* R± necessary signature is possible.

    This condition is preserved by the currently implemented F1/F2 and
    strict F3' moves, all of which conserve individual vertex totals.
    """
    return all(abs(w) != 1 for w in vertex_weights(obj))
