"""A4.4: F1/F2-equivariance of §12.5 restricted ring abbreviation.

Theorem 12.21 (Bonahon--Siebenmann §12.7.3) uses V4-labelled free
bonds and F1, F2, F3' and modified ring moves on *canonical* trees.
This module addresses ONLY source-defined F1 and F2, which are already
coded in tangle_labels.py.

Crucial scope: F1/F2 on a ring-numbered *abbreviated* representation
reverse the SAME cyclic rotations as on a weighted tree and modify
free-bond labels using the source's V4 law. Ring numbers remain
attached to their original surviving vertices.

For disjoint signed caps in the strict zero-base-sector domain:
    abbreviate(Fi(T)) == Fi(abbreviate(T))  for i in {1,2}
as explicit (indexed) Abbreviated objects when the F2 pivot is identified
by the SAME surviving free bond on both sides.

This is an exact *formal presentation* equivariance test, not an
independent 3D proof of arborescent tangle isotopy.
"""
from __future__ import annotations

from .ring_abbreviation import Abbreviated, abbreviate_disjoint_rings
from .tangle_labels import (
    IDENTITY, V4LabeledTree, tangle_flip_f1, tangle_flip_f2,
)


def _free_owner(labeled: V4LabeledTree, free_bond: str) -> int:
    """Return the unique vertex incident to an explicitly named free bond."""
    if not isinstance(labeled, V4LabeledTree):
        raise ValueError("expected V4-labelled presentation")
    if not isinstance(free_bond, str) or not free_bond.startswith("f"):
        raise ValueError("pivot must be an existing free bond")
    owners = [
        i for i, vertex in enumerate(labeled.tree.vertices)
        if free_bond in vertex.slots
    ]
    if len(owners) != 1:
        raise ValueError("pivot must be a present uniquely owned free bond")
    return owners[0]


def abbreviated_flip_f1(source: Abbreviated) -> Abbreviated:
    """Source §12.7.1 F1: mirror every vertex, XOR ζ on every free label."""
    if not isinstance(source, Abbreviated):
        raise ValueError("expected ring-numbered abbreviated presentation")
    return Abbreviated(
        tangle_flip_f1(source.labeled), source.rings, source.sign
    )


def abbreviated_flip_f2(
    source: Abbreviated, pivot_free_bond: str = "f0"
) -> Abbreviated:
    """Source §12.7.1 F2 pivot identified by a surviving marked free bond.

    F2 parity is computed on the *abbreviated* residual tree. Removing
    terminal ring caps does not change distances between surviving
    vertices, so parity agrees with the full source tree in the strict
    disjoint-cap abbreviation domain.
    """
    if not isinstance(source, Abbreviated):
        raise ValueError("expected ring-numbered abbreviated presentation")
    pivot = _free_owner(source.labeled, pivot_free_bond)
    return Abbreviated(
        tangle_flip_f2(source.labeled, pivot), source.rings, source.sign
    )


def abbreviated_flip_group(
    source: Abbreviated, root: str = "f0"
) -> tuple[Abbreviated, Abbreviated, Abbreviated, Abbreviated]:
    """Exactly the four source F1/F2 actions, preserving ring numbers."""
    _free_owner(source.labeled, root)
    f2 = abbreviated_flip_f2(source, root)
    return (
        source,
        abbreviated_flip_f1(source),
        f2,
        abbreviated_flip_f1(f2),
    )


def abbreviated_flip_gauge(
    source: Abbreviated, root: str = "f0"
) -> Abbreviated:
    """Unique subgroup transform whose V4 root label is identity.

    This is a gauge *within F1/F2 only*, NOT a complete isotopy class.
    """
    group = abbreviated_flip_group(source, root)
    selected = tuple(x for x in group
                     if x.labeled.label(root) == IDENTITY)
    if len(selected) != 1:
        raise AssertionError("root V4 labels must distinguish subgroup actions")
    return selected[0]


def abbreviated_flip_orbit_key(source: Abbreviated, root: str = "f0") -> tuple:
    """Rooted plane-presentation isomorphism modulo only F1/F2."""
    return abbreviated_flip_gauge(source, root).rooted_code(root)


def flip_abbreviation_square(
    source: V4LabeledTree, sign: int, root: str = "f0"
) -> tuple[bool, bool]:
    """Check BOTH source F1/F2 equivariance squares on a full presentation.

    F2 pivot is identified by the same surviving free bond, NOT by
    changing vertex indices after ring caps have been erased.
    Invalid/unsupported abbreviations or missing roots fail closed.

    Returns (F1 commutes, F2 commutes). False means a real mismatch,
    not automatic evidence about topology.
    """
    pivot = _free_owner(source, root)
    ordinary = abbreviate_disjoint_rings(source, sign)
    _free_owner(ordinary.labeled, root)
    lhs_f1 = abbreviate_disjoint_rings(tangle_flip_f1(source), sign)
    rhs_f1 = abbreviated_flip_f1(ordinary)
    lhs_f2 = abbreviate_disjoint_rings(tangle_flip_f2(source, pivot), sign)
    rhs_f2 = abbreviated_flip_f2(ordinary, root)
    return lhs_f1 == rhs_f1, lhs_f2 == rhs_f2
