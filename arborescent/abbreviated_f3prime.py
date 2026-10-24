"""A4.5: restricted F3'-abbreviation commuting square on surviving edges.

This is a *presentation-level* naturality check for source §12.7.3
modified flip F3' on V4-labelled weighted trees, after §12.5's
separate, lossy ring-numbered abbreviation.

Only an F3' whose pivot and branch edge BOTH survive strict signed
abbreviation is considered. The pivot is identified by an explicit
surviving free bond, never a stale vertex index. Every ring cap must
remain strictly recognized after the full move; otherwise the function
fails closed. All surviving ring counts are transported unchanged.

No F3' crossing erased ring stems, full F3' calculus, canonicality
exception, or independent isotopy theorem is asserted.
"""
from __future__ import annotations

from .ring_abbreviation import Abbreviated, abbreviate_disjoint_rings, strict_ring_caps
from .ring_signatures import ring_signatures
from .abbreviated_flips import _free_owner
from .tangle_labels import V4LabeledTree, tangle_flip_f3prime


def abbreviated_f3prime(
    source: Abbreviated,
    pivot_free_bond: str,
    branch_edge: str,
    *,
    from_side: str = "before",
) -> Abbreviated:
    """Apply already implemented strict F3' on residual weighted tree.

    The number and attachment vertex of every abbreviated ring is
    unchanged. The F3' pivot's free bond must survive abbreviation.
    """
    if not isinstance(source, Abbreviated):
        raise ValueError("expected ring-numbered abbreviated presentation")
    pivot = _free_owner(source.labeled, pivot_free_bond)
    next_labeled = tangle_flip_f3prime(
        source.labeled, pivot, branch_edge, from_side=from_side
    )
    return Abbreviated(next_labeled, source.rings, source.sign)


def f3prime_abbreviation_square(
    source: V4LabeledTree,
    sign: int,
    pivot_free_bond: str,
    branch_edge: str,
    *,
    from_side: str = "before",
) -> bool:
    """Check `abbreviate(F3'(T)) = F3'*(abbreviate(T))`.

    Fail closed if a cap is absent/unrecognized after the full move,
    if a pivot or chosen edge is erased by abbreviation, or if the
    relevant F3' strict preconditions fail.

    The square concerns the *existing partial* F3' implementation and
    the strict abstraction; it does not certify topological isotopy.
    """
    if not isinstance(source, V4LabeledTree):
        raise ValueError("expected original V4-labelled presentation")
    # Ensure that every necessary ring motif is in the strict
    # supported cap representation *before* any reduction.
    original_caps = strict_ring_caps(source, sign)
    if len(original_caps) != len(ring_signatures(source)):
        raise ValueError("incomplete source ring-signature recognition")
    abbreviation = abbreviate_disjoint_rings(source, sign)
    if sum(abbreviation.rings) != len(original_caps):
        raise ValueError("incomplete source cap recognition")
    pivot = _free_owner(source, pivot_free_bond)
    _free_owner(abbreviation.labeled, pivot_free_bond)
    # The branch must be an *interior surviving* edge, not a ring stem
    # whose deletion is part of the separate abbreviation operation.
    remaining = {bond for v in abbreviation.labeled.tree.vertices
                 for bond in v.slots}
    if branch_edge not in remaining or not branch_edge.startswith("e"):
        raise ValueError("F3' branch must survive ring abbreviation")

    full_target = tangle_flip_f3prime(
        source, pivot, branch_edge, from_side=from_side
    )
    destination_caps = strict_ring_caps(full_target, sign)
    if (len(destination_caps) != len(original_caps)
            or {cap.center for cap in destination_caps}
               != {cap.center for cap in original_caps}):
        raise ValueError("F3' exits the restricted ring recognition domain")
    from_full = abbreviate_disjoint_rings(full_target, sign)
    if from_full.rings != abbreviation.rings:
        raise ValueError("F3' changed ring attachment counts")
    from_abbreviation = abbreviated_f3prime(
        abbreviation, pivot_free_bond, branch_edge, from_side=from_side
    )
    return from_full == from_abbreviation
