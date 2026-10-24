"""A4.8: necessary rooted invariant and honest answers for the partial calculus.

A4.7's typed F1/F2, strict surviving-edge F3', and source-*motivated*
formal ring slides preserve the underlying rooted incidence tree,
TOTAL weight at each surviving vertex, its surviving free-leg count,
and the per-vertex strict §12.5 ring multiplicity. They need not
preserve sector LOCATION, cyclic order, or V4 labels individually.

A rooted UNORDERED ring-numbered skeleton thus gives a computable
NECESSARY invariant. Unequal fingerprints prove nonreachability under
the explicitly implemented PARTIAL FORMAL moves (not nonisotopy).

Equal fingerprints are NOT sufficient. Even exhaustive A4.7 BFS is
only exhaustive on the fixed indexed presentation graph; arbitrary
vertex/bond renaming and the full topological calculus are excluded.
A truncated search is explicitly UNKNOWN, never a negative answer.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .partial_joint_f3prime import (
    MixedTrace, _guard, explore_mixed,
)
from .ring_abbreviation import abbreviate_disjoint_rings
from .tangle_labels import V4LabeledTree

DecisionStatus = Literal[
    "reachable_formally",
    "excluded_by_necessary_invariant",
    "not_reached_in_complete_indexed_graph",
    "undetermined_cutoff",
]


def unordered_abbreviated_fingerprint(
    source: V4LabeledTree, sign: int, root: str = "f0"
) -> tuple:
    """Necessary isomorphism invariant of the restricted A4.7 move graph.

    Uses complete signed strict-cap recognition. Free V4 labels,
    incident cyclic order and angular-sector positions are discarded
    on purpose: all may change under the declared F1/F2/F3' moves.
    The root free bond remains distinguished; other free bonds are
    counted, not identified by their incidental storage names.

    Do NOT interpret equality as isotopy or even formal reachability.
    """
    aligned = _guard(source, sign, root)
    abbr = abbreviate_disjoint_rings(aligned, sign)
    tree = abbr.labeled.tree
    owners: dict[str, list[int]] = {}
    for i, vertex in enumerate(tree.vertices):
        for bond in vertex.slots:
            owners.setdefault(bond, []).append(i)
    if len(owners.get(root, ())) != 1:
        raise ValueError("root must survive abbreviation exactly once")

    def visit(vertex_index: int, incoming: str) -> tuple:
        vertex = tree.vertices[vertex_index]
        branches = []
        for bond in vertex.slots:
            if bond == incoming:
                continue
            if bond.startswith("f"):
                branches.append(("free",))
                continue
            ends = owners.get(bond, [])
            if len(ends) != 2:
                raise ValueError("interior bond must have two owners")
            next_index = next(i for i in ends if i != vertex_index)
            branches.append(("branch", visit(next_index, bond)))
        return (
            vertex.weight,
            abbr.rings[vertex_index],
            tuple(sorted(branches, key=repr)),
        )

    return (sign, visit(owners[root][0], root))


@dataclass(frozen=True)
class FormalReachabilityDecision:
    """A decision scoped to A4.7's indexed *formal presentation* graph."""
    status: DecisionStatus
    fingerprint_start: tuple
    fingerprint_target: tuple
    examined_states: int
    complete: bool
    trace: MixedTrace | None

    def verify_positive(self) -> bool:
        """Only a positive typed formal path is independently replayable."""
        return (
            self.status == "reachable_formally"
            and self.trace is not None
            and self.trace.verify()
        )


def decide_partial_formal_reachability(
    start: V4LabeledTree, target: V4LabeledTree, sign: int,
    root: str = "f0", max_states: int = 256,
) -> FormalReachabilityDecision:
    """Four-way result; never conflate index-complete BFS with isotopy.

    Unequal invariants: mathematically excluded under these moves,
    irrespective of search depth. Equal invariants: run the bounded
    indexed exact-state BFS. A reached target gives a replayable
    positive witness. An exhausted, complete graph only excludes
    the exact *indexed* target, not its renamings and not isotopy.
    A cutoff leaves reachability undetermined.
    """
    if type(max_states) is not int or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    canonical_start = _guard(start, sign, root)
    canonical_target = _guard(target, sign, root)
    start_key = unordered_abbreviated_fingerprint(canonical_start, sign, root)
    target_key = unordered_abbreviated_fingerprint(canonical_target, sign, root)
    if start_key != target_key:
        return FormalReachabilityDecision(
            "excluded_by_necessary_invariant", start_key, target_key,
            0, False, None,
        )
    orbit = explore_mixed(canonical_start, sign, root, max_states)
    if canonical_target in orbit.incoming:
        witness = orbit.trace(canonical_target)
        return FormalReachabilityDecision(
            "reachable_formally", start_key, target_key,
            len(orbit.incoming), orbit.complete, witness,
        )
    return FormalReachabilityDecision(
        ("not_reached_in_complete_indexed_graph" if orbit.complete
         else "undetermined_cutoff"),
        start_key, target_key, len(orbit.incoming),
        orbit.complete, None,
    )
