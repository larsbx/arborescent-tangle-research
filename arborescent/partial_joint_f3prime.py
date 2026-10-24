"""A4.7: replayable mixed F1/F2/F3'-strict/formal-ring reachability.

A4.4/A4.5 established restricted commuting squares with Bonahon--
Siebenmann §12.5 *abbreviation*; A4.6 proved a normal form only for
F1/F2 and the strictly formal ring-position swaps. This module does
NOT extend that normal-form completeness theorem to F3'.

Here we search a finite, size-bounded *indexed presentation graph*:
- F1 and root-pivot F2 are the already source-defined V4 flips.
- F3' is the existing STRICT guarded surviving-edge operator, and is
  accepted only when its abbreviation commuting square is verified.
- R-formal is A4.1's adjacent ring-stem presentation interchange,
  NEVER represented as a certified topology-preserving R± move.

Every intermediate state must admit complete strict ring recognition
and reversible rooted cyclic storage alignment. The BFS produces
replayable *positive formal reachability witnesses*. A cutoff is
never negative evidence. Even complete search is only complete for
this declared partial finite-index move graph, not for tangle isotopy.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from .abbreviated_f3prime import f3prime_abbreviation_square
from .abbreviated_flips import _free_owner, flip_abbreviation_square
from .ring_abbreviation import abbreviate_disjoint_rings, strict_ring_caps
from .ring_normalizer import align_rooted_cyclic_slots, normalize_rooted_ring_positions
from .ring_relocation import RingSlide, apply_ring_slide, ring_slide_successors
from .ring_signatures import ring_signatures
from .tangle_labels import (
    V4LabeledTree, tangle_flip_f1, tangle_flip_f2, tangle_flip_f3prime,
)


@dataclass(frozen=True)
class MixedMove:
    """A typed, replayable *formal* step. Not a 3D isotopy certificate."""
    kind: str
    pivot_free_bond: str = ""
    branch_edge: str = ""
    from_side: str = ""
    ring_slide: RingSlide | None = None

    def __post_init__(self) -> None:
        if self.kind in ("F1", "F2"):
            if any((self.pivot_free_bond, self.branch_edge,
                    self.from_side, self.ring_slide)):
                raise ValueError("global/root flip takes no per-move arguments")
        elif self.kind == "F3prime":
            if (not isinstance(self.pivot_free_bond, str)
                    or not self.pivot_free_bond.startswith("f")
                    or not isinstance(self.branch_edge, str)
                    or not self.branch_edge.startswith("e")
                    or self.from_side not in ("before", "after")
                    or self.ring_slide is not None):
                raise ValueError("F3prime needs a free pivot, branch edge and side")
        elif self.kind == "R-formal":
            if (not isinstance(self.ring_slide, RingSlide)
                    or any((self.pivot_free_bond,
                            self.branch_edge, self.from_side))):
                raise ValueError("R-formal needs exactly one typed RingSlide")
        else:
            raise ValueError("unknown mixed formal move")


def _guard(source: V4LabeledTree, sign: int, root: str) -> V4LabeledTree:
    """Reject incomplete ring data, bad anchors and unsupported signs."""
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled weighted tree required")
    aligned = align_rooted_cyclic_slots(source, root)
    _free_owner(aligned, root)
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("ring sign must be ±1")
    signatures = ring_signatures(aligned)
    caps = strict_ring_caps(aligned, sign)
    if len(signatures) != len(caps):
        raise ValueError("incomplete signed ring-cap recognition")
    normalized = normalize_rooted_ring_positions(aligned, sign, root)
    if not normalized.verify():
        raise AssertionError("unverifiable strict root alignment")
    abbreviated = abbreviate_disjoint_rings(aligned, sign)
    if sum(abbreviated.rings) != len(caps):
        raise AssertionError("unmatched abbreviation multiplicity")
    return aligned


def apply_mixed(
    source: V4LabeledTree, sign: int, move: MixedMove, root: str = "f0"
) -> V4LabeledTree:
    """Apply one guarded source flip or formal ring slide, then root-align.

    For F3' and F1/F2, checks the appropriate actual before/after
    §12.5 commuting square independently. For R-formal, ensures
    abbreviation is unchanged. In all cases surviving ring counts
    are preserved, but F3' may change the F1/F2 abbreviated key.
    """
    if not isinstance(move, MixedMove):
        raise ValueError("typed MixedMove required")
    state = _guard(source, sign, root)
    before = abbreviate_disjoint_rings(state, sign)
    if move.kind == "F1":
        if flip_abbreviation_square(state, sign, root) != (True, True):
            raise ValueError("invalid F1/F2 abbreviation square")
        target = tangle_flip_f1(state)
    elif move.kind == "F2":
        if flip_abbreviation_square(state, sign, root) != (True, True):
            raise ValueError("invalid F1/F2 abbreviation square")
        target = tangle_flip_f2(state, _free_owner(state, root))
    elif move.kind == "F3prime":
        if not f3prime_abbreviation_square(
                state, sign, move.pivot_free_bond, move.branch_edge,
                from_side=move.from_side):
            raise ValueError("F3prime exits source-compatible abbreviation domain")
        target = tangle_flip_f3prime(
            state, _free_owner(state, move.pivot_free_bond),
            move.branch_edge, from_side=move.from_side
        )
    elif move.kind == "R-formal":
        target = apply_ring_slide(state, move.ring_slide)
    else:
        raise AssertionError("MixedMove constructor guarantees type")
    target = _guard(target, sign, root)
    after = abbreviate_disjoint_rings(target, sign)
    if before.rings != after.rings:
        raise ValueError("move changed the per-base ring-number vector")
    if move.kind == "R-formal" and before.rooted_code(root) != after.rooted_code(root):
        raise AssertionError("formal ring slide changed abbreviation")
    return target


def mixed_successors(
    source: V4LabeledTree, sign: int, root: str = "f0",
) -> tuple[tuple[MixedMove, V4LabeledTree], ...]:
    """Exact currently implemented generators; no guessed R± or full F3'."""
    state = _guard(source, sign, root)
    options = [MixedMove("F1"), MixedMove("F2")]
    residual = abbreviate_disjoint_rings(state, sign)
    surviving = {s for v in residual.labeled.tree.vertices
                 for s in v.slots if s.startswith("e")}
    free = tuple(sorted(b for b, _ in state.labels))
    for pivot_bond in free:
        vertex = state.tree.vertices[_free_owner(state, pivot_bond)]
        for edge in vertex.slots:
            if edge in surviving:
                for side in ("before", "after"):
                    options.append(MixedMove(
                        "F3prime", pivot_bond, edge, side
                    ))
    options.extend(MixedMove("R-formal", ring_slide=move)
                   for move, _ in ring_slide_successors(state, sign))
    seen: set[V4LabeledTree] = set()
    results = []
    for move in options:
        try:
            child = apply_mixed(state, sign, move, root)
        except ValueError:
            continue
        if child not in seen:
            seen.add(child)
            results.append((move, child))
    return tuple(results)


@dataclass(frozen=True)
class MixedTrace:
    """A replayable path in the explicitly partial formal move graph."""
    start: V4LabeledTree
    target: V4LabeledTree
    sign: int
    root: str
    moves: tuple[MixedMove, ...]

    def verify(self) -> bool:
        current = _guard(self.start, self.sign, self.root)
        for move in self.moves:
            current = apply_mixed(current, self.sign, move, self.root)
        return current == self.target


@dataclass
class MixedOrbit:
    start: V4LabeledTree
    sign: int
    root: str
    incoming: dict[V4LabeledTree, tuple[V4LabeledTree, MixedMove] | None]
    complete: bool

    def trace(self, target: V4LabeledTree) -> MixedTrace:
        if target not in self.incoming:
            raise ValueError("target was not reached")
        reverse = []
        cur = target
        while (row := self.incoming[cur]) is not None:
            parent, move = row
            reverse.append(move)
            cur = parent
        if cur != self.start:
            raise AssertionError("corrupt parent map")
        result = MixedTrace(
            self.start, target, self.sign, self.root,
            tuple(reversed(reverse)),
        )
        if not result.verify():
            raise AssertionError("formal move trace replay failed")
        return result


def explore_mixed(
    source: V4LabeledTree, sign: int, root: str = "f0",
    max_states: int = 256,
) -> MixedOrbit:
    """Bounded exact-state BFS, with explicit *nonexhaustive* cutoff."""
    if type(max_states) is not int or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    start = _guard(source, sign, root)
    incoming: dict[V4LabeledTree, tuple[V4LabeledTree, MixedMove] | None] = {
        start: None
    }
    queue = deque([start])
    while queue:
        state = queue.popleft()
        for move, child in mixed_successors(state, sign, root):
            if child in incoming:
                continue
            if len(incoming) >= max_states:
                return MixedOrbit(start, sign, root, incoming, False)
            incoming[child] = (state, move)
            queue.append(child)
    return MixedOrbit(start, sign, root, incoming, True)
