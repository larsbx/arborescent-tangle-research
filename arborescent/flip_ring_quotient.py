"""A4.6: certified joint *formal presentation* quotient of F1/F2 + ring slides.

Bonahon--Siebenmann §12.7 tangle flips F1/F2 act on all cyclic
rotations and every free-bond V4 label. A4.1 strict *formal* ring
slides rearrange a full ±2--±1--±2 cap at an all-zero-sector base.
A4.2/A4.3 canonically stabilize only this source-motivated restricted
ring-slide subsystem; §12.5 ring abbreviation is a DIFFERENT
representation change.

For completely recognized disjoint uniform-sign caps, with all base
sectors zero, define a two-layer formal certificate:
    (1) unique F1/F2 variant carrying identity root V4;
    (2) rooted cyclic alignment followed by inversion-decreasing
        strict ring-slide normalization.

The resulting marked presentation code is invariant under the
implemented F1/F2 + strict-slide subsystem; it can also be checked
against the F1/F2 orbit key of the restricted abbreviated notation.

Neither this formal code nor the BFS below certifies 3D tangle isotopy,
full R± / F3', or exception-complete canonicality.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from .abbreviated_flips import abbreviated_flip_orbit_key, _free_owner
from .ring_abbreviation import abbreviate_disjoint_rings
from .ring_normalizer import (
    RootAlignedRingCertificate, normalize_rooted_ring_positions,
    align_rooted_cyclic_slots,
)
from .ring_relocation import RingSlide, apply_ring_slide, ring_slide_successors
from .tangle_labels import (
    IDENTITY, V4LabeledTree, tangle_flip_f1, tangle_flip_f2,
    tangle_flip_group,
)


def _domain(source: V4LabeledTree, sign: int, root: str) -> None:
    """Require *complete* strict cap recognition and rooted anchors."""
    _free_owner(source, root)
    aligned = normalize_rooted_ring_positions(source, sign, root)
    if not aligned.verify():
        raise AssertionError("invalid strict formal ring normalization")
    abbr = abbreviate_disjoint_rings(source, sign)
    _free_owner(abbr.labeled, root)
    if not isinstance(aligned.normal, V4LabeledTree):
        raise AssertionError("unexpected normal form type")


@dataclass(frozen=True)
class FlipRingNormalCertificate:
    """One sourced finite flip, one cyclic isomorphism, typed formal slides."""
    source: V4LabeledTree
    sign: int
    root: str
    flip_index: int
    flip_gauge: V4LabeledTree
    slide: RootAlignedRingCertificate

    def verify(self) -> bool:
        if type(self.flip_index) is not int or self.flip_index not in range(4):
            return False
        if type(self.sign) is not int or self.sign not in (-1, 1):
            return False
        if (self.slide.root != self.root
                or self.slide.slide_certificate.sign != self.sign):
            return False
        group = tangle_flip_group(self.source, self.root)
        if self.flip_gauge != group[self.flip_index]:
            return False
        if self.flip_gauge.label(self.root) != IDENTITY:
            return False
        if self.slide.source != self.flip_gauge or not self.slide.verify():
            return False
        return True

    @property
    def normal(self) -> V4LabeledTree:
        return self.slide.normal


def normalize_flip_ring(
    source: V4LabeledTree, sign: int, root: str = "f0",
    max_slides: int | None = None,
) -> FlipRingNormalCertificate:
    """Produce the unique rooted-V4-gauge, formal-ring-stable representative.

    Does not apply F3', exceptional ring moves, or isotopy reduction.
    In particular *no ring subtree is ever deleted* in the returned
    representative: abbreviation is used solely as a consistency
    check and as an alternative quotient invariant.
    """
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled full presentation required")
    _domain(source, sign, root)
    choices = tangle_flip_group(source, root)
    selected = [(i, q) for i, q in enumerate(choices)
                if q.label(root) == IDENTITY]
    if len(selected) != 1:
        raise AssertionError("F1/F2 root gauge must be unique")
    i, gauge = selected[0]
    step = normalize_rooted_ring_positions(gauge, sign, root, max_slides)
    certificate = FlipRingNormalCertificate(source, sign, root, i, gauge, step)
    if not certificate.verify():
        raise AssertionError("joint formal certificate did not replay")
    if (abbreviated_flip_orbit_key(abbreviate_disjoint_rings(source, sign), root)
            != abbreviated_flip_orbit_key(
                abbreviate_disjoint_rings(certificate.normal, sign), root)):
        raise AssertionError("flip+slide changed F1/F2 abbreviated key")
    return certificate


def flip_ring_normal_key(
    source: V4LabeledTree, sign: int, root: str = "f0"
) -> tuple:
    """Stable presentation-isomorphism key under this limited joint action."""
    return normalize_flip_ring(source, sign, root).normal.labelled_rooted_code(root)


@dataclass(frozen=True)
class JointMove:
    """F1, rooted-pivot F2, or an A4.1 typed *formal* ring slide."""
    kind: str
    slide: RingSlide | None = None

    def __post_init__(self) -> None:
        if self.kind in ("F1", "F2"):
            if self.slide is not None:
                raise ValueError("flip move must not carry ring slide arguments")
        elif self.kind == "R-formal":
            if not isinstance(self.slide, RingSlide):
                raise ValueError("formal ring move requires a typed RingSlide")
        else:
            raise ValueError("unsupported formal move")


def apply_joint(
    source: V4LabeledTree, move: JointMove, root: str = "f0"
) -> V4LabeledTree:
    """Apply exactly one implemented move with explicit type/precondition guards."""
    if not isinstance(move, JointMove) or not isinstance(source, V4LabeledTree):
        raise ValueError("expected V4 presentation and typed joint move")
    pivot = _free_owner(source, root)
    if move.kind == "F1":
        return tangle_flip_f1(source)
    if move.kind == "F2":
        return tangle_flip_f2(source, pivot)
    if move.kind == "R-formal":
        # Slot zero is a rooted anchor, never an arbitrary storage start.
        aligned = align_rooted_cyclic_slots(source, root)
        target = apply_ring_slide(aligned, move.slide)
        # Restore the caller's cyclic storage convention at EVERY vertex.
        # This is only a presentation isomorphism; it preserves exact replay
        # and avoids counting cyclic starts as extra indexed orbit states.
        from .weighted_trees import Vertex, WeightedPlaneTree
        vertices = []
        for before, after in zip(source.tree.vertices, target.tree.vertices):
            offset = after.slots.index(before.slots[0])
            vertices.append(Vertex(after.slots[offset:] + after.slots[:offset],
                                   after.sectors[offset:] + after.sectors[:offset]))
        return V4LabeledTree(WeightedPlaneTree(tuple(vertices)), target.labels)
    raise AssertionError("invalid move passed constructor")


def joint_successors(
    source: V4LabeledTree, sign: int, root: str = "f0"
) -> tuple[tuple[JointMove, V4LabeledTree], ...]:
    """Finite exact-state successors for the *restricted formal* action."""
    _domain(source, sign, root)
    moves = [(JointMove("F1"), apply_joint(source, JointMove("F1"), root)),
             (JointMove("F2"), apply_joint(source, JointMove("F2"), root))]
    moves.extend(
        (JointMove("R-formal", slide),
         apply_joint(source, JointMove("R-formal", slide), root))
        for slide, _ in ring_slide_successors(
            align_rooted_cyclic_slots(source, root), sign)
    )
    seen = set()
    distinct = []
    for move, child in moves:
        if child not in seen:
            seen.add(child)
            distinct.append((move, child))
    return tuple(distinct)


@dataclass(frozen=True)
class JointTrace:
    start: V4LabeledTree
    target: V4LabeledTree
    moves: tuple[JointMove, ...]
    root: str

    def verify(self) -> bool:
        current = self.start
        for move in self.moves:
            current = apply_joint(current, move, self.root)
        return current == self.target


@dataclass
class JointOrbit:
    start: V4LabeledTree
    root: str
    incoming: dict[V4LabeledTree, tuple[V4LabeledTree, JointMove] | None]
    complete: bool

    def trace(self, target: V4LabeledTree) -> JointTrace:
        if target not in self.incoming:
            raise ValueError("target is not in the explored formal joint graph")
        moves = []
        position = target
        while (prior := self.incoming[position]) is not None:
            parent, move = prior
            moves.append(move)
            position = parent
        if position != self.start:
            raise AssertionError("corrupt predecessor graph")
        witness = JointTrace(self.start, target, tuple(reversed(moves)), self.root)
        if not witness.verify():
            raise AssertionError("joint trace did not replay")
        return witness


def explore_joint(
    source: V4LabeledTree, sign: int, root: str = "f0",
    max_states: int = 4096
) -> JointOrbit:
    """Breadth-first search; cutoff marked incomplete, never a false orbit count."""
    if type(max_states) is not int or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    _domain(source, sign, root)
    incoming: dict[V4LabeledTree, tuple[V4LabeledTree, JointMove] | None] = {
        source: None
    }
    todo = deque((source,))
    while todo:
        state = todo.popleft()
        for move, child in joint_successors(state, sign, root):
            if child in incoming:
                continue
            if len(incoming) >= max_states:
                return JointOrbit(source, root, incoming, False)
            incoming[child] = (state, move)
            todo.append(child)
    return JointOrbit(source, root, incoming, True)

