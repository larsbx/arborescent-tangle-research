"""A4.1: a conservative, reversible *formal* ring-position slide.

Bonahon--Siebenmann §12.4, printed p.155, has modified ring
moves R± that RELOCATE rather than delete a ±2--±1--±2 ring subtree.
§12.5, p.157, says the ring's planar location is irrelevant within
the source's canonical setting and separately erases rings into ring
numbers as a *representation change*.

This module swaps the cap's attachment stem with ONE adjacent ordinary
bond at the SAME base. Only zero-sector bases, strict picture-matching
caps, a fixed nonwrap anchor bond and no adjacent/mixed ring-center
exception are supported. It checks that the restricted §12.5
abbreviation does not change.

It is a source-motivated formal move on weighted-plane *presentations*,
NOT independently established topological isotopy, NOT a complete
Bonahon--Siebenmann R± implementation and NOT a canonicality test.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from .ring_abbreviation import abbreviate_disjoint_rings, strict_ring_caps
from .ring_signatures import ring_signatures
from .tangle_labels import V4LabeledTree
from .weighted_trees import Vertex, WeightedPlaneTree


@dataclass(frozen=True)
class RingSlide:
    """A replayable adjacent interchange at one cap attachment base."""
    sign: int
    center: int
    past_bond: str

    def __post_init__(self) -> None:
        if type(self.sign) is not int or self.sign not in (-1, 1):
            raise ValueError("sign must be ±1")
        if type(self.center) is not int or self.center < 0:
            raise ValueError("center must be a nonnegative vertex index")
        if not isinstance(self.past_bond, str) or not self.past_bond:
            raise ValueError("past_bond must be an incident ordinary bond")


def _eligible_context(source: V4LabeledTree, sign: int) -> tuple:
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled presentation required")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("ring sign must be ±1")
    motifs = ring_signatures(source)
    if any(m.sign != sign for m in motifs):
        raise ValueError("mixed-sign ring centers require separate source analysis")
    centers = {m.center for m in motifs}
    caps = strict_ring_caps(source, sign)
    if any(cap.base in centers for cap in caps):
        raise ValueError("ring attached to ring requires exceptional analysis")
    return caps


def apply_ring_slide(source: V4LabeledTree, move: RingSlide) -> V4LabeledTree:
    """Swap one cap stem and its adjacent noncap bond, without deleting.

    Base incident slots are a *linear* representation of a cyclic order.
    We disallow crossing index 0, retaining its anchor (normally the
    incoming or marked free root bond), and require all its sector
    weights zero. This is deliberately stricter than source R±.
    """
    if not isinstance(move, RingSlide):
        raise ValueError("expected a typed RingSlide")
    caps = _eligible_context(source, move.sign)
    matches = [cap for cap in caps if cap.center == move.center]
    if len(matches) != 1:
        raise ValueError("center has no unique strict ring cap")
    cap = matches[0]
    base = source.tree.vertices[cap.base]
    if any(base.sectors):
        raise ValueError("nonzero base sectors unsupported")
    slots = list(base.slots)
    old = slots.index(cap.stem)
    if move.past_bond not in slots:
        raise ValueError("past bond is not attached to this base")
    neighbor = slots.index(move.past_bond)
    if old == 0 or neighbor == 0 or abs(neighbor-old) != 1:
        raise ValueError("interchange must be adjacent and not cross anchor")
    # In particular, do not move one ring through a second ring; the
    # source has special adjacent-ring exceptions we do not implement.
    cap_stems = {c.stem for c in caps if c.base == cap.base}
    if move.past_bond in cap_stems:
        raise ValueError("adjacent cap interchange unsupported")
    slots[old], slots[neighbor] = slots[neighbor], slots[old]
    vertices = list(source.tree.vertices)
    vertices[cap.base] = Vertex(tuple(slots), base.sectors)
    target = V4LabeledTree(WeightedPlaneTree(tuple(vertices)), source.labels)
    if (abbreviate_disjoint_rings(source, move.sign).rooted_code()
            != abbreviate_disjoint_rings(target, move.sign).rooted_code()):
        raise AssertionError("formal slide changed abbreviated presentation")
    return target


def ring_slide_successors(
    source: V4LabeledTree, sign: int
) -> tuple[tuple[RingSlide, V4LabeledTree], ...]:
    """Deterministic, finite set of directly admissible strict slides."""
    caps = _eligible_context(source, sign)
    candidates = []
    cap_by_base = {}
    for c in caps:
        cap_by_base.setdefault(c.base, set()).add(c.stem)
    for cap in caps:
        base = source.tree.vertices[cap.base]
        if any(base.sectors):
            continue
        stem = base.slots.index(cap.stem)
        for neighbor in (stem-1, stem+1):
            if neighbor <= 0 or neighbor >= len(base.slots):
                continue
            other = base.slots[neighbor]
            if other in cap_by_base[cap.base]:
                continue
            move = RingSlide(sign, cap.center, other)
            candidates.append((move, apply_ring_slide(source, move)))
    return tuple(candidates)


@dataclass(frozen=True)
class RingTrace:
    start: V4LabeledTree
    target: V4LabeledTree
    moves: tuple[RingSlide, ...]

    def verify(self) -> bool:
        current = self.start
        for move in self.moves:
            current = apply_ring_slide(current, move)
        return current == self.target


@dataclass
class RingOrbit:
    start: V4LabeledTree
    incoming: dict[V4LabeledTree, tuple[V4LabeledTree, RingSlide] | None]
    complete: bool

    def trace(self, target: V4LabeledTree) -> RingTrace:
        if target not in self.incoming:
            raise ValueError("target not in explored formal orbit")
        sequence = []
        state = target
        while (previous := self.incoming[state]) is not None:
            parent, move = previous
            sequence.append(move)
            state = parent
        if state != self.start:
            raise AssertionError("malformed predecessor graph")
        witness = RingTrace(self.start, target, tuple(reversed(sequence)))
        if not witness.verify():
            raise AssertionError("failed to replay move certificate")
        return witness


def explore_ring_slides(
    start: V4LabeledTree, sign: int, max_states: int = 4096
) -> RingOrbit:
    """Exhaustive iff complete; an incomplete result is never a quotient."""
    if type(max_states) is not int or max_states < 1:
        raise ValueError("max_states must be positive")
    _eligible_context(start, sign)
    incoming: dict[V4LabeledTree, tuple[V4LabeledTree, RingSlide] | None] = {
        start: None
    }
    q = deque((start,))
    while q:
        state = q.popleft()
        for move, child in ring_slide_successors(state, sign):
            if child in incoming:
                continue
            if len(incoming) >= max_states:
                return RingOrbit(start, incoming, complete=False)
            incoming[child] = (state, move)
            q.append(child)
    return RingOrbit(start, incoming, complete=True)
