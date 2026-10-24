"""A4.9: bounded A4.7 formal reachability modulo *rooted presentation isomorphism*.

A4.8 showed that a complete fixed-index BFS can miss a differently
numbered storage representation of exactly the SAME rooted weighted
plane tree. This module quotients the explored graph by the existing
faithful V4-labelled rooted presentation code.

Every positive certificate consists of:
    (1) an exact typed A4.7 formal move trace to a stored representative;
    (2) equality of that representative's complete marked rooted
        presentation-isomorphism code with the requested target.

A complete BFS covers the declared A4.7 *partial formal move graph*
modulo rooted plane-tree presentation isomorphism, conditional on
the explicitly implemented moves being equivariant under renaming.
It does NOT cover unimplemented R±/F3', exceptional canonicity,
Conway/plumbing topology, or actual 3D tangle isotopy.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Literal

from .mixed_invariants import unordered_abbreviated_fingerprint
from .partial_joint_f3prime import (
    MixedMove, MixedTrace, _guard, mixed_successors,
)
from .tangle_labels import V4LabeledTree

ClassStatus = Literal[
    "reachable_up_to_presentation_isomorphism",
    "excluded_by_necessary_invariant",
    "not_reached_in_complete_formal_class_graph",
    "undetermined_cutoff",
]


def rooted_presentation_key(source: V4LabeledTree, root: str = "f0") -> tuple:
    """Faithful rooted-V4 presentation code, NOT topology or isotopy."""
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled weighted tree required")
    return source.labelled_rooted_code(root)


@dataclass(frozen=True)
class PresentationClassCertificate:
    """Typed formal path ending at a proven rooted-presentation isomorph."""
    target: V4LabeledTree
    trace: MixedTrace
    root: str

    def verify(self) -> bool:
        return (
            self.trace.verify()
            and rooted_presentation_key(self.trace.target, self.root)
                == rooted_presentation_key(self.target, self.root)
        )


@dataclass
class ClassOrbit:
    """Exact discovered rooted presentation classes for partial A4.7 moves."""
    start: V4LabeledTree
    sign: int
    root: str
    representatives: dict[tuple, V4LabeledTree]
    incoming: dict[tuple, tuple[tuple, MixedMove] | None]
    complete: bool

    def certificate(self, target: V4LabeledTree) -> PresentationClassCertificate:
        code = rooted_presentation_key(target, self.root)
        if code not in self.representatives:
            raise ValueError("target presentation class has not been reached")
        sequence: list[MixedMove] = []
        position = code
        while (prior := self.incoming[position]) is not None:
            parent_code, move = prior
            sequence.append(move)
            position = parent_code
        if position != rooted_presentation_key(self.start, self.root):
            raise AssertionError("corrupt class-predecessor graph")
        representative = self.representatives[code]
        witness = MixedTrace(
            self.start, representative, self.sign, self.root,
            tuple(reversed(sequence)),
        )
        cert = PresentationClassCertificate(target, witness, self.root)
        if not cert.verify():
            raise AssertionError("formal path/isomorphism certificate failed replay")
        return cert


def explore_presentation_classes(
    source: V4LabeledTree, sign: int, root: str = "f0",
    max_classes: int = 256,
) -> ClassOrbit:
    """Quotient A4.7 successor graph by rooted presentation code.

    Every discovered *class* is expanded once using a representative.
    Rooted-code congruence of the declared moves is tested
    independently on reordered/renamed representatives.
    `complete=False` means the requested class budget was hit.
    """
    if type(max_classes) is not int or max_classes < 1:
        raise ValueError("max_classes must be a positive integer")
    start = _guard(source, sign, root)
    key = rooted_presentation_key(start, root)
    reps: dict[tuple, V4LabeledTree] = {key: start}
    incoming: dict[tuple, tuple[tuple, MixedMove] | None] = {key: None}
    todo = deque([key])
    while todo:
        parent = todo.popleft()
        representative = reps[parent]
        for move, child in mixed_successors(representative, sign, root):
            code = rooted_presentation_key(child, root)
            if code in reps:
                continue
            if len(reps) >= max_classes:
                return ClassOrbit(start, sign, root, reps, incoming, False)
            reps[code] = child
            incoming[code] = (parent, move)
            todo.append(code)
    return ClassOrbit(start, sign, root, reps, incoming, True)


@dataclass(frozen=True)
class ClassReachabilityDecision:
    status: ClassStatus
    start_invariant: tuple
    target_invariant: tuple
    examined_classes: int
    complete: bool
    certificate: PresentationClassCertificate | None

    def verify_positive(self) -> bool:
        return (
            self.status == "reachable_up_to_presentation_isomorphism"
            and self.certificate is not None
            and self.certificate.verify()
        )


def decide_formal_class_reachability(
    source: V4LabeledTree, target: V4LabeledTree, sign: int,
    root: str = "f0", max_classes: int = 256,
) -> ClassReachabilityDecision:
    """Four-way partial *formal* class decision; NEVER an isotopy oracle."""
    if type(max_classes) is not int or max_classes < 1:
        raise ValueError("max_classes must be a positive integer")
    start = _guard(source, sign, root)
    goal = _guard(target, sign, root)
    left = unordered_abbreviated_fingerprint(start, sign, root)
    right = unordered_abbreviated_fingerprint(goal, sign, root)
    if left != right:
        return ClassReachabilityDecision(
            "excluded_by_necessary_invariant", left, right,
            0, False, None,
        )
    orbit = explore_presentation_classes(start, sign, root, max_classes)
    target_code = rooted_presentation_key(goal, root)
    if target_code in orbit.representatives:
        certificate = orbit.certificate(target)
        return ClassReachabilityDecision(
            "reachable_up_to_presentation_isomorphism",
            left, right, len(orbit.representatives),
            orbit.complete, certificate,
        )
    return ClassReachabilityDecision(
        ("not_reached_in_complete_formal_class_graph" if orbit.complete
         else "undetermined_cutoff"),
        left, right, len(orbit.representatives),
        orbit.complete, None,
    )
