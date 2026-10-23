"""A2.2: V4 boundary label bookkeeping for *tangle* F1/F2 moves.

Bonahon--Siebenmann sec. 12.7.1, printed p. 166 (PDF p. 177):
F1 applies zeta to every free bond; F2 applies xi at even distance
and eta at odd distance from its pivot. The three nonidentity elements
of the Klein four-group obey xi*eta=zeta. We encode V4 as F2^2:

    identity=0, xi=1, eta=2, zeta=3; composition is XOR.

Unlike the earlier unlabelled pair flips, these actions are tangle
presentation moves only when labels are updated. "Boundary fixed"
here means the labelled Conway parametrizations in the sense of
the cited theorem; it is NOT a pointwise-fixed ambient isotopy oracle.

This is solely the F1/F2 subgroup, NOT F3', rings, canonical
normalization, or a complete isotopy decision procedure.
"""
from __future__ import annotations

from dataclasses import dataclass

from .flips import bipartition, flip_f1, flip_f2
from .weighted_trees import (
    Spec, Vertex, WeightedPlaneTree, enumerate_specs, to_weighted_tree,
)

IDENTITY, XI, ETA, ZETA = 0, 1, 2, 3
V4 = frozenset((IDENTITY, XI, ETA, ZETA))


@dataclass(frozen=True)
class V4LabeledTree:
    """Weighted plane presentation with *every* free bond V4-labelled.

    Labels are stored as sorted (bond_name, element) tuples.
    The free-bond names themselves are not topological invariants.
    """

    tree: WeightedPlaneTree
    labels: tuple[tuple[str, int], ...]

    def __post_init__(self) -> None:
        if not isinstance(self.tree, WeightedPlaneTree):
            raise ValueError("expected a weighted plane tree")
        free = {
            slot
            for vertex in self.tree.vertices
            for slot in vertex.slots if slot.startswith("f")
        }
        if (not isinstance(self.labels, tuple)
                or any(not isinstance(item, tuple) or len(item) != 2
                       for item in self.labels)):
            raise ValueError("labels must be a tuple of pairs")
        names: list[str] = []
        for name, value in self.labels:
            if (not isinstance(name, str) or type(value) is not int
                    or value not in V4):
                raise ValueError("invalid V4 label")
            names.append(name)
        if len(set(names)) != len(names) or set(names) != free:
            raise ValueError("label exactly every free bond, once")

    def label(self, bond: str) -> int:
        for name, value in self.labels:
            if name == bond:
                return value
        raise ValueError("missing free bond label")

    def labelled_rooted_code(self, root: str = "f0") -> tuple:
        """Faithful marked *presentation*-isomorphism code with free labels.

        Includes the root V4 label and other free V4 labels, unlike the
        unlabelled pair code in weighted_trees.py.
        """
        incidence: dict[str, list[int]] = {}
        for i, v in enumerate(self.tree.vertices):
            for s in v.slots:
                incidence.setdefault(s, []).append(i)
        if root not in incidence or not root.startswith("f"):
            raise ValueError("root must be a present free bond")

        def code(i: int, incoming: str) -> tuple:
            v = self.tree.vertices[i]
            shift = v.slots.index(incoming)
            slots = v.slots[shift:] + v.slots[:shift]
            sectors = v.sectors[shift:] + v.sectors[:shift]
            children = []
            for slot in slots[1:]:
                if slot[0] == "f":
                    children.append(("free", self.label(slot)))
                else:
                    other = next(j for j in incidence[slot] if j != i)
                    children.append(code(other, slot))
            return (sectors, tuple(children))

        return (self.label(root), code(incidence[root][0], root))


def label_free_ports(tree: WeightedPlaneTree, default: int = IDENTITY) -> V4LabeledTree:
    """Create the explicit all-equal free-label presentation."""
    if type(default) is not int or default not in V4:
        raise ValueError("invalid V4 element")
    free = sorted({
        s for v in tree.vertices for s in v.slots if s.startswith("f")
    })
    return V4LabeledTree(tree, tuple((s, default) for s in free))


def _relabel(
    source: V4LabeledTree, target: WeightedPlaneTree, group_by_vertex: dict[int, int]
) -> V4LabeledTree:
    owner = {
        s: i for i, v in enumerate(source.tree.vertices)
        for s in v.slots if s.startswith("f")
    }
    return V4LabeledTree(target, tuple(
        (name, value ^ group_by_vertex[owner[name]])
        for name, value in source.labels
    ))


def tangle_flip_f1(source: V4LabeledTree) -> V4LabeledTree:
    """Source §12.7.1 F1: reverse all rotations; zeta on every free label."""
    return _relabel(
        source, flip_f1(source.tree),
        {i: ZETA for i in range(len(source.tree.vertices))},
    )


def tangle_flip_f2(source: V4LabeledTree, pivot: int) -> V4LabeledTree:
    """Source §12.7.1 F2: xi at even distance and eta at odd distance."""
    even, odd = bipartition(source.tree, pivot)
    return _relabel(
        source, flip_f2(source.tree, pivot),
        {**{i: XI for i in even}, **{i: ETA for i in odd}},
    )


def tangle_flip_group(source: V4LabeledTree, root: str = "f0") -> tuple[V4LabeledTree, ...]:
    """The four actions I,F1,F2(root-vertex),F1∘F2(root-vertex)."""
    owners = [
        i for i, v in enumerate(source.tree.vertices) if root in v.slots
    ]
    if not root.startswith("f") or len(owners) != 1:
        raise ValueError("root must be a present free bond")
    p = owners[0]
    f2 = tangle_flip_f2(source, p)
    return (source, tangle_flip_f1(source), f2, tangle_flip_f1(f2))


def normalize_root_label(source: V4LabeledTree, root: str = "f0") -> V4LabeledTree:
    """Unique transform in the F1/F2 *action subgroup* with root label 0.

    This is a gauge slice for a labelled presentation, not a
    canonical form for all topologically equivalent tangles.
    """
    matches = [s for s in tangle_flip_group(source, root)
               if s.label(root) == IDENTITY]
    if len(matches) != 1:
        raise AssertionError("root label must distinguish all 4 subgroup actions")
    return matches[0]


def tangle_flip_orbit_key(source: V4LabeledTree, root: str = "f0") -> tuple:
    """Group quotient, then marked plane-tree isomorphism; nothing else."""
    return normalize_root_label(source, root).labelled_rooted_code(root)


def tangle_grammar_orbit_count(n: int, labels: tuple[int, ...] = (2, 3)) -> int:
    """Enumerate the labelled F1/F2 quotient; small n only.

    One free root port has all 4 labels, so the action is free and
    quotient count equals the number of underlying ordered syntax
    presentations (NOT the earlier unlabelled pair-flip quotient).
    """
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive number of forks")
    keys = set()
    for spec in enumerate_specs(n, labels):
        tree = to_weighted_tree(spec, labels)
        for v4 in V4:
            keys.add(tangle_flip_orbit_key(
                V4LabeledTree(tree, (("f0", v4),))
            ))
    return len(keys)
