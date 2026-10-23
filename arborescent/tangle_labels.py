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


def tangle_flip_f3prime(
    source: V4LabeledTree, pivot: int, branch_edge: str, *,
    from_side: str = "before",
) -> V4LabeledTree:
    """Restricted W-preserving F3' move across an internal edge.

    Bonahon--Siebenmann §§12.4, 12.7.3, printed pp.155,167:
    relocate the single nonzero sector a across an incident branch
    edge. If a is odd, reverse rotations at vertices *within that
    branch* at odd distance from the pivot. For tangles, multiply
    V4 free labels in the branch by xi (even distance) or eta
    (odd distance). No boundary label outside that branch changes.

    This implementation requires that one of the two angular sectors
    adjacent to branch_edge has nonzero weight a and the other zero.
    It is a deliberately *strict* subset of the depicted general
    move; F3'', ring, and generalized weight transport are absent.
    """
    if type(pivot) is not int or pivot < 0 or pivot >= len(source.tree.vertices):
        raise ValueError("invalid pivot")
    if from_side not in ("before", "after"):
        raise ValueError("from_side must be before or after")
    v = source.tree.vertices[pivot]
    if branch_edge not in v.slots or not branch_edge.startswith("e"):
        raise ValueError("choose an incident internal edge")
    idx = v.slots.index(branch_edge)
    before, after = (idx - 1) % len(v.slots), idx
    old_slot = before if from_side == "before" else after
    new_slot = after if from_side == "before" else before
    a = v.sectors[old_slot]
    if a == 0 or v.sectors[new_slot] != 0:
        raise ValueError("weight transport requires one nonzero adjacent sector")
    if before == after:
        raise ValueError("no distinct adjacent sectors at pivot")

    adjacency = [set() for _ in source.tree.vertices]
    owner: dict[str, int] = {}
    for u, vertex in enumerate(source.tree.vertices):
        for s in vertex.slots:
            if s.startswith("e"):
                if s in owner:
                    other = owner[s]
                    adjacency[u].add(other)
                    adjacency[other].add(u)
                else:
                    owner[s] = u
    # The branch is exactly the component on the other side of edge
    # branch_edge after deleting its incidence to the pivot.
    other = next(j for j, vertex in enumerate(source.tree.vertices)
                 if j != pivot and branch_edge in vertex.slots)
    distance = {other: 1}
    pending = [other]
    for u in pending:
        for w in sorted(adjacency[u]):
            if w == pivot or w in distance:
                continue
            distance[w] = distance[u] + 1
            pending.append(w)

    sectors = list(v.sectors)
    sectors[old_slot] = 0
    sectors[new_slot] = a
    new_vertices = list(source.tree.vertices)
    new_vertices[pivot] = Vertex(v.slots, tuple(sectors))
    if a % 2:
        from .flips import reverse_vertex
        for u, depth in distance.items():
            if depth % 2:
                new_vertices[u] = reverse_vertex(new_vertices[u])
    transformed = WeightedPlaneTree(tuple(new_vertices))
    port_owners = {
        s: u for u, vertex in enumerate(source.tree.vertices)
        for s in vertex.slots if s.startswith("f")
    }
    return V4LabeledTree(transformed, tuple(
        (name, label ^ (ETA if distance[port_owners[name]] % 2 else XI)
         if a % 2 and port_owners[name] in distance else label)
        for name, label in source.labels
    ))
