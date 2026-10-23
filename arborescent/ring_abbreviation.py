"""Restricted Bonahon--Siebenmann ring *abbreviation*, NOT an R± move.

§12.4 p.155 R± moves relocate a ring subtree in the plane, without
deleting it. §12.5 pp.157-158 *separately defines* abbreviated trees,
where disjoint ring subtrees are erased and their multiplicities are
recorded as integers at attachment vertices. This module implements
only that lossy-to-ring-number *notation change*, with conservative
guards and explicit refusal on overlapping caps or an empty base.

It is NOT a knot/tangle isotopy classifier, and the source's special
two-adjacent-ring configuration is rejected rather than guessed.
"""
from __future__ import annotations
from dataclasses import dataclass

from .ring_signatures import ring_signatures
from .tangle_labels import V4LabeledTree
from .weighted_trees import Vertex, WeightedPlaneTree


@dataclass(frozen=True)
class RingCap:
    sign: int
    base: int
    center: int
    tips: tuple[int, int]
    stem: str


@dataclass(frozen=True)
class Abbreviated:
    """Explicit ring-numbered tree: weights and V4 labels still survive."""
    labeled: V4LabeledTree
    rings: tuple[int, ...]
    sign: int

    def __post_init__(self):
        if self.sign not in (-1, 1) or type(self.sign) is not int:
            raise ValueError("sign must be ±1")
        if len(self.rings) != len(self.labeled.tree.vertices):
            raise ValueError("one ring count for each surviving vertex")
        if any(type(x) is not int or x < 0 for x in self.rings):
            raise ValueError("ring multiplicities must be nonnegative integers")

    def rooted_code(self, root: str = "f0") -> tuple:
        """Marked ring-numbered *presentation isomorphism* code, not isotopy."""
        tree = self.labeled.tree
        occur: dict[str, list[int]] = {}
        for i, vertex in enumerate(tree.vertices):
            for bond in vertex.slots:
                occur.setdefault(bond, []).append(i)
        if root not in occur or not root.startswith("f") or len(occur[root]) != 1:
            raise ValueError("root must be an existing free bond")

        def walk(i: int, incoming: str):
            vertex = tree.vertices[i]
            ix = vertex.slots.index(incoming)
            slots = vertex.slots[ix:] + vertex.slots[:ix]
            sectors = vertex.sectors[ix:] + vertex.sectors[:ix]
            descendants = []
            for bond in slots[1:]:
                if bond.startswith("f"):
                    descendants.append(("free", self.labeled.label(bond)))
                else:
                    other = next(j for j in occur[bond] if j != i)
                    descendants.append(walk(other, bond))
            return (self.rings[i], sectors, tuple(descendants))

        return (self.sign, self.labeled.label(root),
                walk(occur[root][0], root))


def strict_ring_caps(source: V4LabeledTree, sign: int) -> tuple[RingCap, ...]:
    """Conservatively recognize diagram-shaped ±2/±1/±2 ring subtrees.

    Two terminal ±2 tips meet a valence-three center of weight ±1,
    with center weight *in the wedge between the tips*. The third
    edge must attach to a distinct base vertex of valence >=3. Both
    base sectors adjacent to the stem must be zero so sector fusion
    requires no unsupported arithmetic move.

    Presence does NOT establish all hypotheses of R±/canonicality.
    """
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labeled tree required")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("sign must be ±1")
    tree = source.tree
    occ: dict[str, list[int]] = {}
    for i, vertex in enumerate(tree.vertices):
        for bond in vertex.slots:
            occ.setdefault(bond, []).append(i)
    out = []
    for motif in ring_signatures(source):
        if motif.sign != sign or not motif.stem_bond.startswith("e"):
            continue
        center = tree.vertices[motif.center]
        owners = occ[motif.stem_bond]
        if len(owners) != 2:
            continue
        base_idx = next(j for j in owners if j != motif.center)
        base = tree.vertices[base_idx]
        if len(base.slots) < 3:
            continue
        j = base.slots.index(motif.stem_bond)
        if base.sectors[j] or base.sectors[(j - 1) % len(base.slots)]:
            continue
        i = center.slots.index(motif.stem_bond)
        sectors = center.sectors[i:] + center.sectors[:i]
        if sectors != (0, sign, 0):
            continue
        out.append(RingCap(sign, base_idx, motif.center,
                           motif.leaf_vertices, motif.stem_bond))
    return tuple(sorted(out, key=lambda c: (c.base, c.center, c.tips)))


def abbreviate_disjoint_rings(source: V4LabeledTree, sign: int) -> Abbreviated:
    """Encode disjoint, separately attached ring subtrees by ring numbers.

    No ring *move* or arbitrary deletion equivalence is being applied.
    Reject overlapping ring candidates, empty surviving bases, and
    angular configurations outside the strict source-picture subset.
    """
    caps = strict_ring_caps(source, sign)
    tree = source.tree
    deleted = set()
    for cap in caps:
        region = {cap.center, *cap.tips}
        if region & deleted:
            raise ValueError("overlapping ring caps require exceptional handling")
        deleted.update(region)
    if any(cap.base in deleted for cap in caps):
        raise ValueError("rings attached to rings need exceptional handling")

    slots = [list(v.slots) for v in tree.vertices]
    sectors = [list(v.sectors) for v in tree.vertices]
    ring_counts = [0] * len(slots)
    for cap in caps:
        i = cap.base
        index = slots[i].index(cap.stem)
        if len(slots[i]) < 2:
            raise ValueError("abbreviation cannot leave an empty base")
        before = (index - 1) % len(slots[i])
        if sectors[i][before] != 0 or sectors[i][index] != 0:
            raise ValueError("sector split not supported")
        sectors[i][before] += sectors[i][index]
        del sectors[i][index]
        del slots[i][index]
        ring_counts[i] += 1
    remaining = [i for i in range(len(slots)) if i not in deleted]
    if any(not slots[i] for i in remaining):
        raise ValueError("abbreviation would leave a zero-valence vertex")
    weighted = WeightedPlaneTree(tuple(
        Vertex(tuple(slots[i]), tuple(sectors[i])) for i in remaining
    ))
    return Abbreviated(
        V4LabeledTree(weighted, source.labels),
        tuple(ring_counts[i] for i in remaining), sign,
    )
