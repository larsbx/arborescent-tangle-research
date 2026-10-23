"""A2.6: conservative transcription of the *pictured* R+/R- ring cap.

Bonahon--Siebenmann, The Plumbing Calculus, §12.4, printed p.155,
PDF zero-based page 166: the ring diagram deletes a THREE-VERTEX
cap (tip ±2 -- fork ±1 -- tip ±2) attached via its stem to a
SEPARATE base vertex. A stem that is merely a free bond on the
weighted cap vertex is NOT that pictured configuration.

The move implemented here is a strict, locally checkable subset:
- center has exactly 3 bonds, two to valence-one ±2 tips,
  one to a separate base vertex;
- the only nonzero center *sector* is the ±1 between its tip bonds;
- base has >=3 bonds and both sectors adjacent to the cap stem
  are zero, so removing the stem needs no unknown weight split;
- removed cap has no free-bond labels; every other V4 label persists.

The inverse inserts this explicitly pictured cap into a ZERO sector
between bonds at a base of bond-valence >=2. Bond IDs are chosen fresh.
No arbitrary ring move, exceptional canonicality or isotopy decision
procedure is claimed. The manuscript is the topology source, not an
independent diagrammatic validation of this Python implementation.
"""
from __future__ import annotations

from dataclasses import dataclass

from .ring_signatures import RingSignature, ring_signatures
from .tangle_labels import V4LabeledTree
from .weighted_trees import Vertex, WeightedPlaneTree


@dataclass(frozen=True)
class StrictRingCap:
    """A picture-matching ±(2,1,2) cap plus its *separate* base vertex."""
    sign: int
    base: int
    center: int
    leaves: tuple[int, int]
    stem_edge: str


def _owner_map(tree: WeightedPlaneTree) -> dict[str, tuple[int, ...]]:
    occ: dict[str, list[int]] = {}
    for i, v in enumerate(tree.vertices):
        for s in v.slots:
            occ.setdefault(s, []).append(i)
    return {s: tuple(v) for s, v in occ.items()}


def strict_ring_caps(source: V4LabeledTree) -> tuple[StrictRingCap, ...]:
    """Conservatively recognize the full pictured *cap attachment*.

    The permissive detector in ring_signatures.py does not test
    attachment, angular-sector placement, or base sectors.
    """
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled weighted tree required")
    tree = source.tree
    owners = _owner_map(tree)
    out: set[StrictRingCap] = set()
    for rough in ring_signatures(source):
        center = tree.vertices[rough.center]
        stem = rough.stem_bond
        if not stem.startswith("e"):
            continue  # a free stem is NOT a cap attached to a base
        other = owners[stem]
        if len(other) != 2 or rough.center not in other:
            continue
        base_index = next(i for i in other if i != rough.center)
        base = tree.vertices[base_index]
        if len(base.slots) < 3:
            continue
        idx = base.slots.index(stem)
        if base.sectors[(idx - 1) % len(base.slots)] != 0 or base.sectors[idx] != 0:
            continue
        # Rotate the center to put its base stem first. The tip weights
        # must then occur in the wedge between the two tip bonds.
        i = center.slots.index(stem)
        sectors = center.sectors[i:] + center.sectors[:i]
        slots = center.slots[i:] + center.slots[:i]
        if sectors != (0, rough.sign, 0):
            continue
        tips = []
        for tipbond in slots[1:]:
            endpoints = owners[tipbond]
            if len(endpoints) != 2:
                break
            tip = next(j for j in endpoints if j != rough.center)
            if (len(tree.vertices[tip].slots) != 1
                    or tree.vertices[tip].weight != 2 * rough.sign):
                break
            tips.append(tip)
        if len(tips) != 2 or tuple(sorted(tips)) != rough.leaf_vertices:
            continue
        out.add(StrictRingCap(
            rough.sign, base_index, rough.center,
            tuple(sorted(tips)), stem,
        ))
    return tuple(sorted(out, key=lambda s: (s.base, s.center, s.sign)))


def contract_ring_cap(source: V4LabeledTree, cap: StrictRingCap) -> V4LabeledTree:
    """Formal R± cap deletion, preserving every surviving V4 free label.

    Deletes the two weighted ±2 tips, the weighted ±1 cap center and
    its attachment edge. At the separate base vertex the two ZERO
    angular sectors adjacent to the attachment fuse into one zero
    sector. All other local data are unchanged.
    """
    if cap not in strict_ring_caps(source):
        raise ValueError("cap does not satisfy strict pictured R± preconditions")
    vertices = list(source.tree.vertices)
    base = vertices[cap.base]
    i = base.slots.index(cap.stem_edge)
    previous = (i - 1) % len(base.slots)
    slots = tuple(s for j, s in enumerate(base.slots) if j != i)
    sectors = tuple(
        w + (base.sectors[i] if j == previous else 0)
        for j, w in enumerate(base.sectors) if j != i
    )
    # Special case: when the cap stem is at index zero, the previous
    # slot is the last original slot; merging zero sectors leaves
    # sectors on surviving slots (all adjacent values are zero).
    vertices[cap.base] = Vertex(slots, sectors)
    deleted = {cap.center, *cap.leaves}
    survivors = tuple(v for j, v in enumerate(vertices) if j not in deleted)
    return V4LabeledTree(WeightedPlaneTree(survivors), source.labels)


def expand_ring_cap(
    source: V4LabeledTree,
    base_index: int,
    after_bond: str,
    sign: int,
) -> V4LabeledTree:
    """Inverse formal diagram: insert ±2--±1--±2 on a zero sector.

    Exact inverse of strict cap contraction *up to bond renaming and
    marked plane-presentation isomorphism*; no boundary labels change.
    """
    if not isinstance(source, V4LabeledTree):
        raise ValueError("V4-labelled tree required")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("ring sign must be +1 or -1")
    if type(base_index) is not int or not 0 <= base_index < len(source.tree.vertices):
        raise ValueError("invalid base vertex")
    vertices = list(source.tree.vertices)
    base = vertices[base_index]
    if len(base.slots) < 2 or after_bond not in base.slots:
        raise ValueError("base requires two prior incident bonds and named zero sector")
    k = base.slots.index(after_bond)
    if base.sectors[k] != 0:
        raise ValueError("new cap must occupy zero angular sector")
    occupied = {s for v in vertices for s in v.slots}
    used = {int(s[1:]) for s in occupied if s.startswith("e")}
    candidate = max(used, default=-1) + 1
    news = []
    while len(news) < 3:
        bond = f"e{candidate}"
        candidate += 1
        if bond not in occupied:
            news.append(bond)
    stem, left, right = news
    vertices[base_index] = Vertex(
        base.slots[:k + 1] + (stem,) + base.slots[k + 1:],
        base.sectors[:k + 1] + (0,) + base.sectors[k + 1:],
    )
    vertices.extend((
        Vertex((stem, left, right), (0, sign, 0)),
        Vertex((left,), (2 * sign,)),
        Vertex((right,), (2 * sign,)),
    ))
    return V4LabeledTree(WeightedPlaneTree(tuple(vertices)), source.labels)
