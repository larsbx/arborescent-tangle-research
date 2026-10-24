"""Exact placement fibers of §12.5 *abbreviation* (NOT R± isotopy).

A star has a base with d>=3 free bonds, one marked f0, and r identical
three-vertex ring caps attached in any of the d gaps among free bonds.
All cyclic sectors at the base are zero; cap signs and wedge weights
follow the restricted source-picture recognizer.

After cyclic rooting at f0, cap placements are words with d-1 identical
free ports and r identical ring caps; the number is binom(d+r-1,r).
Every such presentation has the *same ring-numbered abbreviation*:
one zero-sector base with d free bonds, ring number r, and V4 identities.

This is syntactic ambiguity. It neither proves a topological R±
equivalence nor applies any ambient isotopy/classification theorem.
"""
from __future__ import annotations

from itertools import combinations
from math import comb

from .ring_abbreviation import abbreviate_disjoint_rings, strict_ring_caps
from .tangle_labels import V4LabeledTree
from .weighted_trees import Vertex, WeightedPlaneTree


def placement_fiber_formula(free_bonds: int, rings: int) -> int:
    """Exact rooted cyclic-interleaving count on the corolla family."""
    if (type(free_bonds) is not int or free_bonds < 3
            or type(rings) is not int or rings < 0):
        raise ValueError("free_bonds>=3 and rings>=0 integers required")
    return comb(free_bonds + rings - 1, rings)


def placement_presentations(
    free_bonds: int, rings: int, sign: int = 1
) -> tuple[V4LabeledTree, ...]:
    """Enumerate indexed source-pattern cap insertions, small sizes only.

    f0 anchors the base rotation. All other free legs have V4 label 0;
    ring caps are combinatorially indistinguishable after bond renaming.
    The root code is the exact object counted by the formula.
    """
    placement_fiber_formula(free_bonds, rings)
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("sign must be ±1")
    size = free_bonds - 1 + rings
    presentations = []
    for ring_positions in combinations(range(size), rings):
        chosen = set(ring_positions)
        slots = ["f0"]
        fidx, ridx = 1, 0
        for j in range(size):
            if j in chosen:
                slots.append(f"e{3*ridx}")
                ridx += 1
            else:
                slots.append(f"f{fidx}")
                fidx += 1
        verts = [Vertex(tuple(slots), (0,) * len(slots))]
        for j in range(rings):
            stem, tip1, tip2 = (f"e{3*j+k}" for k in range(3))
            verts.extend((
                Vertex((stem, tip1, tip2), (0, sign, 0)),
                Vertex((tip1,), (2 * sign,)),
                Vertex((tip2,), (2 * sign,)),
            ))
        labels = tuple((f"f{j}", 0) for j in range(free_bonds))
        presentations.append(V4LabeledTree(WeightedPlaneTree(tuple(verts)), labels))
    return tuple(presentations)


def placement_fiber_experiment(
    free_bonds: int, rings: int, sign: int = 1
) -> tuple[int, int, int]:
    """Return (indexed presentations, marked tree codes, abbreviations).

    This is exponential; use only small inputs. No topology quotient.
    """
    objs = placement_presentations(free_bonds, rings, sign)
    before = {obj.labelled_rooted_code() for obj in objs}
    after = set()
    for obj in objs:
        if len(strict_ring_caps(obj, sign)) != rings:
            raise AssertionError("fixture is outside strict source recognizer")
        abbr = abbreviate_disjoint_rings(obj, sign)
        if abbr.rings != (rings,):
            raise AssertionError("abbreviation multiplicity was lost")
        after.add(abbr.rooted_code())
    return len(objs), len(before), len(after)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=5)
    opts = parser.parse_args()
    if opts.max_n < 1:
        parser.error("max-n >= 1 required")
    print("n,free_bonds,ring_caps,placement_codes,abbreviated_codes,formula")
    for n in range(1, opts.max_n + 1):
        count, distinct, abbreviated = placement_fiber_experiment(n + 2, n)
        print(f"{n},{n+2},{n},{distinct},{abbreviated},"
              f"{placement_fiber_formula(n+2, n)}")


if __name__ == "__main__":
    main()



def bounded_valence_chain_presentations(
    bases: int, sign: int = 1
) -> tuple[V4LabeledTree, ...]:
    """2^k distinct placements; max weighted-tree bond valence is four.

    For k>=2 connect k zero-sector base vertices in a path. Each base
    has one marked/V4-zero free port and one identical signed ring cap.
    The ring cap attaches in either of two local angular positions,
    while the chain and free-port order is kept fixed. This produces
    two independent source-pattern placements at *each* base.

    Every abbreviation gives the SAME ring-numbered path with one
    ring recorded at each base. This shows that even bounded valence
    does not give polynomial fibers for UNREDUCED placements.
    """
    from itertools import product

    if type(bases) is not int or bases < 2:
        raise ValueError("at least two path bases required")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("sign must be ±1")
    def spine(j: int) -> str:
        return f"e{3*bases+j}"
    presentations = []
    for bits in product((0, 1), repeat=bases):
        vertices: list[Vertex] = []
        for i in range(bases):
            cap = f"e{3*i}"
            free = f"f{i}"
            if i == 0:
                other = spine(0)
                slots = (free, other, cap) if bits[i] == 0 else (free, cap, other)
            elif i == bases-1:
                incoming = spine(i-1)
                slots = ((incoming, free, cap) if bits[i] == 0
                         else (incoming, cap, free))
            else:
                incoming, outgoing = spine(i-1), spine(i)
                slots = ((incoming, outgoing, free, cap) if bits[i] == 0
                         else (incoming, outgoing, cap, free))
            vertices.append(Vertex(slots, (0,)*len(slots)))
        for i in range(bases):
            stem, left, right = (f"e{3*i+j}" for j in range(3))
            vertices.extend((
                Vertex((stem, left, right), (0, sign, 0)),
                Vertex((left,), (2*sign,)),
                Vertex((right,), (2*sign,)),
            ))
        labels = tuple((f"f{i}", 0) for i in range(bases))
        presentations.append(V4LabeledTree(WeightedPlaneTree(tuple(vertices)), labels))
    return tuple(presentations)


def bounded_valence_fiber_experiment(bases: int, sign: int = 1
                                   ) -> tuple[int, int, int, int]:
    """Return (source count, unique sources, abbreviation codes, max valence)."""
    source = bounded_valence_chain_presentations(bases, sign)
    originals = {t.labelled_rooted_code() for t in source}
    reduced = set()
    highest_valence = 0
    for t in source:
        highest_valence = max(highest_valence, *(len(v.slots) for v in t.tree.vertices))
        if len(strict_ring_caps(t, sign)) != bases:
            raise AssertionError("all k ring caps should be eligible")
        abbr = abbreviate_disjoint_rings(t, sign)
        if abbr.rings != (1,)*bases:
            raise AssertionError("lost a ring number or attachment base")
        reduced.add(abbr.rooted_code())
    return len(source), len(originals), len(reduced), highest_valence
