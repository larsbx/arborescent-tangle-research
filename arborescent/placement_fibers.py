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
