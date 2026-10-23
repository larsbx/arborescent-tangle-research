"""A2.4: first non-vacuous grammar for the *strict implemented* F3-prime.

Presentations have a single nonzero-weight root fork (arity 2 or 3),
A1 zero-weight internal descendants and weighted terminal leaves.
A lone free root bond f0 is labelled by identity in V4.

This is a formal, bounded-local grammar of *presentations*, not a
Bonahon--Siebenmann canonical tangle grammar. In particular its W/S/P/N
canonicity and any ring-move normal forms are NOT established.

Strict F3-prime is applied ONLY across a root-to-child edge. It moves
the root's nonzero weight to an adjacent angular sector, and for odd
weight, reflects the affected subtree as specified in §12.7.3.
"""
from __future__ import annotations

from dataclasses import dataclass

from .tangle_labels import V4LabeledTree, tangle_flip_f3prime
from .weighted_trees import (
    Fork, Leaf, Spec, Vertex, WeightedPlaneTree,
    enumerate_specs, syntax_counts, to_weighted_tree,
)


@dataclass(frozen=True)
class PivotSpec:
    """A1 fork with one nonzero root sector and a fixed root twist."""
    root: Fork
    root_weight: int
    root_sector: int

    def __post_init__(self) -> None:
        if not isinstance(self.root, Fork) or len(self.root.children) not in (2, 3):
            raise ValueError("root must be a binary or ternary Fork")
        if type(self.root_weight) is not int or abs(self.root_weight) < 2:
            raise ValueError("root weight must be integer of absolute value >= 2")
        if (type(self.root_sector) is not int
                or not 0 <= self.root_sector <= len(self.root.children)):
            raise ValueError("invalid root sector")


def compile_pivot(spec: PivotSpec, leaf_weights: tuple[int, ...] = (2, 3)) -> V4LabeledTree:
    """Compile a PivotSpec to the preexisting sector/free-port model."""
    if not isinstance(spec, PivotSpec):
        raise ValueError("expected PivotSpec")
    base = to_weighted_tree(spec.root, leaf_weights)
    root = base.vertices[0]
    sector_weights = [0] * len(root.sectors)
    sector_weights[spec.root_sector] = spec.root_weight
    weighted = WeightedPlaneTree((
        Vertex(root.slots, tuple(sector_weights)),
        *base.vertices[1:],
    ))
    return V4LabeledTree(weighted, (("f0", 0),))


def root_sector(source: V4LabeledTree) -> int:
    """Return position of the unique nonzero root angular sector.

    Requires the marked root f0 to occur at the *first slot* of vertex
    zero and the other root incident bonds to be internal. Both
    requirements hold for all compiled PivotSpec presentations and
    are preserved by their strict root-F3-prime moves.
    """
    if not isinstance(source, V4LabeledTree):
        raise ValueError("expected V4-labelled tree")
    v = source.tree.vertices[0]
    if (v.slots[0] != "f0"
            or any(not s.startswith("e") for s in v.slots[1:])
            or len(v.slots) not in (3, 4)
            or set(s for vv in source.tree.vertices
                   for s in vv.slots if s.startswith("f")) != {"f0"}):
        raise ValueError("not a one-root-port pivot presentation")
    nonzero = [i for i, w in enumerate(v.sectors) if w != 0]
    if len(nonzero) != 1 or abs(v.sectors[nonzero[0]]) < 2:
        raise ValueError("expected one nonzero root sector of magnitude >= 2")
    return nonzero[0]


def slide_root(source: V4LabeledTree, step: int) -> V4LabeledTree:
    """Move root weight one sector left (-1) or right (+1).

    Raises on attempts to pass the free boundary bond f0.
    """
    if type(step) is not int or step not in (-1, 1):
        raise ValueError("step must be -1 or +1")
    pos = root_sector(source)
    last = len(source.tree.vertices[0].sectors) - 1
    if pos + step < 0 or pos + step > last:
        raise ValueError("cannot cross the free root bond")
    edge = source.tree.vertices[0].slots[pos if step == -1 else pos + 1]
    return tangle_flip_f3prime(
        source, 0, edge, from_side="after" if step == -1 else "before"
    )


def normalized_pivot(source: V4LabeledTree) -> V4LabeledTree:
    """Unique sector-zero representative under *root F3-prime only*.

    Other possible topological moves (F1/F2, R±, arithmetic moves) are
    expressly NOT quotiented by this function.
    """
    while root_sector(source) != 0:
        source = slide_root(source, -1)
    return source


def pivot_orbit_key(source: V4LabeledTree) -> tuple:
    """Presentation-isomorphism key after root-sector normalizing."""
    return normalized_pivot(source).labelled_rooted_code()


def pivot_specs(
    n: int,
    leaf_weights: tuple[int, ...] = (2, 3),
    root_weights: tuple[int, ...] = (3, 4),
) -> tuple[PivotSpec, ...]:
    """Generate all root-weight-position presentations for small n."""
    if type(n) is not int or n < 1:
        raise ValueError("n >= 1 required")
    if (not root_weights or len(set(root_weights)) != len(root_weights)
            or any(type(w) is not int or abs(w) < 2 for w in root_weights)):
        raise ValueError("root weights must be distinct integers |w| >= 2")
    return tuple(
        PivotSpec(root, a, sector)
        for root in enumerate_specs(n, leaf_weights) if isinstance(root, Fork)
        for a in root_weights
        for sector in range(len(root.children) + 1)
    )


def pivot_syntax_counts(max_n: int, leaf_alphabet_size: int, root_alphabet_size: int
                        ) -> tuple[tuple[int, int], ...]:
    """Exact (raw, root-F3-orbits) counts from algebraic F=m+x(F²+F³).

    H=p*x*(3F²+4F³), J=p*x*(F²+F³) = p*(F-m).
    The n=0 term is (0,0): the weighted root must branch.
    """
    if (any(type(v) is not int for v in (max_n, leaf_alphabet_size, root_alphabet_size))
            or max_n < 0 or leaf_alphabet_size < 1 or root_alphabet_size < 1):
        raise ValueError("invalid sizes")
    a = syntax_counts(max_n, leaf_alphabet_size)
    rows = [(0, 0)]
    for n in range(1, max_n + 1):
        pair = sum(a[i] * a[n - 1 - i] for i in range(n))
        triple = sum(a[i] * a[j] * a[n - 1 - i - j]
                     for i in range(n) for j in range(n - i))
        p = root_alphabet_size
        rows.append((p * (3 * pair + 4 * triple), p * (pair + triple)))
    return tuple(rows)


def brute_pivot_counts(n: int, leaf_weights: tuple[int, ...] = (2, 3),
                       root_weights: tuple[int, ...] = (3, 4)) -> tuple[int, int]:
    """Independent exhaustive enumeration + computed move-orbit keys."""
    specs = pivot_specs(n, leaf_weights, root_weights)
    raw_codes = [compile_pivot(s, leaf_weights) for s in specs]
    return (
        len({s.labelled_rooted_code() for s in raw_codes}),
        len({pivot_orbit_key(s) for s in raw_codes}),
    )


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=6)
    parser.add_argument("--leaf-labels", type=int, default=2)
    parser.add_argument("--root-labels", type=int, default=2)
    opts = parser.parse_args()
    print("n,raw_pivot_presentations,strict_root_F3_orbits")
    for n, (raw, orbits) in enumerate(
        pivot_syntax_counts(opts.max_n, opts.leaf_labels, opts.root_labels)
    ):
        print(f"{n},{raw},{orbits}")


if __name__ == "__main__":
    main()
