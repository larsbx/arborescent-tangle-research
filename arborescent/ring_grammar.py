"""A3: exact grammar of *abbreviated* weighted plane-tree syntax with rings.

Bonahon--Siebenmann §12.5 associates a nonnegative ring number to
vertices of a weighted plane tree. This module puts bounded ring numbers
only on 2/3-ary fork vertices of the prior A1 grammar. It counts those
syntactic decorations, NOT canonical arborescent tangles or isotopy
classes; topology-valid ring moves and exceptional canonicity remain
outside the model.

For finite R, F=m+x(1+y+...+y^R)(F^2+F^3) where x counts forks
and y counts abbreviated rings. At fixed fork size n, ring-number
decoration is independent on the n forks, so the bivariate
coefficient is a_n(m)*[y^q](1+...+y^R)^n.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterator

from .ring_abbreviation import Abbreviated
from .tangle_labels import label_free_ports
from .weighted_trees import Fork, Leaf, Spec, enumerate_specs, syntax_counts, to_weighted_tree


def _nonnegative(value: int, label: str) -> None:
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")


def ring_number_poly(n: int, maximum_per_fork: int) -> tuple[int, ...]:
    """Coefficients of (1+y+...+y^R)^n, exact and independent of A1."""
    _nonnegative(n, "number of forks")
    _nonnegative(maximum_per_fork, "maximum ring number")
    coefficients = [1]
    for _ in range(n):
        updated = [0] * (len(coefficients) + maximum_per_fork)
        for i, old in enumerate(coefficients):
            for q in range(maximum_per_fork + 1):
                updated[i + q] += old
        coefficients = updated
    return tuple(coefficients)


def bounded_ring_counts(n: int, m: int, maximum_per_fork: int) -> tuple[int, ...]:
    """Exact b_(n,q) for each q=0..nR, where n is the fork count."""
    _nonnegative(n, "number of forks")
    _nonnegative(maximum_per_fork, "maximum ring number")
    if type(m) is not int or m < 1:
        raise ValueError("leaf alphabet size must be positive")
    a_n = syntax_counts(n, m)[n]
    return tuple(a_n * count for count in ring_number_poly(n, maximum_per_fork))


def total_bounded_rings(n: int, m: int, maximum_per_fork: int) -> int:
    """Sum over ring counts: a_n(m)*(R+1)^n."""
    _nonnegative(n, "number of forks")
    _nonnegative(maximum_per_fork, "maximum ring number")
    if type(m) is not int or m < 1:
        raise ValueError("leaf alphabet size must be positive")
    return syntax_counts(n, m)[n] * (maximum_per_fork + 1) ** n


def unbounded_fixed_ring_count(n: int, q: int, m: int) -> int:
    """a_n(m)*C(n+q-1,q) for n>=1; formal leaf base handled separately."""
    from math import comb
    _nonnegative(n, "number of forks")
    _nonnegative(q, "number of rings")
    if type(m) is not int or m < 1:
        raise ValueError("leaf alphabet size must be positive")
    if n == 0:
        return m if q == 0 else 0
    return syntax_counts(n, m)[n] * comb(n + q - 1, q)


def _fork_count(spec: Spec) -> int:
    if isinstance(spec, Leaf):
        return 0
    if not isinstance(spec, Fork):
        raise ValueError("expected A1 Leaf/Fork")
    return 1 + sum(_fork_count(child) for child in spec.children)


def _decorate(spec: Spec, numbers: Iterator[int], output: list[int]) -> None:
    """Match the preorder vertex numbering in to_weighted_tree."""
    if isinstance(spec, Leaf):
        output.append(0)
        return
    if not isinstance(spec, Fork):
        raise ValueError("expected A1 syntax node")
    try:
        count = next(numbers)
    except StopIteration as exc:
        raise ValueError("insufficient fork ring numbers") from exc
    _nonnegative(count, "fork ring number")
    output.append(count)
    for child in spec.children:
        _decorate(child, numbers, output)


def annotate_abbreviated(
    spec: Fork,
    fork_ring_numbers: tuple[int, ...],
    leaf_weights: tuple[int, ...] = (2, 3),
    sign: int = 1,
) -> Abbreviated:
    """Compile ordered fork ring numbers into a labelled abbreviated tree.

    Fork numbers are in preorder, all leaf numbers are zero. The marked
    root free bond is f0 and has identity V4 label. This is a *syntax*
    map; Bonahon--Siebenmann's canonicality hypotheses are not checked.
    """
    if not isinstance(spec, Fork):
        raise ValueError("root must be a fork")
    if not isinstance(fork_ring_numbers, tuple):
        raise ValueError("ring numbers must be a tuple")
    if len(fork_ring_numbers) != _fork_count(spec):
        raise ValueError("one ring number per fork required")
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError("sign must be ±1")
    numbers = iter(fork_ring_numbers)
    ring_list: list[int] = []
    _decorate(spec, numbers, ring_list)
    labeled = label_free_ports(to_weighted_tree(spec, leaf_weights))
    return Abbreviated(labeled, tuple(ring_list), sign)


def enumerate_abbreviated(
    n: int, leaf_weights: tuple[int, ...] = (2, 3),
    maximum_per_fork: int = 1, sign: int = 1,
) -> tuple[Abbreviated, ...]:
    """Exhaustive, intentionally exponential small-n presentation generator."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be at least one fork")
    _nonnegative(maximum_per_fork, "maximum ring number")
    if (not isinstance(leaf_weights, tuple) or len(set(leaf_weights)) != len(leaf_weights)
            or any(type(w) is not int or abs(w) < 2 for w in leaf_weights)):
        raise ValueError("distinct finite leaf weights with absolute value >=2 required")
    return tuple(
        annotate_abbreviated(spec, tuple(rings), leaf_weights, sign)
        for spec in enumerate_specs(n, leaf_weights)
        for rings in product(range(maximum_per_fork + 1), repeat=n)
    )


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=6)
    parser.add_argument("--leaf-labels", type=int, default=2)
    parser.add_argument("--max-rings-per-fork", type=int, default=1)
    args = parser.parse_args()
    print("forks,ring-count-polynomial-coefficients,total")
    for n in range(args.max_n + 1):
        vector = bounded_ring_counts(n, args.leaf_labels, args.max_rings_per_fork)
        print(f"{n},{':'.join(map(str, vector))},{sum(vector)}")


if __name__ == "__main__":
    main()
