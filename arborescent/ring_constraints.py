"""A3.1: nontrivial bounded-local ring-number *syntax* grammar.

Experimental combinatorial restriction: no edge joins two fork vertices
whose ring numbers are both positive. This is a declared finite-state
guard, NOT Bonahon--Siebenmann's full adjacent-ring exception and
NOT a theorem about allowable canonical topological presentations.

Let Z(x,y) enumerate leaves and zero-ring-root forks; let P(x,y)
enumerate positive-ring-root forks. With F=Z+P and
D_+(y)=y+...+y^R:

    Z = m + x(F^2 + F^3)
    P = x D_+(y) (Z^2 + Z^3)
    F = Z + P.

This coupled polynomial grammar is algebraic over Q(x,y); its exact
coefficients are computed by integer convolution below.
"""
from __future__ import annotations

from collections import defaultdict

from .ring_grammar import enumerate_abbreviated, total_bounded_rings
from .ring_abbreviation import Abbreviated


def _convolve(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    out: dict[int, int] = defaultdict(int)
    for i, av in a.items():
        for j, bv in b.items():
            out[i+j] += av*bv
    return dict(out)


def _add(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    out: dict[int, int] = defaultdict(int, a)
    for j, bv in b.items():
        out[j] += bv
    return dict(out)


def guarded_ring_counts(max_n: int, m: int, max_rings: int
                        ) -> tuple[tuple[int, ...], ...]:
    """Exact [x^n y^q] counts, q=0..n*max_rings.

    Every ring-bearing *fork* is separated from ring-bearing child
    forks by at least one non-ring fork. Leaves carry no ring number.
    """
    if (any(type(x) is not int for x in (max_n, m, max_rings))
            or max_n < 0 or m < 1 or max_rings < 0):
        raise ValueError("max_n>=0, m>=1, max_rings>=0 required")
    z = [{0: m}]
    p = [{}]
    out = [(m,)]
    for n in range(1, max_n + 1):
        f = [_add(a, b) for a, b in zip(z, p)]
        bare_root: dict[int, int] = defaultdict(int)
        occupied_childless: dict[int, int] = defaultdict(int)
        for i in range(n):
            for q, val in _convolve(f[i], f[n-1-i]).items():
                bare_root[q] += val
            for q, val in _convolve(z[i], z[n-1-i]).items():
                occupied_childless[q] += val
        for i in range(n):
            for j in range(n-i):
                for q, val in _convolve(
                    _convolve(f[i], f[j]), f[n-1-i-j]
                ).items():
                    bare_root[q] += val
                for q, val in _convolve(
                    _convolve(z[i], z[j]), z[n-1-i-j]
                ).items():
                    occupied_childless[q] += val
        occupied: dict[int, int] = defaultdict(int)
        for q, val in occupied_childless.items():
            for ring_number in range(1, max_rings+1):
                occupied[q+ring_number] += val
        z.append(dict(bare_root))
        p.append(dict(occupied))
        total = _add(z[-1], p[-1])
        out.append(tuple(total.get(q, 0) for q in range(n*max_rings+1)))
    return tuple(out)


def no_adjacent_ring_forks(source: Abbreviated) -> bool:
    """Syntactic guard on explicit ring-numbered weighted trees."""
    if not isinstance(source, Abbreviated):
        raise ValueError("expected ring-numbered abbreviated presentation")
    return all(
        not (source.rings[i] > 0 and source.rings[j] > 0)
        for i, j in source.labeled.tree.edges()
    )


def brute_guarded_counts(n: int, leaf_weights: tuple[int, ...] = (2, 3),
                         max_rings: int = 1) -> tuple[int, ...]:
    """Independent exhaustive filtering for small n; exponential cost."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be positive")
    objects = enumerate_abbreviated(n, leaf_weights, max_rings)
    histogram = [0]*(n*max_rings+1)
    for obj in objects:
        if no_adjacent_ring_forks(obj):
            histogram[sum(obj.rings)] += 1
    return tuple(histogram)


def control_ratio(n: int, m: int, max_rings: int) -> tuple[int, int]:
    """(restricted total, independent-label total) at fixed fork count."""
    if type(n) is not int or n < 0:
        raise ValueError("n must be nonnegative")
    restricted = sum(guarded_ring_counts(n, m, max_rings)[n])
    unrestricted = total_bounded_rings(n, m, max_rings)
    return (restricted, unrestricted)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=5)
    parser.add_argument("--m", type=int, default=2)
    parser.add_argument("--max-rings", type=int, default=1)
    a = parser.parse_args()
    print("forks,admissible_ring_polynomial,total,independent_label_total")
    rows = guarded_ring_counts(a.max_n, a.m, a.max_rings)
    for n, row in enumerate(rows):
        print(f"{n},{':'.join(map(str,row))},{sum(row)},"
              f"{total_bounded_rings(n,a.m,a.max_rings)}")


if __name__ == "__main__":
    main()
