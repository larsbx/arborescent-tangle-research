"""Presentation combinatorics, NOT a general tangle-isotopy classifier.

Full ordered binary bracketings -> ordered leaf words ->
rotation/reversal orbits of cyclic *link-closure* presentations.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
from math import comb


def catalan(n: int) -> int:
    """Number of plane full binary trees with n internal vertices."""
    if n < 0:
        raise ValueError("n must be nonnegative")
    return comb(2 * n, n) // (n + 1)


def raw_binary_syntax_count(k: int, alphabet_size: int) -> int:
    """Full binary parenthesizations of k leaves with independent labels."""
    if k < 1 or alphabet_size < 1:
        raise ValueError("k and alphabet_size must be positive")
    return catalan(k - 1) * alphabet_size**k


def validate_slopes(slopes: tuple[str, ...]) -> tuple[Fraction, ...]:
    """Canonical fractions in (0,1); exclude repeated rational-tangle slopes."""
    if not slopes:
        raise ValueError("provide at least one rational slope")
    try:
        fractions = tuple(Fraction(s) for s in slopes)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("slopes must be valid rational numbers") from exc
    if any(not 0 < x < 1 for x in fractions):
        raise ValueError("slopes must be strictly between 0 and 1")
    if len(set(fractions)) != len(fractions):
        raise ValueError("duplicate slopes after rational normalization")
    return fractions


def dihedral_normal_form(word: tuple[int, ...]) -> tuple[int, ...]:
    """Lexicographically minimal cyclic rotation or reversed rotation."""
    if not word:
        raise ValueError("word must be nonempty")
    k = len(word)
    forward = [word[i:] + word[:i] for i in range(k)]
    backwards = word[::-1]
    reverse_rotations = [backwards[i:] + backwards[:i] for i in range(k)]
    return min(forward + reverse_rotations)


def orbit_count(k: int, alphabet_size: int) -> int:
    """Independent exhaustive count for small sizes; exponential work."""
    if k < 1 or alphabet_size < 1:
        raise ValueError("k and alphabet_size must be positive")
    return len({dihedral_normal_form(word)
                for word in product(range(alphabet_size), repeat=k)})


def _totient(n: int) -> int:
    result = n
    factor = 2
    while factor * factor <= n:
        if n % factor == 0:
            while n % factor == 0:
                n //= factor
            result -= result // factor
        factor += 1
    if n > 1:
        result -= result // n
    return result


def bracelet_count(k: int, alphabet_size: int) -> int:
    """Burnside count of k-bead m-color bracelets (dihedral orbits)."""
    if k < 1 or alphabet_size < 1:
        raise ValueError("k and alphabet_size must be positive")
    rotations = sum(_totient(d) * alphabet_size ** (k // d)
                    for d in range(1, k + 1) if k % d == 0)
    if k % 2:
        reflections = k * alphabet_size ** ((k + 1) // 2)
    else:
        reflections = (k // 2) * (alphabet_size ** (k // 2 + 1)
                                  + alphabet_size ** (k // 2))
    numerator = rotations + reflections
    assert numerator % (2 * k) == 0
    return numerator // (2 * k)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slopes", nargs="+", default=["1/2", "1/3", "2/5"],
                        help="distinct rational slopes strictly between 0 and 1")
    parser.add_argument("--max-k", type=int, default=9)
    args = parser.parse_args()
    labels = validate_slopes(tuple(args.slopes))
    if args.max_k < 3:
        parser.error("max-k must be >= 3 for the Montesinos closure test regime")
    m = len(labels)
    print(f"rational leaf alphabet: {', '.join(str(s) for s in labels)}")
    print("k  Catalan  raw-binary  ordered-words  dihedral-orbits")
    for k in range(3, args.max_k + 1):
        print(f"{k:>2} {catalan(k-1):>8} {raw_binary_syntax_count(k,m):>11}"
              f" {m**k:>14} {bracelet_count(k,m):>16}")
    print("NOTE: dihedral-orbits count cyclic closure presentations, not"
          " boundary-fixed tangle isotopy classes.")


if __name__ == "__main__":
    main()
