"""Catalan syntax versus dihedral classes of Montesinos *link parameters*.

No diagram/isotopy computation is performed here. Interpretation of the
parameter-orbit count as a link count requires the hypotheses of the
Montesinos classification theorem (see docs/montesinos-prototype.md).
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
from math import comb, gcd

# Proper slopes; denominators >= 5 put k>=3 in the usual non-exceptional
# Montesinos classification range for this fixed experiment.
SLOPES = (Fraction(1, 5), Fraction(2, 5), Fraction(1, 7))


def _check(k: int, alphabet_size: int) -> None:
    if type(k) is not int or k < 1:
        raise ValueError("k must be a positive integer (not a bool)")
    if type(alphabet_size) is not int or alphabet_size < 1:
        raise ValueError("alphabet_size must be a positive integer (not a bool)")


def catalan(n: int) -> int:
    if type(n) is not int or n < 0:
        raise ValueError("n must be a nonnegative integer")
    return comb(2 * n, n) // (n + 1)


def raw_syntax_count(k: int, alphabet_size: int) -> int:
    """Binary parenthesizations of k ordered, independently labeled leaves."""
    _check(k, alphabet_size)
    return catalan(k - 1) * alphabet_size**k


def ordered_word_count(k: int, alphabet_size: int) -> int:
    """Associativity quotient retaining all ordered leaf labels."""
    _check(k, alphabet_size)
    return alphabet_size**k


def bracelet_count(k: int, alphabet_size: int) -> int:
    """Burnside count under rotations and reversal, not a tangle isotopy test."""
    _check(k, alphabet_size)
    rotations = sum(alphabet_size ** gcd(k, j) for j in range(k))
    if k % 2:
        reflections = k * alphabet_size ** ((k + 1) // 2)
    else:
        reflections = (k // 2) * (
            alphabet_size ** (k // 2 + 1) + alphabet_size ** (k // 2)
        )
    numerator = rotations + reflections
    assert numerator % (2 * k) == 0
    return numerator // (2 * k)


def dihedral_key(word: tuple[int, ...]) -> tuple[int, ...]:
    """Canonical color word up to D_k, NOT a tangle normal form."""
    if not word:
        raise ValueError("word must be nonempty")
    seq = tuple(word)
    reversed_seq = seq[::-1]
    return min(
        candidate[i:] + candidate[:i]
        for candidate in (seq, reversed_seq)
        for i in range(len(seq))
    )


def enumerated_bracelet_count(k: int, alphabet_size: int) -> int:
    """Slow independent enumeration for small test instances."""
    _check(k, alphabet_size)
    return len({dihedral_key(word) for word in product(range(alphabet_size), repeat=k)})


def rows(max_k: int, alphabet_size: int) -> list[tuple[int, int, int, int]]:
    _check(max_k, alphabet_size)
    return [
        (k, raw_syntax_count(k, alphabet_size),
         ordered_word_count(k, alphabet_size), bracelet_count(k, alphabet_size))
        for k in range(1, max_k + 1)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-k", type=int, default=10)
    parser.add_argument("--alphabet-size", type=int, default=3)
    args = parser.parse_args()
    _check(args.max_k, args.alphabet_size)
    print("k,raw_binary_syntax,ordered_words,dihedral_orbits")
    for line in rows(args.max_k, args.alphabet_size):
        print(",".join(map(str, line)))


if __name__ == "__main__":
    main()
