"""Deterministic regression tests; exhaustive orbits independently check Burnside."""
import unittest
from fractions import Fraction

from arborescent.prototype import (
    bracelet_count, catalan, dihedral_normal_form, orbit_count,
    raw_binary_syntax_count, validate_slopes,
)


class PrototypeTests(unittest.TestCase):
    def test_catalan_numbers(self):
        self.assertEqual([catalan(n) for n in range(9)],
                         [1, 1, 2, 5, 14, 42, 132, 429, 1430])

    def test_four_leaf_associativity_collapses_five_trees(self):
        self.assertEqual(catalan(3), 5)
        self.assertEqual(raw_binary_syntax_count(4, 3), 5 * 3**4)

    def test_burnside_matches_independent_exhaustive_orbits(self):
        for m in (1, 2, 3):
            for k in range(1, 9):
                with self.subTest(m=m, k=k):
                    self.assertEqual(bracelet_count(k, m), orbit_count(k, m))

    def test_cyclic_and_reversal_normalization(self):
        word = (0, 1, 1, 2, 0)
        self.assertEqual(dihedral_normal_form(word),
                         dihedral_normal_form(word[2:] + word[:2]))
        self.assertEqual(dihedral_normal_form(word),
                         dihedral_normal_form(word[::-1]))

    def test_nontrivial_counts(self):
        self.assertEqual([bracelet_count(k, 2) for k in range(3, 9)],
                         [4, 6, 8, 13, 18, 30])

    def test_slopes(self):
        self.assertEqual(validate_slopes(("2/4", "2/3")),
                         (Fraction(1, 2), Fraction(2, 3)))
        for slopes in (("1/2", "2/4"), ("-1/3",), ("1",), (), ("x",), ("1/0",)):
            with self.subTest(slopes=slopes), self.assertRaises(ValueError):
                validate_slopes(slopes)

    def test_invalid_sizes(self):
        for func, args in ((catalan, (-1,)),
                           (bracelet_count, (0, 2)),
                           (bracelet_count, (3, 0)),
                           (orbit_count, (0, 3)),
                           (raw_binary_syntax_count, (0, 2)),
                           (raw_binary_syntax_count, (3, 0)),
                           (dihedral_normal_form, ((),))):
            with self.subTest(func=func.__name__), self.assertRaises(ValueError):
                func(*args)


if __name__ == "__main__":
    unittest.main()
