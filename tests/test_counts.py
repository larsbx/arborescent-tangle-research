"""Independent exact and brute-force checks; no knot isotopy claims."""

import unittest

from arborescent_tangles.counts import (
    SLOPES, bracelet_count, catalan, dihedral_key,
    enumerated_bracelet_count, ordered_word_count, raw_syntax_count,
)


class CombinatorialBaselineTests(unittest.TestCase):
    def test_catalan_numbers(self):
        self.assertEqual([catalan(n) for n in range(9)],
                         [1, 1, 2, 5, 14, 42, 132, 429, 1430])

    def test_parenthesization_collapse(self):
        for k in range(1, 10):
            self.assertEqual(raw_syntax_count(k, 3),
                             catalan(k - 1) * ordered_word_count(k, 3))
        self.assertEqual(raw_syntax_count(4, 1), 5)
        self.assertEqual(ordered_word_count(4, 1), 1)

    def test_burnside_against_direct_orbit_enumeration(self):
        for alphabet_size in (1, 2, 3, 4):
            for k in range(1, 8):
                with self.subTest(k=k, m=alphabet_size):
                    self.assertEqual(bracelet_count(k, alphabet_size),
                                     enumerated_bracelet_count(k, alphabet_size))

    def test_dihedral_invariance_and_non_equivalence(self):
        word = (0, 1, 2, 1, 2)
        key = dihedral_key(word)
        for shift in range(len(word)):
            self.assertEqual(dihedral_key(word[shift:] + word[:shift]), key)
            reversed_word = word[::-1]
            self.assertEqual(dihedral_key(reversed_word[shift:] + reversed_word[:shift]), key)
        self.assertNotEqual(dihedral_key((0, 0, 1)), dihedral_key((0, 1, 1)))

    def test_slopes_restricted_and_distinct(self):
        self.assertEqual(len(set(SLOPES)), 3)
        self.assertTrue(all(0 < q < 1 and q.denominator >= 5 for q in SLOPES))
        self.assertLessEqual(sum(1 / q.denominator for q in SLOPES), 1)

    def test_invalid_inputs(self):
        for fn in (raw_syntax_count, ordered_word_count, bracelet_count):
            for args in ((0, 3), (3, 0), (True, 3), (3, True)):
                with self.subTest(fn=fn.__name__, args=args):
                    with self.assertRaises(ValueError):
                        fn(*args)
        with self.assertRaises(ValueError):
            catalan(-1)
        with self.assertRaises(ValueError):
            dihedral_key(())


if __name__ == "__main__":
    unittest.main()
