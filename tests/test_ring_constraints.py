"""A3.1 tests: finite local compatibility is not independent decoration."""
import unittest

from arborescent.ring_constraints import (
    brute_guarded_counts, control_ratio, guarded_ring_counts,
    no_adjacent_ring_forks,
)
from arborescent.ring_grammar import (
    annotate_abbreviated, bounded_ring_counts, total_bounded_rings,
)
from arborescent.weighted_trees import Fork, Leaf, syntax_counts


class GuardedRingGrammarTests(unittest.TestCase):
    def test_known_small_vectors(self):
        self.assertEqual(
            guarded_ring_counts(4, 2, 1),
            ((2,), (12, 12), (192, 384, 0),
             (4080, 12240, 4080, 0),
             (99264, 397056, 297792, 17856, 0)),
        )
        self.assertEqual(
            guarded_ring_counts(3, 2, 2),
            ((2,), (12, 12, 12), (192, 384, 384, 0, 0),
             (4080, 12240, 16320, 8160, 4080, 0, 0)),
        )

    def test_exact_brute_force_independent_from_dp(self):
        for m, bound, upto in ((1, 2, 3), (2, 1, 3),
                               (3, 0, 2), (2, 2, 2)):
            w = tuple(range(2, 2+m))
            counts = guarded_ring_counts(upto, m, bound)
            for n in range(1, upto+1):
                with self.subTest(m=m, R=bound, n=n):
                    self.assertEqual(
                        brute_guarded_counts(n, w, bound), counts[n]
                    )

    def test_no_rings_reproduces_a1_counting(self):
        for m in (1, 2, 3):
            a = syntax_counts(6, m)
            self.assertEqual(
                guarded_ring_counts(6, m, 0),
                tuple((v,) for v in a),
            )
            for n in range(7):
                self.assertEqual(control_ratio(n, m, 0), (a[n], a[n]))

    def test_guard_strictly_reduces_independent_grammar(self):
        for m in (1, 2, 3):
            for R in (1, 2):
                for n in (2, 3, 4):
                    with self.subTest(m=m, R=R, n=n):
                        restricted, raw = control_ratio(n, m, R)
                        self.assertLess(restricted, raw)
                        self.assertEqual(raw, sum(bounded_ring_counts(n, m, R)))
                        self.assertEqual(guarded_ring_counts(n, m, R)[n][0],
                                         syntax_counts(n, m)[n])
        self.assertEqual(control_ratio(2, 2, 1), (576, 768))

    def test_local_edge_guard_separates_parent_child(self):
        chain = Fork((Fork((Leaf(2), Leaf(3))), Leaf(2)))
        # Rings assigned to outer then inner fork in preorder.
        violation = annotate_abbreviated(chain, (1, 1))
        allowed_a = annotate_abbreviated(chain, (1, 0))
        allowed_b = annotate_abbreviated(chain, (0, 1))
        self.assertFalse(no_adjacent_ring_forks(violation))
        self.assertTrue(no_adjacent_ring_forks(allowed_a))
        self.assertTrue(no_adjacent_ring_forks(allowed_b))
        with self.assertRaises(ValueError):
            no_adjacent_ring_forks(None)

    def test_fixed_one_ring_count_unchanged_and_two_ring_count_reduced(self):
        # One ring cannot violate an adjacency restriction.
        for m in (1, 2, 3):
            n = 4
            guarded = guarded_ring_counts(n, m, 1)[n]
            independent = bounded_ring_counts(n, m, 1)
            self.assertEqual(guarded[0], independent[0])
            self.assertEqual(guarded[1], independent[1])
            self.assertLess(guarded[2], independent[2])

    def test_invalid_parameters(self):
        cases = ((-1, 2, 1), (1, 0, 1), (1, 2, -1),
                 (True, 2, 1), (1, True, 1), (1, 2, False))
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                guarded_ring_counts(*args)
        with self.assertRaises(ValueError):
            brute_guarded_counts(0)
        with self.assertRaises(ValueError):
            control_ratio(True, 2, 1)


if __name__ == "__main__":
    unittest.main()
