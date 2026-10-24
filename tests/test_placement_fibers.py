"""A3.2 exact fibers are presentation fibers, not isotopy fibers."""
import unittest
from math import comb

from arborescent.placement_fibers import (
    placement_fiber_experiment, placement_fiber_formula, placement_presentations,
)
from arborescent.ring_abbreviation import abbreviate_disjoint_rings
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import WeightedPlaneTree, Vertex


class RingPlacementFiberTests(unittest.TestCase):
    def test_formula_independent_exhaustive_marked_codes(self):
        for d in range(3, 7):
            for n in range(4):
                for sign in (1, -1):
                    with self.subTest(d=d, n=n, sign=sign):
                        expected = comb(d + n - 1, n)
                        self.assertEqual(placement_fiber_formula(d, n), expected)
                        raw, marked, abbreviations = placement_fiber_experiment(
                            d, n, sign
                        )
                        self.assertEqual((raw, marked, abbreviations),
                                         (expected, expected, 1))

    def test_first_nontrivial_examples(self):
        self.assertEqual(placement_fiber_experiment(3, 1), (3, 3, 1))
        self.assertEqual(placement_fiber_experiment(3, 2), (6, 6, 1))
        self.assertEqual(placement_fiber_experiment(4, 2), (10, 10, 1))
        self.assertEqual(placement_fiber_experiment(5, 3), (35, 35, 1))

    def test_multiring_counts_and_boundary_labels_persist(self):
        for d, n in ((3, 2), (5, 3)):
            for obj in placement_presentations(d, n, -1):
                reduced = abbreviate_disjoint_rings(obj, -1)
                self.assertEqual(reduced.rings, (n,))
                self.assertEqual(dict(reduced.labeled.labels),
                                 {f"f{i}": 0 for i in range(d)})
                self.assertEqual(len(obj.tree.vertices), 1 + 3*n)
                self.assertEqual(len(reduced.labeled.tree.vertices), 1)
                self.assertEqual(reduced.labeled.tree.vertices[0].weight, 0)

    def test_root_marking_matters_even_with_identical_free_labels(self):
        source = placement_presentations(3, 1)[0]
        code = source.labelled_rooted_code()
        # Different V4 on the *root* leg makes a different code
        # while keeping the source ring count unchanged.
        relabeled = V4LabeledTree(
            source.tree, tuple((b, 1 if b == "f0" else 0)
                               for b, _ in source.labels)
        )
        self.assertNotEqual(code, relabeled.labelled_rooted_code())
        self.assertNotEqual(
            abbreviate_disjoint_rings(source, 1).rooted_code(),
            abbreviate_disjoint_rings(relabeled, 1).rooted_code(),
        )

    def test_negative_controls_invalid_inputs(self):
        for d, r in ((2, 1), (0, 0), (3, -1), (True, 1),
                     (3, True), (3, 2.0)):
            with self.subTest(d=d, r=r), self.assertRaises(ValueError):
                placement_fiber_formula(d, r)
        for sign in (0, 2, True):
            with self.subTest(sign=sign), self.assertRaises(ValueError):
                placement_presentations(3, 1, sign)

    def test_exponential_ambiguity_under_growing_valence(self):
        # Exact identities: with d=n>=3, the preimage cardinality is
        # binom(2n-1,n) = 1/2*binom(2n,n); ~4^n/(2 sqrt(pi*n)).
        for n in range(3, 9):
            self.assertEqual(
                placement_fiber_formula(n, n), comb(2*n, n)//2
            )
        self.assertGreater(placement_fiber_formula(8, 8), 5000)


if __name__ == "__main__":
    unittest.main()
