"""Exact A3 ring-number grammar tests; all claims are syntactic."""
import unittest
from math import comb
from itertools import product

from arborescent.ring_abbreviation import Abbreviated, abbreviate_disjoint_rings
from arborescent.ring_grammar import (
    annotate_abbreviated, bounded_ring_counts, enumerate_abbreviated,
    ring_number_poly, total_bounded_rings, unbounded_fixed_ring_count,
)
from arborescent.weighted_trees import (
    Leaf, Fork, Vertex, WeightedPlaneTree, syntax_counts,
)
from arborescent.tangle_labels import V4LabeledTree


class RingGrammarTests(unittest.TestCase):
    def test_exact_coefficient_vectors_and_totals(self):
        self.assertEqual(ring_number_poly(0, 2), (1,))
        self.assertEqual(ring_number_poly(3, 1), (1, 3, 3, 1))
        self.assertEqual(ring_number_poly(2, 2), (1, 2, 3, 2, 1))
        self.assertEqual(bounded_ring_counts(0, 2, 1), (2,))
        self.assertEqual(bounded_ring_counts(1, 2, 1), (12, 12))
        self.assertEqual(bounded_ring_counts(2, 2, 1), (192, 384, 192))
        self.assertEqual(
            bounded_ring_counts(3, 2, 1),
            (4080, 12240, 12240, 4080),
        )
        self.assertEqual(bounded_ring_counts(2, 2, 2),
                         (192, 384, 576, 384, 192))
        self.assertEqual(total_bounded_rings(2, 2, 2), 1728)

    def test_summation_and_fixed_ring_coefficients_independent(self):
        for m in (1, 2, 3):
            a = syntax_counts(6, m)
            for R in range(4):
                for n in range(7):
                    with self.subTest(m=m, R=R, n=n):
                        b = bounded_ring_counts(n, m, R)
                        self.assertEqual(sum(b), a[n] * (R + 1)**n)
                        self.assertEqual(sum(b), total_bounded_rings(n, m, R))
                        # For 1-per-fork bound, coefficient = choose(n,q).
                        if R == 1:
                            self.assertEqual(
                                b, tuple(a[n]*comb(n, q) for q in range(n+1))
                            )
            for n in range(1, 7):
                for q in range(5):
                    self.assertEqual(
                        unbounded_fixed_ring_count(n, q, m),
                        a[n] * comb(n+q-1, q),
                    )
        self.assertEqual(unbounded_fixed_ring_count(0, 0, 2), 2)
        self.assertEqual(unbounded_fixed_ring_count(0, 1, 2), 0)

    def test_small_exhaustive_decorations_match_closed_count(self):
        for m, R, n in ((1, 2, 2), (2, 1, 1), (2, 1, 2), (3, 0, 2)):
            labels = tuple(range(2, m+2))
            arr = enumerate_abbreviated(n, labels, R)
            self.assertEqual(len(arr), sum(bounded_ring_counts(n, m, R)))
            keys = [a.rooted_code() for a in arr]
            self.assertEqual(len(set(keys)), len(keys))
            histogram = [sum(sum(a.rings) == q for a in arr)
                         for q in range(n*R+1)]
            self.assertEqual(tuple(histogram), bounded_ring_counts(n, m, R))

    def test_ring_numbers_really_attached_to_forks_not_leaves(self):
        spec = Fork((Leaf(2), Fork((Leaf(3), Leaf(2), Leaf(3)))))
        abbr = annotate_abbreviated(spec, (2, 1))
        self.assertEqual(abbr.rings, (2, 0, 1, 0, 0, 0))
        self.assertEqual(abbr.labeled.label("f0"), 0)
        self.assertEqual(abbr.sign, 1)
        self.assertNotEqual(abbr.rooted_code(),
                            annotate_abbreviated(spec, (1, 2)).rooted_code())

    def test_rooted_code_invariant_under_bond_vertex_renumbering(self):
        spec = Fork((Leaf(2), Fork((Leaf(3), Leaf(2)))))
        abbr = annotate_abbreviated(spec, (1, 2))
        old = abbr.labeled.tree
        vertices = tuple(
            Vertex(
                tuple("e" + str(100 + int(s[1:]))
                      if s.startswith("e") else s
                      for s in v.slots[1:] + v.slots[:1]),
                v.sectors[1:] + v.sectors[:1]
            )
            for v in reversed(old.vertices)
        )
        shifted = Abbreviated(
            V4LabeledTree(WeightedPlaneTree(vertices), abbr.labeled.labels),
            abbr.rings[::-1], +1,
        )
        self.assertEqual(abbr.rooted_code(), shifted.rooted_code())

    def test_two_distinct_ring_placements_have_identical_abbreviation(self):
        # Distinct *marked* weighted-plane presentations. The ring cap
        # attaches at a different wedge on a valence-4 base. After
        # §12.5 representation change, their ring counts and residual
        # cyclic child orders agree. This does NOT independently prove
        # the associated tangles are isotopic.
        rest = (
            Vertex(("e0", "e1", "e2"), (0, 1, 0)),
            Vertex(("e1",), (2,)),
            Vertex(("e2",), (2,)),
            Vertex(("e3",), (3,)),
            Vertex(("e4",), (4,)),
        )
        one = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e0", "e3", "e4"), (0, 0, 0, 0)),
            *rest,
        )), (("f0", 0),))
        two = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e3", "e0", "e4"), (0, 0, 0, 0)),
            *rest,
        )), (("f0", 0),))
        self.assertNotEqual(one.labelled_rooted_code(),
                            two.labelled_rooted_code())
        abbr_one = abbreviate_disjoint_rings(one, +1)
        abbr_two = abbreviate_disjoint_rings(two, +1)
        self.assertEqual(abbr_one.rings, (1, 0, 0))
        self.assertEqual(abbr_one.rooted_code(), abbr_two.rooted_code())

    def test_invalid_inputs_fail_closed(self):
        for fn, args in (
            (ring_number_poly, (-1, 1)),
            (ring_number_poly, (1, True)),
            (bounded_ring_counts, (1, 0, 1)),
            (bounded_ring_counts, (1, 1, -1)),
            (total_bounded_rings, (True, 2, 1)),
            (unbounded_fixed_ring_count, (1, -1, 2)),
            (unbounded_fixed_ring_count, (0, 0, 0)),
            (enumerate_abbreviated, (0,)),
            (annotate_abbreviated, (Leaf(2), ())),
            (annotate_abbreviated, (Fork((Leaf(2), Leaf(3))), (1, 2))),
            (annotate_abbreviated, (Fork((Leaf(2), Leaf(3))), (True,))),
        ):
            with self.subTest(fn=fn.__name__, args=args):
                with self.assertRaises(ValueError):
                    fn(*args)


if __name__ == "__main__":
    unittest.main()
