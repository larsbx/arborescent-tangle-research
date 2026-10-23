"""A2.4 exact pivot grammar, reversible strict F3, and orbit certificates."""
import unittest

from arborescent.root_pivot import (
    PivotSpec, brute_pivot_counts, compile_pivot, normalized_pivot,
    pivot_orbit_key, pivot_specs, pivot_syntax_counts, root_sector, slide_root,
)
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Fork, Leaf, Vertex, WeightedPlaneTree, syntax_counts


class WeightedRootPivotTests(unittest.TestCase):
    def test_formula_exact_coefficients(self):
        self.assertEqual(
            pivot_syntax_counts(6, 2, 2),
            (
                (0, 0), (88, 24), (1440, 384), (30816, 8160),
                (752256, 198528), (19853568, 5229312),
                (551725056, 145135104),
            ),
        )
        for m in (1, 2, 3):
            for p in (1, 2):
                f = syntax_counts(6, m)
                rows = pivot_syntax_counts(6, m, p)
                for n in range(1, 7):
                    # J = p*(F-m): coefficients agree identically.
                    self.assertEqual(rows[n][1], p * f[n])
                    self.assertGreaterEqual(rows[n][0], 3 * rows[n][1])
                    self.assertLessEqual(rows[n][0], 4 * rows[n][1])

    def test_brute_force_pivot_orbits_against_closed_formula(self):
        for leaves, roots in (((2,), (3,)), ((2, 3), (3, 4))):
            rows = pivot_syntax_counts(2, len(leaves), len(roots))
            for n in (1, 2):
                with self.subTest(leaves=leaves, roots=roots, n=n):
                    self.assertEqual(brute_pivot_counts(n, leaves, roots), rows[n])

    def test_root_weight_transports_and_inverse_odd_even(self):
        # A child with internal structure makes odd-weight F3 reflection
        # visibly change the bond rotation system in that subtree.
        rooted = Fork((
            Fork((Leaf(2), Leaf(3), Leaf(2))), Leaf(3), Leaf(2),
        ))
        for a in (3, 4):
            with self.subTest(a=a):
                baseline = compile_pivot(PivotSpec(rooted, a, 0))
                self.assertEqual(root_sector(baseline), 0)
                nxt = slide_root(baseline, +1)
                self.assertEqual(root_sector(nxt), 1)
                self.assertEqual(nxt.tree.vertices[0].weight, a)
                self.assertEqual(slide_root(nxt, -1), baseline)
                self.assertEqual(pivot_orbit_key(nxt), pivot_orbit_key(baseline))
                if a % 2:
                    self.assertNotEqual(nxt.tree.vertices[1].slots,
                                        baseline.tree.vertices[1].slots)
                else:
                    self.assertEqual(nxt.tree.vertices[1].slots,
                                     baseline.tree.vertices[1].slots)
                last = nxt
                while root_sector(last) < len(rooted.children):
                    last = slide_root(last, +1)
                self.assertEqual(root_sector(last), len(rooted.children))
                self.assertEqual(normalized_pivot(last), baseline)
                with self.assertRaises(ValueError):
                    slide_root(last, +1)
                with self.assertRaises(ValueError):
                    slide_root(baseline, -1)

    def test_every_sector_in_single_orbit_and_unique_root_zero_slice(self):
        for root in (
            Fork((Leaf(2), Leaf(3))),
            Fork((Leaf(2), Fork((Leaf(3), Leaf(2))), Leaf(3))),
        ):
            for a in (3, 4):
                for sector in range(len(root.children) + 1):
                    with self.subTest(root=root, weight=a, sector=sector):
                        state = compile_pivot(PivotSpec(root, a, sector))
                        norm = normalized_pivot(state)
                        self.assertEqual(root_sector(norm), 0)
                        # Every possible root-sector state is reachable
                        # reversibly from its unique zero-sector gauge.
                        s = norm
                        for _ in range(sector):
                            s = slide_root(s, +1)
                        self.assertEqual(s, state)
                        self.assertEqual(pivot_orbit_key(s), pivot_orbit_key(state))

    def test_w_f_contrast_and_out_of_a1_nonvacuity(self):
        spec = PivotSpec(Fork((Leaf(2), Leaf(3))), 3, 0)
        pivot = compile_pivot(spec)
        self.assertTrue(all(sum(w != 0 for w in v.sectors) <= 1
                            for v in pivot.tree.vertices))
        self.assertIn(("F", 0), pivot.tree.wfs_failures())
        # Pair free-bond (F) does NOT certify a tangle: it is not among
        # the arborescent-tangle canonical conditions in §12.7.2.
        self.assertNotEqual(slide_root(pivot, +1), pivot)

    def test_invalid_pivots_and_guarded_inputs(self):
        root = Fork((Leaf(2), Leaf(3)))
        for a, sector in ((0, 0), (1, 1), (True, 0),
                          (3, -1), (3, 3), (3, True)):
            with self.subTest(a=a, sector=sector), self.assertRaises(ValueError):
                PivotSpec(root, a, sector)
        for bad in (PivotSpec, None, 42):
            with self.assertRaises(ValueError):
                compile_pivot(bad)
        with self.assertRaises(ValueError):
            PivotSpec(Fork((Leaf(2),)), 3, 0)
        with self.assertRaises(ValueError):
            pivot_specs(0)
        with self.assertRaises(ValueError):
            pivot_specs(1, root_weights=(3, 3))
        with self.assertRaises(ValueError):
            pivot_syntax_counts(-1, 2, 2)
        with self.assertRaises(ValueError):
            pivot_syntax_counts(2, 2, True)
        t = compile_pivot(PivotSpec(root, 3, 0))
        for invalid in (0, 2, True):
            with self.assertRaises(ValueError):
                slide_root(t, invalid)

    def test_presentation_normal_form_refuses_invalid_root(self):
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (0, 0, 0)),
            Vertex(("e0",), (2,)),
            Vertex(("e1",), (3,)),
        ))
        source = V4LabeledTree(tree, (("f0", 0),))
        with self.assertRaises(ValueError):
            root_sector(source)


if __name__ == "__main__":
    unittest.main()
