"""Independent regressions for §12.7 V4-labelled tangle F1/F2 subgroup."""
import unittest

from arborescent.flips import flip_quotient_count
from arborescent.tangle_labels import (
    IDENTITY, XI, ETA, ZETA, V4LabeledTree,
    label_free_ports, normalize_root_label, tangle_flip_f1,
    tangle_flip_f2, tangle_flip_group, tangle_flip_orbit_key,
    tangle_grammar_orbit_count, tangle_flip_f3prime,
)
from arborescent.weighted_trees import (
    Fork, Leaf, Vertex, WeightedPlaneTree, syntax_counts, to_weighted_tree,
)


class TangleLabelTests(unittest.TestCase):
    def test_klein_group_and_root_labels_all_four(self):
        self.assertEqual(XI ^ ETA, ZETA)
        spec = Fork((Leaf(2), Fork((Leaf(3), Leaf(2))), Leaf(2)))
        original = label_free_ports(to_weighted_tree(spec))
        orbit = tangle_flip_group(original)
        self.assertEqual({v.label("f0") for v in orbit}, {0, 1, 2, 3})
        self.assertEqual([v.label("f0") for v in orbit],
                         [IDENTITY, ZETA, XI, ETA])
        self.assertTrue(all(
            tangle_flip_orbit_key(v) == tangle_flip_orbit_key(original)
            for v in orbit
        ))

    def test_tangle_flip_involutions_and_composition(self):
        spec = Fork((Leaf(2), Fork((Leaf(2), Leaf(3)), Leaf(3))))
        obj = label_free_ports(to_weighted_tree(spec), ETA)
        self.assertEqual(tangle_flip_f1(tangle_flip_f1(obj)), obj)
        self.assertEqual(tangle_flip_f2(tangle_flip_f2(obj, 0), 0), obj)
        self.assertEqual(tangle_flip_f1(tangle_flip_f2(obj, 0)),
                         tangle_flip_f2(tangle_flip_f1(obj), 0))
        even = tangle_flip_f2(obj, 0)
        odd = tangle_flip_f2(obj, 1)
        self.assertEqual(tangle_flip_f2(even, 1), tangle_flip_f1(obj))
        self.assertEqual(tangle_flip_f2(odd, 0), tangle_flip_f1(obj))

    def test_nontrivial_multiple_free_bond_labels(self):
        # General weighted plane tree with two free boundary components.
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "f1"), (0, 0, 0)),
            Vertex(("e0", "e1", "e2"), (0, 0, 0)),
            Vertex(("e1",), (2,)), Vertex(("e2",), (3,)),
        ))
        t = V4LabeledTree(tree, (("f0", XI), ("f1", ETA)))
        for transformed in tangle_flip_group(t):
            self.assertEqual(tangle_flip_orbit_key(transformed),
                             tangle_flip_orbit_key(t))
        # Both free ports attached to the same vertex: F1 applies zeta
        # to both, F2-even applies xi to both.
        self.assertEqual(tangle_flip_f1(t).label("f1"), ETA ^ ZETA)
        self.assertEqual(tangle_flip_f2(t, 0).label("f1"), ETA ^ XI)
        # Transforming all labels independently is NOT an allowed F1/F2 move.
        wrong = V4LabeledTree(tree, (("f0", XI), ("f1", XI)))
        self.assertNotEqual(tangle_flip_orbit_key(t),
                            tangle_flip_orbit_key(wrong))

    def test_invariance_under_bond_and_vertex_renumbering(self):
        old = to_weighted_tree(Fork((Leaf(2), Leaf(3), Fork((Leaf(3), Leaf(2))))))
        source = V4LabeledTree(old, (("f0", ZETA),))
        new = WeightedPlaneTree(tuple(
            Vertex(
                tuple("e" + str(int(s[1:]) + 101) if s.startswith("e") else s
                      for s in v.slots[1:] + v.slots[:1]),
                v.sectors[1:] + v.sectors[:1],
            ) for v in reversed(old.vertices)
        ))
        shifted = V4LabeledTree(new, (("f0", ZETA),))
        self.assertEqual(source.labelled_rooted_code(),
                         shifted.labelled_rooted_code())
        self.assertEqual(tangle_flip_orbit_key(source),
                         tangle_flip_orbit_key(shifted))

    def test_gauge_slice_for_single_vertex_with_free_port(self):
        tree = WeightedPlaneTree((Vertex(("f0",), (2,)),))
        source = V4LabeledTree(tree, (("f0", ETA),))
        self.assertEqual(len(tangle_flip_group(source)), 4)
        self.assertEqual(len({
            g.label("f0") for g in tangle_flip_group(source)
        }), 4)
        fixed = normalize_root_label(source)
        self.assertEqual(fixed.label("f0"), 0)
        self.assertEqual(tangle_flip_orbit_key(fixed),
                         tangle_flip_orbit_key(source))

    def test_restricted_f3prime_parity_and_inverse(self):
        # Two free ports in the moved branch witness respectively:
        # distance 1 => eta, distance 2 => xi. f0 is outside.
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (3, 0, 0)),
            Vertex(("e0", "e2", "f1"), (0, 0, 0)),
            Vertex(("e1",), (2,)),
            Vertex(("e2", "f2"), (0, 0)),
        ))
        t = V4LabeledTree(
            tree, (("f0", XI), ("f1", ETA), ("f2", ZETA))
        )
        u = tangle_flip_f3prime(t, 0, "e0", from_side="before")
        self.assertEqual(u.tree.vertices[0].sectors, (0, 3, 0))
        self.assertEqual(u.label("f0"), XI)
        self.assertEqual(u.label("f1"), ETA ^ ETA)
        self.assertEqual(u.label("f2"), ZETA ^ XI)
        self.assertEqual(u.tree.vertices[1].slots,
                         ("e0", "f1", "e2"))
        self.assertEqual(u.tree.vertices[3].slots,
                         tree.vertices[3].slots)
        self.assertEqual(
            tangle_flip_f3prime(u, 0, "e0", from_side="after"), t
        )
        self.assertEqual(
            tuple(sorted(v.weight for v in t.tree.vertices)),
            tuple(sorted(v.weight for v in u.tree.vertices)),
        )

    def test_restricted_f3prime_even_no_flip(self):
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (4, 0, 0)),
            Vertex(("e0", "f1"), (0, 0)),
            Vertex(("e1",), (2,)),
        ))
        t = V4LabeledTree(tree, (("f0", 0), ("f1", 1)))
        u = tangle_flip_f3prime(t, 0, "e0")
        self.assertEqual(u.tree.vertices[0].sectors, (0, 4, 0))
        self.assertEqual(u.tree.vertices[1], t.tree.vertices[1])
        self.assertEqual(u.labels, t.labels)
        self.assertEqual(tangle_flip_f3prime(u, 0, "e0", from_side="after"), t)

    def test_restricted_f3prime_preconditions(self):
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (0, 3, 2)),
            Vertex(("e0",), (2,)), Vertex(("e1",), (3,)),
        ))
        t = label_free_ports(tree)
        for pivot, edge, side in (
            (-1, "e0", "before"), (True, "e0", "before"),
            (0, "e1", "before"), (0, "e0", "before"),
            (0, "f0", "before"), (0, "e99", "before"),
            (1, "e0", "middle"),
        ):
            with self.subTest(pivot=pivot, edge=edge, side=side):
                with self.assertRaises(ValueError):
                    tangle_flip_f3prime(t, pivot, edge, from_side=side)

    def test_label_and_root_invalid_guards(self):
        tree = to_weighted_tree(Fork((Leaf(2), Leaf(3))))
        for labels in (
            (), (("f0", True),), (("f0", -1),), (("f0", 4),),
            (("f1", 0),), (("f0", 0), ("f0", 1)),
        ):
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                V4LabeledTree(tree, labels)
        with self.assertRaises(ValueError):
            label_free_ports(tree, True)
        source = label_free_ports(tree)
        with self.assertRaises(ValueError):
            tangle_flip_group(source, "e0")
        with self.assertRaises(ValueError):
            normalize_root_label(source, "f123")
        with self.assertRaises(ValueError):
            source.labelled_rooted_code("f1")
        with self.assertRaises(ValueError):
            tangle_grammar_orbit_count(0)

    def test_exact_tangle_quotient_differs_from_unlabelled_pair_quotient(self):
        # For *one labelled root port* the four group elements act freely
        # on its boundary label. Orbit size = 4, so 4*raw/4 = raw.
        for m in (1, 2, 3):
            labels = tuple(range(2, 2 + m))
            for n in (1, 2):
                with self.subTest(n=n, m=m):
                    raw = syntax_counts(n, m)[n]
                    self.assertEqual(tangle_grammar_orbit_count(n, labels), raw)
                    if m > 1:
                        self.assertLess(flip_quotient_count(n, labels), raw)
        # The missing boundary label in A2.1 caused non-applicability of
        # the unlabelled quotient to the tangle classification theorem.
        self.assertEqual(flip_quotient_count(2, (2, 3)), 81)
        self.assertEqual(tangle_grammar_orbit_count(2, (2, 3)), 192)


if __name__ == "__main__":
    unittest.main()
