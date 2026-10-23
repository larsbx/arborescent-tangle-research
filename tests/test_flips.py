"""A2.1: independently check precisely sourced F1/F2 subgroup actions.

The mirror-on-spec computation below is independent of the weighted
plane-rotation-system implementation and tests complete small orbits.
"""
import unittest

from arborescent.flips import (
    bipartition, flip_f1, flip_f2, flip_orbit_codes, flip_orbit_key,
    flip_quotient_count, flip_vertices, parity_flip_representatives,
    reverse_vertex,
)
from arborescent.weighted_trees import (
    Fork, Leaf, Vertex, WeightedPlaneTree, enumerate_specs, to_weighted_tree,
)


def mirror_spec(node, parity: int, depth: int = 0):
    """Independent ordered-tree flip of all nodes at one depth parity."""
    if isinstance(node, Leaf):
        return node
    children = tuple(mirror_spec(child, parity, depth + 1)
                     for child in node.children)
    if depth % 2 == parity:
        children = children[::-1]
    return Fork(children)


def spec_flip_key(spec):
    even = mirror_spec(spec, 0)
    odd = mirror_spec(spec, 1)
    both = mirror_spec(even, 1)
    return min(map(repr, (spec, even, odd, both)))


class FlipSubgroupTests(unittest.TestCase):
    def test_sector_convention_not_naive_slot_reversal(self):
        v = Vertex(("e0", "e1", "e2"), (2, 0, -3))
        mirrored = reverse_vertex(v)
        self.assertEqual(mirrored.slots, ("e0", "e2", "e1"))
        self.assertEqual(mirrored.sectors, (-3, 0, 2))
        self.assertEqual(reverse_vertex(mirrored), v)
        self.assertEqual(mirrored.weight, v.weight)

    def test_f1_f2_involutive_and_bipartition_identity(self):
        specimen = Fork((Leaf(2), Fork((Leaf(3), Leaf(2), Leaf(2))), Leaf(3)))
        tree = to_weighted_tree(specimen)
        even, odd = bipartition(tree, 0)
        self.assertFalse(even & odd)
        self.assertEqual(even | odd, frozenset(range(len(tree.vertices))))
        self.assertEqual(flip_f1(flip_f1(tree)), tree)
        self.assertEqual(flip_f2(flip_f2(tree, 0), 0), tree)
        self.assertEqual(flip_f2(flip_f2(tree, 2), 2), tree)
        # F1 is the composition of the independent two parity flips.
        self.assertEqual(flip_f2(flip_f2(tree, 0), 2), flip_f1(tree))
        self.assertEqual(flip_f2(flip_f2(tree, 2), 0), flip_f1(tree))
        self.assertEqual(len(parity_flip_representatives(tree)), 4)
        # This asymmetric specimen witnesses four *different* marked codes.
        self.assertEqual(len(flip_orbit_codes(tree)), 4)

    def test_weights_graph_wfs_preserved_by_both_flips(self):
        # W fails deliberately, to verify that failures are transported,
        # rather than spuriously cured by cyclic reflection.
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (0, 1, 2)),
            Vertex(("e0",), (2,)),
            Vertex(("e1",), (-3,)),
        ))
        before = tuple(v.weight for v in tree.vertices)
        failures = tree.wfs_failures()
        self.assertTrue(failures)
        for target in parity_flip_representatives(tree):
            self.assertEqual(tuple(v.weight for v in target.vertices), before)
            self.assertEqual(target.edges(), tree.edges())
            self.assertEqual(target.wfs_failures(), failures)

    def test_wrong_pivot_and_wrong_vertex_set_rejected(self):
        tree = to_weighted_tree(Fork((Leaf(2), Leaf(3))))
        for p in (-1, len(tree.vertices), True, "0"):
            with self.subTest(p=p), self.assertRaises(ValueError):
                bipartition(tree, p)
        for bad in (frozenset((-1,)), frozenset((len(tree.vertices),)), frozenset((True,))):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                flip_vertices(tree, bad)

    def test_normal_form_same_on_move_orbit(self):
        specimen = Fork((Leaf(2), Fork((Leaf(3), Leaf(2), Leaf(2))), Leaf(3)))
        tree = to_weighted_tree(specimen)
        key = flip_orbit_key(tree)
        for candidate in parity_flip_representatives(tree):
            self.assertEqual(flip_orbit_key(candidate), key)
        different = to_weighted_tree(Fork((Leaf(3), Fork((Leaf(3), Leaf(2), Leaf(2))), Leaf(3))))
        self.assertNotEqual(key, flip_orbit_key(different))

    def test_conjugacy_by_presentation_isomorphism(self):
        specimen = Fork((Leaf(2), Fork((Leaf(3), Leaf(2))), Leaf(3)))
        tree = to_weighted_tree(specimen)
        relabeled = WeightedPlaneTree(
            tuple(Vertex(
                tuple("e" + str(int(s[1:]) + 100) if s[0] == "e" else s
                      for s in v.slots[1:] + v.slots[:1]),
                v.sectors[1:] + v.sectors[:1],
            ) for v in tree.vertices[::-1])
        )
        self.assertEqual(flip_orbit_key(tree), flip_orbit_key(relabeled))

    def test_exact_small_orbits_vs_independent_syntax_reflection(self):
        known = {
            2: ((12, 9), (192, 81), (4080, 1251)),
            3: ((36, 24), (1188, 432)),  # m=3, n=3 is relatively large
        }
        for m, rows in known.items():
            labels = tuple(range(2, m + 2))
            for n, (raw, expected) in enumerate(rows, 1):
                with self.subTest(m=m, n=n):
                    specs = enumerate_specs(n, labels)
                    independent = len({spec_flip_key(s) for s in specs})
                    self.assertEqual(len(specs), raw)
                    self.assertEqual(independent, expected)
                    self.assertEqual(flip_quotient_count(n, labels), independent)

    def test_lower_bound_ninety_proportion_and_single_vertex(self):
        # Each orbit has at most four elements. A bound, not an
        # asymptotic claim for isotopy classes.
        for n in range(1, 4):
            raw = len(enumerate_specs(n))
            quotient = flip_quotient_count(n)
            self.assertLessEqual(quotient, raw)
            self.assertGreaterEqual(4 * quotient, raw)
        single = WeightedPlaneTree((Vertex(("f0",), (2,)),))
        self.assertLessEqual(len(flip_orbit_codes(single)), 2)
        with self.assertRaises(ValueError):
            flip_quotient_count(0)


if __name__ == "__main__":
    unittest.main()
