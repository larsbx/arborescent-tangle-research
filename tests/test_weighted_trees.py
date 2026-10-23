"""Independent exact counts and structural negative controls."""
import unittest
from arborescent.weighted_trees import (
    Fork, Leaf, Vertex, WeightedPlaneTree, count_by_arity,
    enumerate_specs, syntax_count_closed, syntax_counts, to_weighted_tree,
)


class WeightedTreeTests(unittest.TestCase):
    def test_exact_grammar_closed_formula_and_exhaustive(self):
        for m in (1, 2, 3):
            counts = syntax_counts(6, m)
            self.assertEqual(counts[0], m)
            for n in range(7):
                self.assertEqual(counts[n], syntax_count_closed(n, m))
            for n in range(3):
                labels = tuple(range(2, m + 2))
                specs = enumerate_specs(n, labels)
                self.assertEqual(len(specs), counts[n])
                self.assertEqual(len(set(specs)), len(specs))

    def test_nonbinary_branches_and_actual_sector_data(self):
        spec = Fork((Leaf(2), Fork((Leaf(3), Leaf(2), Leaf(3))), Leaf(2)))
        tree = to_weighted_tree(spec)
        self.assertEqual(len(tree.vertices), 7)
        self.assertEqual(len(tree.edges()), 6)
        self.assertEqual(tree.vertices[0].slots[0], 'f0')
        self.assertEqual(tree.vertices[0].sectors, (0, 0, 0, 0))
        self.assertEqual(tree.wfs_failures(), ())
        self.assertEqual(sum(1 for v in tree.vertices if len(v.slots) >= 3), 2)

    def test_presentation_code_ignores_renumbering_and_slot_offset(self):
        spec = Fork((Leaf(2), Fork((Leaf(3), Leaf(2))), Leaf(3)))
        tree = to_weighted_tree(spec)
        remapped = []
        for vertex in reversed(tree.vertices):
            slots = tuple('e' + str(int(s[1:]) + 27) if s[0] == 'e' else s
                          for s in vertex.slots)
            remapped.append(Vertex(slots[1:] + slots[:1],
                                   vertex.sectors[1:] + vertex.sectors[:1]))
        equivalent = WeightedPlaneTree(tuple(remapped))
        self.assertEqual(tree.rooted_plane_code(), equivalent.rooted_plane_code())
        distinct = to_weighted_tree(Fork((Leaf(3), Fork((Leaf(3), Leaf(2))), Leaf(2))))
        self.assertNotEqual(tree.rooted_plane_code(), distinct.rooted_plane_code())
        with self.assertRaises(ValueError):
            tree.rooted_plane_code('f9')

    def test_all_bounded_small_trees_satisfy_strict_wfs(self):
        for n in range(1, 4):
            specs = enumerate_specs(n, (2, 3))
            self.assertEqual(len(specs), syntax_counts(n, 2)[-1])
            for spec in specs:
                self.assertEqual(to_weighted_tree(spec).wfs_failures(), ())

    def test_w_f_s_detect_counterexamples(self):
        bad_w = WeightedPlaneTree((
            Vertex(('f0', 'e0', 'e1'), (1, 1, 0)),
            Vertex(('e0',), (2,)), Vertex(('e1',), (3,)),
        ))
        self.assertIn(('W', 0), bad_w.wfs_failures())
        self.assertIn(('F', 0), bad_w.wfs_failures())
        bad_f = WeightedPlaneTree((Vertex(('f0',), (2,)),))
        self.assertIn(('F', 0), bad_f.wfs_failures())
        bad_s = WeightedPlaneTree((
            Vertex(('e0',), (2,)),
            Vertex(('e0', 'e1'), (3, 0)),
            Vertex(('e1',), (-1,)),
        ))
        self.assertIn(('S-sign', 0), bad_s.wfs_failures())
        self.assertIn(('S-end', 0), bad_s.wfs_failures())
        bad_zero = WeightedPlaneTree((
            Vertex(('e0',), (2,)), Vertex(('e0',), (0,)),
        ))
        self.assertIn(('S-zero', 0), bad_zero.wfs_failures())

    def test_malformed_graphs_rejected(self):
        bad = (
            (),
            (Vertex(('e0',), (2,)),),
            (Vertex(('f0',), (2,)), Vertex(('f1',), (3,))),
            (Vertex(('e0', 'e0'), (0, 0)), Vertex(('e0',), (2,))),
            (Vertex(('e0',), (2,)), Vertex(('e0',), (3,)), Vertex(('f0',), (0,))),
            (Vertex(('f0',), (True,)),),
        )
        for vertices in bad:
            with self.subTest(vertices=vertices), self.assertRaises(ValueError):
                WeightedPlaneTree(vertices)

    def test_invalid_compiler_inputs(self):
        for spec in (Leaf(2), Fork((Leaf(2),)), Fork((Leaf(2), Leaf(5)))):
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                to_weighted_tree(spec)

    def test_arity_refinement(self):
        self.assertEqual(count_by_arity(2, 0, 2), 16)
        self.assertEqual(count_by_arity(2, 1, 2), 80)
        self.assertEqual(count_by_arity(2, 2, 2), 96)
        self.assertEqual(syntax_counts(2, 2), (2, 12, 192))


if __name__ == '__main__':
    unittest.main()
