"""Exact, replayable *partial-move* closure tests and A1 vacuity control."""
import unittest

from arborescent.partial_moves import (
    Move, Trace, a1_has_applicable_f3prime, apply, explore, successors,
)
from arborescent.tangle_labels import V4LabeledTree, label_free_ports
from arborescent.weighted_trees import Fork, Leaf, Vertex, WeightedPlaneTree, to_weighted_tree


class MoveExplorerTests(unittest.TestCase):
    def test_a1_f3prime_is_vacuous_exhaustive_small_instances(self):
        # All nonzero weights belong to valence-one leaves: F3' can
        # never transport one across the pair of distinct sectors
        # adjacent to a branch edge. This is a theorem for all n.
        for n in range(1, 4):
            self.assertFalse(a1_has_applicable_f3prime(n, (2, 3)))
        with self.assertRaises(ValueError):
            a1_has_applicable_f3prime(0)

    def test_all_reachable_states_and_certificates_single_free_port(self):
        tree = to_weighted_tree(Fork((Leaf(2), Leaf(3))))
        initial = label_free_ports(tree)
        orbit = explore(initial)
        self.assertTrue(orbit.complete)
        self.assertEqual(orbit.size, 4)
        for state in orbit.predecessors:
            self.assertTrue(orbit.trace(state).verify())
            self.assertEqual(state.tree.wfs_failures(), ())
        self.assertEqual(set(
            s.label("f0") for s in orbit.predecessors
        ), {0, 1, 2, 3})

    def test_f3_witness_and_exact_trace_from_move_search(self):
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (3, 0, 0)),
            Vertex(("e0", "e2", "f1"), (0, 0, 0)),
            Vertex(("e1",), (2,)),
            Vertex(("e2",), (3,)),
        ))
        start = V4LabeledTree(tree, (("f0", 0), ("f1", 0)))
        f3 = Move("F3prime", 0, "e0", "before")
        target = apply(start, f3)
        self.assertNotEqual(target, start)
        self.assertNotEqual(target.label("f1"), start.label("f1"))
        self.assertIn(f3, [m for m, _ in successors(start)])
        # Search is exact on the fixed *indexed* weighted tree, with no
        # assumption that it is a topological normal form.
        orbit = explore(start, max_states=32768)
        self.assertTrue(orbit.complete)
        self.assertGreater(orbit.size, 4)
        self.assertLess(orbit.size, 32768)
        path = orbit.trace(target)
        self.assertTrue(path.verify())
        self.assertEqual(path.moves, (f3,))
        # Every stored predecessor chain is rechecked against the rule.
        for reachable in orbit.predecessors:
            self.assertTrue(orbit.trace(reachable).verify())

    def test_cutoff_is_flagged_never_claimed_complete(self):
        tree = to_weighted_tree(Fork((Leaf(2), Leaf(3))))
        obj = label_free_ports(tree)
        partial = explore(obj, max_states=1)
        self.assertFalse(partial.complete)
        self.assertEqual(partial.size, 1)
        self.assertTrue(partial.trace(obj).verify())
        with self.assertRaises(ValueError):
            explore(obj, max_states=0)
        with self.assertRaises(ValueError):
            explore(obj, max_states=True)
        with self.assertRaises(ValueError):
            explore(None)
        with self.assertRaises(ValueError):
            partial.trace(label_free_ports(tree, 1))

    def test_invalid_move_types_are_explicit(self):
        invalid = (
            ("Rplus",),
            ("F1", 0),
            ("F2", True),
            ("F2", -1),
            ("F3prime", 0, "", "before"),
            ("F3prime", 0, "e0", ""),
        )
        for args in invalid:
            with self.subTest(args=args), self.assertRaises(ValueError):
                Move(*args)

    def test_trace_rejects_fake_end_state(self):
        obj = label_free_ports(to_weighted_tree(Fork((Leaf(2), Leaf(3)))))
        flipped = apply(obj, Move("F1"))
        true = Trace(obj, flipped, (Move("F1"),))
        false = Trace(obj, obj, (Move("F1"),))
        self.assertTrue(true.verify())
        self.assertFalse(false.verify())


if __name__ == "__main__":
    unittest.main()
