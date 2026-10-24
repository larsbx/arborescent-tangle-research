"""A4.9 class-BFS tests: no fabricated negative isotopy conclusions."""
import unittest

from arborescent.presentation_class_bfs import (
    decide_formal_class_reachability, explore_presentation_classes,
    rooted_presentation_key,
)
from arborescent.mixed_invariants import (
    decide_partial_formal_reachability,
    unordered_abbreviated_fingerprint,
)
from arborescent.partial_joint_f3prime import (
    mixed_successors, explore_mixed, MixedMove, apply_mixed,
)
from arborescent.placement_fibers import (
    placement_presentations, bounded_valence_chain_presentations,
)
from arborescent.tangle_labels import V4LabeledTree, label_free_ports
from arborescent.weighted_trees import (
    Fork, Leaf, Vertex, WeightedPlaneTree, to_weighted_tree,
)


def isomorphic_renaming(source: V4LabeledTree) -> V4LabeledTree:
    """Reverse stored vertex order, rename internal bonds, rotate cyclic starts."""
    vertices = []
    for vertex in reversed(source.tree.vertices):
        names = tuple("e"+str(800+int(b[1:])) if b.startswith("e") else b
                      for b in vertex.slots)
        vertices.append(Vertex(names[1:] + names[:1],
                               vertex.sectors[1:] + vertex.sectors[:1]))
    return V4LabeledTree(WeightedPlaneTree(tuple(vertices)), source.labels)


def f3_witness(a, sign):
    return V4LabeledTree(WeightedPlaneTree((
        Vertex(("f0", "e0", "e1"), (a, 0, 0)),
        Vertex(("e0", "e2", "e3", "f1"), (0, 0, 0, 0)),
        Vertex(("e2", "e4", "e5"), (0, sign, 0)),
        Vertex(("e4",), (2*sign,)),
        Vertex(("e5",), (2*sign,)),
        Vertex(("e3",), (3*sign,)),
        Vertex(("e1",), (2*sign,)),
    )), (("f0", 1), ("f1", 2)))


class FormalPresentationClassTests(unittest.TestCase):
    def test_rooted_isomorphism_reaches_renumbered_target_with_empty_trace(self):
        for sign in (-1, 1):
            start = placement_presentations(3, 1, sign)[0]
            target = isomorphic_renaming(start)
            self.assertEqual(rooted_presentation_key(start),
                             rooted_presentation_key(target))
            fixed_index = decide_partial_formal_reachability(
                start, target, sign, max_states=64
            )
            self.assertEqual(fixed_index.status,
                             "not_reached_in_complete_indexed_graph")
            result = decide_formal_class_reachability(
                start, target, sign, max_classes=64
            )
            self.assertEqual(result.status,
                             "reachable_up_to_presentation_isomorphism")
            self.assertTrue(result.verify_positive())
            self.assertEqual(result.certificate.trace.moves, ())
            self.assertEqual(result.certificate.target, target)

    def test_nontrivial_ring_move_to_renumbered_target(self):
        for sign in (-1, 1):
            first, second = placement_presentations(3, 1, sign)[:2]
            goal = isomorphic_renaming(second)
            decision = decide_formal_class_reachability(
                first, goal, sign, max_classes=64
            )
            self.assertEqual(decision.status,
                             "reachable_up_to_presentation_isomorphism")
            self.assertTrue(decision.verify_positive())
            self.assertEqual(len(decision.certificate.trace.moves), 1)
            self.assertEqual(
                rooted_presentation_key(decision.certificate.trace.target),
                rooted_presentation_key(goal),
            )

    def test_all_complete_small_corolla_class_codes_equal_indexed_quotient(self):
        for sign in (-1, 1):
            for d, caps in ((3, 1), (3, 2), (4, 1)):
                with self.subTest(sign=sign, d=d, caps=caps):
                    start = placement_presentations(d, caps, sign)[0]
                    raw = explore_mixed(start, sign, max_states=256)
                    classes = explore_presentation_classes(
                        start, sign, max_classes=256
                    )
                    self.assertTrue(raw.complete)
                    self.assertTrue(classes.complete)
                    self.assertEqual(
                        set(classes.representatives),
                        {rooted_presentation_key(x) for x in raw.incoming},
                    )
                    for obj in raw.incoming:
                        self.assertTrue(classes.certificate(obj).verify())

    def test_successor_class_codes_are_equivariant_under_storage_renaming(self):
        for sign in (-1, 1):
            for state in (
                placement_presentations(3, 1, sign)[0],
                bounded_valence_chain_presentations(2, sign)[0],
                f3_witness(3, sign),
            ):
                renamed = isomorphic_renaming(state)
                self.assertEqual(rooted_presentation_key(state),
                                 rooted_presentation_key(renamed))
                before = {rooted_presentation_key(child)
                          for _, child in mixed_successors(state, sign)}
                after = {rooted_presentation_key(child)
                         for _, child in mixed_successors(renamed, sign)}
                self.assertEqual(before, after)

    def test_complete_formal_class_negative_for_ordinary_branch_reordering(self):
        labels = (2, 3, 4)
        a = label_free_ports(to_weighted_tree(
            Fork((Leaf(2), Leaf(3), Leaf(4))), labels))
        b = label_free_ports(to_weighted_tree(
            Fork((Leaf(3), Leaf(2), Leaf(4))), labels))
        self.assertEqual(unordered_abbreviated_fingerprint(a, 1),
                         unordered_abbreviated_fingerprint(b, 1))
        self.assertNotEqual(rooted_presentation_key(a),
                            rooted_presentation_key(b))
        decision = decide_formal_class_reachability(
            a, b, 1, max_classes=64
        )
        self.assertEqual(decision.status,
                         "not_reached_in_complete_formal_class_graph")
        self.assertTrue(decision.complete)
        self.assertIsNone(decision.certificate)
        self.assertFalse(decision.verify_positive())

    def test_nonroot_v4_label_not_discarded_in_full_class_key(self):
        a = bounded_valence_chain_presentations(2, 1)[0]
        b = V4LabeledTree(
            a.tree, tuple((bond, 2 if bond == "f1" else value)
                          for bond, value in a.labels)
        )
        self.assertEqual(unordered_abbreviated_fingerprint(a, 1),
                         unordered_abbreviated_fingerprint(b, 1))
        result = decide_formal_class_reachability(
            a, b, 1, max_classes=64
        )
        self.assertEqual(result.status,
                         "not_reached_in_complete_formal_class_graph")
        self.assertTrue(result.complete)

    def test_cutoff_and_necessary_exclusion_are_distinct(self):
        start, target = placement_presentations(3, 1)[:2]
        budget = decide_formal_class_reachability(
            start, target, 1, max_classes=1
        )
        self.assertEqual(budget.status, "undetermined_cutoff")
        self.assertEqual(budget.examined_classes, 1)
        self.assertFalse(budget.complete)
        changed = label_free_ports(to_weighted_tree(
            Fork((Leaf(2), Leaf(3))), (2, 3)))
        distinct = label_free_ports(to_weighted_tree(
            Fork((Leaf(2), Leaf(4))), (2, 4)))
        excluded = decide_formal_class_reachability(
            changed, distinct, 1, max_classes=1
        )
        self.assertEqual(excluded.status,
                         "excluded_by_necessary_invariant")
        self.assertEqual(excluded.examined_classes, 0)
        self.assertFalse(excluded.verify_positive())

    def test_invalid_root_type_sign_and_budget(self):
        source = placement_presentations(3, 1)[0]
        for cap in (0, -1, 1.2, True):
            with self.subTest(cap=cap), self.assertRaises(ValueError):
                explore_presentation_classes(source, 1, max_classes=cap)
        for root in ("missing", "e0", "", None):
            with self.subTest(root=root), self.assertRaises(ValueError):
                decide_formal_class_reachability(
                    source, source, 1, root=root
                )
        with self.assertRaises(ValueError):
            decide_formal_class_reachability(None, source, 1)
        with self.assertRaises(ValueError):
            decide_formal_class_reachability(source, source, True)


if __name__ == "__main__":
    unittest.main()
