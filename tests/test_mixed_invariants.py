"""A4.8 independent necessary-invariant and decision-status regressions."""
import unittest

from arborescent.mixed_invariants import (
    decide_partial_formal_reachability,
    unordered_abbreviated_fingerprint,
)
from arborescent.partial_joint_f3prime import (
    MixedMove, apply_mixed, mixed_successors,
)
from arborescent.placement_fibers import (
    bounded_valence_chain_presentations, placement_presentations,
)
from arborescent.ring_relocation import RingSlide
from arborescent.tangle_labels import (
    V4LabeledTree, label_free_ports, tangle_flip_f1,
)
from arborescent.weighted_trees import (
    Fork, Leaf, Vertex, WeightedPlaneTree, to_weighted_tree,
)


def f3_witness(a: int, sign: int) -> V4LabeledTree:
    return V4LabeledTree(WeightedPlaneTree((
        Vertex(("f0", "e0", "e1"), (a, 0, 0)),
        Vertex(("e0", "e2", "e3", "f1"), (0, 0, 0, 0)),
        Vertex(("e2", "e4", "e5"), (0, sign, 0)),
        Vertex(("e4",), (2*sign,)),
        Vertex(("e5",), (2*sign,)),
        Vertex(("e3",), (3*sign,)),
        Vertex(("e1",), (2*sign,)),
    )), (("f0", 1), ("f1", 2)))


def renamed_presentation(source: V4LabeledTree) -> V4LabeledTree:
    """Same marked plane presentation; deliberately new vertex/edge IDs."""
    vertices = tuple(
        Vertex(
            tuple("e"+str(1000+int(bond[1:])) if bond.startswith("e")
                  else bond for bond in vertex.slots),
            vertex.sectors,
        )
        for vertex in reversed(source.tree.vertices)
    )
    return V4LabeledTree(WeightedPlaneTree(vertices), source.labels)


class FormalNecessaryInvariantTests(unittest.TestCase):
    def test_every_one_step_source_or_formal_move_preserves_fingerprint(self):
        for sign in (-1, 1):
            fixtures = (
                placement_presentations(3, 1, sign)[0],
                bounded_valence_chain_presentations(2, sign)[0],
                f3_witness(3, sign),
                f3_witness(4, sign),
            )
            for source in fixtures:
                before = unordered_abbreviated_fingerprint(source, sign)
                moves = mixed_successors(source, sign)
                self.assertTrue(moves)
                for move, child in moves:
                    with self.subTest(sign=sign, move=move):
                        self.assertEqual(
                            unordered_abbreviated_fingerprint(child, sign),
                            before,
                        )

    def test_explicit_strict_f3_move_preserves_total_weights(self):
        for sign in (-1, 1):
            for a in (-4, -3, 3, 4):
                source = f3_witness(a, sign)
                move = MixedMove("F3prime", "f0", "e0", "before")
                dest = apply_mixed(source, sign, move)
                self.assertNotEqual(source.labelled_rooted_code(),
                                    dest.labelled_rooted_code())
                self.assertEqual(unordered_abbreviated_fingerprint(source, sign),
                                 unordered_abbreviated_fingerprint(dest, sign))

    def test_mismatching_weighted_leaf_is_provably_unreachable_formally(self):
        a = label_free_ports(
            to_weighted_tree(Fork((Leaf(2), Leaf(3))), (2, 3)))
        b = label_free_ports(
            to_weighted_tree(Fork((Leaf(2), Leaf(4))), (2, 4)))
        decision = decide_partial_formal_reachability(a, b, 1, max_states=1)
        self.assertEqual(decision.status, "excluded_by_necessary_invariant")
        self.assertNotEqual(decision.fingerprint_start,
                            decision.fingerprint_target)
        self.assertEqual(decision.examined_states, 0)
        self.assertFalse(decision.complete)
        self.assertIsNone(decision.trace)
        self.assertFalse(decision.verify_positive())

    def test_reachable_formal_pair_has_shortest_replay_certificate(self):
        for sign in (-1, 1):
            source, target = placement_presentations(3, 1, sign)[:2]
            result = decide_partial_formal_reachability(
                source, target, sign, max_states=64
            )
            self.assertEqual(result.status, "reachable_formally")
            self.assertIsNotNone(result.trace)
            self.assertTrue(result.verify_positive())
            self.assertEqual(len(result.trace.moves), 1)
            self.assertEqual(result.fingerprint_start,
                             result.fingerprint_target)

    def test_cutoff_never_claims_inequivalence(self):
        source, target = placement_presentations(3, 1)[:2]
        decision = decide_partial_formal_reachability(
            source, target, 1, max_states=1
        )
        self.assertEqual(decision.status, "undetermined_cutoff")
        self.assertEqual(decision.examined_states, 1)
        self.assertFalse(decision.complete)
        self.assertIsNone(decision.trace)
        self.assertFalse(decision.verify_positive())

    def test_complete_indexed_search_is_not_isomorphism_quotient(self):
        source = placement_presentations(3, 1)[0]
        target = renamed_presentation(source)
        self.assertEqual(source.labelled_rooted_code(),
                         target.labelled_rooted_code())
        self.assertEqual(unordered_abbreviated_fingerprint(source, 1),
                         unordered_abbreviated_fingerprint(target, 1))
        answer = decide_partial_formal_reachability(
            source, target, 1, max_states=64
        )
        # Despite presentation isomorphism, no indexed move renames IDs.
        # This is why even a complete indexed search is NOT an
        # isotopy or (unindexed) presentation-inequivalence oracle.
        self.assertEqual(
            answer.status, "not_reached_in_complete_indexed_graph"
        )
        self.assertTrue(answer.complete)
        self.assertIsNone(answer.trace)
        self.assertFalse(answer.verify_positive())

    def test_storage_and_v4_changes_do_not_break_necessary_invariant(self):
        source = bounded_valence_chain_presentations(3)[0]
        reflected = tangle_flip_f1(source)
        self.assertEqual(unordered_abbreviated_fingerprint(source, 1),
                         unordered_abbreviated_fingerprint(reflected, 1))
        renamed = renamed_presentation(source)
        self.assertEqual(unordered_abbreviated_fingerprint(source, 1),
                         unordered_abbreviated_fingerprint(renamed, 1))
        changed_labels = V4LabeledTree(
            source.tree,
            tuple((bond, label ^ (2 if bond == "f2" else 0))
                  for bond, label in source.labels)
        )
        self.assertEqual(unordered_abbreviated_fingerprint(source, 1),
                         unordered_abbreviated_fingerprint(changed_labels, 1))

    def test_invalid_sign_root_type_budget_fail_closed(self):
        source = placement_presentations(3, 1)[0]
        for sign in (0, True, -7):
            with self.subTest(sign=sign), self.assertRaises(ValueError):
                unordered_abbreviated_fingerprint(source, sign)
        for root in ("f99", "e0", "", None):
            with self.subTest(root=root), self.assertRaises(ValueError):
                unordered_abbreviated_fingerprint(source, 1, root)
        for budget in (0, -1, True, 1.2):
            with self.subTest(budget=budget), self.assertRaises(ValueError):
                decide_partial_formal_reachability(source, source, 1,
                                                    max_states=budget)
        with self.assertRaises(ValueError):
            decide_partial_formal_reachability(None, source, 1)


if __name__ == "__main__":
    unittest.main()
