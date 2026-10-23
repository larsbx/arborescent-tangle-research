"""Diagram-transcribed strict ring caps: positive/negative and controls.

No independent knot-isotopy claim; tests verify formal local moves,
rooted-presentation isomorphisms, all surviving V4 labels and source
picture preconditions.
"""
import unittest

from arborescent.internal_forks import InternalFork, Tip, Port, compile_fork, tangle_screen
from arborescent.ring_caps import (
    StrictRingCap, contract_ring_cap, expand_ring_cap, strict_ring_caps,
)
from arborescent.ring_signatures import ring_signatures
from arborescent.tangle_labels import (
    V4LabeledTree, label_free_ports, tangle_flip_f1, tangle_flip_f2,
)
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def witness(sign: int) -> V4LabeledTree:
    """Separate base with two branches, attached to ±2--±1--±2 cap."""
    return compile_fork(InternalFork((
        InternalFork((Tip(2 * sign), Tip(2 * sign)), sign, 1),
        Tip(-3 * sign),
    ), 3 * sign, 2))


class StrictRingCapTests(unittest.TestCase):
    def test_actual_positive_and_negative_picture_subgraphs(self):
        for sign in (1, -1):
            with self.subTest(sign=sign):
                source = witness(sign)
                self.assertEqual(len(ring_signatures(source)), 1)
                caps = strict_ring_caps(source)
                self.assertEqual(len(caps), 1)
                cap = caps[0]
                self.assertEqual((cap.sign, cap.base, cap.center), (sign, 0, 1))
                self.assertEqual(cap.leaves, (2, 3))
                self.assertEqual(cap.stem_edge, "e0")
                contracted = contract_ring_cap(source, cap)
                self.assertEqual(len(source.tree.vertices) - len(contracted.tree.vertices), 3)
                self.assertEqual(contracted.tree.vertices[0].slots, ("f0", "e3"))
                self.assertEqual(contracted.tree.vertices[0].sectors, (0, 3 * sign))
                self.assertEqual(contracted.tree.vertices[1].weight, -3 * sign)
                self.assertEqual(contracted.labels, source.labels)
                # This small example also passes current nonexceptional
                # tangle screens; their completeness is NOT claimed.
                polarity = "+" if sign > 0 else "-"
                self.assertTrue(tangle_screen(source, polarity).passed)
                self.assertTrue(tangle_screen(contracted, polarity).passed)

    def test_ring_inverse_replays_up_to_marked_plane_isomorphism(self):
        for sign in (1, -1):
            original = witness(sign)
            contracted = contract_ring_cap(original, strict_ring_caps(original)[0])
            expanded = expand_ring_cap(contracted, 0, "f0", sign)
            self.assertEqual(original.labelled_rooted_code(),
                             expanded.labelled_rooted_code())
            self.assertEqual(len(strict_ring_caps(expanded)), 1)
            self.assertEqual(
                contract_ring_cap(expanded, strict_ring_caps(expanded)[0]),
                contracted,
            )
            self.assertEqual(expanded.labels, contracted.labels)

    def test_conservative_picture_vs_prior_permissive_signature(self):
        # A lone +1 fork with two +2 tips and root free-bond f0 matches
        # the earlier NECESSARY center motif, but is NOT the ring
        # diagram (which has a separate base vertex under its stem).
        false_candidate = compile_fork(InternalFork((Tip(2), Tip(2)), 1, 1))
        self.assertEqual(len(ring_signatures(false_candidate)), 1)
        self.assertEqual(strict_ring_caps(false_candidate), ())
        # Even a separate base fails if center weight is not in the
        # angular sector *between* its two ±2 leaves.
        wrong_wedge = compile_fork(InternalFork((
            InternalFork((Tip(2), Tip(2)), 1, 0), Tip(-3),
        ), 3, 2))
        self.assertTrue(ring_signatures(wrong_wedge))
        self.assertEqual(strict_ring_caps(wrong_wedge), ())
        # Nonzero angular sectors adjacent to the cap stem are beyond
        # this conservative zero-split local rule.
        wrong_base = compile_fork(InternalFork((
            InternalFork((Tip(2), Tip(2)), 1, 1), Tip(-3),
        ), 3, 0))
        self.assertEqual(strict_ring_caps(wrong_base), ())

    def test_preserves_nontrivial_v4_labels_and_flip_orientations(self):
        src = witness(1)
        marked = V4LabeledTree(src.tree, (("f0", 2),))
        reduced = contract_ring_cap(marked, strict_ring_caps(marked)[0])
        self.assertEqual(reduced.label("f0"), 2)
        for action in (tangle_flip_f1(marked), tangle_flip_f2(marked, 0)):
            self.assertEqual(len(strict_ring_caps(action)), 1)
            cap = strict_ring_caps(action)[0]
            self.assertEqual(contract_ring_cap(action, cap).label("f0"),
                             action.label("f0"))

    def test_two_free_ports_on_base_survive(self):
        # A second free bond can be part of the base's other arms.
        s = compile_fork(InternalFork((
            InternalFork((Tip(2), Tip(2)), 1, 1), Port(),
        ), 3, 2))
        self.assertEqual(len(s.labels), 2)
        self.assertEqual(len(strict_ring_caps(s)), 1)
        labeled = V4LabeledTree(s.tree, (("f0", 1), ("f1", 3)))
        reduced = contract_ring_cap(labeled, strict_ring_caps(labeled)[0])
        self.assertEqual(dict(reduced.labels), {"f0": 1, "f1": 3})
        self.assertEqual(reduced.tree.vertices[0].slots, ("f0", "f1"))

    def test_guards_for_unsourced_replacements(self):
        s = witness(1)
        wrong = StrictRingCap(1, 0, 1, (2, 4), "e0")
        with self.assertRaises(ValueError):
            contract_ring_cap(s, wrong)
        with self.assertRaises(ValueError):
            strict_ring_caps(None)
        reduced = contract_ring_cap(s, strict_ring_caps(s)[0])
        for base, bond, sign in ((0, "e3", 1), (0, "missing", 1),
                                 (-1, "f0", 1), (0, "f0", 0),
                                 (0, "f0", True), (True, "f0", 1)):
            with self.subTest(base=base, bond=bond, sign=sign):
                with self.assertRaises(ValueError):
                    expand_ring_cap(reduced, base, bond, sign)


if __name__ == "__main__":
    unittest.main()
