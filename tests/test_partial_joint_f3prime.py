"""A4.7 finite formal mixed move graph, with no topological isotopy claim."""
import unittest

from arborescent.partial_joint_f3prime import (
    MixedMove, apply_mixed, mixed_successors, explore_mixed,
)
from arborescent.ring_relocation import RingSlide
from arborescent.placement_fibers import placement_presentations
from arborescent.ring_abbreviation import abbreviate_disjoint_rings
from arborescent.abbreviated_f3prime import f3prime_abbreviation_square
from arborescent.abbreviated_flips import abbreviated_flip_orbit_key
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def decorated_witness(weight, sign):
    """Full surviving-root-edge strict F3' witness plus movable ring cap."""
    return V4LabeledTree(WeightedPlaneTree((
        Vertex(("f0", "e0", "e1"), (weight, 0, 0)),
        Vertex(("e0", "e2", "e3", "f1"), (0, 0, 0, 0)),
        Vertex(("e2", "e4", "e5"), (0, sign, 0)),
        Vertex(("e4",), (2*sign,)),
        Vertex(("e5",), (2*sign,)),
        Vertex(("e3",), (3*sign,)),
        Vertex(("e1",), (2*sign,)),
    )), (("f0", 1), ("f1", 2)))


class MixedPartialFormalTests(unittest.TestCase):
    def test_known_source_f3_commuting_square_and_one_step_certificate(self):
        for sign in (-1, 1):
            for weight in (-4, -3, 3, 4):
                with self.subTest(sign=sign, weight=weight):
                    start = decorated_witness(weight, sign)
                    move = MixedMove("F3prime", "f0", "e0", "before")
                    self.assertTrue(f3prime_abbreviation_square(
                        start, sign, "f0", "e0", from_side="before"
                    ))
                    target = apply_mixed(start, sign, move)
                    self.assertNotEqual(target, start)
                    self.assertEqual(
                        abbreviate_disjoint_rings(start, sign).rings,
                        abbreviate_disjoint_rings(target, sign).rings,
                    )
                    self.assertEqual(
                        apply_mixed(
                            target, sign,
                            MixedMove("F3prime", "f0", "e0", "after")
                        ),
                        start,
                    )
                    possibilities = dict(mixed_successors(start, sign))
                    self.assertIn(move, possibilities)
                    self.assertEqual(possibilities[move], target)
                    orbit = explore_mixed(start, sign, max_states=128)
                    self.assertIn(target, orbit.incoming)
                    self.assertTrue(orbit.trace(target).verify())
                    self.assertEqual(len(orbit.trace(target).moves), 1)
                    # Bound exhaustion is not a no-equivalence decision:
                    self.assertTrue(all(len(s.tree.vertices) == 7
                                        for s in orbit.incoming))

    def test_explicit_f3_vs_formal_ring_slide_diamond(self):
        # F3' on e0 does not erase the cap. With odd pivot weight, it
        # mirrors the cap-bearing base, reversing slide direction;
        # nevertheless the two formal paths reach the identical
        # indexed V4-labelled presentation.
        for sign in (-1, 1):
            for weight in (-4, -3, 3, 4):
                with self.subTest(sign=sign, weight=weight):
                    start = decorated_witness(weight, sign)
                    f3 = MixedMove("F3prime", "f0", "e0", "before")
                    ring = MixedMove("R-formal",
                                     ring_slide=RingSlide(sign, 2, "e3"))
                    left = apply_mixed(
                        apply_mixed(start, sign, ring), sign, f3
                    )
                    right = apply_mixed(
                        apply_mixed(start, sign, f3), sign, ring
                    )
                    self.assertEqual(left, right)
                    self.assertEqual(left.labels, right.labels)
                    self.assertEqual(left.tree, right.tree)
                    self.assertEqual(abbreviate_disjoint_rings(left, sign),
                                     abbreviate_disjoint_rings(right, sign))

    def test_no_f3_candidate_on_all_zero_corolla_matches_joint_syntax_orbit(self):
        # The root/cap bases are zero; nonzero ring tips have no free
        # pivot, so no new strict F3' move is allowed. Thus the joint
        # F1/F2+formal ring orbit is complete and still four-sheeted.
        for sign in (-1, 1):
            source = placement_presentations(3, 1, sign)[0]
            choices = mixed_successors(source, sign)
            self.assertTrue(choices)
            self.assertFalse(any(m.kind == "F3prime" for m, _ in choices))
            orbit = explore_mixed(source, sign, max_states=64)
            self.assertTrue(orbit.complete)
            self.assertEqual(len(orbit.incoming), 12)
            self.assertTrue(all(orbit.trace(s).verify()
                                for s in orbit.incoming))
            self.assertEqual(
                len({abbreviated_flip_orbit_key(
                    abbreviate_disjoint_rings(s, sign)
                ) for s in orbit.incoming}),
                1,
            )

    def test_cutoff_not_negative_evidence(self):
        start = decorated_witness(3, 1)
        partial = explore_mixed(start, 1, max_states=1)
        self.assertFalse(partial.complete)
        self.assertEqual(len(partial.incoming), 1)
        self.assertTrue(partial.trace(partial.start).verify())
        with self.assertRaises(ValueError):
            explore_mixed(start, 1, max_states=0)
        with self.assertRaises(ValueError):
            explore_mixed(start, 1, max_states=True)

    def test_erased_stem_and_unrecognized_ring_fail_closed(self):
        source = decorated_witness(3, 1)
        with self.assertRaisesRegex(ValueError, "survive"):
            apply_mixed(
                source, 1, MixedMove("F3prime", "f1", "e2", "before")
            )
        with self.assertRaises(ValueError):
            apply_mixed(
                source, 1, MixedMove("F3prime", "f99", "e0", "before")
            )
        original = source.tree.vertices
        corrupt_base = V4LabeledTree(WeightedPlaneTree((
            original[0], Vertex(original[1].slots, (0, 7, 0, 0)),
            *original[2:]
        )), source.labels)
        with self.assertRaisesRegex(ValueError, "incomplete"):
            explore_mixed(corrupt_base, 1)
        with self.assertRaises(ValueError):
            explore_mixed(source, -1)

    def test_invalid_move_types_and_guards(self):
        examples = (
            lambda: MixedMove("invalid"),
            lambda: MixedMove("F1", "f0"),
            lambda: MixedMove("F2", ring_slide=RingSlide(1, 2, "e3")),
            lambda: MixedMove("F3prime", "f0", "e0", ""),
            lambda: MixedMove("F3prime", "f0", "f1", "before"),
            lambda: MixedMove("R-formal"),
            lambda: MixedMove("R-formal", "f0", ring_slide=RingSlide(1, 2, "e3")),
        )
        for constructor in examples:
            with self.assertRaises(ValueError):
                constructor()
        obj = decorated_witness(3, 1)
        with self.assertRaises(ValueError):
            apply_mixed(obj, 1, object())
        with self.assertRaises(ValueError):
            apply_mixed(None, 1, MixedMove("F1"))
        with self.assertRaises(ValueError):
            explore_mixed(obj, 1, root="f99")


if __name__ == "__main__":
    unittest.main()
