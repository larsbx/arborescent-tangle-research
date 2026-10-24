"""A4.5 formal restricted F3'-abbreviation naturality, not knot isotopy."""
import unittest

from arborescent.abbreviated_f3prime import (
    abbreviated_f3prime, f3prime_abbreviation_square,
)
from arborescent.abbreviated_flips import flip_abbreviation_square
from arborescent.internal_forks import tangle_screen
from arborescent.ring_abbreviation import abbreviate_disjoint_rings
from arborescent.tangle_labels import (
    V4LabeledTree, tangle_flip_f3prime,
)
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def witness(weight: int, sign: int) -> V4LabeledTree:
    """Surviving root weighted pivot and downstream attached ring cap."""
    return V4LabeledTree(WeightedPlaneTree((
        Vertex(("f0", "e0", "e1"), (weight, 0, 0)),
        Vertex(("e0", "e2", "e3", "f1"), (0, 0, 0, 0)),
        Vertex(("e2", "e4", "e5"), (0, sign, 0)),
        Vertex(("e4",), (2*sign,)),
        Vertex(("e5",), (2*sign,)),
        Vertex(("e3",), (3*sign,)),
        Vertex(("e1",), (2*sign,)),
    )), (("f0", 1), ("f1", 2)))


class AbbreviatedF3PrimeTests(unittest.TestCase):
    def test_odd_even_weights_and_both_ring_signs(self):
        for weight in (3, 4, -3, -4):
            for sign in (-1, 1):
                with self.subTest(weight=weight, sign=sign):
                    source = witness(weight, sign)
                    polarity = "+" if sign == 1 else "-"
                    # These examples pass the *implemented branched*
                    # tangle W/S/P-or-N structural screen. This is not
                    # a proof that every source exception is modeled.
                    self.assertTrue(tangle_screen(source, polarity).passed)
                    # The older *pair* Free-Bond condition F excludes
                    # a free port on a weighted root; §12.7 *tangles*
                    # do not require F. Do not conflate the screens.
                    self.assertIn(("F", 0), source.tree.wfs_failures())
                    self.assertTrue(f3prime_abbreviation_square(
                        source, sign, "f0", "e0", from_side="before"
                    ))
                    direct = tangle_flip_f3prime(
                        source, 0, "e0", from_side="before"
                    )
                    reduced_before = abbreviate_disjoint_rings(source, sign)
                    reduced_after = abbreviate_disjoint_rings(direct, sign)
                    self.assertEqual(reduced_before.rings, reduced_after.rings)
                    self.assertEqual(sum(reduced_after.rings), 1)
                    self.assertEqual(
                        abbreviated_f3prime(
                            reduced_before, "f0", "e0", from_side="before"
                        ), reduced_after,
                    )
                    self.assertTrue(tangle_screen(direct, polarity).passed)
                    self.assertEqual(direct.label("f0"), source.label("f0"))
                    if weight % 2:
                        self.assertEqual(direct.label("f1"), source.label("f1") ^ 2)
                    else:
                        self.assertEqual(direct.label("f1"), source.label("f1"))

    def test_inverse_square_and_ring_counts(self):
        for weight in (3, 4):
            for sign in (-1, 1):
                original = witness(weight, sign)
                changed = tangle_flip_f3prime(
                    original, 0, "e0", from_side="before"
                )
                self.assertTrue(f3prime_abbreviation_square(
                    changed, sign, "f0", "e0", from_side="after"
                ))
                self.assertEqual(
                    tangle_flip_f3prime(
                        changed, 0, "e0", from_side="after"
                    ), original,
                )
                ab = abbreviate_disjoint_rings(changed, sign)
                restored = abbreviated_f3prime(
                    ab, "f0", "e0", from_side="after"
                )
                self.assertEqual(
                    restored, abbreviate_disjoint_rings(original, sign)
                )

    def test_f1_f2_and_f3prime_square_are_separate_obligations(self):
        for sign in (-1, 1):
            src = witness(3, sign)
            self.assertEqual(flip_abbreviation_square(src, sign), (True, True))
            self.assertTrue(f3prime_abbreviation_square(
                src, sign, "f0", "e0"
            ))

    def test_unsupported_erased_stem_and_missing_pivot_rejected(self):
        src = witness(3, 1)
        with self.assertRaisesRegex(ValueError, "survive"):
            f3prime_abbreviation_square(src, 1, "f1", "e2")
        with self.assertRaises(ValueError):
            f3prime_abbreviation_square(src, 1, "f99", "e0")
        with self.assertRaises(ValueError):
            f3prime_abbreviation_square(None, 1, "f0", "e0")
        with self.assertRaises(ValueError):
            abbreviated_f3prime(None, "f0", "e0")
        reduced = abbreviate_disjoint_rings(src, 1)
        with self.assertRaises(ValueError):
            abbreviated_f3prime(reduced, "f99", "e0")
        with self.assertRaises(ValueError):
            abbreviated_f3prime(reduced, "f0", "e2")

    def test_zero_weight_or_two_nonzero_pivot_sectors_fail_closed(self):
        for sectors in ((0, 0, 0), (3, 0, 5)):
            src = witness(3, 1)
            old = src.tree.vertices[0]
            changed_tree = WeightedPlaneTree((
                Vertex(old.slots, sectors),
                *src.tree.vertices[1:],
            ))
            changed = V4LabeledTree(changed_tree, src.labels)
            with self.assertRaises(ValueError):
                f3prime_abbreviation_square(changed, 1, "f0", "e0")

    def test_unrecognized_candidate_and_nonzero_cap_base_fail_closed(self):
        src = witness(3, 1)
        old = src.tree.vertices[1]
        # One nonzero sector at base adjacent to e2 prevents strict
        # abbreviation. F3' must not label this an equivalence square.
        bad_tree = WeightedPlaneTree((
            src.tree.vertices[0],
            Vertex(old.slots, (0, 7, 0, 0)),
            *src.tree.vertices[2:],
        ))
        bad = V4LabeledTree(bad_tree, src.labels)
        with self.assertRaisesRegex(ValueError, "incomplete"):
            f3prime_abbreviation_square(bad, 1, "f0", "e0")


if __name__ == "__main__":
    unittest.main()
