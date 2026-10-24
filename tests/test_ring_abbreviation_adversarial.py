"""Adversarial §12.4–12.5 controls: ring notation is not an R± deletion.

Compare the source diagrams (printed pp. 155, 157–158) before changing
the restricted implementation. These are presentation-level assertions,
not independently established isotopies of marked tangles.
"""
import unittest

from arborescent.ring_abbreviation import abbreviate_disjoint_rings, strict_ring_caps
from arborescent.ring_signatures import ring_signatures
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def two_rings_with_one_other_branch(base_slots):
    """Two ring caps and a distinguishable weight-3 branch at one base."""
    return V4LabeledTree(WeightedPlaneTree((
        Vertex(base_slots, (0, 0, 0, 0)),
        Vertex(("e0", "e1", "e2"), (0, 1, 0)),
        Vertex(("e1",), (2,)),
        Vertex(("e2",), (2,)),
        Vertex(("e3", "e4", "e5"), (0, 1, 0)),
        Vertex(("e4",), (2,)),
        Vertex(("e5",), (2,)),
        Vertex(("e6",), (3,)),
    )), (("f0", 0),))


class AdversarialRingAbbreviationTests(unittest.TestCase):
    def test_abbreviation_is_many_to_one_not_graph_deleting_R_move(self):
        # A change in placement of ring branches relative to the weight-3
        # branch changes the marked full tree. Erasing both ring subtrees
        # intentionally forgets that placement. This equality is NOT itself
        # a proof of an allowed R+ move in a marked topological category.
        a = two_rings_with_one_other_branch(("f0", "e0", "e3", "e6"))
        b = two_rings_with_one_other_branch(("f0", "e0", "e6", "e3"))
        self.assertNotEqual(a.tree.rooted_plane_code(), b.tree.rooted_plane_code())
        abbrev_a = abbreviate_disjoint_rings(a, +1)
        abbrev_b = abbreviate_disjoint_rings(b, +1)
        self.assertEqual(abbrev_a.rooted_code(), abbrev_b.rooted_code())
        self.assertEqual(abbrev_a.rings, (2, 0))
        self.assertEqual(len(a.tree.vertices), 8)
        self.assertEqual(len(abbrev_a.labeled.tree.vertices), 2)

    def test_nonroot_V4_boundary_marker_survives_and_distinguishes(self):
        tree = WeightedPlaneTree((
            Vertex(("f0", "e0", "f1", "e3"), (0, 0, 0, 0)),
            Vertex(("e0", "e1", "e2"), (0, 1, 0)),
            Vertex(("e1",), (2,)),
            Vertex(("e2",), (2,)),
            Vertex(("e3",), (3,)),
        ))
        a = V4LabeledTree(tree, (("f0", 0), ("f1", 1)))
        b = V4LabeledTree(tree, (("f0", 0), ("f1", 2)))
        aa = abbreviate_disjoint_rings(a, +1)
        bb = abbreviate_disjoint_rings(b, +1)
        self.assertEqual(aa.rings, (1, 0))
        self.assertEqual(aa.labeled.labels, a.labels)
        self.assertEqual(bb.labeled.labels, b.labels)
        self.assertNotEqual(aa.rooted_code(), bb.rooted_code())

    def test_adjacent_center_exception_rejected_for_both_signs(self):
        # The printed p.157 special two-center tree may not be
        # represented by independently incrementing two ring numbers.
        for sign in (+1, -1):
            with self.subTest(sign=sign):
                source = V4LabeledTree(WeightedPlaneTree((
                    Vertex(("e0", "e1", "e2"), (0, sign, 0)),
                    Vertex(("e0", "e3", "e4"), (0, sign, 0)),
                    Vertex(("e1",), (2 * sign,)),
                    Vertex(("e2",), (2 * sign,)),
                    Vertex(("e3",), (2 * sign,)),
                    Vertex(("e4",), (2 * sign,)),
                )), ())
                self.assertEqual(len(strict_ring_caps(source, sign)), 2)
                with self.assertRaisesRegex(ValueError, "rings attached to rings"):
                    abbreviate_disjoint_rings(source, sign)

    def test_source_ring_variants_not_silently_normalized(self):
        # §12.5 p.157 draws more than one angular placement of the ±1
        # central weight. The strict matcher currently recognizes only
        # the between-tips variant; it does not invoke F3 to normalize.
        for sign in (+1, -1):
            with self.subTest(sign=sign):
                source = V4LabeledTree(WeightedPlaneTree((
                    Vertex(("f0", "e0", "e3"), (0, 0, 0)),
                    Vertex(("e0", "e1", "e2"), (sign, 0, 0)),
                    Vertex(("e1",), (2 * sign,)),
                    Vertex(("e2",), (2 * sign,)),
                    Vertex(("e3",), (3 * sign,)),
                )), (("f0", 3),))
                self.assertEqual(len(ring_signatures(source)), 1)
                self.assertEqual(strict_ring_caps(source, sign), ())
                untouched = abbreviate_disjoint_rings(source, sign)
                self.assertEqual(untouched.rings, (0,) * 5)
                self.assertEqual(untouched.labeled, source)


if __name__ == "__main__":
    unittest.main()
