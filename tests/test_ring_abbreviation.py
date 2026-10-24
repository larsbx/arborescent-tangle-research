"""§12.5 abbreviation tests: NO false R± graph-deleting move."""
import unittest

from arborescent.internal_forks import InternalFork, Tip, compile_fork
from arborescent.ring_abbreviation import (
    Abbreviated, abbreviate_disjoint_rings, strict_ring_caps,
)
from arborescent.ring_signatures import ring_signatures
from arborescent.tangle_labels import V4LabeledTree, tangle_flip_f1
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def ring_fixture(sign=1):
    # Ring is attached to a separate base; root port f0 survives.
    return compile_fork(InternalFork((
        InternalFork((Tip(2 * sign), Tip(2 * sign)), sign, 1),
        Tip(3 * sign),
    ), 3 * sign, 2))


class RingAbbreviationTests(unittest.TestCase):
    def test_both_signs_source_picture_and_ring_number(self):
        for sign in (1, -1):
            original = ring_fixture(sign)
            caps = strict_ring_caps(original, sign)
            self.assertEqual(len(caps), 1)
            cap = caps[0]
            self.assertEqual((cap.base, cap.center, cap.tips),
                             (0, 1, (2, 3)))
            result = abbreviate_disjoint_rings(original, sign)
            self.assertIsInstance(result, Abbreviated)
            self.assertEqual(len(result.labeled.tree.vertices), 2)
            self.assertEqual(result.rings, (1, 0))
            self.assertEqual(result.labeled.tree.vertices[0].slots,
                             ("f0", "e3"))
            self.assertEqual(result.labeled.tree.vertices[0].sectors,
                             (0, 3 * sign))
            self.assertEqual(result.labeled.labels, original.labels)
            self.assertEqual(sum(result.rings), 1)

    def test_necessity_not_sufficiency_false_one_fork_pattern(self):
        source = compile_fork(InternalFork((Tip(2), Tip(2)), 1, 1))
        self.assertEqual(len(ring_signatures(source)), 1)
        self.assertEqual(strict_ring_caps(source, +1), ())
        unchanged = abbreviate_disjoint_rings(source, +1)
        self.assertEqual(unchanged.labeled, source)
        self.assertEqual(unchanged.rings, (0, 0, 0))

    def test_multiple_disjoint_rings_at_same_base(self):
        source = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e0", "e3", "e6"), (0, 0, 0, 0)),
            Vertex(("e0", "e1", "e2"), (0, 1, 0)),
            Vertex(("e1",), (2,)), Vertex(("e2",), (2,)),
            Vertex(("e3", "e4", "e5"), (0, 1, 0)),
            Vertex(("e4",), (2,)), Vertex(("e5",), (2,)),
            Vertex(("e6",), (3,)),
        )), (("f0", 3),))
        self.assertEqual(len(strict_ring_caps(source, +1)), 2)
        result = abbreviate_disjoint_rings(source, +1)
        self.assertEqual(len(result.labeled.tree.vertices), 2)
        self.assertEqual(result.rings, (2, 0))
        self.assertEqual(result.labeled.label("f0"), 3)
        self.assertEqual(result.labeled.tree.vertices[0].slots,
                         ("f0", "e6"))
        self.assertEqual(result.labeled.tree.vertices[0].sectors,
                         (0, 0))

    def test_adjacent_ring_centers_are_rejected(self):
        # The source explicitly treats the adjacent-ring exception
        # separately; two central ±1 nodes attached to each other
        # must not be abbreviated as ordinary disjoint rings.
        source = V4LabeledTree(WeightedPlaneTree((
            Vertex(("e0", "e1", "e2"), (0, 1, 0)),
            Vertex(("e0", "e3", "e4"), (0, 1, 0)),
            Vertex(("e1",), (2,)), Vertex(("e2",), (2,)),
            Vertex(("e3",), (2,)), Vertex(("e4",), (2,)),
        )), ())
        self.assertEqual(len(strict_ring_caps(source, +1)), 2)
        with self.assertRaisesRegex(ValueError, "rings attached to rings"):
            abbreviate_disjoint_rings(source, +1)

    def test_sector_placement_and_base_guard(self):
        wrong_center = compile_fork(InternalFork((
            InternalFork((Tip(2), Tip(2)), 1, 0), Tip(3)
        ), 3, 2))
        self.assertEqual(strict_ring_caps(wrong_center, +1), ())
        wrong_base = compile_fork(InternalFork((
            InternalFork((Tip(2), Tip(2)), 1, 1), Tip(3)
        ), 3, 0))
        self.assertEqual(strict_ring_caps(wrong_base, +1), ())
        self.assertEqual(abbreviate_disjoint_rings(wrong_base, +1).rings,
                         (0, 0, 0, 0, 0))

    def test_rooted_abbreviation_code_under_vertex_bond_renaming(self):
        source = ring_fixture()
        abbr = abbreviate_disjoint_rings(source, +1)
        old = abbr.labeled.tree
        relabeled = WeightedPlaneTree(tuple(
            Vertex(
                tuple("e" + str(100 + int(s[1:]))
                      if s.startswith("e") else s
                      for s in v.slots[1:] + v.slots[:1]),
                v.sectors[1:] + v.sectors[:1],
            ) for v in old.vertices[::-1]
        ))
        new = Abbreviated(
            V4LabeledTree(relabeled, abbr.labeled.labels),
            abbr.rings[::-1], +1,
        )
        self.assertEqual(new.rooted_code(), abbr.rooted_code())
        with self.assertRaises(ValueError):
            abbr.rooted_code("f999")

    def test_mirrors_and_free_boundary_labels(self):
        source = ring_fixture()
        mirror = tangle_flip_f1(source)
        # A strict single-oriented wedge may still be recognized after F1
        # because (0,1,0) is symmetric when the cap-tip order is reversed.
        self.assertEqual(len(strict_ring_caps(mirror, +1)), 1)
        self.assertEqual(abbreviate_disjoint_rings(mirror, +1).labeled.label("f0"), 3)

    def test_sign_and_type_guards(self):
        source = ring_fixture()
        for bad in (0, True, -3):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    strict_ring_caps(source, bad)
        with self.assertRaises(ValueError):
            strict_ring_caps(None, +1)
        with self.assertRaises(ValueError):
            Abbreviated(source, (0, -1, 0, 0, 0), +1)
        with self.assertRaises(ValueError):
            Abbreviated(source, (), +1)


if __name__ == "__main__":
    unittest.main()
