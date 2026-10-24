"""A4.4 source V4 F1/F2 equivariance; NOT an isotopy-class computation."""
import unittest

from arborescent.abbreviated_flips import (
    abbreviated_flip_f1, abbreviated_flip_f2,
    abbreviated_flip_gauge, abbreviated_flip_group,
    abbreviated_flip_orbit_key, flip_abbreviation_square,
)
from arborescent.placement_fibers import (
    bounded_valence_chain_presentations, placement_presentations,
)
from arborescent.ring_abbreviation import (
    Abbreviated, abbreviate_disjoint_rings,
)
from arborescent.tangle_labels import V4LabeledTree, tangle_flip_f1
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


class AbbreviatedFlipTests(unittest.TestCase):
    def test_both_commuting_squares_exhaustive_small_corollas(self):
        # Direct source F1/F2 on a weighted presentation, *then*
        # restrictively abbreviate must exactly equal abbreviation
        # followed by corresponding V4-labelled ring-numbered flip.
        for sign in (-1, 1):
            for d, rings in ((3, 0), (3, 1), (3, 2),
                             (4, 2), (4, 3), (5, 2)):
                for source in placement_presentations(d, rings, sign):
                    with self.subTest(sign=sign, d=d, rings=rings):
                        self.assertEqual(
                            flip_abbreviation_square(source, sign),
                            (True, True),
                        )
                        self.assertEqual(
                            flip_abbreviation_square(tangle_flip_f1(source), sign),
                            (True, True),
                        )

    def test_independent_bounded_valence_nonroot_f2_pivots(self):
        for sign in (-1, 1):
            for k in range(2, 6):
                for source in bounded_valence_chain_presentations(k, sign):
                    for bond in ("f0", f"f{k-1}"):
                        with self.subTest(sign=sign, k=k, root=bond):
                            self.assertEqual(
                                flip_abbreviation_square(source, sign, bond),
                                (True, True),
                            )

    def test_four_action_group_and_root_gauge(self):
        for sign in (-1, 1):
            obj = abbreviate_disjoint_rings(
                bounded_valence_chain_presentations(3, sign)[0], sign
            )
            group = abbreviated_flip_group(obj)
            self.assertEqual(len(group), 4)
            self.assertEqual(
                {x.labeled.label("f0") for x in group}, {0, 1, 2, 3}
            )
            self.assertTrue(all(x.rings == (1, 1, 1) for x in group))
            self.assertEqual(
                {abbreviated_flip_orbit_key(x) for x in group},
                {abbreviated_flip_orbit_key(obj)},
            )
            for element in group:
                self.assertEqual(abbreviated_flip_gauge(element).labeled.label("f0"), 0)
                self.assertEqual(
                    abbreviated_flip_f1(abbreviated_flip_f1(element)),
                    element,
                )
                self.assertEqual(
                    abbreviated_flip_f2(abbreviated_flip_f2(element, "f0"), "f0"),
                    element,
                )
                self.assertEqual(
                    abbreviated_flip_f1(abbreviated_flip_f2(element)),
                    abbreviated_flip_f2(abbreviated_flip_f1(element)),
                )

    def test_two_stage_abbreviation_retains_all_nonroot_v4(self):
        source = bounded_valence_chain_presentations(3, 1)[0]
        marked = V4LabeledTree(
            source.tree,
            tuple((bond, (i * 3) % 4) for i, (bond, _) in enumerate(source.labels))
        )
        self.assertEqual(flip_abbreviation_square(marked, 1), (True, True))
        obj = abbreviate_disjoint_rings(marked, 1)
        self.assertEqual(dict(obj.labeled.labels),
                         {"f0": 0, "f1": 3, "f2": 2})
        for transformed in abbreviated_flip_group(obj):
            self.assertEqual(transformed.rings, obj.rings)
            self.assertEqual(transformed.sign, obj.sign)
        altered = Abbreviated(
            V4LabeledTree(obj.labeled.tree,
                          (("f0", 0), ("f1", 1), ("f2", 2))),
            obj.rings, obj.sign,
        )
        self.assertNotEqual(abbreviated_flip_orbit_key(altered),
                            abbreviated_flip_orbit_key(obj))

    def test_ring_number_and_sign_are_not_flipped(self):
        obj = abbreviate_disjoint_rings(placement_presentations(3, 2)[0], 1)
        first = abbreviated_flip_f1(obj)
        second = abbreviated_flip_f2(obj)
        self.assertEqual(first.rings, obj.rings)
        self.assertEqual(second.rings, obj.rings)
        self.assertEqual(first.sign, 1)
        self.assertEqual(second.sign, 1)
        fabricated = Abbreviated(obj.labeled, (1,), 1)
        self.assertNotEqual(abbreviated_flip_orbit_key(fabricated),
                            abbreviated_flip_orbit_key(obj))

    def test_rooted_code_invariance_under_storage_rotation(self):
        full = placement_presentations(4, 1)[0]
        old = abbreviate_disjoint_rings(full, 1)
        vertices = tuple(
            Vertex(v.slots[1:] + v.slots[:1],
                   v.sectors[1:] + v.sectors[:1])
            for v in reversed(old.labeled.tree.vertices)
        )
        changed = Abbreviated(
            V4LabeledTree(WeightedPlaneTree(vertices), old.labeled.labels),
            old.rings[::-1], old.sign,
        )
        self.assertEqual(old.rooted_code(), changed.rooted_code())
        self.assertEqual(abbreviated_flip_orbit_key(old),
                         abbreviated_flip_orbit_key(changed))

    def test_invalid_pivots_and_wrong_input_classes(self):
        example = abbreviate_disjoint_rings(placement_presentations(3, 1)[0], 1)
        source = placement_presentations(3, 1)[0]
        for fn in (abbreviated_flip_f1, abbreviated_flip_f2,
                   abbreviated_flip_gauge, abbreviated_flip_group,
                   abbreviated_flip_orbit_key):
            with self.subTest(fn=fn.__name__), self.assertRaises(ValueError):
                fn(None)
        for root in ("e0", "f100", "", 123, True):
            with self.subTest(root=root), self.assertRaises(ValueError):
                abbreviated_flip_f2(example, root)
            with self.subTest(root=root), self.assertRaises(ValueError):
                flip_abbreviation_square(source, 1, root)
        with self.assertRaises(ValueError):
            flip_abbreviation_square(None, 1)

    def test_unsupported_adjacent_centers_are_not_declared_equivalent(self):
        # The sign-specific source exception remains unresolved and
        # must not be 'fixed' by ordinary flip/abbreviation bookkeeping.
        source = V4LabeledTree(WeightedPlaneTree((
            Vertex(("e0", "e1", "e2"), (0, 1, 0)),
            Vertex(("e0", "e3", "e4"), (0, -1, 0)),
            Vertex(("e1",), (2,)), Vertex(("e2",), (2,)),
            Vertex(("e3",), (-2,)), Vertex(("e4",), (-2,)),
        )), ())
        with self.assertRaises(ValueError):
            flip_abbreviation_square(source, 1)


if __name__ == "__main__":
    unittest.main()
