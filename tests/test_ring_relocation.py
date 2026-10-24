"""A4.1 strict picture-motivated *formal* R-position slide tests.

These are graph/sector/V4 replay proofs only, not 3D isotopy proofs.
"""
import unittest

from arborescent.placement_fibers import (
    bounded_valence_chain_presentations, placement_fiber_formula,
    placement_presentations,
)
from arborescent.ring_abbreviation import abbreviate_disjoint_rings, strict_ring_caps
from arborescent.ring_relocation import (
    RingSlide, apply_ring_slide, explore_ring_slides, ring_slide_successors,
)
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


class RingRelocationTests(unittest.TestCase):
    def test_single_cap_moves_without_deleting_anything(self):
        for sign in (-1, 1):
            originals = placement_presentations(3, 1, sign)
            source = originals[0]
            before_code = source.labelled_rooted_code()
            cap = strict_ring_caps(source, sign)[0]
            move = RingSlide(sign, cap.center, "f1")
            target = apply_ring_slide(source, move)
            self.assertEqual(target, originals[1])
            self.assertNotEqual(target.labelled_rooted_code(), before_code)
            self.assertEqual(source.tree.vertices[cap.center],
                             target.tree.vertices[cap.center])
            self.assertEqual(len(source.tree.vertices), len(target.tree.vertices))
            self.assertEqual(source.labels, target.labels)
            self.assertEqual(
                tuple(v.weight for v in source.tree.vertices),
                tuple(v.weight for v in target.tree.vertices),
            )
            self.assertEqual(
                abbreviate_disjoint_rings(source, sign).rooted_code(),
                abbreviate_disjoint_rings(target, sign).rooted_code(),
            )
            self.assertEqual(apply_ring_slide(target, move), source)

    def test_corolla_interleavings_exhaustively_connected(self):
        # Exactly as many *reachable indexed presentations* as
        # the exact binomial marked-placement family, for each sign.
        for sign in (-1, 1):
            for d, rings in ((3, 1), (3, 2), (4, 2), (4, 3), (5, 2)):
                with self.subTest(sign=sign, d=d, rings=rings):
                    originals = placement_presentations(d, rings, sign)
                    start = originals[0]
                    orbit = explore_ring_slides(start, sign)
                    self.assertTrue(orbit.complete)
                    expected = placement_fiber_formula(d, rings)
                    self.assertEqual(len(originals), expected)
                    self.assertEqual(len(orbit.incoming), expected)
                    self.assertEqual(set(orbit.incoming), set(originals))
                    for target in originals:
                        self.assertTrue(orbit.trace(target).verify())

    def test_bounded_valence_chain_orbits_have_size_two_power_k(self):
        for sign in (-1, 1):
            for k in range(2, 6):
                with self.subTest(sign=sign, k=k):
                    originals = bounded_valence_chain_presentations(k, sign)
                    orbit = explore_ring_slides(originals[0], sign)
                    self.assertTrue(orbit.complete)
                    self.assertTrue(set(originals).issubset(set(orbit.incoming)))
                    # At each interior degree-four base the cap has
                    # THREE formal positions relative to the other
                    # two nonanchor bonds (not merely the two in
                    # the old A3.2 binary test subfamily).
                    self.assertEqual(len(orbit.incoming), 4*3**(k-2))
                    self.assertEqual(
                        len({abbreviate_disjoint_rings(s, sign).rooted_code()
                             for s in orbit.incoming}), 1,
                    )
                    self.assertTrue(all(orbit.trace(s).verify() for s in originals))

    def test_cutoff_is_reported_and_not_a_closed_orbit(self):
        source = placement_presentations(3, 2)[0]
        orbit = explore_ring_slides(source, 1, max_states=2)
        self.assertFalse(orbit.complete)
        self.assertEqual(len(orbit.incoming), 2)
        self.assertTrue(orbit.trace(source).verify())
        other = placement_presentations(3, 2)[-1]
        if other not in orbit.incoming:
            with self.assertRaises(ValueError):
                orbit.trace(other)
        with self.assertRaises(ValueError):
            explore_ring_slides(source, 1, max_states=0)
        with self.assertRaises(ValueError):
            explore_ring_slides(source, 1, max_states=True)

    def test_source_guards_refuse_unspecified_local_configurations(self):
        source = placement_presentations(3, 1)[0]
        c = strict_ring_caps(source, 1)[0]
        for move in (RingSlide(1, c.center, "f0"),
                     RingSlide(1, c.center, "f2"),
                     RingSlide(1, c.center + 99, "f1"),
                     RingSlide(-1, c.center, "f1")):
            with self.subTest(move=move), self.assertRaises(ValueError):
                apply_ring_slide(source, move)
        with self.assertRaises(ValueError):
            apply_ring_slide(source, object())
        for sign, center, other in ((0, 1, "f1"), (True, 1, "f1"),
                                    (1, -1, "f1"), (1, True, "f1"),
                                    (1, 1, "")):
            with self.subTest(sign=sign, center=center, other=other):
                with self.assertRaises(ValueError):
                    RingSlide(sign, center, other)

    def test_mixed_sign_ring_centers_and_nonzero_base_sectors_rejected(self):
        # Mixed signs are outside a single positively/negatively
        # canonical source regime and are deliberately fail-closed.
        plus = placement_presentations(3, 1, +1)[0]
        src = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e0", "e3", "f1"), (0, 0, 0, 0)),
            Vertex(("e0", "e1", "e2"), (0, 1, 0)),
            Vertex(("e1",), (2,)), Vertex(("e2",), (2,)),
            Vertex(("e3", "e4", "e5"), (0, -1, 0)),
            Vertex(("e4",), (-2,)), Vertex(("e5",), (-2,)),
        )), (("f0", 0), ("f1", 0)))
        with self.assertRaisesRegex(ValueError, "mixed-sign"):
            explore_ring_slides(src, +1)
        # All base sectors must be zero (stricter than abbreviation).
        base = plus.tree.vertices[0]
        nonzero = V4LabeledTree(WeightedPlaneTree((
            Vertex(base.slots, (0, 0, 4, 0)), *plus.tree.vertices[1:]
        )), plus.labels)
        with self.assertRaises(ValueError):
            apply_ring_slide(nonzero, RingSlide(1, 1, "f1"))

    def test_ordinary_bond_only_cap_cap_neighbor_refused(self):
        source = placement_presentations(3, 2, +1)[0]
        self.assertEqual(source.tree.vertices[0].slots[:3], ("f0", "e0", "e3"))
        with self.assertRaisesRegex(ValueError, "adjacent cap"):
            apply_ring_slide(source, RingSlide(1, 1, "e3"))
        allowed = ring_slide_successors(source, +1)
        self.assertTrue(allowed)
        self.assertTrue(all(m.past_bond.startswith("f") for m, _ in allowed))


if __name__ == "__main__":
    unittest.main()
