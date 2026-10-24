"""A4.2 formal normal form: exact confluence on the restricted slide graph."""
import unittest

from arborescent.placement_fibers import (
    placement_fiber_formula, placement_presentations,
    bounded_valence_chain_presentations,
)
from arborescent.ring_abbreviation import abbreviate_disjoint_rings
from arborescent.ring_normalizer import (
    RingNormalForm, normalize_ring_positions, ring_position_inversions,
    strict_ring_normal_key,
)
from arborescent.ring_relocation import explore_ring_slides
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree
from arborescent.internal_forks import InternalFork, Tip, compile_fork


class RingNormalizerTests(unittest.TestCase):
    def test_all_corolla_positions_have_one_normal_code(self):
        for sign in (-1, 1):
            for d, r in ((3, 1), (3, 2), (4, 2), (4, 3), (5, 2)):
                with self.subTest(sign=sign, d=d, rings=r):
                    family = placement_presentations(d, r, sign)
                    self.assertEqual(len(family), placement_fiber_formula(d, r))
                    keys = set()
                    for source in family:
                        inv = ring_position_inversions(source, sign)
                        cert = normalize_ring_positions(source, sign)
                        self.assertTrue(cert.verify())
                        self.assertEqual(len(cert.moves), inv)
                        self.assertEqual(ring_position_inversions(cert.normal, sign), 0)
                        self.assertEqual(
                            normalize_ring_positions(cert.normal, sign).moves, ()
                        )
                        self.assertEqual(
                            abbreviate_disjoint_rings(source, sign).rooted_code(),
                            abbreviate_disjoint_rings(cert.normal, sign).rooted_code(),
                        )
                        keys.add(strict_ring_normal_key(source, sign))
                    self.assertEqual(len(keys), 1)

    def test_every_reachable_state_has_same_normal_code(self):
        for sign in (-1, 1):
            initial = placement_presentations(4, 2, sign)[0]
            orbit = explore_ring_slides(initial, sign)
            self.assertTrue(orbit.complete)
            keys = {strict_ring_normal_key(s, sign) for s in orbit.incoming}
            self.assertEqual(len(keys), 1)
            for target in orbit.incoming:
                self.assertTrue(orbit.trace(target).verify())
                self.assertTrue(normalize_ring_positions(target, sign).verify())

    def test_bounded_valence_path_full_orbits(self):
        for sign in (-1, 1):
            for k in range(2, 5):
                with self.subTest(sign=sign, k=k):
                    original = bounded_valence_chain_presentations(k, sign)
                    full = explore_ring_slides(original[0], sign)
                    self.assertTrue(full.complete)
                    self.assertEqual(len(full.incoming), 4*3**(k-2))
                    self.assertEqual(
                        len({strict_ring_normal_key(s, sign)
                             for s in full.incoming}), 1
                    )
                    self.assertTrue(all(
                        normalize_ring_positions(s, sign).verify()
                        for s in original
                    ))

    def test_step_budget_is_fail_closed_and_counts_inversions(self):
        obj = placement_presentations(3, 2)[0]
        n = ring_position_inversions(obj, 1)
        self.assertGreater(n, 0)
        with self.assertRaisesRegex(ValueError, "budget"):
            normalize_ring_positions(obj, 1, max_moves=n-1)
        cert = normalize_ring_positions(obj, 1, max_moves=n)
        self.assertEqual(len(cert.moves), n)
        self.assertTrue(cert.verify())
        self.assertFalse(
            RingNormalForm(obj, cert.normal, 1, cert.moves[:-1]).verify()
        )

    def test_root_and_nonroot_v4_markings_not_dropped(self):
        obj = placement_presentations(3, 1)[0]
        altered = V4LabeledTree(
            obj.tree, tuple((b, 3 if b == "f1" else v)
                            for b, v in obj.labels)
        )
        self.assertNotEqual(strict_ring_normal_key(obj, 1),
                            strict_ring_normal_key(altered, 1))
        self.assertEqual(normalize_ring_positions(altered, 1).normal.labels,
                         altered.labels)

    def test_unsupported_ring_candidates_are_explicitly_rejected(self):
        # One-fork ±2/±1/±2 is a permissive A2.5 necessary motif,
        # but has no separate base for the source's strict cap.
        unsupported = compile_fork(InternalFork((Tip(2), Tip(2)), 1, 1))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            normalize_ring_positions(unsupported, 1)
        valid = placement_presentations(3, 1)[0]
        vertices = list(valid.tree.vertices)
        base = vertices[0]
        # Rotate storage cyclically to place a cap in anchor slot zero.
        vertices[0] = Vertex(base.slots[1:]+base.slots[:1],
                             base.sectors[1:]+base.sectors[:1])
        anchored_cap = V4LabeledTree(WeightedPlaneTree(tuple(vertices)),
                                     valid.labels)
        with self.assertRaisesRegex(ValueError, "anchored"):
            normalize_ring_positions(anchored_cap, 1)

    def test_mixed_sign_sector_and_parameter_guards(self):
        plus = placement_presentations(3, 1)[0]
        for value in (0, True, -2):
            with self.subTest(sign=value), self.assertRaises(ValueError):
                normalize_ring_positions(plus, value)
        for bound in (-1, True, 0.25):
            with self.subTest(bound=bound), self.assertRaises(ValueError):
                normalize_ring_positions(plus, 1, max_moves=bound)
        base = plus.tree.vertices[0]
        vertices = list(plus.tree.vertices)
        vertices[0] = Vertex(base.slots, (0, 0, 5, 0))
        weighted = V4LabeledTree(WeightedPlaneTree(tuple(vertices)), plus.labels)
        with self.assertRaises(ValueError):
            normalize_ring_positions(weighted, 1)

    def test_empty_ring_family_normalizes_to_identity(self):
        obj = compile_fork(InternalFork((Tip(2), Tip(3)), 0))
        cert = normalize_ring_positions(obj, 1)
        self.assertEqual(cert.normal, obj)
        self.assertEqual(cert.moves, ())
        self.assertTrue(cert.verify())
        self.assertEqual(ring_position_inversions(obj, 1), 0)


if __name__ == "__main__":
    unittest.main()
