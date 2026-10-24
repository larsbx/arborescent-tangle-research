"""A4.6 joint formal F1/F2 + ring-slide equivalence, never isotopy."""
import unittest

from arborescent.abbreviated_flips import (
    abbreviated_flip_orbit_key, flip_abbreviation_square,
)
from arborescent.flip_ring_quotient import (
    JointMove, apply_joint, explore_joint, flip_ring_normal_key,
    joint_successors, normalize_flip_ring,
)
from arborescent.placement_fibers import (
    bounded_valence_chain_presentations, placement_fiber_formula,
    placement_presentations,
)
from arborescent.ring_abbreviation import abbreviate_disjoint_rings
from arborescent.ring_relocation import RingSlide
from arborescent.tangle_labels import (
    V4LabeledTree, tangle_flip_f1, tangle_flip_f2, tangle_flip_group,
)
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


class JointFormalQuotientTests(unittest.TestCase):
    def test_all_small_corolla_joint_orbits_have_four_flip_sheets(self):
        for sign in (-1, 1):
            for d, rings in ((3, 1), (3, 2), (4, 2)):
                with self.subTest(sign=sign, d=d, rings=rings):
                    originals = placement_presentations(d, rings, sign)
                    orbit = explore_joint(originals[0], sign)
                    self.assertTrue(orbit.complete)
                    expected = 4 * placement_fiber_formula(d, rings)
                    self.assertEqual(len(orbit.incoming), expected)
                    self.assertTrue(set(originals).issubset(orbit.incoming))
                    self.assertEqual(
                        {t.label("f0") for t in orbit.incoming}, {0, 1, 2, 3}
                    )
                    keys = {flip_ring_normal_key(t, sign)
                            for t in orbit.incoming}
                    abbreviations = {
                        abbreviated_flip_orbit_key(
                            abbreviate_disjoint_rings(t, sign)
                        ) for t in orbit.incoming
                    }
                    self.assertEqual(len(keys), 1)
                    self.assertEqual(len(abbreviations), 1)
                    self.assertTrue(all(orbit.trace(t).verify()
                                        for t in orbit.incoming))

    def test_bounded_valence_joint_orbit_four_times_ring_orbit(self):
        for sign in (-1, 1):
            for k in (2, 3):
                with self.subTest(sign=sign, k=k):
                    originals = bounded_valence_chain_presentations(k, sign)
                    orbit = explore_joint(originals[0], sign)
                    self.assertTrue(orbit.complete)
                    self.assertEqual(len(orbit.incoming), 16*3**(k-2))
                    self.assertTrue(set(originals).issubset(orbit.incoming))
                    self.assertEqual(len({
                        flip_ring_normal_key(t, sign)
                        for t in orbit.incoming
                    }), 1)
                    self.assertEqual(len({
                        abbreviated_flip_orbit_key(abbreviate_disjoint_rings(t, sign))
                        for t in orbit.incoming
                    }), 1)
                    for t in originals:
                        self.assertTrue(orbit.trace(t).verify())

    def test_all_source_flip_variants_gauge_to_one_ring_normal(self):
        for sign in (-1, 1):
            for source in placement_presentations(3, 2, sign):
                codes = set()
                abbr_codes = set()
                for q in tangle_flip_group(source):
                    cert = normalize_flip_ring(q, sign)
                    self.assertTrue(cert.verify())
                    self.assertEqual(cert.flip_gauge.label("f0"), 0)
                    self.assertEqual(
                        len(cert.slide.slide_certificate.moves),
                        len(cert.slide.slide_certificate.moves)
                    )
                    self.assertEqual(len(q.tree.vertices),
                                     len(cert.normal.tree.vertices))
                    self.assertEqual(
                        tuple(sorted(v.weight for v in q.tree.vertices)),
                        tuple(sorted(v.weight for v in cert.normal.tree.vertices))
                    )
                    codes.add(cert.normal.labelled_rooted_code())
                    abbr_codes.add(abbreviated_flip_orbit_key(
                        abbreviate_disjoint_rings(q, sign)
                    ))
                self.assertEqual(len(codes), 1)
                self.assertEqual(len(abbr_codes), 1)

    def test_one_step_generators_preserve_full_formal_key(self):
        for sign in (-1, 1):
            source = placement_presentations(4, 2, sign)[0]
            target = flip_ring_normal_key(source, sign)
            for move, child in joint_successors(source, sign):
                with self.subTest(sign=sign, move=move):
                    self.assertEqual(flip_ring_normal_key(child, sign), target)
                    self.assertEqual(flip_abbreviation_square(child, sign),
                                     (True, True))
                    self.assertEqual(
                        abbreviated_flip_orbit_key(
                            abbreviate_disjoint_rings(child, sign)),
                        abbreviated_flip_orbit_key(
                            abbreviate_disjoint_rings(source, sign)),
                    )

    def test_root_storage_rotation_and_nonroot_v4_label_sensitivity(self):
        source = bounded_valence_chain_presentations(3)[0]
        vertices = tuple(
            Vertex(v.slots[1:]+v.slots[:1],
                   v.sectors[1:]+v.sectors[:1])
            for v in reversed(source.tree.vertices)
        )
        rotated = V4LabeledTree(WeightedPlaneTree(vertices), source.labels)
        self.assertEqual(flip_ring_normal_key(rotated, 1),
                         flip_ring_normal_key(source, 1))
        changed = V4LabeledTree(
            source.tree, tuple((b, v ^ 2 if b == "f2" else v)
                               for b, v in source.labels)
        )
        self.assertNotEqual(flip_ring_normal_key(changed, 1),
                            flip_ring_normal_key(source, 1))
        self.assertNotEqual(
            abbreviated_flip_orbit_key(abbreviate_disjoint_rings(changed, 1)),
            abbreviated_flip_orbit_key(abbreviate_disjoint_rings(source, 1)),
        )

    def test_bounded_graph_does_not_claim_completeness_on_cutoff(self):
        source = placement_presentations(3, 2)[0]
        partial = explore_joint(source, 1, max_states=2)
        self.assertFalse(partial.complete)
        self.assertEqual(len(partial.incoming), 2)
        self.assertTrue(partial.trace(source).verify())
        with self.assertRaises(ValueError):
            explore_joint(source, 1, max_states=0)
        with self.assertRaises(ValueError):
            explore_joint(source, 1, max_states=True)
        with self.assertRaises(ValueError):
            partial.trace(tangle_flip_f1(tangle_flip_f2(source, 0)))

    def test_unsupported_types_mixed_signs_and_overlapping_centers(self):
        source = placement_presentations(3, 1)[0]
        for maker in (
            lambda: JointMove("other"),
            lambda: JointMove("F1", RingSlide(1, 1, "f1")),
            lambda: JointMove("R-formal"),
        ):
            with self.assertRaises(ValueError):
                maker()
        with self.assertRaises(ValueError):
            apply_joint(source, object())
        with self.assertRaises(ValueError):
            normalize_flip_ring(source, -1)
        for invalid in (None, "not a tangle"):
            with self.assertRaises(ValueError):
                flip_ring_normal_key(invalid, 1)
        with self.assertRaises(ValueError):
            flip_ring_normal_key(source, 1, root="f99")
        with self.assertRaises(ValueError):
            flip_ring_normal_key(source, 1, root="e0")
        # Valid no-ring presentations remain admissible: the same
        # gauge comprises a four-action F1/F2 orbit.
        plain = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e0", "e1"), (0, 0, 0)),
            Vertex(("e0",), (2,)),
            Vertex(("e1",), (3,)),
        )), (("f0", 0),))
        self.assertTrue(normalize_flip_ring(plain, 1).verify())


if __name__ == "__main__":
    unittest.main()
