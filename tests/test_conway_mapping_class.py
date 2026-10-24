"""A5.2b: exact 4-puncture mapping class obstruction, not tangle isotopy."""
import unittest

from arborescent.conway_mapping_class import (
    PureConwayMappingClass, pure_twist, algebraic_slope_det,
    seam_ambiguity_witness, ConwaySeamAmbiguityWitness,
)
from arborescent.conway_braid_diagrams import (
    PuncturePermutation, proposed_seam_puncture_matching,
)
from arborescent.atomic_plumbing_contract import Port, atomic_plumbing_plan
from arborescent.placement_fibers import bounded_valence_chain_presentations
from arborescent.ring_abbreviation import abbreviate_disjoint_rings


class PureConwayMappingClassTests(unittest.TestCase):
    def test_infinite_pairwise_distinct_classes_fix_all_four_punctures(self):
        seen_keys = set()
        seen_slopes = set()
        for n in range(-24, 25):
            g = pure_twist(n)
            self.assertEqual(g.matrix, ((1, 2*n), (0, 1)))
            self.assertTrue(g.fixes_four_punctures())
            self.assertEqual(g.image_slope(0, 1),
                             (2*n, 1) if n > 0
                             else ((-2*n, -1) if n < 0 else (0, 1)))
            self.assertEqual(g.image_slope(1, 0), (1, 0))
            self.assertEqual(g.compose(g.inverse()).projective_key(),
                             pure_twist(0).projective_key())
            seen_keys.add(g.projective_key())
            seen_slopes.add(g.image_slope(0, 1))
        self.assertEqual(len(seen_keys), 49)
        self.assertEqual(len(seen_slopes), 49)

    def test_sl2_integer_group_composition_and_projective_sign(self):
        for n in range(-5, 6):
            for m in range(-5, 6):
                t = pure_twist(n).compose(pure_twist(m))
                self.assertEqual(t.projective_key(),
                                 pure_twist(n+m).projective_key())
                self.assertEqual(t.image_slope(0, 1),
                                 pure_twist(n+m).image_slope(0, 1))
                self.assertEqual(t.compose(t.inverse()).projective_key(),
                                 pure_twist(0).projective_key())
                for a, b in ((1, 0), (0, 1), (1, 1), (3, -2)):
                    self.assertEqual(
                        algebraic_slope_det((a, b), (0, 1)),
                        algebraic_slope_det(
                            t.image_slope(a, b), t.image_slope(0, 1)
                        ),
                    )
        minus_identity = PureConwayMappingClass(-1, 0, 0, -1)
        self.assertEqual(minus_identity.projective_key(),
                         pure_twist(0).projective_key())
        self.assertEqual(minus_identity.image_slope(1, 1), (1, 1))

    def test_sl2_and_puncture_action_are_distinct_data(self):
        for v4 in range(4):
            fixed_punctures = PuncturePermutation(v4)
            seen = set()
            for n in range(13):
                pure = pure_twist(n)
                self.assertTrue(pure.fixes_four_punctures())
                seen.add(pure.projective_key())
                self.assertEqual(
                    pure.branch_point_action(),
                    ((0, 0), (1, 0), (0, 1), (1, 1)),
                )
                # The non-pure, separately chosen V4 permutation is
                # identical independently of the hidden pure twist.
                self.assertEqual(
                    fixed_punctures.image_indices,
                    PuncturePermutation(v4).image_indices,
                )
            self.assertEqual(len(seen), 13)

    def test_real_a5_ports_produce_identical_puncture_bijections_but_distinct_classes(self):
        for sign in (-1, 1):
            for k in (2, 3, 4):
                full = bounded_valence_chain_presentations(k, sign)[0]
                plan = atomic_plumbing_plan(
                    abbreviate_disjoint_rings(full, sign)
                )
                self.assertEqual(plan.topology_status, "UNVERIFIED_GEOMETRY")
                for seam in plan.seams:
                    for v4 in range(4):
                        matching = proposed_seam_puncture_matching(
                            plan, seam.parent, seam.child, v4
                        )
                        for power in (-3, -1, 1, 3):
                            cert = seam_ambiguity_witness(
                                plan, seam.parent, seam.child, v4, power
                            )
                            self.assertTrue(cert.verify())
                            self.assertEqual(cert.matching, matching)
                            self.assertEqual(cert.first.image_slope(0, 1), (0, 1))
                            self.assertEqual(
                                algebraic_slope_det(
                                    cert.first.image_slope(0, 1),
                                    cert.second.image_slope(0, 1),
                                ),
                                2*abs(power),
                            )
                            with self.assertRaisesRegex(
                                NotImplementedError, "UNVERIFIED_GEOMETRY"
                            ):
                                cert.second.require_verified_conway_gluing()
                            with self.assertRaises(NotImplementedError):
                                cert.matching.certify_seam_isotopy()

    def test_generator_distinguishes_slopes_even_though_punctures_fixed(self):
        identity = pure_twist(0)
        twist = pure_twist(1)
        self.assertEqual(identity.image_slope(0, 1), (0, 1))
        self.assertEqual(twist.image_slope(0, 1), (2, 1))
        self.assertEqual(algebraic_slope_det((0, 1), (2, 1)), 2)
        self.assertNotEqual(identity.projective_key(), twist.projective_key())
        self.assertTrue(identity.fixes_four_punctures())
        self.assertTrue(twist.fixes_four_punctures())

    def test_bad_matrices_and_slope_inputs_fail_closed(self):
        for coeff in ((1, 1, 1, 1), (2, 0, 0, 1), (1, 1, 0, 1),
                      (1, 0, 1, 1), (True, 0, 0, 1),
                      (1, 0, 0, 1.0)):
            with self.subTest(coeff=coeff), self.assertRaises(ValueError):
                PureConwayMappingClass(*coeff)
        for power in (True, None, 1.4):
            with self.subTest(power=power), self.assertRaises(ValueError):
                pure_twist(power)
        for slope in ((0, 0), (2, 2), (True, 1), (1, 1.0)):
            with self.subTest(slope=slope), self.assertRaises(ValueError):
                pure_twist(1).image_slope(*slope)
        for pair in (((2, 0), (0, 1)), ((0, 0), (1, 1)),
                     (("1", 0), (0, 1))):
            with self.subTest(pair=pair), self.assertRaises(ValueError):
                algebraic_slope_det(*pair)
        with self.assertRaises(ValueError):
            pure_twist(0).compose("not a matrix")

    def test_invalid_a5_seams_and_zero_twist_fail_closed(self):
        plan = atomic_plumbing_plan(abbreviate_disjoint_rings(
            bounded_valence_chain_presentations(3)[0], 1
        ))
        seam = plan.seams[0]
        for n in (0, True, None, 1.2):
            with self.subTest(n=n), self.assertRaises(ValueError):
                seam_ambiguity_witness(plan, seam.parent, seam.child, 0, n)
        with self.assertRaises(ValueError):
            seam_ambiguity_witness(
                plan, seam.child, seam.parent, 0, 1
            )
        with self.assertRaises(ValueError):
            seam_ambiguity_witness(
                plan, Port((), 100), seam.child, 0, 1
            )
        with self.assertRaises(ValueError):
            seam_ambiguity_witness(
                None, seam.parent, seam.child, 0, 1
            )


if __name__ == "__main__":
    unittest.main()
