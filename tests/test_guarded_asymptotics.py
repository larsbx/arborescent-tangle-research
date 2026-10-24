"""Independent exact and critical-system checks for A3.3 (syntax only)."""
import unittest
from fractions import Fraction

from arborescent.guarded_asymptotics import (
    P_COEFFICIENTS_2_1, discriminant_factor_2_1,
    h_derivative_g, h_value, rho_bracket_2_1,
    total_guarded_coefficients,
)
from arborescent.ring_constraints import guarded_ring_counts
from arborescent.ring_grammar import total_bounded_rings
from arborescent.weighted_trees import syntax_counts


class GuardedAsymptoticTests(unittest.TestCase):
    def test_univariate_independent_from_bivariate(self):
        for m in (1, 2, 3):
            for R in (0, 1, 2, 3):
                actual = total_guarded_coefficients(7, m, R)
                expected = tuple(sum(row) for row in guarded_ring_counts(7, m, R))
                self.assertEqual(actual, expected)
                self.assertEqual(actual[0], m)
                self.assertTrue(all(v > 0 for v in actual))
                if R == 0:
                    self.assertEqual(actual, syntax_counts(7, m))
                if R:
                    self.assertLess(actual[2], total_bounded_rings(2, m, R))

    def test_known_coefficients_for_two_leaf_one_ring_label(self):
        self.assertEqual(
            total_guarded_coefficients(10, 2, 1),
            (
                2, 24, 576, 20400, 811968, 35127552, 1601049600,
                75745871616, 3684334669824, 183099176982528,
                9256526405812224,
            ),
        )

    def test_exact_discriminant_root_bracket(self):
        lo, hi = rho_bracket_2_1(iterations=70)
        self.assertLess(discriminant_factor_2_1(lo), 0)
        self.assertGreater(discriminant_factor_2_1(hi), 0)
        self.assertLess(hi-lo, Fraction(1, 10**25))
        self.assertGreater(lo, Fraction(169943, 10**7))
        self.assertLess(hi, Fraction(169944, 10**7))
        mu_lo, mu_hi = 1/hi, 1/lo
        self.assertLess(mu_lo, Fraction(58844, 1000))
        self.assertGreater(mu_lo, Fraction(58843, 1000))
        self.assertGreater(mu_hi, mu_lo)
        self.assertEqual(len(P_COEFFICIENTS_2_1), 12)

    def test_critical_polynomial_characteristic_residuals(self):
        # Published numerical illustration of an exact algebraic system:
        # residuals are an approximation check, not the formal proof.
        rho = 0.01699430960763456
        g = 3.304691838094537
        self.assertLess(abs(h_value(rho,g,2,1)-g), 1e-12)
        self.assertLess(abs(h_derivative_g(rho,g,2,1)-1), 1e-12)
        self.assertGreater(h_derivative_g(rho,g+0.1,2,1), 1)
        self.assertLess(h_derivative_g(rho,g-0.1,2,1), 1)

    def test_guard_invalid_inputs(self):
        for args in ((-1, 2, 1), (True, 2, 1),
                     (1, 0, 1), (1, 2, -1), (1, 2, False)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                total_guarded_coefficients(*args)
        for bad in (-1, 2049, True, 1.2):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                rho_bracket_2_1(bad)
        with self.assertRaises(ValueError):
            h_value(1, 1, 0, 1)
        with self.assertRaises(ValueError):
            h_derivative_g(1, 1, 1, -1)


if __name__ == "__main__":
    unittest.main()
