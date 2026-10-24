"""A3.3: exact coefficient recurrence and certified algebraic radius bracket.

For the *synthetic* ring-adjacency grammar of ring_constraints.py,
sum over all grades by setting y=1. Let G=Z+P, with:

    Z = m + x(G^2+G^3)
    P = R*x(Z^2+Z^3)
    G = Z+P.

Eliminating Z gives G=H(x,G), where
    A=m+x*(G^2+G^3)
    H=A+R*x*(A^2+A^3).

The smooth implicit-function schema yields a square-root asymptotic
for fixed m>=1,R>=1. This is *combinatorial syntax only*.
"""
from __future__ import annotations

from fractions import Fraction


def total_guarded_coefficients(max_n: int, m: int, max_rings: int
                               ) -> tuple[int, ...]:
    """Independent exact univariate recurrence (no bivariate DP imports)."""
    if (any(type(v) is not int for v in (max_n, m, max_rings))
            or max_n < 0 or m < 1 or max_rings < 0):
        raise ValueError("invalid grammar parameters")
    z = [m]
    p = [0]
    result = [m]
    for n in range(1, max_n + 1):
        g = [z[k]+p[k] for k in range(n)]
        size = n-1
        def pair(seq: list[int]) -> int:
            return sum(seq[i]*seq[size-i] for i in range(size+1))
        def triple(seq: list[int]) -> int:
            return sum(seq[i]*seq[j]*seq[size-i-j]
                       for i in range(size+1) for j in range(size+1-i))
        zn = pair(g)+triple(g)
        pn = max_rings*(pair(z)+triple(z))
        z.append(zn)
        p.append(pn)
        result.append(zn+pn)
    return tuple(result)


def h_value(x, g, m: int, max_rings: int):
    """Evaluate the algebraic positive polynomial H(x,g)."""
    if (type(m) is not int or m < 1
            or type(max_rings) is not int or max_rings < 0):
        raise ValueError("invalid alphabets")
    z = m + x * (g*g+g*g*g)
    return z+max_rings*x*(z*z+z*z*z)


def h_derivative_g(x, g, m: int, max_rings: int):
    """H_g; characteristic singularity has H_g=1."""
    if (type(m) is not int or m < 1
            or type(max_rings) is not int or max_rings < 0):
        raise ValueError("invalid alphabets")
    z = m+x*(g*g+g*g*g)
    return x*(2*g+3*g*g)*(1+max_rings*x*(2*z+3*z*z))


# Exact discriminant factor of the degree-nine polynomial
#   g - H(x,g), with m=2 and R=1:
#   discriminant_g(g-H) = x^28 * sum P_COEFFICIENTS[k]*x^k.
# Coefficients from symbolic polynomial elimination (see A3.3 note).
P_COEFFICIENTS_2_1 = (
    -20502784, 818039424, 20330611824, 143077927472,
    321532769813, 70236192676, 6070016512, -1866365412,
    -70907536, 42625536, 4025344, 98304,
)


def discriminant_factor_2_1(x):
    """Evaluate exact degree-11 polynomial (supports Fraction)."""
    result = 0
    for c in reversed(P_COEFFICIENTS_2_1):
        result = result*x+c
    return result


def rho_bracket_2_1(iterations: int = 60) -> tuple[Fraction, Fraction]:
    """Bisection bracket isolating a positive discriminant root.

    An exact root bracket of the *polynomial*, not by itself a
    certified topology or analytic uniqueness proof. The attached
    mathematical note verifies its correspondence to the positive
    characteristic solution and explains all analytic hypotheses.
    """
    if type(iterations) is not int or iterations < 0 or iterations > 2048:
        raise ValueError("invalid bisection iterations")
    lo = Fraction(169943, 10_000_000)
    hi = Fraction(169944, 10_000_000)
    if not (discriminant_factor_2_1(lo) < 0
            < discriminant_factor_2_1(hi)):
        raise AssertionError("critical-root bracket must have opposite signs")
    for _ in range(iterations):
        mid = (lo+hi)/2
        if discriminant_factor_2_1(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def main() -> None:
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--max-n", type=int, default=10)
    p.add_argument("--m", type=int, default=2)
    p.add_argument("--R", type=int, default=1)
    args = p.parse_args()
    print("n,guarded_total")
    for n, count in enumerate(total_guarded_coefficients(args.max_n,args.m,args.R)):
        print(f"{n},{count}")
    if args.m == 2 and args.R == 1:
        lo, hi = rho_bracket_2_1()
        print(f"rho in ({float(lo):.13f},{float(hi):.13f})")
        print(f"mu in ({float(1/hi):.11f},{float(1/lo):.11f})")


if __name__ == "__main__":
    main()
