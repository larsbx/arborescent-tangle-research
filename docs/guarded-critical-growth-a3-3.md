# A3.3 — Square-root law for the *locally constrained* ring grammar

**Theorem class:** rigorous asymptotic enumeration of **finite-state
weighted plane-tree syntax**. It is **not** an asymptotic for canonical
arborescent tangles or for the topological isotopy quotient.

## Context and fixed local constraint

[A3.1](local-ring-constraints-a3-1.md) assigns each 2/3-ary ordered
fork an integer ring number 0,...,R and excludes adjacent fork
vertices simultaneously carrying positive ring numbers. This is a
**synthetic** experimental constraint; it is **not** claimed to be
the Bonahon–Siebenmann §12.5 adjacent-ring exception.

At y=1, let `m>=1` be the finite leaf alphabet size and `R>=1`
the number of available positive ring labels. The exact two-state
series are

```text
G=Z+P
Z=m+x(G²+G³)
P=R*x(Z²+Z³).
```

Set

```text
A(x,g)=m+x(g²+g³)
H(x,g)=A(x,g)+R*x(A(x,g)²+A(x,g)³).
```

Then `G(x)=H(x,G(x))`. The polynomial equation is degree nine
in g when R>0, and has a unique formal power-series solution with
constant term m.

## Theorem A3.3.1 — exact smooth implicit singularity

For every fixed pair of integers `m>=1,R>=1`, there exist
strictly positive algebraic constants `rho=ρ_(m,R)`,
`mu=1/rho`, `gamma`, and `K=gamma/(2 sqrt(pi))` such that

```text
[x^n] G(x) ~ K*mu^n*n^(-3/2).
```

The critical algebraic system defining `rho,z,g` is

```text
z = m + rho*(g²+g³)
g = z + R*rho*(z²+z³)
1 = rho*(2g+3g²)*(1+R*rho*(2z+3z²)).
```

Choose its positive branch reached by the coefficientwise-minimal
solution as x increases from zero. The square-root amplitude is

```text
gamma = sqrt(2*rho*H_x(rho,g)/H_gg(rho,g)),
H_x = (g²+g³)*(1+R*rho*(2z+3z²)) + R*(z²+z³),
H_gg = rho*(2+6g)*(1+R*rho*(2z+3z²))
     + R*rho*(2+6z)*(rho*(2g+3g²))².
```

### Proof

Shift the unknown `u=G-m` and put
`Phi(x,u)=H(x,m+u)-m`. This is a polynomial with
**nonnegative** coefficients, `Phi(0,u)=0`, and genuinely nonlinear
dependence on u at positive x. Its coefficientwise-minimal
solution has `u(0)=0). Every coefficient at degree n≥1 is positive:
a tree with n binary forks, zero ring numbers and m-labelled leaves
is always admitted. Thus the solution is **aperiodic**.

A positive characteristic point exists: for sufficiently small
positive x, the scalar map `g -> H(x,g)-g` is negative at a
suitable fixed g>m (because H(0,g)=m); for any sufficiently large
positive x and all g>=m, H(x,g)>g because its positive
quadratic/cubic terms dominate. Each x>0 has H strictly convex in
g>=m and H increases strictly in x; as x rises, the minimum crosses
zero exactly once. At the crossing, the minimum is zero and
`H_g=1`. There `H_x>0` and `H_gg>0`, so the characteristic
singularity is nondegenerate.

Apply the smooth implicit-function schema after the shift by m
(Flajolet–Sedgewick, *Analytic Combinatorics*, Theorem VII.3,
Chapter VII.4). The square-root singularity is

```text
G(x)=g-gamma*sqrt(1-x/rho)+O(1-x/rho),
```

and aperiodicity makes it the unique dominant singularity.
Transfer of coefficients yields the asymptotic above. □

Reference:
https://sedgewick.io/books/analytic-combinatorics/
(Theorem VII.3 and the Meir–Moon smooth implicit-function schema).

**Consequence:** A3.1's exact matching upper bound showed that
the constraint reduces the exponential rate strictly below
`(R+1) mu_m` from independent ring decoration.
A3.3 proves the **constrained** total nevertheless belongs
to a square-root `n^(-3/2)` class. The resulting growth constant
is an algebraic solution to the displayed critical system, not a
chosen arbitrary sequence.

## Exact m=2, R=1 instance and algebraic elimination

Write `z=2+x(g²+g³)`. Then the univariate equation is

```text
P(x,g)=g-z-x(z²+z³)=0.
```

Taking the discriminant in g yields an exact factor
`x^28 D(x)` with

```text
D(x) =
  98304*x^11 + 4025344*x^10 + 42625536*x^9
  - 70907536*x^8 - 1866365412*x^7 + 6070016512*x^6
  + 70236192676*x^5 + 321532769813*x^4
  + 143077927472*x^3 + 20330611824*x^2
  + 818039424*x - 20502784.
```

The positive characteristic root has the numerical values

```text
rho ≈ 0.01699430960763456
z   ≈ 2.798927791206695
g   ≈ 3.304691838094537
mu  ≈ 58.8432259437452
K   ≈ 0.62282901645.
```

These decimal values are **illustrative**; the exact mathematical
objects are algebraic solutions. The standard-library module
`guarded_asymptotics.py` proves by exact `Fraction` arithmetic
that `D(x)` changes sign over a small positive rational interval,
and performs rational bisection for an arbitrary fixed number of
iterations. This certifies *a polynomial root bracket*, while the
positive characteristic branch is selected by the theorem above.

## Exact coefficient cross-check

For `m=2,R=1`, the totals for fork counts n=0..10 are

```text
2
24
576
20400
811968
35127552
1601049600
75745871616
3684334669824
183099176982528
9256526405812224
```

A new **independent univariate integer convolution recurrence**
uses the two equations for Z and P without referring to the A3.1
bivariate program; the tests compare the two implementations for
multiple m and R, check the exact low-order sequence and the
discriminant-root sign bracket, and verify the critical equation's
small numerical residual at the specified illustrative approximation.

```sh
python -m unittest discover -s tests -v
python -m arborescent.guarded_asymptotics --max-n 10 --m 2 --R 1
```

## Topological boundary

This result proves the Catalan square-root exponent for a **specified
bounded-local nonindependent ring-decoration grammar**. It does not
prove that any actual knot/tangle canonical subclass has this grammar.
The general R± move, exceptional ring arrangements, W/S/P/N,
independent 3-dimensional plumbing soundness, and the full
tangle-isotopy quotient remain open under issue #16.

A final topology-level theorem requires first a correct **source-
canonical normalization** and only then an independent fiber or
isotopy-class count.
