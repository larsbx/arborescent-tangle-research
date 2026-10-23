# A3 — Exact bivariate grammar of ring-numbered abbreviated syntax

**Status:** exact presentation enumeration with rigorous generating functions;
**not** a classification or isotopy count. It is built on [A2.5/A2.6
ring-numbered abbreviation](correct-ring-abbreviation-a2-6.md), which is
proposed in PR #12, not yet an R± transformation engine.

## Link to the source (what is and is not being counted)

Bonahon–Siebenmann §12.5, printed pp.157–158, introduces a new
weighted planar tree carrying a **nonnegative ring number at each
vertex**. Their abbreviated trees result, under stated canonicality
hypotheses and with a special adjacent-ring exception, from
*erasing* certain ring subtrees as a **representation change**. This is
not the same as the §12.4 modified R± location-changing move.

Source (diagrams checked): https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

To get a provable finite combinatorial model **without asserting a new
topological theorem**, we decorate only the fork vertices of the
previous A1 2/3-ary ordered plane-tree grammar by ring numbers
`r_v in {0,...,R}`. All terminal vertices have ring number zero,
the root is marked by one identity-V4 free bond, and all twist
weights are the same finite alphabet of A1 integer leaf weights.
This is a **restricted language of candidate abbreviated presentations**.
It is not established that every object meets the full positive/negative
canonicality criteria, is the image of a canonical source tree, or
represents a unique arborescent tangle.

## Theorem A3.1 — exact bivariate generating series

Fix `m>=1` distinct leaf weights and `R>=0`. Let `a_n(m)` be
the A1 number of syntactic 2/3-ary trees with `n` fork vertices,
given by

```
F_m(x) = m + x(F_m(x)^2 + F_m(x)^3).
```

Write `b_(n,q)^(R)(m)` for the number of ring-numbered trees with
`n` fork vertices and **total abbreviated ring number q**. Then

```
D_R(y) = 1 + y + ... + y^R
F_(m,R)(x,y) = m + x D_R(y) [F_(m,R)(x,y)^2+F_(m,R)(x,y)^3]
b_(n,q)^(R)(m) = a_n(m) * [y^q] D_R(y)^n
sum(q) b_(n,q)^(R)(m) = (R+1)^n a_n(m).
```

In particular `F_(m,R)` is algebraic over `Q(x,y)`.

**Proof.** Forks are the only ring-decorated vertices. There are
exactly `n` of them, and their ring numbers are chosen independently
from `0,...,R`, so for any fixed undecorated tree the total ring
number distribution is `D_R(y)^n`. The underlying tree has
`a_n(m)` choices. At a fork, the selected ring number contributes
`x y^r`, and its ordered child list contains either two or three
recursive subtrees. This proves the implicit polynomial equation and
the coefficient identities. Because `D_R(y)` is a polynomial, the
generating series is algebraic by construction. □

This is a genuine structural factorization, **not** the vacuous
`b_n=C_n a_n` decomposition: the decoration distribution and its
validity are fixed before counting.

## Theorem A3.2 — graded asymptotic effects

Assume `m>=2`, and let the established A1 tree singularity be

```
a_n(m) ~ K_m * mu_m^n * n^(-3/2), K_m,mu_m > 0.
```

**Sum over all allowed ring counts.** By the exact identity,

```
sum(q) b_(n,q)^(R)(m)
  ~ K_m * ((R+1)*mu_m)^n * n^(-3/2).
```

Thus finite independent ring numbers **change the exponential
growth constant** but **preserve the square-root exponent**.

**Fix total ring count `q` first.** For every fixed integer `q>=0`
and every fixed `R>=1`,

```
[y^q] D_R(y)^n = n^q/q! + O(n^(q-1))   (n -> infinity)
b_(n,q)^(R)(m)
  ~ (K_m/q!) * mu_m^n * n^(q-3/2).
```

**Proof.** At leading order, `q` distinct forks receive one ring
each, producing `binom(n,q)=n^q/q!+O(n^(q-1))` assignments. A repeated
ring at a fork uses at most `q-1` distinct occupied forks and
contributes only `O(n^(q-1))` assignments. Multiply by the A1
asymptotic. □

The prefactor **does not universally stay `n^(-3/2)` across different
gradings**. That is a concrete, falsifiable warning for future
Catalan–arborescent enumerations.

For unrestricted (finite **total**, but unbounded per-fork) ring
numbers, the exact fixed-q count is

```
b_(n,q)^(∞)(m) = a_n(m) * binom(n+q-1,q), for n>=1.
```

**Proof:** weak compositions of `q` into `n` ordered fork counts. At
`n=0`, only the formal `m` leaf bases contribute at `q=0`.

**Different size gradings.** Each *abbreviated* ring number stands
for a deleted **three-vertex ring subtree** in the source's
restricted abbreviation setting. If a virtual expanded-size
parameter is desired, substitute `y=x^3` in the bivariate
series. This is a bookkeeping statistic, not a proven crossing
number or a canonical geometric complexity.

## Reproducible exact data

For `m=2,R=1`:

| forks n | total q=0,1,...,n coefficients | total |
|---:|---|---:|
| 0 | 2 | 2 |
| 1 | 12, 12 | 24 |
| 2 | 192, 384, 192 | 768 |
| 3 | 4,080, 12,240, 12,240, 4,080 | 32,640 |

For `m=2,R=2,n=2`, the vector is
`192,384,576,384,192` (total `1,728`).

The module `arborescent/ring_grammar.py` implements exact integer
dynamic programming for the ring-number polynomial, the closed
counting identity, a stars-and-bars fixed-q formula, and an
independent exhaustive generator that compiles **every small-n
decoration** to a marked `Abbreviated` presentation containing
explicit ring numbers and V4 boundary labels. Regression tests
check uniqueness of rooted presentation codes under the defined
grammar, invariance under graph/edge renumbering and cyclic local
rotation, and explicit invalid-input guards.

```sh
python -m unittest discover -s tests -v
python -m arborescent.ring_grammar --max-n 6 --leaf-labels 2 --max-rings-per-fork 1
```

## Topological obligations before using this to count tangles

- Establish that the selected labelled abbreviated trees lie inside
  the source's correct positive/negative **canonical** classes,
  including ring-adjacent exceptions.
- Implement the actual R± location-changing operation; prove its
  compatibility with abbreviation under the relevant boundary/V4
  conventions, without confusing ring abbreviation with isotopy.
- Prove whether the grammar is closed under *all* other canonical
  equivalence moves, or explicitly restrict to a certified subfamily.
- Only then compare algebraic counts with isotopy classes and ask
  whether the shape exponent survives that *full* quotient.

The result is an exact taxonomy of **syntactic growth laws**, not
an unconditional Catalan theorem about all non-rational tangles.
