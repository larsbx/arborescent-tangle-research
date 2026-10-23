# A2.1 — Certified F1/F2 move subsystem and exact finite quotient

**Scope:** Bonahon–Siebenmann's **two elementary flips** (F1) and
(F2) acting on the weighted planar-tree presentation structure.
Not a complete topological normal form, nor a general isotopy decision
procedure. The manuscript's statement is **isomorphism of associated
arborescent pairs**, not pointwise boundary-fixing of tangles.

## Provenance and move definitions

Francis Bonahon and Laurence C. Siebenmann,
*New Geometric Splittings of Classical Knots and the Classification and
Symmetries of Arborescent Knots*, version 26 January 2016, **§12.3,
printed pp. 146–149** (PDF pp. 157–160), Proposition 12.6 and its
proof on printed p. 148. Source:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The source specifies:

- **(F1):** reverse the cyclic order at every vertex;
- **(F2):** reverse it at a chosen vertex and at every vertex at even
  graph distance from that vertex.

A vertex with cyclic bond slots `(s0,s1,...,sd-1)` has sector weights
`(a0,a1,...,ad-1)`, where `ai` lies between `si` and `s(i+1)`.
A literal reflection fixing the anchor bond `s0` is

```text
slots    -> (s0,sd-1,...,s1)
sectors  -> (ad-1,...,a0)
```

**Both** lists must transform; independently reversing just the bond
slots changes the weights' placement and is not the stated move.

The implementation (`arborescent/flips.py`) preserves the integer
sector multiset and total weight at each vertex, undirected edge
incidence, number of free bonds, and the W/F/S structural screens.
Its comparisons employ `rooted_plane_code(f0)` *after* applying moves:
this is a **marked presentation-isomorphism** code, not a proof of
boundary-fixed isotopy.

## Theorem A2.1 — the elementary flip subgroup has order at most four

Let `Γ` be a connected weighted planar tree and bipartition its
vertices as `V0 ∪ V1` according to distance parity from a chosen root.
Write `R0` and `R1` for the operations reflecting exactly the
vertices in `V0` and `V1`, respectively. Since the operations act
on disjoint vertices and reverse each cyclic order, we have

```text
R0^2 = R1^2 = id,
R0 R1 = R1 R0,
F1 = R0 R1.
F2(pivot) = R0 or R1 according to the pivot parity.
```

**Proof.** In a connected tree, the vertices at even distance from a
pivot are exactly the pivot's bipartition class. Reflecting at any
vertex twice restores its bond and sector orders. Reflections at
different vertices commute. The claims follow, and the group generated
by (F1) and all (F2) is a quotient of `C2 × C2`. Its orbits have at
most four presentations (possibly fewer after marked plane-tree
isomorphisms).

**Consequences:** quotienting by **only F1/F2** cannot erase an
exponential number of distinct weighted-tree presentations. This does
*not* bound fibers after (F3'), (F3''), modified rings, other plumbing
moves, exceptional splitting, or global link equivalence.

## Theorem A2.2 — exact fixed-point and Burnside series

Apply the A1 **ordered 2/3-fork grammar** with `m≥1` distinct allowed
leaf weights; the zero-fork leaf terms are formal bases. Let

```text
F(x) = m + x(F(x)^2 + F(x)^3)
S(x) = F(x^2)
```

and let `E(x)`, `O(x)`, `G(x)` count trees fixed under reflecting
forks at even depth, at odd depth, and at both depth parities,
respectively. Then

```text
E = m + x S(1 + O)
O = m + x(E^2 + E^3)
G = m + x S(1 + G)
Q = (F + E + O + G)/4
```

where `Q(x)` is the exact generating function for rooted **syntax**
orbits under the elementary flip subgroup.

**Proof.** Under an involution reversing a node's children, a binary
fork is fixed precisely when its second child is the transformed first
child (enumerated by `F(x²)`). A ternary fork additionally has a
middle child fixed under the child-depth-parity transform. A fork that
is *not* reversed has independent ordered children, each separately
fixed under the child transform. The equations follow recursively,
with base `m` fixed leaves. Burnside's lemma averages fixed-point
counts of the four group elements. For one-vertex formal bases, all
four operations act trivially, so the formula still applies.

The system is algebraic: `F(x²)` is algebraic, `G` is rational
in `F(x²)`, and eliminating `O` gives a polynomial equation in
`E` over `Q(x,F(x²))`. Thus **Q is algebraic over Q(x)**.

### Exact integer results

| forks n | raw m=2 | F1/F2 orbits m=2 | raw m=3 | F1/F2 orbits m=3 |
|--:|--:|--:|--:|--:|
| 0 | 2 | 2 | 3 | 3 |
| 1 | 12 | 9 | 36 | 24 |
| 2 | 192 | 81 | 1,188 | 432 |
| 3 | 4,080 | 1,251 | 52,164 | 14,688 |
| 4 | 99,264 | 26,688 | 2,623,428 | 678,132 |
| 5 | 2,614,656 | 669,420 | 142,863,588 | 36,032,040 |

The fixed-point equations give **exact arbitrary-precision integer
recurrences**; separately, the implementation generates every small
ordered syntax tree, compiles to weighted planar rotations, applies
the four moves and quotients marked presentation codes. A third
independent check mirrors nested ordered syntax trees directly.

## Theorem A2.3 — preserved square-root growth for the limited quotient

For fixed integer `m≥2`, let `τ>m` be the positive solution of

```text
2τ² + (1-3m)τ - 2m = 0,
ρ = 1/(2τ+3τ²),  μ = 1/ρ,
c = sqrt(2(τ²+τ³)/(2+6τ)).
```

Then

```text
[x^n]F(x) ~ c/(2sqrt(π)) * μ^n * n^(-3/2),
[x^n]Q(x) ~ c/(8sqrt(π)) * μ^n * n^(-3/2).
```

**Proof sketch.** The univariate positive implicit-tree equation for
`F` has its first critical pair `(ρ,τ)` satisfying
`τ=m+ρ(τ²+τ³)` and `1=ρ(2τ+3τ²)`; the grammar is aperiodic.
The smooth implicit function/transfer theorem yields a square-root
singularity with coefficient `c`. At `x=ρ`, the series
`S=F(x²)` is analytic (`0<ρ<1`). The fixed-point series satisfy
coefficientwise `E,O,G≤F`. The Jacobian spectral-radius square for
the coupled E/O system is

```text
J² = ρ² S(ρ) (2E(ρ)+3E(ρ)²)
   ≤ ρ² τ (2τ+3τ²)
   = ρτ = 1/(2+3τ) < 1.
```

Also `ρS(ρ)<1`, so `G` is regular. Standard positive-majorant
implicit-function arguments extend these fixed-point functions
holomorphically across the `|x|=ρ` circle, giving a strictly larger
radius of convergence than `F`. Thus their coefficients are
exponentially smaller than `[x^n]F`. Since
`Q=(F+E+O+G)/4`, the asserted `1/4` leading coefficient follows.

**Interpretation:** the square-root universality class is preserved
by this particular finite, explicitly sourced subgroup; it is **not**
thereby proved to survive the full topological equivalence calculus.
In particular the earlier Montesinos parenthesization negative control
remains intact.

## Reproduction and remaining obligations

```sh
python -m unittest discover -s tests -v
python -m arborescent.flips --max-n 6 --alphabet-size 2
```

Still open in [A2 issue #4](https://github.com/larsbx/arborescent-tangle-research/issues/4):
modified (F3'), (F3''), (R±), boundary/arithmetic moves; (P)/(N)
and exceptional canonicity; the underlying plumbing pair and marked
boundary requirements; independent tests of the pair-isomorphism
implication; topological completeness and true isotopy-class growth.

## Addendum: do not misidentify the unlabelled quotient with a tangle quotient

The A2.1 system Q counts **unlabelled weighted-tree presentation orbits**.
Bonahon–Siebenmann §12.7.1 prescribes additional Klein-four V4
free-bond labels for tangles, with explicit relabeling under F1/F2.
The correct labelled counting experiment is derived in
[A2.2](v4-tangle-moves-a2-2.md). With one free root bond taking all
four V4 labels, the F1/F2 action is free, so the labelled presentation
quotient counts a_n rather than Q_n. This does *not* alter the A2.1
theorem about its originally specified unlabelled quotient.
