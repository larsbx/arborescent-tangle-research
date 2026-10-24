# Arborescent tangle research

Research question: **which Catalan-type enumeration laws survive reduction and isotopy in explicitly specified arborescent tangle/link families?**

This repo begins with a deterministic **negative control**, not a proof of a general enumeration theorem.

## Three different counts

For an alphabet of m rational-tangle slopes and k ordered leaves:

1. Full ordered binary presentation trees: **C_(k-1) m^k**.
2. Presentations modulo formal associativity of horizontal sum: **m^k** ordered words. Catalan parenthesizations disappear.
3. Words modulo cyclic rotation and reversal: **B_m(k)** (bracelets), modeling cyclic **Montesinos link-closure presentations**, NOT isotopy classes of boundary-marked tangles.

The third counting series is logarithmic and non-algebraic for m>=2; Catalan syntax by itself does not guarantee a Catalan or algebraic series after quotienting.

## Run

Standard-library-only Python 3.11+:

```sh
python -m unittest discover -s tests -v
python -m arborescent.prototype --slopes 1/2 1/3 2/5 --max-k 8
```

The slopes must be distinct rational numbers in (0,1). Their values distinguish the atomic labels; the resulting presentation-orbit counts depend on the alphabet size, not the numerical values.

See [research program](docs/research-program.md), [exact prototype counts](docs/prototype.md), and [literature](docs/references.md). The program does not decide tangle isotopy, test non-rationality, calculate crossing numbers, or implement weighted arborescent-tree normal forms.

Next steps: establish the precise link-classification applicability and exceptional cases; specify and implement genuinely branching reduced weighted planar trees, including sector weights and the complete Bonahon–Siebenmann equivalence moves; test normal-form completeness; only then conjecture and attempt to prove a tree-type asymptotic for a specific nondegenerate grammar.

## A1 implemented: genuine branching weighted plane-tree syntax

[The A1 note](docs/weighted-tree-a1.md) specifies an explicit sector-weighted
planar rotation system with free bonds, checks **strict, nonexceptional**
Bonahon–Siebenmann (W), (F), (S) conditions, and proves exact counts for a
finite binary/ternary weighted-fork syntax grammar:

```text
F_m(x) = m + x*(F_m(x)^2 + F_m(x)^3)
a_n(m) = sum(j=0..n) (2n+j)! / ((n-j)! j! (n+j+1)!) * m^(n+j+1)
```

For m=2 and n=1..4, the syntax counts are 12, 192, 4080, 99264.
Exhaustive small-size enumeration, closed formulas, and W/F/S negative
controls have tests. This is still **not** the full Bonahon–Siebenmann
calculus, a topological canonical form, or an isotopy-class enumeration.

Run all tests: `python -m unittest discover -s tests -v`.

## A2.1: sourced elementary plumbing flips and exact restricted quotient

The [A2.1 flip-subgroup proof and experiment](docs/flip-subgroup-a2-1.md)
implements the Bonahon–Siebenmann **F1/F2** moves on weighted planar-tree
rotation systems, respecting angular sectors. On a connected tree these
moves generate a subgroup of size **at most four**: independent reflection
of the two vertex-bipartition classes. The quotient of the A1 finite
2/3-fork grammar is **exactly algebraic**, with fixed-point series
`E=m+x F(x²)(1+O)`, `O=m+x(E²+E³)`, `G=m+x F(x²)(1+G)`,
and `Q=(F+E+O+G)/4`. For m=2, the first positive-size F1/F2
presentation orbit counts are **9, 81, 1251, 26688, 669420**.
A square-root singularity proof shows the quotient's leading
asymptotic is one quarter of the raw syntax count.

These are **pair-isomorphism** moves from §12.3, not general
boundary-fixed-tangle isotopies. **A2 is not complete:** F3'/F3'', ring,
boundary and arithmetic moves, P/N canonicality and exceptional trees
are still outstanding. Reproduce with:

```sh
python -m arborescent.flips --max-n 6 --alphabet-size 2
python -m unittest discover -s tests -v
```

## A2.2: Klein-four boundary labels and a strict F3′ move

[The A2.2 tangle-calculus note](docs/v4-tangle-moves-a2-2.md) records a
critical distinction from A2.1: Bonahon–Siebenmann §12.7 requires **V4
labels on every free bond**. Its tangle F1 and F2 moves alter those labels
(ζ globally; ξ at even distance and η at odd distance). A restricted
W-preserving **F3′** weight transport is now implemented, including the
odd/even parity effect on the attached subtree and its boundary labels.

**New exact correction:** adding all four root-boundary labels makes
the F1/F2 action free on the A1 one-root-port grammar. Its correctly
V4-labelled presentation quotient has coefficients **a_n**, not the
unlabelled pair-presentation quotient **Q_n**. For m=2,n=2: raw=192,
unlabelled pair quotient=81, labelled tangle presentation quotient=192.
Neither is a full arborescent isotopy enumeration.

The §12.4 Free-Bond canonicity condition **F** pertains to the *pair*
calculus; §12.7.2 omits it for *tangles*. Tangle W/S/P/N and ring
equivalence remain open; the existing `wfs_failures()` must not be
interpreted as a tangle-canonicality certificate.

## A2.3: bounded replayable move search and new obstruction

[The A2.3 exact audit](docs/partial-move-closure-a2-3.md) fixes a
previously merged A2.2 test-fixture bug that caused hosted CI to fail.
It adds bounded exhaustive search under only implemented F1/F2/strict
F3′ moves on fixed indexed graphs, with replayable shortest move traces
and explicit incomplete status when the budget is reached.

**New exact obstruction:** The A1 grammar places nonzero weights only
at valence-one leaves. The strict F3′ weight-transport move requires two
distinct adjacent sectors at a weighted vertex. Therefore **F3′ is
inapplicable on every A1 syntax tree, even after F1/F2**. Its earlier
implementation must not be used to infer quotient behavior on A1.
Future nontrivial enumeration requires a different, predeclared grammar
with weighted internal forks and topology-valid tangle W/S/P/N conditions.

## A2.4: non-vacuous strict F3′ root-pivot quotient

[The exact root-pivot grammar and proof](docs/root-pivot-a2-4.md)
adds **one weighted internal root fork** to the A1 zero-weight-fork
language. Its strict root-F3′ transport is genuinely applicable,
unlike the old leaf-weight-only A1 grammar. With `m` distinct leaf
weights and `p` nonzero root weights, we have the *presentation*
generating functions:

```text
F = m + x(F² + F³)
H = p*x*(3F² + 4F³)        # all positions of the root weight
J = p*x*(F² + F³) = p*(F-m)  # quotient by root F3′ only
```

The quotient has an exact sector-zero normal form because the allowed
root-sector moves form a path. It preserves the same square-root
syntax asymptotic but does **not** establish a complete arborescent
isotopy-class enumeration or Bonahon–Siebenmann canonicality.
Reproduce with `python -m arborescent.root_pivot --max-n 6`.

## A2.5 — ring-signature obstruction

[The source-scoped necessary-pattern detector](docs/ring-signature-a2-5.md)
locates ±1-weight valence-three vertices with adjacent matching ±2
terminal weights. It does **not** implement R±. A1 and default A2.4
grammars have no ±1 vertex; since F1/F2/strict-F3′ preserve individual
vertex totals, they cannot acquire this necessary modified-ring motif.
The new unit-weight/±2-tip test alphabet offers a concrete next input
for a figure-exact ring rewrite with V4 boundary bookkeeping.

## A2.6 source correction: ring number is *not* a ring move

[The A2.6 ring-abbreviation audit](docs/correct-ring-abbreviation-a2-6.md)
corrects closed, **unmerged** PR #11. The source's modified R± moves
**relocate** an attached ±2–±1–±2 ring subtree in the plane; they do
**not delete its three vertices**. The source *separately* introduces
ring-numbered abbreviated trees (§12.5), erasing suitable **disjoint**
ring subtrees and storing their multiplicity at attachment vertices.

`ring_abbreviation.py` implements only a conservative, lossful
**notation change**, preserving weighted sectors, surviving V4 labels,
and explicit ring multiplicities. The exceptional case of one ring
attached to another is rejected, not silently normalized. No full R±
move, tangle isotopy classifier, or general Catalan quotient is claimed.

## A3: exact ring-numbered abbreviated tree grammar

[The A3 bivariate enumeration theorem](docs/ring-number-grammar-a3.md)
counts *syntactic* abbreviated presentations with an independent
ring number 0..R at each internal 2/3-ary fork, retaining explicit
V4-marked weighted-plane-tree codes. With `m` finite leaf weights:

```text
F_(m,R)(x,y) = m+x*(1+y+...+y^R)*(F_(m,R)^2+F_(m,R)^3)
[x^n y^q]F_(m,R) = a_n(m)*[y^q](1+y+...+y^R)^n
sum_q [x^n y^q]F_(m,R) = (R+1)^n*a_n(m)
```

Summing over ring counts retains the tree square-root exponent
`n^(-3/2)` while multiplying exponential growth by `R+1`.
Holding total ring number q fixed instead yields a
`n^(q-3/2)` polynomial prefactor. These are **combinatorial
presentation results**, not classification of canonical tangles or
a proven complete R± isotopy quotient.

## A3.1: finite-local ring compatibility alters growth

[The exact constrained grammar theorem](docs/local-ring-constraints-a3-1.md)
adds a **synthetic** (not source-canonical) local exclusion rule:
adjacent fork vertices may not both have positive ring numbers.
Its finite two-state algebraic equations are
`Z=m+x((Z+P)^2+(Z+P)^3)`,
`P=x(y+...+y^R)(Z^2+Z^3)`.
For `m=2,R=1,n=2` the guarded syntax count is **576**, versus
**768** under independent ring numbers. At *fixed total ring grade* q,
the leading factor remains `n^q/q!` (up to lower order);
when all grades are counted, a matching argument proves a **strict
exponential deficit** relative to the unrestricted grammar.
This is a falsifiable demonstration of finite-local sequence
classification, not a Bonahon–Siebenmann canonicality theorem.

## A3.2: exact exponential *presentation* fibers of ring abbreviation

[The explicit ring-location fiber theorem](docs/ring-placement-fibers-a3-2.md)
constructs a marked zero-sector base with `d` V4-labelled free legs
and `r` attached disjoint ±2–±1–±2 ring caps. There are exactly
`binom(d+r-1,r)` different rooted cyclic placement codes, but all
have **one** identical ring-numbered abbreviated code. For `d=r=n`,
this fiber is `binom(2n-1,n)~4^n/(2sqrt(pi n))`—exponential.
This proves **no polynomial-fiber bound for unreduced placements at
unbounded base valence**. It says nothing by itself about fibers of
already-reduced bounded-valence grammars or topological isotopy.

**Stronger A3.2 negative control:** an exact path-of-bases family has
**maximum bond valence four** and one ring cap at each of `k` bases.
Its two local cyclic cap positions per base produce **2^k** different
marked presentation codes, yet §12.5 ring-number abbreviation gives
**one** identical code with one ring at each base. Thus bounded
valence **alone** cannot give polynomial unreduced-placement fibers;
ring-location normalization is essential. This is not an isotopy count.

## A3.3: square-root theorem for the constrained ring grammar

[The A3.3 analytic and exact-algebraic proof](docs/guarded-critical-growth-a3-3.md)
shows the A3.1 **synthetic finite-local** ring exclusion retains
`K*mu^n*n^(-3/2)` total presentation growth for all fixed finite
leaf alphabet sizes `m>=1` and positive ring alphabet sizes `R>=1`.
Eliminating the two grammar states gives the scalar smooth implicit
function `G=H(x,G)`, with a positive nondegenerate critical pair.
For `m=2,R=1`, the exact critical discriminant contains a
specified degree-11 polynomial whose positive characteristic root
has `mu≈58.8432259437452`. Exact integer recurrence and rational
root-bracketing tests are included. This does **not** establish a
canonical tangle growth law or an isotopy quotient.

## A4.1 — source-motivated formal ring *relocation* (no cap deletion)

[The A4.1 proof and tested finite certificates](docs/ring-relocation-a4-1.md)
add a conservative, typed **adjacent ring-position slide** on
zero-sector weighted plane-tree presentations. It swaps a full
±2–±1–±2 ring subtree's stem with an adjacent ordinary bond **at
the same base**; it never deletes vertices or changes V4 labels,
and preserves the §12.5 abbreviated code. Overlapping rings, mixed
signs, nonzero base sectors, root-anchor crossings and ring-ring
swaps are rejected.

For the previously established A3.2 families, the finite formal
slide graph is connected on all binomially many corolla
placements and all `2^k` bounded-valence path placements, with
replayable shortest move witnesses. These are **presentation graph
moves**, not independently certified 3D R± tangle isotopies.
Complete Bonahon–Siebenmann Theorem 12.21 compatibility remains
open in issue #16.

## A4.2 — stable normalization of the strict formal ring slides

[The A4.2 proof and replay certificates](docs/ring-normal-form-a4-2.md)
orient the A4.1 formal adjacent ring/ordinary-bond slide as
`RO→OR`. The inversion count `#(R before O)` decreases exactly by
one on each step. Thus the finite slide system has a **unique
stable normal form** placing every ring after all ordinary
nonanchor bonds at its base, preserving their respective relative
orders. An implementation replays every typed move and checks
termination. For a base with `o` ordinary bonds and `r` identical
ring caps the rooted-presentation orbit contributes exactly
`binom(o+r,r)` positions; factors multiply over independent bases.

This is a mathematically exact **formal syntax normalizer**, not an
independently certified Bonahon–Siebenmann R± topological normal form
or knot/tangle isotopy classifier.

## A4.3 — rooted cyclic-storage invariant ring normal keys

[The A4.3 correction](docs/root-aligned-ring-normal-form-a4-3.md)
fixes a subtle presentation-isomorphism issue in A4.2:
`Vertex.slots[0]` is **not intrinsic** for a cyclic incident-bond
order. A cyclically rotated representation could be wrongly rejected
or normalized differently simply because a ring stem moved into
storage slot zero. A4.3 aligns every local cyclic order on the
unique **incoming bond toward the marked free root**, rotating its
sector weights simultaneously, before running A4.2's explicitly
anchored stable ring sort.

The root-alignment certificate distinguishes the **presentation
isomorphism** from the typed **formal ring slides**; tests verify
invariance under vertex/edge renumbering, cyclic starting slots,
both ring signs, arbitrary path/corolla placements and V4 labels.
Neither layer certifies actual tangle isotopy or a complete R±
classification.
