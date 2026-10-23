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
