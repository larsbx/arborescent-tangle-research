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
