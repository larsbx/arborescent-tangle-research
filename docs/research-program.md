# Research program and conjectures (v0)

## Mathematical scope and definitions

Bonahon–Siebenmann already provides a classification framework for arborescent pairs and knots using almost canonical weighted planar trees. This program does **not** propose that existing classification as new.

Keep distinct: (1) freely generated planar binary syntax; (2) syntax after algebraic rewriting; (3) reduced **weighted arborescent plumbing trees** after the published calculus; and (4) ambient isotopy classes of *closed links* or boundary-fixed classes of *marked tangles*. No two spaces should be identified without a proof.

The current Montesinos M0 experiment is a **negative control**: associativity erases the raw Catalan shape factor. It does not provide evidence that isotopy classes universally have Catalan asymptotics.

The tree-size parameter `n` is always specified along with its model. Raw binary syntax complexity, reduced plumbing-tree size, and geometric crossing number are **different filtrations**.

## Precise targets

**C1 — Fixed reduced plumbing grammar.** Predeclare a finite alphabet of allowed twist/sector data, vertex valences, markings, roots and cyclic orders; explicitly encode the Bonahon–Siebenmann moves. Construct a sound and complete recognizer/normalization algorithm **relative to the established classification**. This is not a new tangle normal-form theorem. Failure cases include additional global moves, exceptional pieces, and boundary-marking dependence.

**C2 — Algebraic grammar generating series, conditional.** For a finite unambiguous context-free **tree** grammar with bounded local compatibility, prove its size generating series is algebraic using explicit functional equations. This claim applies to *syntax*, not automatically to a quotient by ambient isotopy. If the actual reduced grammar is not finite/local/unambiguous, document the obstruction instead.

**C3 — Square-root asymptotic, conditional.** For a fixed nondegenerate, aperiodic, strongly connected positive polynomial tree system satisfying the necessary criticality hypotheses, test/prove

```text
A_n ~ K * mu^n * n^(-3/2)
```

for its **reduced presentations**. Establish the analytic hypotheses and exclude linear, periodic, or singular exceptional grammar components. No universal `-3/2` statement is asserted for all arborescent link types.

**C4 — Isotopy quotient control, open.** For an explicitly fixed C1 grammar, test whether maximum fibers from certified reduced presentations to correctly specified link/tangle isotopy classes grow polynomially in `n`. Polynomial control preserves exponential growth rates, **not necessarily the polynomial prefactor** (`n^(-3/2)` could change). To transfer the full asymptotic, prove a stronger statement on fiber distributions.

**C5 — Boundary and geometry, exploratory.** Compare arborescent/twist-band pieces with polyhedral pieces in characteristic Conway splitting. This is **not** a planar/nonplanar-tree dichotomy: both kinds of pieces can occur in a tree decomposition. Study reduced complexity versus crossing number, and independent fixed-local Kauffman/skein representations only after the topological quotient is validated.

## Explicit negative controls

1. **Associativity:** `C_(k-1)` parenthesizations of the same ordered word collapse to a single sum. The raw-tree polynomial-fiber conjecture is false.
2. **Closure:** cyclic permutations/reversal belong to Montesinos **link** presentations; they cannot be applied blindly to boundary-fixed tangles.
3. **No programmable universality:** arbitrary unbounded sequences are not granted by inventing arbitrary invariants and alphabets after seeing the target. Declare bounded local data beforehand.
4. **Classification exceptions:** carry denominator conditions, integral twists, orientation, mirrors, and marked-boundary conventions explicitly.
5. **Two kinds of trees:** raw binary parenthesization trees are **not** the almost canonical weighted planar plumbing trees.

## Roadmap

| Stage | Deliverable | Objective acceptance | Status |
|---|---|---|---|
| M0 | Catalan/ordered/dihedral counts | Burnside = brute-force orbits, associativity negative control | Implemented |
| M1 | Classification audit | Theorem/section-level check of precise Montesinos link hypotheses and sample diagrams | Open |
| A1 | Weighted plumbing grammar | Vertex/edge/sector types, admissibility, local moves, canonical encoding | Open |
| A2 | Certified reduction | Critical pair/inverse move tests; distinct isotopic and nonisotopic controls | Open |
| A3 | Genuine reduced enumeration | Counts with explicitly proved/checked quotients and reproducible certificates | Open |
| A4 | Analytic growth | Derived functional equations, singularity and periodicity analysis | Open |
| A5 | Geometry/representations | Crossing/skein/Jones comparisons with independent verification | Open |

**Current verified status:** exact combinatorics and a raw-syntax obstruction; **not** a novel knot-theoretic classification, general isotopy counting engine, or proof of an arborescent growth law.
