# A2.5 — modified-ring eligibility obstruction (source-scoped negative result)

## Why a guard is necessary

The arborescent-tangle classification theorem of Bonahon–Siebenmann
(§12.7.3, Theorem 12.21, printed p.167) includes **modified ring
moves R+ or R−** in addition to F1, F2 and F3′.

The actual diagrams for the two ring moves appear in §12.4, printed
p.155 (PDF zero-based page 166):
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The **positive** diagram involves a weighted valence-three junction
with central weight **+1** and two terminal **+2** weights.
The **negative** counterpart reverses those signs.

**This PR does not implement those moves.** An R± move also requires
the precise sector placements, attachment patterns, and boundary
bookkeeping of the source's diagrams. Extracting a graph pattern is
not enough to prove topological soundness. Accordingly the new code
only reports a deliberately permissive **necessary ring signature**:
a weight ±1 valence-three junction with two ±2 weighted valence-one
neighbors and a third stem bond. The signature is an overapproximation;
presence does *not* establish R± applicability or knot/tangle isotopy.

## Proposition A2.5.1 — weight-spectrum obstruction

Let Γ be a weighted plane-tree presentation such that no vertex has
total integer weight +1 or −1. Then:

1. Γ has **no necessary modified-ring signature** of the R± diagrams.
2. No finite sequence of the **currently implemented** F1, F2, and
   strict F3′ moves can create such a signature.

**Proof.** The source's R± diagrams require a central vertex of
weight ±1. This necessary weight is missing, proving (1).
Every current move preserves the graph and the **total weight at each
indexed vertex**: F1 and F2 reverse cyclic sector order, while strict
F3′ moves a weight between sectors of the *same* vertex and optionally
reverses other cyclic sector orders. V4 free-bond relabeling does not
change integer weights. Thus the absence of unit-weight vertices is
inductive under all implemented moves, proving (2). □

This does **not** claim that R± themselves preserve vertex weights:
they are not implemented and may change the local graph presentation.

## Applied to existing finite grammars

| Model | Internal weights | Terminal weights | Necessary ring signature |
|---|---|---|---|
| A1 leaf-weight-only grammar | 0 | 2,3 | Impossible |
| A2.4 root pivot, PR #8 | 3,4 at root; 0 elsewhere | 2,3 | Impossible |
| A2.4 internal-fork default, PR #9 | 0,3 | 3 | Impossible |

The last two grammars successfully exercise strict F3′ in selected
cases, but neither permits R± motifs. Consequently they cannot yet
support an experimental claim about the **full** Theorem 12.21
quotient: adding a correct R± implementation would have no
applicable move *within their existing partial-move orbits*.

For PR #8 this is an algebraic statement about its declared alphabet
and generator; the present PR does not need to import that feature
branch.

## Minimal ring-ready extension — candidate, not certified move

The PR #9 binary-fork generator allows a finite internal-weight alphabet.
Choosing **weights=(0,1)** and **tips=(2,)** gives the first small
syntactic motifs of the required +1/2/2 type. At n=1 internal fork,
with no additional free ports, its grammar generates exactly four
weighted presentations:

- one weight-zero central fork, which is ring-ineligible;
- three presentations carrying the central +1 weight in one of the
  three possible angular sectors, all flagged as **candidates** by
  the permissive local detector.

The analogous negative candidate uses root weight −1 with two tips
−2. These examples pass the existing **branched tangle W/S/P/N**
screen with the matching sign. That screen is incomplete for general
exceptional cases and is **not** an independent topological guarantee.

This is an explicit next experimental input, not an R± orbit count.

## Exact, reproducible tests

Run:

```sh
python -m unittest discover -s tests -v
```

Tests check positive/negative necessary patterns, sign-sensitive tip
requirements, candidate invariance under F1/F2 and weighted-plane-tree
renaming, exhaustive default internal-fork grammar absence at n=1,2,
the unit-weight obstruction on every state in a complete nontrivial
strict-F3′ fixed-graph orbit, and the four-case minimally extended
ring-ready grammar.

## What remains open

Before implementing the **actual** ring move, the next PR must:

1. Transcribe both source diagrams into *typed* local pre/post
   weighted-tree patterns, preserving angular-sector locations and
   V4 labels. Prove each rewrite sound relative to plumbing, not
   merely syntactically invertible.
2. Handle both matching signs and the full source hypotheses,
   including the stem attachment and boundary parametrizations.
3. Supply nontrivial replayable positive/negative examples, inverses,
   graph-isomorphism controls, and changes to the *full* restricted
   orbit count, with no truncation or topological overclaims.
4. Only then investigate whether a **canonical** ring-ready
   grammar retains Catalan-type asymptotic growth after *all*
   permitted moves.

This completes a **negative control and preparation step**, not the
modified-ring implementation or the complete isotopy quotient.
