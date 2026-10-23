# A2.2 — V4 boundary data, tangle flips and strict F3′ weight transport

**State:** A2.2 implements *three carefully delimited moves* on weighted
plane-tree **presentations of tangles**, not the complete plumbing calculus:
labelled (F1), labelled (F2), and the strict nonexceptional internal-edge
(F3′) case described below. The whole arborescent-tangle isotopy classifier
is **not** implemented.

## Source audit, including a correction to the previous milestone

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots and the
Classification and Symmetries of Arborescent Knots* (2016):

- §12.3, printed pp. 146–148 (PDF pages 157–159): pair flips F1/F2,
  weight transport F3, Proposition 12.6 and the boundary parametrisation
  discussion in Complement 12.7.
- §12.4, printed p. 155 (PDF page 166): the canonical-tree modified
  flips F3′ and F3″; **W** is preserved by the modified moves.
- §12.7.1, printed p. 166 (PDF page 177): the calculus of
  *arborescent tangles* includes free-bond labels in the Klein four-group
  V4. **F1** applies ζ at **every** free bond; **F2** applies ξ at
  **even** vertex distance and η at **odd** distance from its pivot.
- §12.7.2, printed pp. 166–167 (PDF pages 177–178):
  *tangle* canonicity uses (W), modified (S) and (P)/(N). The earlier
  pair-calculus **Free Bond condition (F) is ABSENT** for tangles,
  because boundary moves are omitted. The A1 `wfs_failures()`
  filter is a strict screen for the **pair** calculus and must **not**
  be used to certify tangle canonicity.
- §12.7.3, printed p. 167 (PDF page 178): modified
  **F3′**, when the transferred weight `a` is odd, reverses local cyclic
  orders in the selected attached subtree at **odd** graph distance
  from the weighted vertex. It additionally composes the labels of
  free bonds in that subtree with **ξ at even distance**, **η at odd**.
  If `a` is even, those reversals/relabelings do not occur.
- Theorem 12.21, printed p. 167 (PDF page 178) gives classification
  *under its canonical and topological hypotheses*, by F1, F2, F3′
  and R±. **We have not implemented the R± moves, full exceptional
  canonicity or proof that our sample syntax satisfies the theorem.**

Primary source:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

### Why A2.1's unlabelled quotient is not a tangle quotient

The previous [A2.1 calculation](flip-subgroup-a2-1.md) proves a valid
finite-orbit theorem for **unlabelled weighted-tree presentations of
arborescent pairs**. Its series `Q=(F+E+O+G)/4` counts syntax modulo
the graph-theoretic F1/F2 reflections. It is **not** an isotopy count
for boundary-parametrized tangles.

The tangle calculus requires explicit V4 data. Identifying the
Klein four-group with `(Z/2)^2`, write

```text
ι = 0, ξ = 1, η = 2, ζ = 3
composition αβ = α XOR β, so ξη = ζ.
```

The new `V4LabeledTree` pairs each *free bond* with exactly one element
of V4 and validates complete coverage, uniqueness and the integer range.
Its labelled marked-plane-tree code remembers free-bond labels,
unlike A2.1's unlabelled code. It remains a code for **presentations**
up to rooted planar isomorphism, not for all isotopic tangles.

## Theorem A2.2.1 — finite labelled flip group

The two elementary flip operations on the **V4-labelled** trees are
involutions and commute. On a connected tree with root free bond
`f0` at vertex `v0`, the four transformations are

```text
Identity       : no reflections; root label unchanged
F1             : all reflections; root label XOR ζ
F2(v0)         : even-distance reflections; root label XOR ξ
F1 ∘ F2(v0)    : odd-distance reflections; root label XOR η
```

**Proof.** Reversal of a cyclic ordering is an involution. On any vertex,
the reflections combine by parity, as in A2.1. The V4 changes compose
by XOR; F1 applies ζ globally while F2(v0) applies ξ/η on the
even/odd bipartition. Since XOR is commutative, F1 and F2 commute.
On the marked free bond `f0` the four changes are exactly
`{0,ζ,ξ,η}`. Thus, *unlike the unlabelled pair action*, the
labelled action is **free** on any presentation with at least one
marked free bond.

The transformed boundary labels matter: holding `f0` fixed at
`ι` while pretending these four transformations preserve that fixed
value is not legitimate. They generally leave that gauge slice.

### Corollary — an exact, surprising change to the previous count

For the A1 ordered 2/3-fork grammar, add **one free root bond f0 with
any of the four V4 labels** and no other free bonds. If `a_n(m)`
is the raw syntax count, then the number of V4-labelled presentations
is `4a_n(m)`. The F1/F2 group acts freely, so its orbit count is

```text
T_n(m) = (4 a_n(m)) / 4 = a_n(m).
T_m(x) = F_m(x) = m + x(F_m(x)^2+F_m(x)^3).
```

This is an exact **labelled presentation quotient**, *not* an isotopy
count under F3′/R± or a classification of arbitrary tangles.

For example (m=2), at n=2 we now distinguish:

```text
raw unlabelled binary/ternary grammar  : 192
unlabelled pair F1/F2 presentation     :  81
all V4 root-label choices              : 768
correctly labelled F1/F2 quotient      : 192
```

This comparison is important: *forgetting boundary labels changed the
meaning of the quotient and its enumerative predictions.*

The labelled quotient retains the A1 algebraic, square-root syntactic
growth *exactly* for the particular one-free-root grammar. The earlier
A2.1 result `Q ~ F/4` remains true for its **different** unlabelled
pair-syntax quotient. There is no contradiction.

## A2.2.2 — restricted, reversible F3′ local move

`tangle_flip_f3prime(source, pivot, branch_edge, from_side=...)`
requires:

1. An **internal** edge incident to the pivot; the chosen edge
   separates a well-defined attached subtree.
2. A single nonzero angular sector weight `a` at the pivot, adjacent
   to that edge, with the other adjacent sector zero (**strict W**).
3. `from_side` says which adjacent sector currently contains `a`.

The function moves `a` to the other side of the edge. For **odd** a,
it reverses the cyclic orders at vertices of the attached subtree
whose graph distance from the pivot is odd, and transforms each free
bond label in that subtree by `ξ` if its vertex is at even distance,
`η` if odd. For **even** a it moves the weight without reversing
or relabelling. Boundary free bonds outside the selected subtree are
unchanged. The opposite `from_side` undoes the move.

The strict guard prevents generalizing beyond the diagram's W-preserving
configuration. This code is **not** a complete implementation of
the F3′/F3″/R± rewriting calculus, and in particular no move closure
is claimed.

## Exact tests

- An independent V4 XOR-algebra sanity test and full four-root-label orbit.
- F1 and F2 involutions, commutation, independent parity action.
- Marked plane-tree renumbering and cyclic-sector rotation invariance.
- Multiple-free-bond invariants: relative boundary labels are not all
  independently mutable.
- The strict F3′ move with an **odd** weight, and free bonds at graph
  distances 1 and 2, checking η and ξ independently, plus a round trip.
- The **even** F3′ control, with no extra reflection/relabeling.
- Guard failures for invalid roots, pivots, weight placements and labels.
- Exhaustive small-grammar labelled quotient counts; independent
  comparison with `syntax_counts` and the A2.1 pair quotient.

Reproduce:

```sh
python -m unittest discover -s tests -v
```

## Still open

Actual canonicality checks for **tangles** (§12.7.2) including the
modified S and P/N exceptions; F3′ full generality; rings R±; marked
and unmarked free-bond equivalence; independent diagram/isotopy checks;
and a complete, certified topological quotient. A2 issue #4 remains
open. Nothing here proves a Catalan asymptotic for **all** arborescent
tangle or link isotopy classes.

## Corrective addendum (A2.3)

PR #6 was merged **before** hosted CI finished. Hosted Python 3.12
then failed because an A2.2 test fixture constructed a `Fork` with
more than one positional argument. A2.3 corrects the fixture and adds
independent regression coverage.

A more substantial mathematical limitation has also been identified:
**strict F3′ is never applicable to any tree generated by the A1
leaf-weight-only grammar**, even after all F1/F2 moves. Nonzero sector
weights occur exclusively at valence-one vertices, which have no pair
of distinct sectors across which the strict F3′ can transport weight.
The A2.2 implementation remains a useful operator on *other*
weighted-plane-tree presentations, but does not reduce A1's count.
See [A2.3](partial-move-closure-a2-3.md).
