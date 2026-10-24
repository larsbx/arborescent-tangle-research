# A4.4 — Equivariance of ring-numbered abbreviation under V4-labelled flips

**Status:** exact finite group-action statement for the **implemented
restricted presentation domain**. It checks the Bonahon–Siebenmann F1
and F2 bookkeeping on ring-numbered abbreviated trees, but does
**not** independently prove 3D tangle isotopy or implement F3′ on
abbreviated data.

## Primary-source contract

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots
and the Classification and Symmetries of Arborescent Knots*,
§12.7.1 printed p.166 (PDF zero-based p.177) defines tangle flips
with **V4-label changes** at every free bond:

- **F1** reverses all cyclic orders and applies ζ to all free-bond
  labels.
- **F2**, pivoted at a vertex, reverses cyclic orders at even graph
  distance and applies ξ to free labels at even distance and η to
  free labels at odd distance. The three nonzero Klein-group
  elements obey ξη=ζ.

Section §12.7.3 printed p.167 (PDF zero-based p.178) gives the
canonical tangle classification involving F1, F2, F3′ and modified
R±. Its paragraph immediately following Theorem 12.21 explains
that the §12.5 **abbreviated** representation allows R± to be
eliminated, leaving F1/F2/F3′ moves, under the theorem's canonical
hypotheses. This is a source-level theorem, **not** proof that our
arbitrary syntactic grammar is canonical.

Sections §12.5 printed pp.157–158 (PDF zero-based pp.168–169)
define the ring numbers as counts attached to surviving vertices
and include an exceptional adjacent-ring case.

Primary source:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

## The A4.4 commuting-square theorem

Take a V4-labelled weighted plane tree T in the **strict**
A2.6 abbreviation domain, with:

1. Disjoint and uniformly signed ±2—±1—±2 ring caps, with each cap's
   central weight in the recognized sector between terminal tips;
2. A separate attachment base for every cap, originally of valence
   at least three, and two zero angular sectors adjacent to its stem;
3. No ring-center-to-ring-center attachment (including mixed signs);
4. All remaining free bonds explicitly assigned V4 labels; and
5. An **existing marked free bond f** used as the F2 pivot's reference.

Let A_s(T) be the existing restricted, index-preserving
ring-numbered abbreviation: delete caps *only in the representation*
and place their number at the surviving attachment base.

Define F1* and F2* on the abbreviated representation by applying the
same source-defined rotation reflection and V4 relabelling to the
surviving weighted tree **without moving its integer ring numbers**.

**Proposition A4.4.1.** For each T satisfying these source-pattern
guards, and for a free-bond pivot f on a surviving vertex,

```text
 A_s(F1(T)) = F1*(A_s(T)),
 A_s(F2_f(T)) = F2*_f(A_s(T)).
```

Both equalities are as **indexed ring-numbered presentation data**
(so they also hold under rooted presentation-isomorphism coding).

**Proof.** The central ring cap and its ±2 tips occupy a pendant
subtree attached by exactly one internal stem. Removing such a
subtree changes neither the distances nor their parity between
surviving vertices. Hence the even/odd vertex bipartition used by
F2 restricts identically from T to A_s(T), with the pivot specified
by the same surviving free bond rather than a possibly renumbered
integer vertex index.

Reversal of a cyclic slot list commutes with deleting a cap stem
and merging its two **zero** neighboring angular sectors: since
those two values are zero and all other sectors retain their
relative cyclic positions, performing these two operations in
either order yields the same reflected cyclic residual list.
The ring cap's center weight (0,±1,0) is reflected to the same
recognized outer wedge; thus the strict recognizer remains
applicable after either flip. Consequently each retained base
has the same total ring count before and after the flip; that
count stays attached to that surviving vertex.

For every surviving free bond, the parity of its owner's distance
from the F2 pivot is unchanged when pendant caps are removed.
Its V4 label is therefore XOR'd by exactly the same ξ or η
on both routes, while F1 XORs ζ on both routes. The original
integer weights, surviving cyclic sectors, V4 labels and
per-base ring numbers agree. □

The assertion is conditional on the **strict recognized class**;
it is not a theorem that abbreviation is valid on every canonical
arborescent tree, nor that complete topological R± equivalence is
implemented here.

## Exact F1/F2 gauge on abbreviated trees

The abbreviated variants are

```text
{id, F1, F2_f, F1∘F2_f}.
```

These form a four-element action subgroup (on the rooted V4-labelled
presentation). At the marked free root f their labels differ by

```text
{0, ζ, ξ, η},
```

so exactly **one** variant has identity root label. This gives the
explicit F1/F2-only gauge
`abbreviated_flip_gauge` and its rooted
`abbreviated_flip_orbit_key`.

This is not the full canonical classification quotient: F3′
still acts on the actual sector geometry and sometimes on V4
labels, and the underlying tree may have source exceptional
positive/negative canonicality cases. A4.1 formal ring relocation,
A4.2 stable ring normal forms, and A4.3 rooted cyclic alignment
remain separate, more restricted operations.

## Verification

Module: `arborescent/abbreviated_flips.py`.
Tests: `tests/test_abbreviated_flips.py`.

The independent commuting-square checks **apply F1/F2 to full
unabbreviated source trees before erasing any rings** and compare
the result with flips applied after abbreviation. This is not
merely comparing the same normalized code twice.

Test classes cover both signs, corolla ring multiplicity 0 through
3, bounded-valence ring-bearing paths, nonroot F2 pivots, full
four-element group algebra, identity-root V4 gauge, all surviving
per-base ring counts, nonroot V4 label discrimination, cyclic
vertex storage rotation, vertex renumbering, malformed pivots, and
exceptional inputs that must not produce equivalence certificates.

Reproduce:

```sh
python -m unittest discover -s tests -v
```

## Next source-level obligation

The real topological classification in Theorem 12.21 is modulo
F1/F2/**F3′** on the corresponding abbreviated **canonical** trees.
The next work must establish F3′ on ring-numbered atomic tangle
models and verify interactions with the V4-marked boundary
parametrizations, including the source's adjacent-ring exception.
The current theorem certifies the naturality of the two *already
defined source flip transformations* with our restricted notation
change, not actual knot/tangle isotopy by itself.
