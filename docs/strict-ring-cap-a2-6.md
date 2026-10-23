# A2.6 — source-picture R± cap deletion (strict formal subsystem)

## Crucial correction from the original ring-signature search

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots
and the Classification and Symmetries of Arborescent Knots* (2016),
§12.4, **printed page 155, PDF page 167 (one-based)**, contains the
two **modified ring diagrams**. We checked the **actual page image**,
not just extracted text:

https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The diagram is **not** merely a ±1 vertex with two ±2 tips and a
free-bond stem. It depicts a **three-vertex ±2--±1--±2 cap** attached
*by an additional edge* to a **separate lower base vertex** (drawn
with two remaining horizontal attachment directions). The move
deletes the cap and its attachment, leaving that lower base and its
other branches intact. The matching positive/negative signs reverse.

This distinction matters: the one-fork `±1` plus `±2,±2`
examples found by A2.5's necessary-signature detector are **not yet
the full pictured ring move**, because their third center bond is a
free bond, not an attachment edge to a distinct base vertex.

A2.5's mathematical **unit-weight absence obstruction** remains
correct. Its candidate detector is *deliberately permissive* and does
not establish direct ring eligibility.

## Explicitly restricted source-pattern implementation

The new `arborescent/ring_caps.py` implements a strict subset of
this graph-changing ring diagram, as a **formal rewrite on weighted
plane-tree presentations** with complete V4 free-bond label tracking.

A contraction is permitted only if:

1. A valence-three **cap center** is connected by two internal edges
   to distinct valence-one leaves of matching total weight `2s`,
   where `s in {+1,-1}`.
2. The center has exactly one nonzero angular sector of weight `s`,
   in the wedge **between the two cap-tip bonds**; its other two
   sectors are zero. A reflected cyclic bond order is accepted.
3. Its third internal edge (the *stem*) connects to a **different**
   base vertex, which has at least three incident bonds. A **free**
   stem does not match this diagram.
4. Both angular sectors at the base adjacent to that stem are **zero**.
   This intentionally restrictive guard makes sector-weight fusion
   unambiguous and preserves the base's remaining angular weights.
5. None of the three deleted cap vertices has a free bond. All V4
   labels on the surviving tree are exactly unchanged.

The contraction deletes three cap vertices and the three internal
edges comprising their triangle-free *tree cap*: two tip edges plus
the stem to the lower base. It removes the stem's slot from the
base's cyclic order, merging its two adjacent zero sectors. The
result has **three fewer vertices**; the remaining tree is connected.

The inverse `expand_ring_cap(tree, base, after_bond, sign)` inserts
the pictured ±2--±1--±2 cap into a **zero angular sector** at a base
already having at least two incident bonds. It allocates fresh bond
names and does not alter V4 labels. Expansion followed by
contraction is exact; contraction followed by expansion recovers
the input up to **marked weighted-plane-tree presentation isomorphism**
(bond/vertex renaming, cyclic starting-slot rotation).

This does **not** exhaust the modified R± move under arbitrary
plumbing presentations. In particular it does not implement
nonzero-sector splitting at the base or exceptional canonicality.

## Source relation and mathematical status

The monograph says the modified ring moves are consequences of
ordinary ring moves and arithmetic move (2.1); for **canonical
arborescent pairs** they feature in Theorem 12.12 and for **canonical
arborescent tangles** in §12.7.3 / Theorem 12.21.

We have encoded a source-pictured **presentation rewrite**, and
proved its formal graph/sector/V4 consistency in the tested strict
subclass. **No independent 3-dimensional plumbing map was
implemented**, and the current W/S/P/N screen does not certify all
the source's exceptional canonicality hypotheses. Therefore no
unconditional statement identifying arbitrary real tangle isotopy
classes follows from this implementation.

The basic contraction specimen has five vertices: one root base
with free port, the three-vertex cap, and one other terminal
branch of opposite weight. Under strict R+, the cap `+2,+1,+2`
disappears and a two-vertex rooted weighted tree remains; the
analogous R− specimen reverses signs. The reduced two-vertex tree
is deliberately reported **unresolved** by the current screening
code because its *whole unbranched tree exception* has not been
fully formalized. We do not pretend the screen certifies it.

## Exact tests and validation

The new regression family checks:

- Positive and negative figure-specific ±(2,1,2) cap matches;
- Central angular wedge, **separate base**, and zero-sector guards;
- That previously flagged **one-fork** A2.5 motifs remain **false
  positives** for the complete pictured attachment;
- Exactly three deleted vertices, preservation of nondeleted
  sector data and all V4 free-bond labels;
- Inverse expansion/contraction and rooted weighted-tree isomorphism;
- Cyclic orientation changes under the already-sourced F1/F2 flips;
- Nontrivial free labels and multiple free ports at the lower base;
- Rejection of malformed/unsourced moves.

```sh
python -m unittest discover -s tests -v
```

## Next steps

The truly hard A2 requirement is now **not finding candidates** but
proving soundness/completeness of all diagram-transcribed variants
and implementing **full** tangle canonicality with exceptional
unbranched components, plus F1/F2/F3′/R± orbit analysis on *one
fixed canonical, ring-ready finite grammar*.

Before making any Catalan growth claim about tangle **isotopy**,
independently implement the plumbing construction or compare with a
trusted knot/tangle diagram engine. Then prove that the full
reduction system is closed on the chosen grammar (or record all
exits) and its quotient fibers are controlled.
