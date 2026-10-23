# A2.6 CORRECTION: modified ring *move* versus ring-number *abbreviation*

## P0 correction to proposed PR #11

**PR #11 was closed without merging.** It incorrectly equated the
Bonahon–Siebenmann modified ring move (R±) with deletion of a three-vertex
ring cap. This is **not** what the source says. The source **relocates**
the ±2—±1—±2 cap around its attachment in the planar tree. The cap
remains present on both sides of the R± replacement. Accordingly the
graph-shrinking `contract_ring_cap` function in closed PR #11 should
**not** be used or described as R± or as a proof of pair/tangle isotopy.

Primary source visually inspected:

- Bonahon and Siebenmann, *New Geometric Splittings of Classical Knots
  and the Classification and Symmetries of Arborescent Knots*,
  §12.4 **printed page 155 (PDF zero-based page 166)**:
  the paired positive/negative **modified ring diagrams** are
  location-changing replacements, not deletions.
- §12.5 **printed page 157 (PDF zero-based page 168)**: to remove the
  *need to consider* ring moves, the authors separately introduce
  **abbreviated trees**, erasing disjoint ring subtrees while storing
  a **ring number** at each surviving vertex. It is the location of
  the ring subtree in the plane that is immaterial under R±.
- §12.5 **printed page 158 (PDF zero-based page 169)**: the treatment
  requires an exception when erasing one ring subtree destroys another;
  the authors specify a special abbreviated representative for that
  exceptional two-center tree. They primarily abbreviate ring subtrees
  attached at vertices of valence at least three.
- §12.7.3 / Theorem 12.21: tangle canonicality and V4 boundary
  labels still have to be handled in a genuine classification proof.

Source: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

## A2.6 implementation: only the source's abbreviation scheme

`arborescent/ring_abbreviation.py` implements a **restrictive
ring-numbered notation change**, not a topology-preserving R± move:

1. Recognize a signed ±2—±1—±2 ring cap attached by an internal stem
   to a **different** base vertex (base valence at least three), with
   the central ±1 in the wedge between the ±2 leaves.
2. Require zero angular sectors adjacent to the stem at the base to
   avoid any unverified redistribution of nonzero weights. This is
   a restricted subset of the source's abbreviated-tree construction.
3. Identify all compatible caps to abbreviate in **one pass**.
   If one ring center is also the attachment base for another cap,
   reject rather than incorrectly handle the source's exceptional
   overlapping-rings configuration.
4. Erase those *disjoint caps solely for the abbreviation encoding*,
   and store each base's attached-ring count as a nonnegative integer
   alongside the remaining sector-weighted plane-tree data.
5. Preserve all surviving free-bond V4 labels. Provide a marked,
   vertex-renaming-invariant and bond-renaming-invariant
   **abbreviated presentation code** that includes ring numbers.

**No R± rewrite is exposed.** The data type `Abbreviated` prevents a
caller from confusing the resulting diagram with the original weighted
tree; recovering a specific placement of erased rings requires additional
choices and has not been implemented. This is exactly why erasure can be
appropriate as a notation while being false as an original R± move.

## Structural claims and limitations

**Proposition A2.6.1 (correct abstraction, restricted case).**
For a set of pairwise disjoint source-pattern cap subtrees whose bases
are never themselves deleted and whose adjacent base sectors vanish,
the resulting `Abbreviated` data is uniquely determined by the
weighted plane-tree *presentation* and the set of detected caps, up
to marked weighted-tree renaming and cyclic starting-point rotations.
All surviving V4 boundary labels remain unchanged.

**Proof.** Each three-vertex cap has precisely one stem edge to its
distinct base. Its removal, *as an abstract encoding operation*,
deletes exactly the center and the two tip vertices; its stem bond is
removed from the base cyclic list, and the two zero angular sectors
merge to zero. Pairwise disjointness makes removed vertices
unambiguous, and the base-not-deleted condition ensures the resulting
graph remains a connected tree. At each base, the total number of
abbreviated caps defines the ring number. Since zero-sector deletion
commutes between distinct stems, the result is independent of the
order in which the caps are processed. Free bonds cannot be on the
deleted cap vertices under the specified pattern, so every free-bond
label survives. Marked rooted tree-code recursion is invariant under
vertex/edge renaming and cyclic start slots by construction. □

**Important limitations.** This is *not* a canonicality theorem for
every arborescent tangle and not proof that two inputs with the same
abbreviated data are necessarily isotopic in a chosen marked category.
That inference requires all the source's hypotheses, equivalence moves,
ring exceptions, and topology. In particular abbreviation is lossy
because it discards ring placement and should never be applied to an
arbitrary tree without the stated restrictions.

## Regressions and next work

```sh
python -m unittest discover -s tests -v
```

Checks cover both signs, separate attached base, full ±2—±1—±2
wedge, one-fork false positives, multiple disjoint rings at the
same base, explicit rejection of adjacent-ring exceptions, surviving
V4 data, F1-reflected arrangements, and marked ring-number codes.

**A2.7 objective:** introduce a *genuine* R± location-change move
preserving the number of ring vertices (in contrast to deleted-cap
abbreviation), with a vertex/sector/V4-preserving source-backed
local rewrite and positive/negative diagrammatic comparisons. Prove
that its quotient is correctly represented by ring numbers only
within the applicable canonical class; otherwise produce
counterexamples and restrict scope.

A2 issue #4 remains open. No Catalan asymptotics under the full
topological quotient are established here.
