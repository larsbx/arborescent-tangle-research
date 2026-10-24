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

## Independent source-diagram and adversarial audit

The PDF linked above was independently downloaded and rendered for this
review (SHA-256 `c887732f91421ce8bf424f32ab4c3b52b842a5c15f8bc95c62599b8d8baa88a6`).
The inspected diagrams are printed pp.155,157–158,167 (PDF indices
166,168–169,178). The following are observations from those pages:

- **p.155:** each R+ and R− replacement retains the two terminal
  vertices, center and attachment stem. The cap changes angular
  location at its base. Its unit weight is in the wedge between its
  two terminal bonds. No three-vertex deletion occurs.
- **p.157:** ring numbers count erased attached subtrees at each
  surviving base, rather than merely indicating presence. The example
  includes multiplicity two. The source also admits the other two
  central weight sectors, using F3′; our strict recognizer deliberately
  handles only the directly pictured outer wedge.
- **pp.157–158:** the two adjacent positive centers are an exception,
  up to flips; the source assigns a special abbreviated representative.
  We reject this configuration rather than implement that representative.
  The ordinary construction uses bases of original valence at least three.
- **p.167, §12.7.3:** the classification concerns canonical trees with
  V4-labelled free bonds and includes flips as well as ring moves.
  Abbreviated trees retain the labelled tangle setting. Merely comparing
  our rooted presentation codes implements neither that flip quotient
  nor the classification theorem.

The negative implementation is the sign-reversed structural analogue.
The orientation/sign correspondence is stated on printed p.154; §12.5
itself develops the positive notation. Neither sign branch of this API
checks all canonical-tree hypotheses. Arbitrary structurally valid
inputs therefore receive only an encoding, never a topological certificate.
Likewise zero adjacent base sectors are our conservative restriction,
not a claimed general source requirement.

The independent tests in `test_ring_abbreviation_adversarial.py` show:

1. Four insertion positions include distinct original presentations but
   the same ring-numbered code. This is intended loss of placement;
   the source inputs are immutable, and no R± API is thereby implemented.
2. Bare deletion (ring count zero) differs from abbreviation, and counts
   zero through three remain distinct even when the surviving tree agrees.
3. Counts at two different bases remain distinct; untouched nonzero
   sectors and all four free-port V4 labels survive. Changing any of the
   three labels in the angular-placement fixture changes its rooted code.
4. Nonzero sectors beside a stem prevent abbreviation and are preserved.
5. Opposite-sign adjacent centers exposed a guard gap in the original
   head: selecting either sign erased one center and destroyed the other
   cap's attachment. Such a mixed pair is outside the source's uniform
   canonical-sign setting. The API now conservatively rejects a selected
   cap whose base is *any* necessary ring center, regardless of sign or
   sector eligibility. This extends refusal; it does not assert a new
   mixed-sign source equivalence.

The structural invariant remains: each encoded cap deletes three
vertices and three internal bonds only in the new representation,
adds one to exactly its surviving base's ring number, leaves a connected
tree, and preserves every surviving sector total and free-bond label.
Unsupported patterns may be left untouched; abbreviation need not be
exhaustive, idempotent, or canonical under flips. Boundary bond names
are presentation identifiers (the existing model permits their renaming);
the chosen root, cyclic order and associated V4 data enter the code.

Local acceptance after the correction: `python -m unittest discover -s tests -v` passes **72 tests**, and the workflow's prototype
smoke command passes. Hosted Python 3.11 and 3.12 must pass on the updated
head before merge; the earlier green run applies only to `f975ce8`.

**Still open:** genuine R± relocation and its marked boundary bookkeeping;
exceptional canonicality and the special adjacent-ring representative;
independent plumbing soundness; full topological isotopy classification.
