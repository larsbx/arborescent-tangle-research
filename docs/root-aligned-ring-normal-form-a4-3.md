# A4.3 — Rooted cyclic-anchor invariance for formal ring normal forms

**Finding:** The initial A4.2 strict ring normalizer was correct for a
fixed **slot-zero anchored** weighted plane-tree representation. But
`Vertex.slots` is a *cyclic* bond list. Its first stored position is
**not intrinsic**, and an arbitrary cyclic rotation of that list
represents the same planar presentation (provided the corresponding
angular-sector list is rotated by the same amount).

In particular, a rotation could put a ring stem in slot zero, causing
`normalize_ring_positions` to reject a representation that its
unrotated twin could normalize. Applying its result directly as an
unqualified marked-plane-tree **isomorphism key** would therefore
be unjustified. This is an implementation-level normal-form issue,
**not** evidence about genuine Bonahon–Siebenmann isotopy.

## Corrected marked-root rule

Choose an existing marked free root bond `f0` (or another explicit
marked free bond `root`). Traverse the connected weighted tree from
the unique vertex incident to that bond. At every vertex:

1. Select its unique *incoming* bond toward that root.
2. Cyclically rotate **both** incident-bond order and angular-sector
   list to put that incoming bond at index zero.
3. Preserve all vertex/edge identities, weights, and explicit free-bond
   V4 labels.

The result `align_rooted_cyclic_slots(source, root)` is a
**presentation isomorphism** of the original rooted cyclic weighted
plane tree; it is not an F1/F2 flip and not a ring relocation. Its
original and aligned `labelled_rooted_code(root)` are required to
match exactly.

Only **after this non-move alignment** do we invoke the A4.2 strict
formal ring-slide normalizer, which sorts cap stems past ordinary
nonanchor bonds by the proven inversion-decreasing algorithm.

The new `normalize_rooted_ring_positions` returns a
`RootAlignedRingCertificate` explicitly distinguishing its
rooted cyclic-isomorphism stage from the subsequent typed formal
`RingSlide` sequence. Both stages are replayable and audited.
The previously existing `normalize_ring_positions` remains
**deliberately anchored** and preserves its earlier guard behavior.

## Theorem A4.3 — well-definedness under planar-tree presentation isomorphism

For any strict input in A4.2's uniformly signed,
non-overlapping, all-zero-base-sector class, and any two versions
that differ only by:

- cyclic changes to starting slots and equally shifted sector lists;
- renumbering the storage indices of vertices;
- renaming internal edge/bond identifiers while preserving their
  incidence; and
- an identical marked free root bond and V4 boundary labels,

the value of `strict_ring_normal_key(source, sign, root)` is
**identical**.

**Proof.** On each version of the rooted tree, the unique parent
bond at each vertex is intrinsically determined by the marked
root-to-vertex path. Aligning the local slots on these incoming
bonds therefore gives identical cyclic local rotations up to
vertex/edge renaming. The A4.2 sorted output preserves the ordered
ordinary and strict-ring child subsequences, and every ring cap
has the same signed weighted local subtree. Thus the normalized
*rooted code*, which ignores incidental vertex/edge names,
agrees in both versions. □

This theorem has only **presentation-isomorphism plus A4.1 formal
slide** scope. Equality of these keys does not prove actual tangle
isotopy or complete R± canonicality.

## Regressions

`tests/test_ring_root_alignment.py` checks:

- both ring signs and every small corolla location, under two
  independent per-vertex cyclic rotations and complete vertex/edge
  storage renumbering;
- previously rejected slot-zero cap representations, now accepted
  through *rooted alignment* but still rejected by the lower-level
  slot-zero-anchored normalizer;
- full bounded-valence A4.1 formal slide orbits with arbitrary
  cyclic storage starts;
- exact replay and idempotence of the alignment and slide stages;
- preservation and discrimination of nonroot V4 boundary labels;
- invalid or missing roots; and
- retention of nonzero sector weights under cyclic rotation.

Run `python -m unittest discover -s tests -v`.

## Important remaining limits

The root-alignment certificate repairs a **combinatorial data model**
issue. It does not implement a 3D plumbing/Conway tangle diagram,
source-complete W/S/P/N canonicity, exceptional R± cases,
F1/F2/F3′ topological equivalence or a sound isotopy classifier.
[Issue #16](https://github.com/larsbx/arborescent-tangle-research/issues/16)
remains open for those obligations.

Bonahon–Siebenmann primary source: §12.4 printed p.155,
§12.5 printed pp.157–158, and §12.7.3 printed p.167,
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf
