# A2.3: finite, replayable partial-move closure and the F3′ vacuity obstruction

## Trigger and verification debt

On 23 October 2026, the *already merged* A2.2 PR #6 hosted Python
3.12 job failed. The failure was an erroneous test fixture with
`Fork((Leaf(2), Fork((Leaf(2), Leaf(3)), Leaf(3))))`, passing
two positional arguments to a `Fork` constructor that accepts exactly
one tuple of children. This failure was discovered in the completed
[PR #6 Actions workflow](https://github.com/larsbx/arborescent-tangle-research/actions/runs/39986861666).
This patch corrects that fixture and adds independent graph-level
regression tests. **A2.2 was not validated when it was merged.**

## Research finding: the A1 grammar cannot use F3′

**Proposition A2.3.1 (strict-F3′ vacuity on A1).** For every number
`n≥1` of internal forks and every finite alphabet of leaf weights,
the A1 grammar cannot admit the implemented strictly W-preserving F3′
weight-transport move, even after applying any sequence of F1 and F2.

**Proof.** The A1 grammar has zero sector weight at every fork of bond
valence 3 or 4, and a single potentially nonzero sector at each
valence-one leaf. The strict F3′ operation transfers a single nonzero
sector weight between two **distinct** sectors meeting an internal
edge at its weighted pivot. A zero-weight fork fails this precondition.
At a valence-one leaf, the sectors immediately before and after its
sole incident bond coincide, so there are no two distinct sectors.
F1 and F2 merely reverse cyclic orders and rearrange sector positions;
they never transfer a nonzero weight between vertices. Therefore no
F3′ edge can arise under the elementary flip subgroup. QED.

**Consequence:** the A2.2 F3′ implementation is not a genuine new
identification on the A1 enumerated grammar. It can be exercised only
on trees containing nonzero-weight internal vertices, outside A1.
It would be incorrect to use A2.2 alone as evidence that the complete
topological quotient preserves Catalan-type growth.

This is an exact structural obstruction, not a missing-data claim.
An exhaustive check is included for n≤3 and labels `{2,3}`.

## Source-scoped partial rewrite engine

The file `arborescent/partial_moves.py` introduces:

- Typed move data for **F1**, **F2**, and the strict internal-edge
  **F3′** case previously implemented;
- Explicit preconditions, including the F3′ nonzero-sector guard;
- A deterministic finite graph search over *indexed* weighted-tree
  presentations, with redundant F2 pivot choices compressed to the
  two distance-parity classes;
- A caller-selected `max_states` budget, reporting
  `complete=False` on exhaustion rather than misrepresenting a
  truncated search as a full orbit;
- Shortest-path witness traces replayed by the actual move functions.

The underlying graph, number of vertices, sector-weight values and
integer bond names remain fixed. Reversal of rotations, redistribution
of one nonzero sector at its original vertex, and finite V4 changes
make the state set finite. This graph search is exhaustive **for the
implemented source-scoped moves on the fixed indexed graph** when
`complete=True`. It does not quotient arbitrary graph isomorphisms
and is not a tangle-isotopy classifier or a soundness proof of plumbing
moves independent of the published theorem.

## Genuine nontrivial F3′ witness

The regression fixture uses a four-vertex weighted plane tree with
two free ports, a pivot weight `a=3` at a valence-three vertex,
and a branch containing its second free port. Moving the weight to
the other side of edge `e0` changes the subtree's odd-distance
rotations and its V4 free-bond label. The explorer checks a shortest
witness path, complete fixed-graph orbit, all predecessor traces,
and explicit budget truncation.

See Bonahon–Siebenmann, *New Geometric Splittings of Classical
Knots and the Classification and Symmetries of Arborescent Knots*,
§12.7, printed pp. 166–167, and
[the A2.2 source audit](v4-tangle-moves-a2-2.md).
Source: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

## Acceptance and next proof milestones

Run: `python -m unittest discover -s tests -v`.

The next **mathematically consequential** step is to specify an
explicit bounded-local grammar of *weighted internal forks* whose
weighted trees satisfy the appropriate **tangle** W/S/P/N canonicity
conditions. Then test F3′ on that grammar rather than pretending
A1 has such moves. Independently implement and validate modified
ring moves R± with the source's pictures and exceptions. The
Bonahon–Siebenmann classification theorem applies to canonical
objects, **not** to arbitrary output from this exploratory engine.
