# A4.9 — Exact finite class graph modulo rooted plane presentation isomorphism

**Scope:** only the **partial formal** F1/F2, strict surviving-edge
F3′, and strictly guarded *formal ring-slide* system of A4.7.
This is not a knot/tangle isotopy classification or the complete
Bonahon–Siebenmann equivalence calculus.

## The correction to A4.8's indexed search limitation

[A4.8](mixed-necessary-invariants-a4-8.md) established that a
complete *indexed* mixed-move search can miss a target that is
literally the **same marked weighted plane-tree presentation** as
its start, but whose vertex-storage array and internal edge names
are different. Complete fixed-index BFS does not automatically
quotient presentation isomorphism.

The A1/A2.2 data model already has a faithful marked rooted
`V4LabeledTree.labelled_rooted_code(root)`, insensitive to
vertex storage order, internal edge renaming and cyclic
starting-slot choices, but sensitive to signed vertex sector
weights, cyclic **order**, root marking, and **all** V4 free-port
labels.

A4.9 now explores **that mathematical code**, rather than
arbitrary storage object identity.

## Exact quotient construction

Fix a signed, completely strictly recognized input presentation
`T` and a marked root free bond `f`. Write

```text
code(T) = T.labelled_rooted_code(f).
```

A4.7's `mixed_successors` enumerates the expressly declared
partial moves: source-defined F1, root-pivot F2, strict surviving-
edge F3′ with a **free-bond pivot**, and a signed **formal**
adjacent ring-stem interchange. Each move is already guarded
by its independently implemented §12.5-abbreviation commuting
square or fixed ring-number check and by A4.3 root cyclic
alignment.

The class BFS has one stored representative for each code.
For each representative, compute its full finite successor list;
retain a new successor only when its **code** has not previously
appeared. Every stored predecessor records the **actual typed
move** from its stored parent representative.

A positive target certificate consists of:

1. An exact, replayable A4.7 move trace from the aligned start
   to a stored representative U.
2. The exact equality
   `code(U)=code(target)`, certifying a final **rooted
   weighted-plane-tree presentation isomorphism**, not a move.

The code explicitly separates the isomorphism step from
the formal move trace.

## Proposition A4.9.1 — congruence of the declared partial generators

Let `T≅T′` denote a rooted V4-labelled weighted plane-tree
presentation isomorphism preserving the marked free bond, all
integer angular weights under matched cyclic rotations, and all
V4 free-bond labels. Then the sets of successor **rooted codes**
under the A4.7 declared generators agree:

```text
{code(U) : T -> U} = {code(U′) : T′ -> U′}.
```

**Proof.** The isomorphism transports every vertex, internal bond
and nonroot free bond to its counterpart and preserves integer
weights, cyclic rotations and V4 data.

- The global F1 flip applies the same reflection and ζ label
  changes on both sides.
- Root-pivot F2 is characterized intrinsically by tree-distance
  parity from the marked root free bond, which is preserved by
  graph isomorphism, as are ξ/η actions.
- Strict F3′ candidates are enumerated at **each** incident
  surviving internal edge of **each** surviving named free-bond
  pivot, and for both sector sides. Such candidates and their
  parity-sensitive branch reflections/label updates transport
  through the isomorphism. The separate A4.5 commuting-square
  and ring-recognition guards are equivariant since all their
  predicates depend on transported incidences and weights.
- Strict A4.1 formal ring-slide candidates use a recognized
  signed pendant cap and one neighboring **ordinary** incident
  bond. A4.3 first aligns every cyclic list on its intrinsically
  determined parent bond toward the root. Thus the nonanchor
  adjacency, zero-sector and ring-exception guards, and the
  output of a legal formal swap all transport under isomorphism.

Consequently each legal move from T corresponds to one move
with the same resulting rooted code from T′ and conversely. □

The theorem is **about the implemented rules** and does not
supply their independent 3D topological soundness.

## Theorem A4.9.2 — correctness and qualified completeness

For a fixed state limit and signed source T, the quotient BFS
discovers a finite set of rooted presentation-isomorphism
classes reachable under the declared A4.7 partial formal moves.
Every reported positive certificate replays a finite move path
and verifies the final rooted code against the requested target.

If `complete=True`, the BFS has exhausted **all** reachable
rooted presentation-isomorphism classes under those *declared
formal moves*: expanding one representative per class is enough
by Proposition A4.9.1. If `complete=False`, no negative claim
follows.

**Proof.** Every discovered class is reached by a typed move
from a previously reached stored representative, so inductively
every class is reachable. If no budget cutoff occurs, every
representative is expanded. By congruence each legal successor
of *any* representative of a discovered class has the same
class code as a successor from its stored representative, hence
also appears in the discovered set. By induction on path length,
all finite reachable classes are discovered. The replay and
code-equality checks establish a positive witness. □

**Consequence:** Complete absence of a target **code** now
establishes nonreachability up to **marked plane-tree
presentation isomorphism** for the coded *partial formal graph*.
It does **not** imply nonisotopy or nonreachability under any
unimplemented R±, F3′, exceptional, plumbing or boundary moves.

## Four evidence statuses

`decide_formal_class_reachability` returns:

| Status | Restricted mathematical assertion |
|---|---|
| `reachable_up_to_presentation_isomorphism` | Replayable A4.7 formal path to U, plus `code(U)=code(target)`. |
| `excluded_by_necessary_invariant` | A4.8 rooted unordered weighted-ring invariant differs; formal reachability impossible without any BFS. |
| `not_reached_in_complete_formal_class_graph` | The exact finite graph of all **rooted presentation-isomorphism** classes reachable under these coded moves is exhausted; target code absent. No topology claim. |
| `undetermined_cutoff` | Class budget exhausted; no negative assertion. |

The result includes exact fingerprints, explored **class** count,
completion status and the positive two-part certificate where
available.

## Independent regression witnesses

`tests/test_presentation_class_bfs.py` checks:

- A differently indexed and cyclically rotated copy of the
  **same** original presentation, which A4.8's complete
  fixed-index graph correctly does **not** reach by typed
  moves, is recognized by A4.9 **without any formal move**:
  an empty path plus a marked rooted-isomorphism certificate.
- A target reached by one formal ring slide is still positively
  recognized after completely renaming its internal bonds and
  reversing its vertex storage.
- Exact code sets agree with independently enumerated full
  indexed A4.7 orbits on small signed corollas.
- Sets of successor **codes** agree when computed from
  differently numbered/cyclically rotated presentations,
  including strict nonzero-sector F3′ witnesses.
- Two zero-ring trees with the same unordered weighted-ring
  fingerprint but **different cyclic orders of ordinary
  weighted leaves** are not connected by the coded partial
  moves. This produces a genuine complete **formal class**
  nonreachability result without asserting topological
  inequivalence.
- Changing an individually marked nonroot V4 free-port label
  is detected by the *full* rooted class code even though
  A4.8's coarse necessary invariant ignores V4.
- One-class cutoffs remain undetermined; unequal necessary
  fingerprints permit an immediate formal exclusion; invalid
  roots/signs/budgets are rejected.

Reproduce with:

```sh
python -m unittest discover -s tests -v
```

## Remaining obligations

The exact mathematical status is now stronger than the A4.8
fixed-index search, but **still not topological**:

- Source §§12.4,12.5,12.7 do not establish that our arbitrary
  A4.1 **formal** adjacent-stem swaps are genuine R± moves
  for noncanonical weighted tree inputs.
- The implemented strict surviving-edge F3′ is a proper
  subset of the source's complete canonical equivalence rules.
- Exceptional adjacent-ring canonical trees, the tangle-specific
  W/S/P/N hypotheses, full V4 boundary transport and an
  independently verified three-dimensional Conway/plumbing
  realization are still outstanding.
- Thus no test based on these classes can decide ambient isotopy
  or prove Catalan-type growth of genuine isotopy classes.

The source-validation target remains
[Issue #16](https://github.com/larsbx/arborescent-tangle-research/issues/16).
