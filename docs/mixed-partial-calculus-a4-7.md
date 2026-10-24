# A4.7 — Mixed strict F3′ / flip / formal ring-slide reachability

**Status:** bounded, replayable graph search for an explicitly restricted
**formal presentation** calculus. No negative reachability result from
a truncated search, no universal confluence claim, and no theorem about
3D arborescent-tangle isotopy.

## What this stage adds

Earlier stages proved three different statements at three different
scopes:

- A4.4 established that source-defined **V4-labelled F1/F2** flips
  commute with §12.5's disjoint, strict ring-numbered
  representation change.
- A4.5 established a similar **conditional** commuting square for
  the already implemented *strict subset of F3′* when its pivot,
  branch edge, recognized cap centers and per-base ring numbers
  survive abbreviation.
- A4.6 established a complete normal form **only** for the joint
  F1/F2 and A4.1 *formal* ring-stem slides, not when F3′ is also
  admitted.

Therefore it would be mathematically premature to reuse the A4.6
normal-form key as a complete invariant for the joint F1/F2/F3′
quotient. F3′ transports nonzero angular-sector weights, so it can
genuinely change a stable formal ring normal-form key.

This stage instead gives an **honest bounded exact-state automaton**
on a declared partial move graph and verifies each individual
move's compatibility with §12.5 **before accepting it**.

## Declared moves

Given a marked free root bond `f0`, a uniform ring sign s=±1 and a
V4-labelled weighted plane tree:

1. **F1** — the existing source-defined global mirror, with ζ
   applied to every free V4 label; requires the independently
   implemented F1/F2-abbreviation square to pass.
2. **F2_f0** — the existing even-distance mirror at the vertex
   incident to the marked free root, with source ξ/η V4 changes;
   also requires F1/F2-abbreviation squares.
3. **Strict F3′_(f,e,side)** — an already-implemented nonzero-sector
   transport at a vertex incident to a **surviving named free bond
   f**, across a **surviving internal edge e**, with the source
   parity-dependent branch reflection/V4 actions. It is admitted
   only if the A4.5 before/after **commuting square returns true**
   and all strict cap centers and ring counts remain recognized.
4. **R_formal** — an A4.1 strictly guarded adjacent interchange
   of an entire ±2—±1—±2 cap attachment with a neighboring
   ordinary nonanchor bond on an all-zero-sector base. This is
   explicitly a **source-motivated formal presentation rewrite**,
   not an independently verified topology-preserving R± move.

Every input and successor is aligned as a **cyclic plane-tree
presentation isomorphism** toward the marked root (A4.3). The
converter refuses unmatched ring signatures, mixed signs,
adjacent-center exceptional cases, ring stems erased by
abbreviation, and unsupported pivot configurations.

The separate **§12.5 abbreviation** never counts as an action
in this move graph: it is a different representation.

## Proposition A4.7.1 — positive witness correctness

Suppose `MixedOrbit.trace(T')` returns a typed move sequence
`(m1,...,mk)` with `MixedTrace.verify()=True`. Then this
sequence is a valid transformation certificate from the root-aligned
presentation T to T′ **within the exact implemented partial
presentation calculus**.

**Proof.** The predecessor map records only successors computed by
`apply_mixed`. Every application checks the signed strict domain,
applies exactly one existing F1/F2/F3′/formal-R operation, checks
the appropriate before/after §12.5 square or abbreviation
invariant, ensures per-base ring counts remain equal, and
realigns only cyclic storage. Replaying every step without
consulting the predecessor map yields the recorded target.
Therefore the trace is an independently executable witness in
the declared formal graph. □

A failed *complete* finite search means the target is unreachable
**under these coded partial rules** on the explored indexed graph.
It never proves distinct tangle isotopy classes. When a state
budget is exhausted, `complete=False` is reported and no such
negative conclusion is justified.

## Proposition A4.7.2 — explicit F3′/formal ring-slide diamond

Construct the seven-vertex family with:

- a marked root vertex `(f0,e0,e1)` with exactly one nonzero
  angular-sector weight `a∈{±3,±4}`;
- the surviving edge `e0` to a zero-sector ring attachment
  base `(e0,e2,e3,f1)`, with explicitly V4-labelled free
  ports `f0` and `f1`;
- a pendant signed ring cap `e2`: trivalent center sectors
  `(0,s,0)` and two tips each weighted `2s`;
- ordinary weighted terminal branches attached through `e1`
  and `e3`, with weights `2s` and `3s`.

Let `F` be strict F3′ at the root across the surviving edge e0
from its nonzero `before` sector. Let `R` be the **formal**
slide of ring center e2 across ordinary branch e3 at the
attachment base.

For either sign s and all a values above, the two routes
commute as **indexed V4-labelled weighted presentation data**:

```text
                   R_formal
         T -----------------------> R(T)
         |                            |
       F3′                           F3′
         |                            |
         v                            v
       F(T) ----------------------> F(R(T))
                 R_formal
```

**Proof.** F3′ transports the root's sector weight without
changing its base attachment edge e0, and its parity-dependent
reflection affects the ring-bearing base precisely when a is
odd. A formal swap interchanges its neighboring e2 and e3
bond positions; reflection reverses the same local cyclic order,
so the corresponding adjacent interchange remains admissible.
All sectors at the base are zero, and the ring center/tip sectors
remain recognized. Every free V4 label receives the same
parity-dependent XOR along either route. Thus the two resulting
indexed presentations are equal. The tests compare full
`V4LabeledTree` objects for both signs and each a. □

This **specific diamond** is not a proof that every F3′/R±
critical pair resolves: unsupported cases may be noncanonical
or require additional source moves.

## Independently checked negative controls

The test suite checks that the formal graph with a zero-weight
root and recognized rings has **no applicable surviving-edge strict
F3′**: on the corolla with d=3 ordinary free legs and one
ring cap, the complete partial graph equals the previously
computed **12-state** F1/F2 + strict-ring-slide graph.
The new F3′ witness genuinely enlarges the available transitions
when its root carries nonzero angular weight.

Other negative controls require explicit rejection of: an
erased ring stem selected as F3′'s branch; nonexistent pivot
port; unknown or mixed signed cap signatures; nonzero cap-base
sectors; incorrect move constructors; and invalid or exhausted
search budgets. `MixedTrace` reports only replayable positive
witnesses, while a truncated search returns `complete=False`.

## Implementation and validation

- `arborescent/partial_joint_f3prime.py`: typed `MixedMove`,
  audited `apply_mixed`, exact successor enumeration,
  `explore_mixed`, `MixedOrbit` and `MixedTrace`.
- `tests/test_partial_joint_f3prime.py`: signed odd/even
  F3′ witnesses, the explicit F3′/ring-slide diamond,
  zero-F3′ corolla control, bounded reachability, V4
  and ring-number invariants, and fail-closed guards.

Run:

```sh
python -m unittest discover -s tests -v
```

## Explicit topological limitations

Primary source: Bonahon–Siebenmann *New Geometric Splittings of
Classical Knots and the Classification and Symmetries of
Arborescent Knots*, §12.4 printed p.155; §12.5 printed
pp.157–158; §12.7 printed pp.166–167:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

These tests verify **mathematical properties of coded
presentations**, not correctness of a complete source
classification. A4 issue #16 still requires an independently
constructed Conway/plumbing tangle map, actual R± topological
isotopy proof with V4 boundary parametrizations, full F3′,
exceptional adjacent-ring canonical forms, W/S/P/N source
hypotheses and topology-level enumeration.

**No isotopy-class growth or Catalan theorem follows from A4.7.**
