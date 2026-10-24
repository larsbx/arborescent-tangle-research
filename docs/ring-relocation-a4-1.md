# A4.1 — Conservative ring *relocation* certificates, without deletion

**Status:** exact verified rewrite on weighted planar-tree **presentations**.
This implements a precisely limited **formal adjacent-position
interchange** motivated by the modified ring diagrams of
Bonahon–Siebenmann, not an independently certified tangle-isotopy
transformation. It must not be treated as completing issue #16.

## Source-critical distinction

Source: Francis Bonahon and Laurence Siebenmann,
*New Geometric Splittings of Classical Knots and the Classification
and Symmetries of Arborescent Knots* (January 26, 2016),
§12.4, **printed page 155 (PDF zero-based page 166)**;
§12.5 **printed pp.157–158 (PDF pp.168–169)**; and
§12.7.3 **printed page 167 (PDF p.178)**.

PDF: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The modified R± diagram moves an attached signed
`±2—±1—±2` ring subtree to another angular location at its
attachment base. **Both sides contain the entire ring subtree.**
By contrast §12.5 introduces a different **abbreviated ring-number
representation** in which eligible disjoint ring subtrees are erased
from the diagram and counted at surviving bases.

The source states the exact classification in the **positively or
negatively canonical** tree category (Theorem 12.21 for tangles with
V4-labelled free bonds). Its assertion that placement is irrelevant
is not a license to apply arbitrary cyclic reordering to any
weighted tree. We therefore implement a **source-motivated strict
presentation-level model** and record all unproved topology obligations.

## Implemented restricted formal rewrite

`arborescent/ring_relocation.py` exports a typed `RingSlide`,
`apply_ring_slide`, `ring_slide_successors`, and a BFS with
replayable shortest `RingTrace` move certificates.

For a ring center and sign `s in {±1}`, its preconditions are:

1. The center belongs to exactly one **strict signed ring cap** as
   already recognized in `ring_abbreviation.strict_ring_caps`:
   central trivalent ±1 with two ±2 valence-one tips, third internal
   stem attached to a separate base of valence at least three;
   weight located in the wedge between the two tips.
2. The base's angular sector weights are **all zero**. This
   deliberately avoids claiming an arithmetic splitting rule for
   nonzero angular sectors. Existing mixed-sign candidate centers or
   ring-cap attached to another ring center are refused.
3. The named `past_bond` is exactly one position adjacent to the
   cap stem in the base's incident cyclic list; neither bond may
   occupy index zero, which is held fixed as a **linearized cyclic
   anchor**. Swapping with another ring stem is refused.
4. Every free-bond V4 label is preserved *as stored presentation
   data*. This is **not** a separate independent proof of preserved
   boundary-parametrized tangle isotopy.

**Effect:** swap those two incident **bond positions** at the same
base, without altering or deleting any weighted-tree vertices,
internal bonds, sector weights or V4 labels. The operation is its
own inverse when the same neighbor is selected again. Its output
must have exactly the same **restricted §12.5 abbreviated rooted
presentation code** as its input, or the operation raises an
assertion. In contrast to the retracted PR #11, this graph has
**exactly the same vertex count before and after**.

The code does **not** assert that one elementary formal swap is
literally one of the two printed R± figures without additional
topological hypotheses. The source describes ring relocation
equivalences *in the canonical category* and sometimes requires
F3′. To promote these formal witnesses to genuine R± moves would
require an explicit plumbing/Conway-diagram comparison with boundary
parametrizations and a source-valid canonicality certificate.

## Proposition A4.1 — finite formal connectivity of the corolla fibers

Take the exact star family from
[A3.2](ring-placement-fibers-a3-2.md): one base with `d>=3`
V4-identity free bonds `f0,...,f(d-1)` (f0 marked) and `r>=0`
identical signed ring caps, with all base sector weights zero.

For `r>=1`, under the **restricted formal adjacent slides** there
is exactly one orbit containing all

```text
binom(d+r-1,r)
```

different cyclic presentation placements counted in A3.2.
Every presentation in the orbit has the same ring-numbered
abbreviated code.

**Proof:** anchoring the base list at `f0` identifies each
presentation with a word containing `r` symbols `R` and
`d-1` symbols `F`. Each legal adjacent interchange
`RF <-> FR` moves one ring stem past one ordinary free bond
and leaves the two types of symbols, the weighted cap data, and
all zero sectors unchanged. Any two such binary words can be
connected by adjacent swaps (bubble-sorting the `R` symbols).
Conversely no allowed swap changes the multiset or the
invariant represented by the word. Therefore the finite
formal move graph is connected on the claimed family. All
`r` caps remain present, and §12.5 abbreviation remembers only
their number attached to the base. □

**Bounded valence version:** for A3.2's path with `k>=2`
zero-sector base vertices, each carrying exactly one identical ring
cap and one V4-identity free port, the two distinguished positions
at every base are related by one strict adjacent swap. The k
independent commuting choices connect the full set of `2^k`
marked representations, while keeping maximum bond valence ≤4.
In fact the **complete strictly allowed indexed move orbit** is
larger than the 2^k binary seed family. At each endpoint base the
anchored list has two nonanchor slots, so the ring has **two**
reachable positions. At each interior base the anchored list has
three nonanchor slots (forward spine, free leg, cap), so the ring has
**three** reachable positions, connected by adjacent ordinary-bond
interchanges. The positions at distinct bases change independently
and the source's zero-sector and non-overlapping-cap restrictions
remain satisfied. Therefore the exact formal orbit size is

```text
4 * 3^(k-2)   for k >= 2.
```

This is larger than 2^k when k>=3; our first hosted test mistakenly
asserted that the 2^k seed subfamily was the *entire* orbit.
That test failed, revealing a genuine additional set of placements.
The corrected test checks inclusion of all 2^k designated seed
presentations, the exact 4*3^(k-2) full orbit size, and unchanged
abbreviation codes throughout. This distinction is essential when
reporting raw presentation fibers.

Again this is **presentation-graph connectivity**, not
independently established isotopy-class connectivity.

## Verification and fail-closed behavior

Tests include both ring signs; actual vertex/edge/sector preservation;
explicit inverses; all matching corolla placements for
`(d,r)=(3,1),(3,2),(4,2),(4,3),(5,2)`; all 2^k path
placements at k=2..5; unchanged restricted abbreviation codes;
replayable shortest paths; no cap-cap swap; rejected root-anchor
crossings, nonadjacent moves, malformed move types, mixed-sign
centers, nonzero base sectors, and honest incomplete BFS reporting.

```sh
python -m unittest discover -s tests -v
```

The BFS works on **fixed indexed presentations** and marks
`complete=False` if its explicit state budget is exhausted.
There is no unbounded asynchronous search or claim of isotopy
completeness. This code provides finite, replayable combinatorial
evidence to support—but not substitute for—the topological
verification requested by issue #16.

## Next acceptance milestone

A source- and diagram-certified R± tool must:

- construct the underlying *3D* plumbing/Conway-tangle diagram
  for each input and output, preserving the chosen labelled-boundary
  regime;
- establish exact R± correctness for the allowed W/S/P/N
  canonical subtree, including exceptional signs, ring
  adjacency, nonzero base sectors and F3′ substitutions;
- check and transport every V4 boundary label through the
  actual topological local modification;
- prove the restricted formal slides coincide with the published
  equivalence steps wherever claimed, with explicit figure-level
  evidence; and
- only then relate their quotient to isotopy-class enumeration.

**No Catalan asymptotic for actual tangle-isotopy classes is proved
by this stage.**
