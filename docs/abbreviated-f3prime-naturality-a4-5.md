# A4.5 — Restricted F3′-equivariance of V4 ring-numbered abbreviation

**Status:** a checked **formal commuting square** for a nontrivial
class of source-defined V4-labelled F3′ transformations. This does
**not** independently prove topological tangle isotopy or the
full Bonahon–Siebenmann canonicality theorem.

## Source and precise claim

Bonahon–Siebenmann §12.7.3, printed p.167 (PDF page 178
one-based, index 178 zero-based), defines the **modified F3′ flip**
with a single transported nonzero sector, parity-sensitive branch
reflection and associated ξ/η transformations of V4 labels at free
bonds. The paragraph after Theorem 12.21 permits passing to
§12.5's abbreviated ring-numbered presentations, after eliminating
R± through the source's proper canonical-tree equivalences.

The source's separate §12.5 construction erases suitable disjoint
signed ring-cap subtrees **as a representation change**, not as a
graph-shrinking R± move.

Primary monograph:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

[A4.4](abbreviated-flip-equivariance-a4-4.md) established exact
commuting squares for the source F1/F2 V4-labelled flips on
strictly recognized ring-numbered presentations. The present step
addresses the third **modified flip F3′**, but only in a
particularly safe, nontrivial configuration.

## Restricted naturality theorem

Fix a V4-labelled weighted plane-tree presentation T with a
uniform ring sign s=±1 and a **complete strictly recognized**
family of disjoint source-pattern ±2—±1—±2 cap subtrees satisfying
the A2.6 zero-neighbor-sector and nonexceptional attachment guards.

Select a **surviving** free bond f and let v be its vertex on the
full tree. Select an interior edge e that **also survives** in the
abbreviated residual tree, incident to v. Suppose the
already-implemented strict F3′ move is defined at (v,e) on
both T and the residual tree, and the resulting full presentation
has exactly the same family of strictly recognized cap centers and
the same ring counts at their surviving bases.

Write A_s for the restricted ring-numbered abbreviation, and
F3′* for the transformation which applies the identical V4-labelled
strict F3′ rule to the **residual** weighted tree and retains its
per-base ring numbers unchanged.

**Proposition A4.5.1 (restricted formal naturality).**
Under these explicit hypotheses,

```text
  A_s(F3′_(v,e)(T)) = F3′*_(v,e)(A_s(T)).
```

**Proof.** Since e survives abbreviation, the selected
F3′ branch is described by the same bipartition of the
surviving vertex tree before and after deletion of pendant caps.
The parity distance from v of every surviving vertex is unchanged.
Thus parity-dependent cyclic reversals at surviving vertices
and ξ/η V4 changes to their incident free bonds agree on both
routes. If the ring centers and their ring attachment bases remain
strictly recognized as assumed, each erased cap still contributes
exactly one ring count to its original surviving base.
Abbreviation at the base merges two *zero* sectors adjacent to its
cap stem, and this merger commutes with a cyclic reflection of the
base whenever F3′ reverses it. The transported pivot weight is
across the surviving edge e; cap erasure does not touch that
weight or edge. All residual weighted sectors, cyclic bond orders,
free-bond V4 labels and vertex ring numbers therefore agree. □

**Important:** This theorem is asserted *only for cases satisfying
the strict guards* and for our already implemented subset of F3′.
It does not cover moving a weight across an erased ring stem,
the full published F3′ operator, or the source's exceptional
adjacent-ring trees. It does not identify arbitrary isotopic
tangles.

## Fully explicit nonvacuous witness

The test family has a weighted root vertex with V4 free port f0,
root sector weight a=±3 (odd) or a=±4 (even), and a surviving
interior edge e0 to a separate zero-sector base. This base
has another V4-labelled free port f1, terminal branches
weighted 3s and 2s for the chosen ring sign s, and a
signed ±2—±1—±2 ring cap attached through a different edge e2. The cap is fully recognized by A2.6.

The strict F3′ move at the weighted root transports a across
the **surviving edge e0**. For odd a, it reverses the first
surviving base's cyclic order (distance one) and applies the
source's η transformation to f1. For even a the base and f1
are unaffected. In both cases, f0 remains outside the branch.

Both the source witness and its transported image pass the
implemented nonexceptional branched **W/S/P** (s=+1) or **W/S/N**
(s=-1) structural screen. This is a necessary source-canonicity
check, not independent topological validation or complete treatment
of all exceptional trees.

Independent computation along each route verifies equality of
the full `Abbreviated` objects, not merely equality of
their coarser rooted codes. The inverse F3′ direction is also
checked, along with both ring signs and positive/negative a.

## Implementation and fail-closed tests

`arborescent/abbreviated_f3prime.py` provides:

- `abbreviated_f3prime(abbreviated,pivot_free_bond,branch_edge)`;
- `f3prime_abbreviation_square(full,sign,pivot_free_bond,branch_edge)`.

Both resolve the F3′ pivot by a **surviving explicitly marked free
bond**, never a vertex index that might change when pendant ring
caps are erased. The commuting-square checker rejects incomplete
necessary source-pattern recognition, unsupported erased branch
edges, failed strict F3′ preconditions and changed cap-center or
ring-count data rather than returning false topology certificates.

Regressions `tests/test_abbreviated_f3prime.py` include:
both source ring signs; even/odd positive and negative pivot
weights; forward/inverse commuting squares; nonroot V4 η
behaviour; agreement with the previously established F1/F2
commuting squares; missing pivots and erased ring-stem edges;
zero-weight/two-nonzero-sector pivots; and a cap base with a
nonzero adjacent sector, which must fail closed.

Run: `python -m unittest discover -s tests -v`.

## Next proof obligations

The source classification's **topological** content still requires
validating the actual 3D plumbing construction, positive/negative
canonical W/S/P/N hypotheses and their exceptional cases, the full
tangle F3′ calculus, remaining R± location changes, and all V4
boundary-parametrization effects. Combining the A4.4/A4.5
formal equivariance checks does not discharge those requirements.

A4.5 shows a genuinely nontrivial, correctly labelled
F3′/abbreviation interaction **within a predeclared formal
subclass**, not an unconditional enumeration theorem about
arborescent-tangle isotopy classes. Issue #16 remains open.
