# A5.2e2b4c2d1 — Complete vertex/edge-aware PL slit of one connected 3-ball

**A genuine positive PL 3-manifold cut theorem, with limited ambient scope.**

This work continues Issue #58 using the exact connected roof 3-ball from A5.2e2b4c2b (#56), already checked for coherent oriented tetrahedra, real tetrahedral face adjacency and PL vertex links. It upgrades the earlier **dual-face barrier** bookkeeping to an actual **complete split of the vertex and edge stars along an entire separating disk**.

**This is not a Bonahon–Siebenmann source slit:** the product 3-ball covers only the auxiliary roof band region and its source-candidate K is entirely below it. The resulting two true 3-ball boundaries are **unpunctured by K**, not the paper's new four-punctured Conway spheres. No source framed Conway gluing, cyclic band isotopy or ring R± follows.

## Why the earlier dual-face barrier was not a slit

PR #56 constructed one exact 3D rectangular slab

```
  M_roof = [x₀,x_last] × [−2,0] × [4−h,4+h], 0<h≤1/4
```

with n source-candidate chart intervals, n−1 connecting intervals, m=2n−1 longitudinal x-slabs, 12m tetrahedra, and one closed Euler-2 boundary sphere.

For each of the n isolated roof patch disks, removing the incident tetrahedra's dual adjacency still left an alternative route around the marked patch edges through uncut connector material. This did **not** create a globally separating, properly triangulated new boundary: midplane edges/vertices remained shared by neighboring tetrahedra on opposite sides.

The new operation deliberately cuts along **the entire z=4 separating rectangle**, including all n roof patches **and the connecting intervals between them**. This full disk is NOT the Figure 12.10 source plumbing-patch set, but it is a mathematically valid, exactly verifiable separating surface.

## Exact two-sided 3D construction

Classify each original positive-volume tetrahedron by the actual vertex z-heights:

- LOWER: all four z≤4, at least one z<4.
- UPPER: all four z≥4, at least one z>4.
- Reject any tetrahedron crossing the open cut plane.

Copy each tetrahedron to its half with the original **physical Q³ coordinates**. Every physical vertex on the midplane z=4 is replaced by a **distinct typed lower vertex and typed upper vertex**, with no common mesh index across the two cut components. Vertices away from the cut carry an exterior tag; they are physically separated by z anyway. Thus all cut **edges, triangles and vertex links** are split consistently, not just tetrahedral dual-face adjacency.

The independent checker builds each half's own compact vertex/face/tetrahedron index tables from its **actual cut cells**, runs the original coordinate-based tetrahedral manifold verifier on each half (including every vertex link as an interval-boundary disk or interior sphere), and reconstructs each half's complete oriented exterior boundary from tetrahedral face incidences. Only those actually exposed triangular faces whose three vertices have z=4 are accepted as candidate cut faces.

Each actual cut face must be present in **both** halves at identical Q³ coordinates but different typed vertices. Each lower face's outward normal has z>0; each upper face's outward normal has z<0. The two disk boundaries are independently checked as connected PL 2-disks, χ=1, with exactly one simple boundary circuit. Reversing any face, sharing a cut index, omitting a tetrahedron, or forging a cut table is rejected.

## Independently verified exact counts

Let m=2n−1. The original **one** 3-ball has 12m tetrahedra. The complete separating cut yields

```
  lower 3-ball tetrahedra:      6m
  upper 3-ball tetrahedra:      6m
  physical z=4 cut vertices:   2(m+1)   (DUPLICATED on both sides)
  cut triangles per side:      2m
  cut disk V,E,F:              2(m+1), 4m+1, 2m
  cut disk Euler:              1
  boundary Euler, each half:   2
  connected components:       2 (two actual compact 3-balls)
  exact volume, each half:     original volume / 2
```

The independent inverse reglue checks all two-sided cut triangle positions and opposite normals before mapping the two typed midplane copies back to their original unique physical vertex IDs. It verifies that the recovered **complete tetrahedron list** agrees with the original 3D cell complex, with no extra, missing, overlapping or altered simplex.

This is a stronger positive statement than PR #53's single *local* 3-ball slit and PR #56's dual-face barrier: we now have a complete, replayable split of the **connected multi-chart slab**, including the edge and vertex identifications.

## Source-critical negative result

The slab's physical z-range satisfies z≥4−h≥15/4, but the properly embedded source-candidate knot strands have physical z≤1. Therefore **neither side's cut disk has any knot piercing**. The new two 3-ball boundaries are not the source's four-K-puncture Conway boundaries.

Moreover, the cut had to include **extra connector strips** between the candidate plumbing charts to become a separating global disk. The paper's specified source plumbing patches cannot be replaced by that artificial full midplane disk without an independently verified source-equivalence proof. The positive result here is the validity of the **chosen slab's PL cutting operation only**.

The public result is `FULL_ROOF_PL_SLIT_NOT_SOURCE_CONWAY`, and calling `require_source_conway_certificate()` fails with `UNVERIFIED_GEOMETRY`. The earlier flags on the band candidate, slab, and A5.1 plan remain unchanged. The new 3D cut is not an R± isotopy witness.

## Tests and next proof target

`tests/test_full_roof_midplane_slit.py` verifies n=2,3,4 and both atom signs, exact h=1/16,1/8,1/4, source-candidate zero/mixed sectors, the full independent per-half manifold checker, cut disk V/E/F/χ, physical opposite cut-normal orientations, exact volume conservation and inverse cell recovery. Negative tests forge side labels, lose faces/cells, reuse lower triangles as upper faces, and attempt to promote source status. Deleting a boundary tetrahedron can still leave an ordinary valid 3-ball, so the generic 3-manifold verifier is *not* incorrectly required to reject every possible boundary deletion; instead, the **specific canonical-slab and full-cut consistency contracts** require completeness.

Run `python -m unittest discover -s tests -v` under hosted Python 3.11 and 3.12, then check the post-merge `main` workflow.

**Outstanding Issue #58:** triangulate a source-correspondent connected ambient 3-manifold that **includes K**, identify the actual Figure 12.5/12.10 plumbing-patch positions, cut that manifold **along those source patches**, and derive **new** boundary spheres with four geometrically measured K punctures and full arrowed/pure Conway mapping-class framing. Only a separately replayable source-relative ambient isotopy/gluing witness would support certification. The roof-slab cutting infrastructure can inform that implementation but is not its substitute.
