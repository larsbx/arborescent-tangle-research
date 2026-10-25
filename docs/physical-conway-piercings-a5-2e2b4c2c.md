# A5.2e2b4c2c — Compute actual K piercings from physical Conway surface triangles

**Status: EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED.** This stage provides the missing exact *segment–surface intersection* verifier necessary to test four marked knot punctures on a Conway boundary. It does **not** construct the source's *new* 4-punctured Conway spheres by slitting the ambient 3-manifold.

## Mathematical distinction after PR #56

The project now has **one connected, oriented, tetrahedral roof slab** around all n chosen plumbing-patch candidates (PR #56). Its internal two-sided roof faces are physically labelled and face-incidence checked. **But K lies at z≤1 while the slab starts at z≥15/4**, so *these proposed roof slit faces have zero K punctures* and cannot possibly be certified as source four-punctured Conway sphere boundaries.

Separately, PR #38 created **n chosen physical cubical Conway boundaries** before any global ambient slit, and PR #40 created 2n properly embedded actual K arcs meeting all of their 4n named punctures. PR #51 checked that these arcs form a distinguished proper submanifold on the band boundary. These positive facts are necessary but insufficient for source reconstruction: the cubical boundaries were *pre-existing, not induced by cutting*.

The new checker independently **computes the physical K–surface intersections directly from every real cube triangle and every PL K segment** and only *then* compares the resulting coordinates with declared punctures. This is a reusable geometric test for future, genuinely source-cut-generated boundaries.

## Exact rational segment–triangle algorithm

`segment_triangle_hit(a,b,triangle)` takes a nondegenerate directed spatial PL segment and an exact rational 3D oriented triangle. It never samples a continuous surface, uses no floating epsilon and distinguishes:

- `DISJOINT`: segment does not meet the closed triangle.
- `TRANSVERSE`: the segment crosses the triangle's plane at one exact point, with nonzero 3D normal component (the point can be a proper K endpoint).
- `COPLANAR_POINT`: a tangency lying entirely in the triangle plane, touching at only one point.
- `COPLANAR_SEGMENT`: a positive-length coplanar overlap on/in the triangle.

For noncoplanar segments, calculate the exact oriented normal `n=(v₁−v₀)×(v₂−v₀)`, signed plane distances `d_a=n·(a−v₀)`, `d_b=n·(b−v₀)`, and (if 0≤t≤1)

```
    t = d_a/(d_a − d_b),
    p = a + t(b−a).
```

The exact point is inside the oriented triangle iff all three edge half-plane expressions `n·((v_j−v_i)×(p−v_i))` are nonnegative. A point landing exactly on a triangle mesh edge can be reported by multiple faces and is deduplicated at the physical-point level.

For a segment entirely in the triangle plane, use exact rational interval halfspace clipping against all three edge inequalities. This distinguishes a false zero-angle K touch from a truly transverse source pierce and rejects both kinds of coplanar contact from the proper-tangle theorem.

## Whole Conway boundary verification

`audit_existing_cube_piercings(cubes,knot_arcs)` consumes **only the actual 12 oriented PL triangles of each geometric cube and all K arc segments**, not a claimed incidence table. It independently computes all geometric point contacts with each boundary, rejects any coplanar contact, requires **exactly four distinct transverse physical points per cube**, and only afterwards compares them against the four marked Conway puncture coordinates.

It verifies each point belongs to exactly **one** K segment at exactly **one** start/end terminal, and that K approaches the corresponding west/east face from the exterior with nonzero normal component. The checker then separately invokes A5.2e1's full physical PL arc/all-solid-ball collision audit to rule out through-ball interiors or intersections between different K arcs.

The positive result for an existing cubical candidate with n boundary spheres is:

```
  pre-existing physical 4-punctured Conway cubes audited: n
  independently computed transverse K/surface points: 4n
  physical proper K arcs: 2n
  Conway boundaries derived from source ambient slit: ZERO
```

This distinction is **not cosmetic**. A global Fig12.5 plumbing slit must construct its boundary geometry from actual 3D cells, then perform the same kind of independent K triangle-contact audit. It cannot simply reuse pre-existing named cubes or populate four named punctures by looking up A5.1 ports. The resulting type `ExistingConwayPiercingAudit` is explicitly fail-closed under `EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED` and refuses source-boundary promotion.

## Tests and next geometric obligation

`tests/test_physical_conway_piercings.py` verifies exact skew plane intersection, transverse endpoints, rational halfspace clipping, coplanar segment overlap, tangential vertex contact, tiny exact rational gap, zero-area triangle rejection, n=2–4 signed sector candidate integration, computed physical 4n hits, incorrect/missing K arcs, through-ball shortcuts, and false source-cut-generated sphere counts. It also explicitly contrasts PR #56's zero-piercing roof face barrier with the **four positive actual K punctures** on each *already existing* cube.

Both hosted Python 3.11/3.12 and post-merge main must pass.

**Source obligation still OPEN:** A5.2e2b4c2d must triangulate a genuine connected ambient manifold **including K itself as a suitable PL subcomplex**, identify the source's exact plumbing disks, implement complete edge/vertex-aware global slitting, and independently derive newly created Conway boundary components with four *measured* punctures each. It must also transport the source Fig12.6 equatorial arrow and full signed/pure GL₂ Conway framing and construct independent inverse ambient-plumbing isotopy witnesses. No PR56/57 result establishes the source §12.1–12.5 classification or modified R± move.

Original primary source: Bonahon–Siebenmann §12.1–12.5, Figures 12.4–12.10 and 12.29; source PDF provenance in prior geometric docs.
