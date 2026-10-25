# A5.2e2b4c2e2 — Conforming two physical knot rays through one shared three-dimensional corner

**Scope:** a single locally connected rational 3D tetrahedral neighborhood centered at an existing true K vertex Q, containing proper portions of its two adjacent source-candidate K edges P→Q and Q→R. This is an actual **one-corner/two-edge** positive PL certificate, not an all-K or Bonahon–Siebenmann source-ambient construction.

## The essential mathematical obstruction

PR #71 previously proved the existence of a conforming 68-tetrahedron PL ball around one complete original K edge P→Q, but a naive union of two such source-candidate balls across a common Q is invalid. Each independently constructed ball contains almost a full octahedral Q corner: incoming-only and outgoing-only vertex star models share **at least 7/8** of the original corner solid's volume. Matching their knot vertex labels cannot turn that substantial physical overlap into valid simplex incidences.

PR #73 addresses the issue with ONE *doubly truncated* convex K-core at Q, whose two exposed positive-area diamond disks are genuinely disjoint. This milestone is the first attempt to attach an incoming AND outgoing knot-edge tube to that one shared core, rather than superimpose two separate corner 3-balls.

## Actual exact rational 3D construction

Start with the independently audited 15-vertex / 24-tetrahedron `DoublyTruncatedKCorner` at the **physical source-candidate knot vertex Q**.

Its incoming positive-area diamond cap lies on original P→Q at source parameter `1−(3/4)r_Q` measured from P; its outgoing diamond lies on Q→R at `(3/4)r_Q` for genuine bends (and along the correct opposite ray for straight continuations of unequal edge length). The underlying core physically contains its K center and both cap centers as actual tetrahedral 1-simplices.

Construct PR #63's two independent exact-Q³ `RationalKSegmentStar` tube sections on original K segments P→Q and Q→R. The source-centerline tube sections occur at the proper rational parameters 1/4 and 3/4 of their respective original K segments, with four-triangle square fans and real tetrahedral 1-simplex center edges.

Attach the incoming tube's **far** square section to the core's incoming diamond via four source-K-centered triangular prisms subdivided into **three 3D tetrahedra each**. Attach the outgoing diamond to the outgoing tube's **near** square by another four such prisms. All physical vertices are welded by **exact coordinate equality**, including the cap fan triangles that become true interior physical faces. Try the finite cyclic/dihedral mappings of each polygon with independently checked signed 3D orientation and exact self-intersection rejection.

The expected completed local convex-candidate complex contains

```
  15 vertices, 24 tetrahedra — ONE shared dual-cap K-corner core
  10 vertices, 12 tetrahedra — incoming actual P→Q K segment tube
  10 vertices, 12 tetrahedra — outgoing actual Q→R K segment tube
                12 tetrahedra — incoming square-to-diamond transition
                12 tetrahedra — outgoing diamond-to-square transition
  --------------------------
  35 unique Q³ vertices and 72 oriented 3D tetrahedra
```

No incoming-only/outgoing-only overlapping corner solids are added.

## Independent geometric proof checks

`SharedCornerTwoKDirections.audit()` rebuilds and compares the complete source candidate physical vertex table and all 72 tetrahedra against the canonical exact-rational generator. It then independently audits the combined complex with:

1. `verify_connected_tetrahedral_ball`: all oriented positive simplex volumes, face incidences/opposite orientations, every PL vertex-link disk/sphere, tetrahedral dual connectivity and a **single** connected Euler-2 boundary 2-sphere.
2. `audit_embedded_tetrahedra`: ALL `C(72,2)=2,556` unordered closed-tetrahedron pairs, using exact rational convex halfspace intersections (at most 56 exact plane triples per potentially touching pair). Only actual shared mesh vertices, whole edges or whole faces can meet; no nonlocal crossed tetrahedral interiors or T-junctions.
3. `_boundary_faces`: both actual four-triangle square and four-triangle diamond fan interfaces at both ends must become **interior physical 2-simplices**, not exposed or fictitious port labels.
4. The actual original K polyline must remain a strictly parameter-monotone **six-edge 1-skeleton chain** through the physical knot center Q, correctly passing from incoming tube start to its end and incoming diamond, then through Q to outgoing diamond and outgoing tube. The true source K segment indices must be adjacent and share the same unique original vertex.
5. Every nonincident original K segment is clipped against the **closed** transition tetrahedra using rational 3D halfspaces. All new transition cells must lie strictly outside every closed Conway cube solid. Existing corner and tube components independently retain their own source-K and Conway clearance certificates.

Each component is also independently audited before assembly. A guarded finite radius search permits exact rational halves; an unsuccessful cutoff is `UNVERIFIED_GEOMETRY`, never a proof of nonembedding.

## Test fixtures and source limitations

`tests/test_shared_corner_two_k_rays.py` covers real bent K-return geometry and straight K continuation for supported no-ring atomic source-candidate fixtures. It also checks disallowed source status, swapped attachment faces, omitted bridge tetrahedra, forged rotation/source labels, absence of a true proper two-valent corner and forbidden ring realization.

This remains a **LOCAL** geometric K vertex with proper subsegments of two incident knot edges. It does not include adjacent K vertex neighborhoods or Conway-puncture endpoints. Thus it is not yet a *whole proper K arc*, and the disjoint family of neighborhoods of all 2n strands would still not be a single ambient (M,K) without additional complement cells. It does not incorporate the entire twisted ribbon band or Conway boundary surfaces as conforming 2-complexes, identify the source Figure 12.10 plumbing patches, derive NEW four-punctured Conway boundary spheres from Figure 12.5 cutting, or prove any source-isotopy/R± classification.

**Status: `TWO_K_RAYS_LOCAL_NOT_SOURCE_AMBIENT`.** Preserve A5.1 `UNVERIFIED_GEOMETRY`, parent issues #74/#72/#65/#60/#58/#33/#16. Require both hosted Python 3.11 and 3.12 checks passing on final PR head and post-merge main.
