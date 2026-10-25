# A5.2e2b4c2e3 — ONE conforming tetrahedral tube between two dual-cap K corners

**Status: ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT.** This milestone replaces the invalid strategy of superposing individually certified two-ray K neighborhoods at consecutive K vertices. Its positive claim is **one actual source-candidate PL knot edge P→Q**, with P and Q both genuine two-valent source-candidate K vertices, supported by one shared conforming physical tetrahedral tube and two separate doubly truncated corner cores.

It does **not** yet construct an entire proper K arc, its Conway terminal stars, source ribbon framing or an ambient source (M,K), nor a source-generated four-punctured Conway sphere or ring R± isotopy witness.

## Why naive union is invalid

A5.2e2b4c2e2 (#75) constructs one genuine PL 3-ball around a physical two-valent knot vertex Q, connecting incoming and outgoing **proper subsegments** of both incident source K edges. Constructing the corresponding independently certified neighborhood at the previous K vertex P generates a SECOND copy of their shared edge's PR #63 tube. Both tube copies occupy the same positive volume with identical (or overlapping) tetrahedral geometry, not legal adjacent three-dimensional cells. The existing independent convex tetrahedron-contact checker cannot certify this as an embedded simplicial manifold. Matching the shared knot edge's labels or endpoints does not make such an overlap into a valid three-dimensional glue.

**The correct modular unit stores every geometric knot-edge tube exactly once and every geometric shared K-corner core exactly once.**

## Exact rational physical construction

For a supported no-ring `AtomicInstruction`, choose two consecutive **true interior** knot vertices P and Q in the SAME physical source-candidate K arc (original K edge `arc[i] : arc[i+1]`).

1. Build ONE `DoublyTruncatedKCorner` at P and ONE at Q. These are the previously independently certified exact convex 15-vertex / 24-tetrahedron shared-core models. Each has both incoming and outgoing exposed, **disjoint**, positive-area 5-vertex/four-triangle diamond cap disks, with the original K vertex preserved as an internal mesh point.
2. Build ONE `RationalKSegmentStar` covering the proper interior of the original K segment P→Q, with five rational vertices at each of its two physical square sections (source K parameter 1/4 and 3/4), 12 oriented tetrahedra, a real central K mesh edge, and its own independent source-K/closed-Conway-solid clearance checks.
3. Weld P's **outgoing** 4-triangle diamond disk to the tube's near square via FOUR physical triangular prisms triangulated into THREE tetrahedra each (**12 actual bridge tetrahedra**). Weld the far tube square to Q's **incoming** 4-triangle diamond disk by another **12 actual bridge tetrahedra**, with the correct outward-facing physical simplex orientation.
4. Select exact rational cyclic/dihedral section alignments only when the entire assembled 3D geometry passes all independent checks; up to eight exact radius refinements and eight alignments on each side are tried, and failure means UNVERIFIED, not a proof of nonembedding.

The expected **unique physical** mesh has

```
  left shared dual-cap K core:         15 vertices, 24 tetrahedra
  ONE common source K-edge tube:       10 vertices, 12 tetrahedra
  right shared dual-cap K core:        15 vertices, 24 tetrahedra
  left four-prism bridge:              12 tetrahedra
  right four-prism bridge:             12 tetrahedra
  --------------------------------------------------------
  TOTAL:                            40 vertices, 84 tetrahedra
  independently checked closed
  tetrahedron pairs:                 C(84,2) = 3,486
  exterior boundary:                 one connected Euler-2 sphere
```

Mesh vertex identity is established by **exact equality of stored rational Q³ physical coordinates**, never merely by matching A5.1 port identifiers. Every 2D cap/square attachment fan face must be present as an actual internal tetrahedral face. Conversely, P's unused incoming cap and Q's unused outgoing cap must remain genuine exposed 2D disks, ready for connecting the NEXT adjacent original source K edges without duplicating their corner volumes.

## Independent physical and topology audit

`SharedTubeBetweenKVertices.audit()` reconstructs and checks the canonical full mesh from the true K coordinates, rather than trusting a precomputed topology label. It runs the independent `verify_connected_tetrahedral_ball` checker for nondegenerate oriented tetrahedra, coherent face incidence/opposite orientation, all 3D vertex links, connectedness, and one Euler-2 boundary surface. A **separate** `audit_embedded_tetrahedra` computes each of **3,486** physical closed-tetrahedron pair intersections using exact rational halfspaces, permitting only the actual shared complete mesh simplex (empty, common vertex, entire edge or triangular face). Overlap between two ostensibly independent source K models cannot evade this verifier simply because their K edges follow the same line.

The original physical source K segment must be a *real five-edge tetrahedral 1-chain*

```
 P → P outgoing cap center → tube start center
   → tube end center → Q incoming cap center → Q
```

and each successive interior center must appear at a strictly increasing exact rational original K parameter in `(0,1)`. All five source-line subedges must be **literal tetrahedral 1-simplices** with circular interior K edge-link checks; the unused incoming/outgoing K rays in each shared corner remain actual mesh edges. Incorrect source arc, missing vertex, fabricated center, reversed K parameter, lost face or corner duplication must fail closed.

Every source K segment except P→Q is independently clipped against the **closed** NEW transition tetrahedra using rational convex halfspaces. Each new transition must be strictly outside every removed Conway closed solid cube, including braid strands that cross only in a planar projection. The original two convex corners and the tube each retain their own independent source K/cube clearance audits.

## Proof boundary and next steps

`tests/test_shared_tube_between_k_vertices.py` covers a true physically bent-to-bent original source K edge and a bent-to-straight original K edge, complete oriented manifold and convex geometric audits, actual 3D face welding and source K chain, tampered simplex/cap data, forged `SOURCE_VERIFIED` state and unsupported valence-one/ring models.

This gives one **chainable physical unit** with both *opposite* end caps exposed. It is NOT yet itself a full proper K arc: the next implementation must replace successive shared-corner and shared-edge components by a single conforming PL mesh, never duplicating either core or tube. One-valent source Conway puncture terminal stars are still absent. Even disjoint completed K arc neighborhoods would not constitute the single connected source ambient without complement 3D cells, a conforming supporting ribbon and Conway boundary sphere face complexes.

Bonahon–Siebenmann's source Figure 12.10 plumbing patch, Figure 12.5 newly slit **four-K-puncture** boundary, full signed/pure Conway framing and source-relative isotopy, ring R± and classification all remain unproved. Preserve `ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT` and A5.1 `UNVERIFIED_GEOMETRY`; Issue #76 and parent source issues remain open until positive verified acceptances are established.

Validation: GitHub-hosted Python 3.11/3.12 exact-PR-head and post-merge main passing.
