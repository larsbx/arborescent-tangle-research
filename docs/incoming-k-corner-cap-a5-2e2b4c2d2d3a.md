# A5.2e2b4c2d2d3a — Physically exposed incoming K-corner cap

The previously merged PR #68 physically connected a truncated source-candidate knot-corner **outgoing** diamond disk to the **near** square disk of a real tetrahedral segment tube, proving a single local corner→edge 3-ball using independent exact-Q³ tetrahedral collision and PL topology audits. It did not address the far end of the same tube: the *next* knot vertex's **incoming** ray still terminated in one zero-area octahedral apex.

This milestone constructs the **dual positive-area incoming cap** without claiming that the tube is already attached at both ends.

## Exact construction

Take PR #64's independently checked 7-vertex / 8-tetrahedron convex rational octahedron centered on an actual proper source-candidate K vertex P, with generators `u,v,w`, `det(u,v,w)=Δ>0`. By the verified K source convention the **incoming** physical ray P→previous K vertex is along +u, for a signed bend or a straight continuation.

Remove the positive-volume tip cone above the affine section `u-coordinate = 1/2`. Keep the center and the five original octahedron apices other than +u; introduce the cap center at `P+u/2` and four actual ring vertices halfway between the removed +u tip and the neighboring apices `+v,+w,-v,-w`. The cap has **five distinct exact rational vertices, four nondegenerate triangles and a four-edge simple boundary circuit**. Each triangle is an actual exterior face of the same three-dimensional tetrahedral mesh.

Retriangulate the remaining convex polyhedron by coning four untouched facets, eight subdivided clipped facets and four new cap triangles to its original K center P. This yields exactly **11 physical vertices and 16 oriented tetrahedra**. The independently audited volume is

```
V_before = (4/3) Δ,
V_after  = (5/4) Δ = (15/16) V_before.
```

The incoming original K ray P→previous vertex now has a *real 1-simplex segment* P→new cap center on the boundary of the triangulated region, whose edge link is a four-cycle. The untouched outgoing ray toward the next source K vertex remains a real mesh edge: its link is a six-cycle if the two directions form a bend, and a four-cycle if they are opposite along one straight line. Every actual K vertex position and signed crossing height remains unchanged.

## Independent geometric checks

`arborescent/incoming_k_corner_caps.py` builds the 3D simplicial complex and uses the separate `verify_connected_tetrahedral_ball` verifier for signed 3D volumes, proper face incidence, connected ball boundary of Euler χ=2 and all 3D vertex links. It separately reconstructs the **actual exposed** incoming interface disk using tetrahedron face incidences and checks its 5 vertices, 8 edges, 4 triangles and χ=1; no abstract cap label can replace a physical mesh disk.

The original whole octahedron from PR #64 already passed independent exact convex 3D clipping against all unrelated K segments and removed Conway solids. Truncation leaves a subset of that certified solid, so it cannot introduce new spatial collisions.

The K cap center must be at the exact coordinate `P+(radius/2)(previous−P)`; its position is independently recomputed from the true stored source-candidate K polyline, not inferred from a named port.

The test module verifies signed source geometry and two K ray link types, volume 15/16, exposed disk, complete tetrahedral manifold topology and forged mesh/status failure. Python 3.11/3.12 hosted exact PR head and postmerge main must pass.

## Remaining Issue #69 acceptance

This is a **real positive-area 3D incoming attachment interface**, but is still **not welded to the outgoing edge tube** from the preceding corner. The next geometric theorem must conform the far square section of PR #63's tube to this new incoming diamond cap by a true physically embedded transition, with shared tetrahedral faces and correct orientation, independent nonlocal tetrahedron collision checks, no T-junctions, original K centerline 1-skeleton continuity and complete volume/link verification.

After that, the complete two-ended source-candidate edge segment can be audited from P to Q. Entire K arcs, terminal Conway puncture stars, a single whole-K ambient, source Figure 12.10 plumbing-patch positions and any newly source-cut four-punctured Conway spheres remain further OPEN obligations. The source Figure 12.5 band-preserving isotopy, R± moves and full classification are NOT established.

**Status: `INCOMING_K_CAP_NOT_TUBE_GLUED`**, with A5.1 `UNVERIFIED_GEOMETRY` unchanged.
