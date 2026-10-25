# A5.2e2b4c2d2d2 — Exact tetrahedral collision kernel and a conforming K corner-to-edge weld

**Status: ONE_K_CORNER_EDGE_GLUE_NOT_SOURCE_AMBIENT.** The intended finite geometric theorem is one physically embedded 3-ball spanning a **single actual two-valent source-candidate K corner** and a *proper interior segment of its outgoing K edge*. The source Fig12.5/12.10 slit, full K/Conway ambient and R± are **not** implied. This is the first actual attempt to **WELD** the independently verified PR #66 corner 3-ball to the PR #63 outgoing K tube, rather than merely storing nearby geometry.

## Independent exact rational tetrahedron-contact verification

Earlier `verify_connected_tetrahedral_ball` checked only oriented tetrahedron volumes, shared simplex incidence, vertex-link manifoldness, and global boundary topology. A mesh can satisfy these combinatorial checks while unrelated tetrahedra cross in space.

The new `arborescent/exact_tetrahedron_contacts.py` **independently compares every pair of physical closed 3-simplices**. For each tetrahedron, it derives its four exact outward supporting planes `nⱼ·x ≤ cⱼ`, with sign oriented using the opposite tetrahedral vertex. Two intersecting tetrahedra are exactly the intersection of all eight rational halfspaces.

Because this intersection is a bounded convex polytope, each extreme point lies on at least three linearly independent supporting planes. The code enumerates the `C(8,3)=56` possible triples, solves the exact 3×3 rational systems using Cramer's cross-product identity, and retains points satisfying **every** halfspace inequality. No numerical epsilon, random sampling, floating-point collision heuristic or unproved bounding-box-only certificate is used.

For a valid embedded tetrahedral mesh, the **entire physical intersection** between two cells may equal only the convex hull of the vertices that are **literally shared** by their mesh incidence tables: empty, one point, one entire common edge, or one entire common triangular face. Any computed extreme point outside that shared simplex is illegal, even if it corresponds to a tangency, partially intersecting face, or true full-dimensional overlap. Duplicate physical tetrahedra are rejected. Axis-aligned bounding-box rejection is a safe exact broad phase only.

This independent geometric kernel can falsify models that pass the earlier vertex-link and Euler checks. The generic current complexity is quadratic in number of tetrahedra and up to 56 plane-triple solves per non-rejected pair; mesh spatial acceleration is a future engineering task.

## Chosen physical corner-to-edge transition geometry

For an actual independent PR #66 `TruncatedKCorner` and a PR #63 `RationalKSegmentStar` on its **same original outgoing K edge**, their transverse 2D attachment disks are:

- **Near:** 5 physical rational coordinates, a central K point and four diamond corners, with exactly four exposed corner-fan triangles.
- **Far:** 5 physical rational coordinates, a central K point and four square corners, with exactly four exposed tube-start fan triangles.

Both centers are real points on the *same* independently verified original source-candidate K segment. The near point is at the original K-edge parameter `r_corner/2 ≤ 1/8`; the far tube start is at parameter `1/4`. They are **not** automatically in parallel cross-section planes, so the code cannot declare the bridge embedded merely from the existence of compatible cyclic four-point labels.

For each of the four cyclically corresponding triangular fan sectors, build a 3D triangular prism between the near and far triangle. Subdivide it using three actual oriented tetrahedra with a **consistent common diagonal** on shared sector side walls. The total transition contains **12 tetrahedra**. All source coordinates are exact Q³. Try at most four orientation-preserving cyclic square/diamond correspondences and bounded rational tube/corner radius refinements; **accept only** a positive fully audited three-dimensional mesh.

The assembled candidate contains the original corner's **16 tetrahedra**, the transition's **12**, and the outgoing edge tube's **12**, for **40 tetrahedra** total. Its 21 physical vertices are welded by **coordinate equality**, never by a pure A5.1 port alias. The four near/far 2D disk triangles must become **internal actual tetrahedral faces**, leaving one connected closed Euler-2 outer boundary. The complete assembled object is independently passed through both:

1. The existing coordinate-only 3D manifold checker: face orientation/incidence, true vertex-link disks/spheres, connectedness and volume.
2. The new **all-pairs physical tetrahedron intersection** checker: no nonlocal 3D self-overlap, illegal skew crossing or unintended tangency.

The distinguished K subchain **must** be a true ordered 1-skeleton path from the original K corner through the outgoing truncated-cap center, through the transition to the tube start center, and finally to the end of the certified tube core. Every edge of this chain must be a literal 3D tetrahedral **1-simplex** and lie on the correct existing K original source segment, with strictly increasing exact rational edge parameter. Merely drawing a line through the cells or reusing the same two end labels is forbidden.

The transition faces are independently clipped against the **entire original K graph** using exact rational closed segment/convex-tetrahedron clipping; *only* the declared outgoing original K segment is permitted to enter, while the incoming K segment belongs to the pre-existing corner star and no unrelated source strand may intersect the transition. Every transition tetrahedron also needs a strict rational convex-hull-versus-closed-Conway-ball separation certificate.

## Distinct mathematical outcomes

A **positive verified** selected-coordinate PL join is not a proof of Bonahon–Siebenmann's particular global twisted band/framing. These local three-dimensional solids are not yet connected for all K edges and vertices, do not include all ribbon and Conway surfaces as conforming meshes, and have neither a global ambient S³/(M,K) triangulation nor any **new four-puncture sphere reconstructed by a source Figure 12.5 slit**. The framing-compatible equatorial Conway arrow/pure signed GL₂ map, cyclic band-preserving isotopy and r>0 rings/R± remain open.

A **bounded search with no verified join** returns `UNVERIFIED_GEOMETRY`, not proof that no embedding exists. This is important for difficult oblique frames near signed braids: without a replayable orientation and geometric clearance certificate, a plausible picture is insufficient.

## Test and review contract

`tests/test_welded_k_corner_edge.py` exercises closed tetrahedral face/edge/point contacts, true interior overlaps without shared mesh vertices, tiny rational gaps, exact segment/tet clipping and source-candidate bent and straight K corner welds, with tampered source status/face and illegal attachment controls. Require hosted Python 3.11 and 3.12 passing on exact PR head, then post-merge main; any identified geometric overlap must **repair the geometry or fail closed**, not relax the collision kernel.

Source background: Bonahon–Siebenmann §12.1–12.5, Fig12.5–12.10/12.29; full page-image provenance already audited in Issues #52, #60 and #65. Source certificate is **not** granted by this local theorem.
