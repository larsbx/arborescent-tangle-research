# A5.2e2b4c2d2c — Exact PL vertex stars at genuine two-valent K corners

**Positive scope:** a real tetrahedrally triangulated convex 3-ball with the **physical center at an already embedded source-candidate K vertex**, and both incident K subsegments as genuine tetrahedral mesh **1-simplex edges**. This extends the proper-subedge stars of PR #63; it does **not** join those stars into a single conforming whole-K 3-manifold, or source-certify a Conway slit.

## Geometric construction over ℚ

For an independently checked candidate K polyline, choose an *interior* two-valent vertex P with its two adjacent original vertices A,B. Put

`d₁=A−P`, `d₂=B−P`.

For a bent corner (`d₁×d₂≠0`), use the three exact rational octahedral generators

`u=r·d₁`, `v=r·d₂`, `w=r·(d₁×d₂)/max_j |(d₁×d₂)_j|`,

with `0<r≤1/4` rational. Their signed determinant is positive, and both P+u and P+v are strictly inside the corresponding original K edges for sufficiently small r.

For a collinear straight continuation (`d₁×d₂=0`, `d₁·d₂<0`), use `u=r·d₁` and two exact rational orthogonal transverse vectors `v,w` from the existing rational cross-product-frame kernel. The incident K rays correspond to *both* P+u and P−u. Zero-length or geometrically doubled-back original K vertices are refused.

Set the seven Q³ mesh vertices to

`P, P±u, P±v, P±w`.

For each of the eight triples of signs, take the tetrahedron `[P,P±u,P±v,P±w]` with its positive coherent orientation. These are the eight simplicial cones from the octahedron center over its eight triangular faces, forming one physically embedded PL **3-ball**. Its exact occupied 3D volume is

`vol=4/3·det(u,v,w)>0`.

This is an exact rational convex octahedral construction, not a floating-point tube or a drawing in a braid projection. The center P is a real **interior mesh vertex** with a spherical vertex link, and the two adjacent K rays are actual center-to-outer-vertex **mesh edges**. Each such edge has an independent four-cycle PL link in the tetrahedral complex.

## Independent physical collision checks

`RationalKCornerStar.audit()` re-fetches the existing source-candidate K's *actual* three consecutive spatial vertices, recreates the canonical octahedral coordinates and each exact signed tetrahedron independently, verifies both K edge links, and calls the separate full `verify_connected_tetrahedral_ball()` checker. That independently establishes positive volumes, internal-face orientation and incidence, every 3D vertex link as disk/sphere, a connected closed Euler-2 boundary, and the exact octahedron volume.

The physical solid is the **closed affine octahedron**

`P+x·u+y·v+z·w`, with `|x|+|y|+|z|≤1`.

Each original K segment is transformed to the octahedron chart by exact rational Cramer's rule and clipped against **all eight half-spaces** `±x±y±z≤1`. Only the two incident original edges may intersect the closed solid. Both must enter through the physical central vertex and leave the convex 3-ball along a **strict interior subsegment of their original K edge**; no neighboring K segment or unexpected tangent may meet the solid. This includes nonadjacent arcs, signed over/under braid crossings, the outer return rails, and adjacent sharp corners.

For each closed Conway cube, the validator independently computes strict rational axis-aligned separation of the whole octahedron's **actual convex hull vertices** from the cube box. Strict separation on even one axis is sufficient to prove disjointness of the entire two closed solids. This is conservative: if a true oblique separation is present but no axis certificate is found, shrink r or report `UNVERIFIED`, never infer nonembedding.

The deterministic radius search tests `r_max/2ᵏ` for `k=0,…,23`; unsuccessful bounded search is `UNVERIFIED_GEOMETRY` (not an isotopy/nonembedding theorem).

## Tests, source scope and next milestone

`tests/test_rational_k_corner_stars.py` covers rational bent, skew and collinear straight vertices; both incident physical center-to-ray edge links; actual signed braid/return/source-candidate corners for n=2 and n=3; chosen exact volume/3-ball topology; the complete closed-convex K-collision tests; altered radii, vertices/tetrahedra, bogus source status and unsupported r>0/valence-one cases.

Run hosted Python 3.11/3.12 on the exact final PR head and post-merge `main`.

A first hosted test run exposed a **fixture—not a geometry-counterexample—error**: an untwisted (zero-weight) corridor consists of a single K segment with only two vertices, hence it has **no actual interior two-valent knot corner**. The tests were corrected to select physical corners in the neighboring nonzero signed sector, and an explicit negative test now requires fabricated zero-sector corners to fail rather than silently creating a point. The final PR head must independently pass both Python versions.

**This is NOT a conforming global K neighborhood yet:** the previously merged arbitrary segment tubes intentionally cover only proper **edge interiors** and this PR supplies separate **vertex-centered** 3-balls; those solids are not glued to each other. A genuine connected conforming PL knot neighborhood must align each tube's transverse triangulated boundary with its adjacent vertex-star boundary through explicit mesh subdivision, preserving actual K geometry and links without collisions.

Further required: conforming 2D ribbon/conway surfaces within a single tetrahedral ambient, source Fig12.10 patch identification and framed GL₂ boundary transport, an edge/vertex-aware global cut producing **NEW** four-K-puncture Conway spheres, and an independently replayed source-relative isotopy. No ring R± or arborescent classification is established.

All results retain `LOCAL_CORNER_STAR_NOT_WHOLE_K_AMBIENT` / `UNVERIFIED_GEOMETRY`; Issues #62, #60, #58, #54, #52, #49, #33 and #16 remain open. Primary source Bonahon–Siebenmann §12.1–12.5, referenced in the existing source-figure audit documentation.
