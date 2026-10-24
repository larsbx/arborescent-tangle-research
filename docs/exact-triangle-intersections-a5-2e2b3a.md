# A5.2e2b3a — Exact triangle/triangle geometry intersection kernel

**Research target:** Issue #46, the geometric self-intersection checker missing from A5.2e2b2 (#45). This changes the verification *capability* but **does not** certify Bonahon–Siebenmann Figure 12.10 or R±.

The previous `verify_spatial_surface` in `conway_roof_collars.py` reconstructs a valid point-welded triangulation of a 2-manifold with boundary. **Combinatorial manifoldness is not enough**: two physically crossing triangles can have disjoint coordinate vertices/edges and still fail to define an embedded surface. A2--A5.x abstract plans never represented this case. It must be checked geometrically, independently of any returned manifold invariants.

## Exact finite algorithm (no floating tolerances)

`triangle_contact(A,B)` accepts two explicitly nondegenerate triangle triples in Q³ and returns a typed contact class `DISJOINT`, `POINT`, `SEGMENT`, or `AREA`. Coordinates are Python `int`/`Fraction`; even a separation of `1/10¹²` is detected exactly.

1. Reject zero-area input triangles and disjoint axis-aligned bounding boxes.
2. Construct the exact plane normals `n_A = (A₁−A₀)×(A₂−A₀)`, `n_B`; evaluate exact signed vertex distances to the opposite plane. Strict same-side separation proves disjointness.
3. For coplanar triangles, project along an axis with nonzero normal coefficient to an exact two-dimensional affine chart. Clip triangle A against the three inward half-planes of B using rational segment interpolation (Sutherland–Hodgman). Lift the resulting polygon back into the common 3D plane. Any noncollinear triple proves **positive-area overlap**. Zero-area intersections are classified by their farthest two endpoints as a segment, or one point.
4. For nonparallel planes, the triangle-plane intersection for each triangle is the exact segment/point given by signed edge interpolation. Both lie on the same rational intersection line. Intersect their scalar closed intervals along a nonzero coordinate of `n_A×n_B`; reconstruct exact 3D overlap endpoints. The result is empty, a singleton, or a segment.

The implementation does not sample or approximate the intersection. Its scope is closed linear triangles; smooth bands or ambient isotopy are not implied.

`audit_embedded_triangle_mesh(triangles)` runs this independent intersection test over all face pairs, using bounding-box rejection as a cheap exact broad phase. It accepts only:

- No contact between a triangle pair with no shared mesh vertices.
- A **single point** contact exactly equal to the one declared common vertex.
- A **full common-edge segment** exactly equal to the two common mesh vertices.

It rejects shared-vertex penetration with additional interior intersection, coplanar overlapping areas, transverse crossings without shared vertices, partial edge overlap, and unexpected tangency. The existing `verify_spatial_surface` remains responsible for combinatorial link/manifoldness and edge incidence. Neither check substitutes for the other. Pair checking is worst-case quadratic in triangle count; optimizing that broad phase is a later engineering task if warranted.

## Adversarial fixtures

`tests/test_triangle_intersections.py` covers exact coplanar shared vertex, shared edge, containment, area overlap, skew signed-plane crossings and penetrating segments, noncoplanar shared edges, tangencies, tiny rational separation, nondegenerate controls, and a valid two-triangle PL square versus added illegal crossing faces. Contact symmetry is checked under exchanging triangles.

## What this does NOT establish

The currently chosen roof collars of PR #45 lie on portions of their Conway spheres. Even a positive collision audit of those mesh triangles would not establish that they are the source's correct **transverse** Conway attachment, nor determine the pure mapping class/framing that is invisible to a four-puncture permutation. This new checker is a necessary independent geometric validation tool for upcoming off-boundary/source-faithful collar work, not a claim of full r=0 Fig.12.10 realization.

In particular the unresolved research obligations remain:
- Source Fig.12.10 plumbing patch and Conway sphere parametrization.
- Exact geometrically compatible **off-boundary** band attachments.
- Independent positive ambient relative-boundary isotopy witnesses for cyclic rotations and all signed source fixtures.
- Actual r>0 rings and both source R± moves.

Do not promote any A5.1 `AtomicPlumbingPlan` beyond `UNVERIFIED_GEOMETRY`, close #46/#43/#33/#16, or count real isotopy classes using this kernel.

Run `python -m unittest discover -s tests -v`; require hosted Python 3.11 and 3.12 on exact PR head and post-merge `main`.
