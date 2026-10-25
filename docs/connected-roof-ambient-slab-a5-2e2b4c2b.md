# A5.2e2b4c2b — First connected ambient 3D tetrahedral slab around all local slit charts

**Scope:** verified **connected** rational PL 3-ball around the roof plumbing charts, not the entire ambient S³ or boundary-marked source (M,K) manifold. It is a reusable spatial 3D ambient-cell component for Issue #54, not its complete source-isotopy acceptance.

Previously A5.2e2b4c2a (#55) built n **disconnected** product 3-balls with independent local cuts. This stage connects them with **genuine three-dimensional tetrahedral bridge cells sharing actual oriented faces**, not merely an abstract atlas pairing.

## Exact geometric construction

Take the verified A5.2e2b3b chosen-coordinate band candidate with n≥2 source-labelled Conway cubes. Every physically embedded exterior roof patch corresponds to one exact x interval `[Lᵢ,Rᵢ]` with y∈[-2,0], z=4. Adjacent intervals satisfy

```
Lᵢ₊₁ − Rᵢ = 1/2 + 2|aᵢ₊₁| > 0.
```

The previous atlas independently verified these bounds from the *actual* 3D tetrahedral coordinates. Assemble one sorted grid

```
x₀=L₀, x₁=R₀, x₂=L₁, x₃=R₁, ...,
x₂ₙ₋₂=Lₙ₋₁, x₂ₙ₋₁=Rₙ₋₁.
```

Build the **whole connected** exact 3-dimensional rectangular product

```
[x₀,x₂ₙ₋₁] × [-2,0] × [4−h,4+h],     0<h≤1/4, h∈ℚ.
```

Each of the 2n−1 x-slabs contains two z-layer intervals, split into two triangular prisms per layer and then **three positively oriented tetrahedra per triangular prism**. Thus every slab contributes 12 tetrahedra and every adjacent slab **shares a real triangulated yz face** using the *same global coordinate vertex indices*. The result is one connected physical tetrahedral 3-ball, not n independent charts.

Its independently verified counts are:

```
  source site intervals = n
  newly physical connecting intervals = n−1
  total x-slabs = 2n−1
  vertices = 12n
  tetrahedra = 12(2n−1)
  boundary triangles = 12(2n−1)+8
  closed boundary connected components = 1
  boundary Euler characteristic = 2
  exact physical volume = (x_last−x_first)·4h
```

The whole slab lies at z≥4−h≥15/4, strictly above all removed cubical Conway balls (z≤3) and all candidate proper K arcs (z≤1). It does **not** contain the full knot/tangle as an interior subcomplex, and is not yet the source ambient S³.

## Independent three-dimensional topology checker

`verify_connected_tetrahedral_ball(vertices,tetrahedra)` independently reconstructs all actual 3D face incidences from physical exact rational tetrahedra. It checks:

- Unique rational spatial vertices and nonrepeated tetrahedra, each with **strictly positive signed 3D volume**.
- Every face used once on the ambient exterior or twice on an internal manifold face, with exact opposite induced internal face orientation.
- Actual 3D tetrahedral dual adjacency has exactly one connected component (genuine face-connected 3D material).
- The exterior triangle surface is connected, every boundary edge used exactly twice, Euler characteristic **2**.
- **Every 3D vertex link** is separately rebuilt from the incident tetrahedra and checked to be a connected PL **2-disk** at a boundary vertex (χ=1, exactly one link boundary cycle) or PL **2-sphere** at an interior vertex (χ=2); rejects pinched/branched local 3D vertex stars.
- The total absolute oriented simplex volume equals the exact prescribed coordinate volume. A **canonical x-slab geometry** separately proves distinct tetrahedral interiors do not overlap nonlocally; this independent checker is not asserted to solve arbitrary 3D tetrahedron/tetrahedron geometric collision.

## Necessary source boundary obstruction remains explicit

The source-inspired two-arrowed **n physical midplane band rectangles** are now ACTUAL **interior triangular faces** of this connected tetrahedral mesh, not names describing distant boxes.

`ConnectedRoofAmbientSlab.face_barriers()` recomputes, from the actual 3D face incidence table, the **two incident tetrahedra and their oppositely oriented midplane faces** for each of the 2n marked triangles. The source ports are checked against original real collar-face charts.

A necessary distinction: omitting the *dual adjacency* across these marked triangles is only **dual-face detachment**. It does **not** duplicate all cut-edge and cut-vertex links and is therefore not yet a mathematically complete 3D manifold slitting operation. Because the n marked disk faces are separated by uncut connector slabs, the remaining tetrahedral dual graph stays **connected** (one component), unlike cutting the entire global midplane.

Most importantly **the selected roof faces have ZERO intersections with the actual K**: the knot is far below them. Consequently, this new connected ambient slab and its selected cut disks still produce **ZERO** source four-punctured Conway spheres. The source construction requires actual K to cross the newly reconstructed Conway boundary **in four transverse punctures apiece**, with the full GL₂ pure seam and equatorial-arrow framing determined by source geometry. A model that simply renames these unpunctured interior cut faces as source Conway spheres would be incorrect.

The public statuses are `CONNECTED_ROOF_SLAB_NOT_SOURCE_AMBIENT` and `DUAL_FACE_BARRIERS_NOT_PL_MANIFOLD_SLIT`; `require_source_global_slit()` fails closed with `UNVERIFIED_GEOMETRY`. Do not change any earlier A5.1 `UNVERIFIED_GEOMETRY` or A5.2e2b3b `OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED` state.

## Positive and adversarial tests

`tests/test_connected_roof_ambient_slab.py` exercises n=2,3,4 signed/zero sectors, n=3 exhaustive small sector vectors, positive rational h values 1/16, 1/8, 1/4, actual physical connector-face equality, volume, all vertex-link 3-manifold controls, exact source port/marked face constraints and two-sided orientation. Negative fixtures modify vertices, omit tetrahedra, duplicate volumes, forge topology statuses and request r>0 or n=1 unsupported models.

Use hosted Python 3.11 and 3.12 for exact PR head and post-merge main.

## Actual next obligation under Issue #54

Build a **source-identified and globally connected ambient 3D triangulation containing the actual properly embedded marked K as a subcomplex**, not this roof-only slab. Find and independently certify the actual Figure 12.5–12.10 source plumbing patches whose globally slit 3D boundary spheres *physically* have four K-piercings. Implement complete edge/vertex-aware manifold cutting, source equatorial-arrow / signed pure Conway gluing and an exact reverse ambient tangle/plumbing witness. Only then consider band-preserving cyclic rotations, r>0 ring realizations, both modified R± and classification.

The original mathematical source is Bonahon–Siebenmann §12.1–12.5, Figures 12.4–12.10 and 12.29, documented in earlier audited source files. No source-relative isotopy is inferred from this partial connected slab.
