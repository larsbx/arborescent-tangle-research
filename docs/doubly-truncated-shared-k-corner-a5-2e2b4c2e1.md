# A5.2e2b4c2e1 — One shared K corner with TWO disjoint positive-area caps

**Positive claim:** an exact-rational, physically embedded **15-vertex / 24-tetrahedron convex PL 3-ball** centered on an already verified two-valent source-*candidate* K vertex, with **both** original incident knot rays as true tetrahedral 1-simplices and TWO positive-area exposed five-vertex/four-triangle cap disks. The original vertex position is unchanged. Every unordered simplex pair passes independent exact 3D geometric intersection verification.

**Strict nonclaim:** These are *two attachment interfaces*, not yet a conforming edge tube attached at either cap. Source Figure 12.10, new source-cut four-punctured Conway spheres, framed isotopy and R± remain unverified.

## Why neither of the previous two single-sided corner solids can be reused

PR #66 produced the outgoing-only truncated corner 3-ball and PR #70 the incoming-only truncated corner 3-ball. Each individually occupies **15/16** of the same original physical octahedron's 3D volume and shares its true knot vertex. By the volume inclusion–exclusion inequality, the two 3D solids overlap in **at least 7/8 of the original octahedral volume**:

```
Vol(A∩B) ≥ Vol(A)+Vol(B)−Vol(full octahedron)
         = 15/16 V + 15/16 V − V
         = 7/8 V > 0.
```

Therefore naïvely welding two previously verified consecutive edge pieces at that shared K corner would produce **positive-volume tetrahedral overlap** and is mathematically invalid. In particular, joining identical center coordinates Q is not a physical manifold certificate.

## A source-K shared convex core: two *disjoint* tip truncations

Let the physical corner star be the rational octahedron with original K vertex P, generators `u,v,w∈ℚ³`, positive `Δ=det(u,v,w)>0`:

```
  O = {P+xu+yv+zw | |x|+|y|+|z|≤1}.
```

By the existing verified source K convention its incoming ray is `+u`. Its outgoing ray is either `+v` for a genuine bend or `−u` for a straight continuation.

Choose the precise rational truncation parameter **c=3/4**. Clip the convex octahedron with:

```
  incoming:   x≤c
  outgoing:   y≤c  (bend),  OR  −x≤c  (straight).
```

The source-preserving choice **c>1/2** is crucial: in the bent case the two removed tip cones `x>c` and `y>c` are disjoint since `x+y>2c>1` contradicts the original octahedron. In the straight case `x>c` and `−x>c` are also disjoint. Thus the TWO new nonzero-area cap disks neither overlap nor share even a cap vertex. A half-tip c=1/2 would leave two caps touching at the edge midpoint for an actual bent corner, which is inadequate for independently gluing neighboring tubes.

Each new cap is a physical rational diamond at its cut plane with a **true K center** at `P+c r (previous−P)` or `P+c r (next−P)`, where `r` is the independently validated corner radius. Its four diamond ring points lie at the remaining `|transverse_1|+|transverse_2|=1−c=1/4` cross-section. Both cap surfaces are oriented outward and triangulated as four triangles around their physical K center.

## Complete tetrahedralization and independent evidence

For each of the eight triangular facets of the ORIGINAL octahedron, clip the actual rational facet polygon against the two selected halfspaces, using exact rational convex polygon clipping. Triangulate the resulting convex polygon and cone every triangle to the original, still physical K vertex P. Finally add the two oriented four-triangle cap fans and cone these to P as well.

This construction has:

```
    15 unique physical rational vertices
    24 nondegenerate, positively oriented tetrahedra
    2 physically distinct exposed cap disks
    5 vertices, 8 edges, 4 triangular faces per disk; χ=1
    1 connected PL 3-ball; exterior χ=2
    2 genuine source K center-to-cap 1-simplex edges
    1 circular four-cycle tetrahedral link for each K ray
    C(24,2)=276 independently checked physical tetrahedron pairs
```

Each of two disjoint tips removes `(1−c)³·(2Δ/3)=Δ/96` of exact three-dimensional volume. Both tips therefore remove `Δ/48`, and the resulting exact volume is:

```
 V_shared = (4/3−1/48)Δ = (21/16)Δ
          = (63/64)·V_original.
```

The independent full PL-3D manifold checker reconstructs every signed physical tetrahedral volume, internal face adjacency/orientation, connectedness, Euler-2 boundary and all K/macro vertex links. The separately merged exact-Q³ convex simplex-contact kernel independently checks **all 276** tetrahedron pairs and rejects any physical nonlocal overlap, unexpected tangency or T-junction. Every outgoing and incoming disk must consist of actual exposed simplex faces, with positive normal along the genuine corresponding K ray, and must have no shared vertex with the other disk. Every K cap center must be independently recomputed from its actual adjacent original K segment, not from a port label.

The whole new solid is a convex subset of the previously audited original K-corner octahedron. Therefore the previous full exact original-K-segment and Conway closed-cube exclusion checks remain valid for the new shared core, independently of the local knot's projected braid crossings.

## Tests and next acceptance

`tests/test_doubly_truncated_k_corner.py` exercises real signed/zero-sector bent and straight source-candidate K corners, exact 24-tetrahedron topology, both independent cap disks, both K link cycles, outward normals, exact 63/64 volume and 276 physical tetrahedron intersection checks, the **quantified positive-volume overlap** in the naive prior two-solids strategy, rational two-halfspace polygon clipping, and corrupt/fabricated mesh/source status rejection.

**Issue #72 remains open.** Next, physically weld the incoming and outgoing disks of this **SINGLE** shared core to the neighboring K segment tubes, reusing PR #68/#71's conforming transition machinery but *not* superposing their overlapping single-sided corner solids. Rebuild the whole mesh and prove all actual pairs of tetrahedra satisfy geometric intersection constraints. Extend to full knot arcs and eventual one-valent proper Conway endpoint stars. Only afterward can a connected ambient including the whole K, all band triangles and existing Conway sphere surfaces support the actual source Fig. 12.5 slit and NEW four-punctured boundary derivation.

Existing states `DUAL_CAP_K_CORE_NOT_EDGE_GLUED`, `UNVERIFIED_GEOMETRY`, and every source-unverified flag stay unchanged. No source isotopy, ring relocation R± or full classification follows. GitHub-hosted Python 3.11/3.12 exact-head and post-merge main validation required.
