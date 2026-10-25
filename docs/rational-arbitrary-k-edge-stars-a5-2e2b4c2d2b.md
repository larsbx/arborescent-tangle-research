# A5.2e2b4c2d2b — Exact rational K-edge tetrahedral stars in arbitrary spatial directions

**Positive scope:** construct a physically embedded, locally conforming tetrahedral 3-ball around the **proper interior subsegment of any supported existing K edge** in the zero-ring source-candidate geometry, including vertical return lines and diagonal over/under signed braid segments.

**Not claimed:** a conforming tube around an entire K edge including its endpoints, a PL knot **vertex corner star**, a connected whole-K ambient, a global ribbon/Conway face-conforming tetrahedralization, a source-derived four-puncture Conway slit, or a Bonahon–Siebenmann R±/isotopy theorem.

## Exact rational transverse frames without length normalization

For an original independently audited knot segment `a→b`, put `d=b−a≠0`. Select the first standard integer coordinate axis `e_i` for which `d×e_i≠0`, and define

```
 u = (d×e_i) / max_j |(d×e_i)_j|
 v = (d×u) / max_j |(d×u)_j|.
```

All arithmetic is exact in ℚ, using only guarded nonzero rational division; no square roots, angles, trigonometric transcendentals or floating-point geometry. The right-handed frame `(d,u,v)` satisfies

```
 d·u = d·v = u·v = 0,
 det(d,u,v) > 0,
 ||u||∞ = ||v||∞ = 1.
```

For rational radius `r∈(0,1/4]`, use the *proper subsegment* parameter `t∈[1/4,3/4]`. Its enclosing 3D parallelepiped has exact affine coordinates

```
  X(t,s,z) = a + t·d + s·u + z·v
  t∈[1/4,3/4], s,z∈[-r,r].
```

Each of the two transverse sections contains the actual knot core at `(s,z)=(0,0)` plus four rational square-corner points. Four triangles around the core extrude longitudinally into **12 nondegenerate oriented tetrahedra**, with **10 real spatial vertices**. The exact core edge `X(1/4,0,0)→X(3/4,0,0)` is literally one **interior mesh 1-simplex** and is guaranteed to lie inside the previously verified K segment.

## Independent full 3D and geometric checks

The existing independent PL 3D manifold checker derives every tetrahedral face incidence/orientation, vertex link disk/sphere, one Euler-2 connected exterior boundary and exact signed simplex volumes. The separately implemented K core-edge-star validator recomputes the four tetrahedral 2-faces opposite that core and certifies one **four-cycle link**. The physically occupied exact volume is

```
  vol = (1/2)·det(d,u,v)·(2r)².
```

New exact spatial intersection machinery separately rechecks the source K geometry:

- **Segment vs oblique tube:** map both arbitrary physical endpoints into local `(t,s,z)` coordinates by exact 3×3 Cramer determinants, then clip the closed segment against the three rational intervals. The only permitted original K segment that intersects the 3D solid is precisely the declared source segment; its intersection must be the intended exact `[1/4,3/4]` subedge. Adjacent sharp corners and all other K strands must miss the **closed solid**, not just its interior.
- **Tube vs every closed Conway cube:** exact separating-axis test between the parallelotope's eight actual Q³ corners and the cube's eight corners. It uses the three tube face normals, the three cube normals and all nine cross-edge axes; a strict projected separation certificate proves the entire **closed** convex solids are disjoint.
- **Deterministic rational clearance search:** start with r≤1/4 and halve exactly for at most 24 attempts. The first candidate passing all complete audits is returned. If no candidate passes, the result is **UNVERIFIED**, not a nonembedding proof.

This avoids invalid inference from a tube's axis-aligned bounding box, which may overlap a Conway cube even when the actual oriented solid remains disjoint, or from projected planar braid crossings, where the K strands are physically separated by rational z-height.

## Source boundary and validation

`arborescent/rational_k_segment_stars.py` implements the local rational frame, physical 3D tetrahedral star and guarded geometry intersection predicates. `tests/test_rational_k_segment_stars.py` exercises x-, y-, z- and rational-skew directions, *actual signed-braid diagonal K segments* and both return-line directions, positive exact radii, near-contact rational gaps, false source status, changed tetrahedra, altered source arc selection, and ringed/valence-one unsupported cases.

This moves the K-conforming tetrahedral primitive from the two special x-axis return rails of PR #61 to **arbitrary local K segment directions**. It does NOT yet create the **corner vertex stars** that would glue these partial edge tubes into a conforming whole-K ambient. Every tube retains `ARBITRARY_EDGE_STAR_NOT_GLOBAL_K_AMBIENT`; existing `UNVERIFIED_GEOMETRY` remains enforced.

The next Issue #62 stage is to construct exact joint tetrahedral neighborhoods around consecutive 2-valent K vertices **without self-penetration** and a global conforming mesh for all 2n source K arcs, framed band triangles and Conway boundary spheres. Only an actual source-identified Figure 12.5/12.10 global slit generating NEW four-punctured Conway spheres could support source-isotopy certification. No R± or Theorem 12.21 count follows.

Hosted CI: run `python -m unittest discover -s tests -v` on Python 3.11 and 3.12; inspect the exact PR head and post-merge main.
