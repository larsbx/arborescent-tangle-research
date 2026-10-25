# A5.2e2b4c2d2a — A genuine knot edge embedded as a tetrahedral 1-simplex

**Status: LOCAL_K_EDGE_CONFORMING_NOT_GLOBAL_AMBIENT.** This is the first physically conforming **3D tetrahedral star around an actual geometric subedge of the distinguished knot K**, in the existing chosen-coordinate r=0 band candidate. It is not a tetrahedralization of the complete tangle or Bonahon–Siebenmann's source ambient manifold.

## Why PR #59's valid cutting mechanism was insufficient

PR #59 established an exact and fully edge/vertex-aware PL 3D cut and inverse reglue of a connected roof 3-ball. Its z coordinates satisfy z≥15/4, whereas every physical K segment in the source-candidate tangle stays in −1≤z≤1. Consequently **not a single K edge was a simplex of its tetrahedral ambient** and no source four-K-puncture Conway boundary can be derived from that 3D roof cut.

The new work tests the necessary next condition **at actual source K coordinates**: the ambient 3D triangulation must genuinely *contain the knot as an edge subcomplex*, not merely overlay another copy of its coordinates on unrelated tetrahedra.

## Exact constructive local theorem

The existing independently verified proper K arcs from A5.2e1 have two long, straight exterior return rails. The first runs at y=−6,z=0 and the second at y=+6,z=0 from x=−3 to x=R+3 (where R is the last Conway cube's eastern face), in opposite endpoint traversal orientation to our standardized increasing-x mesh convention.

Choose an exact proper subsegment of each existing rail by trimming **one unit at both ends**. Thus x ranges from −2 to R+2. For any rational radius `0<r≤1/4`, construct one actual physical rectangular 3-ball

```
T_± = [−2,R+2] × [±6−r,±6+r] × [−r,r].
```

At each of the two x terminal sections, store **five exact vertices**: a central knot vertex at (x,±6,0), and the four corners of the square width section. Triangulate each section as a disk with four triangles incident at its central K point. Extrude the four triangles along x, dividing each triangular prism into **three actual nondegenerate tetrahedra**, with a consistent canonical staircase face subdivision. The whole tube has:

```
  vertices:                 10
  tetrahedra:               12
  K-internal mesh 1-edge:   [center_left,center_right]
  triangles incident to its edge link: 4 (one 4-cycle)
  boundary Euler:           2
  connected 3-ball:         yes
  exact volume:             (R+4)·(2r)²
```

That core edge is not an invented new knot: `_long_rail_edge` reads the **actual stored independently checked K polyline** and verifies the new segment lies strictly inside its existing long return edge. The 3D simplicial-complex audit checks both core endpoints occur at their actual physical Q³ positions and **the knot core is literally one simplex edge of the 12-tetrahedron mesh**, with a **connected four-cycle edge link**.

The independent `verify_connected_tetrahedral_ball(vertices,tetrahedra)` checker (from PR #56) separately reconstructs every signed tetrahedral volume, interior face incidence and orientation, boundary Euler-2 2-sphere, dual face connectivity and **all vertex-link disk/sphere conditions**. The exact sum of all simplex volumes must equal the real tube's `length×(2r)²`, rejecting overlapping tetrahedra in this chosen cross-section/sweep normal form.

## Independent actual-K/Conway-ball clearance

`_clip_closed_box` independently clips **every actual K segment** against the entire closed 3D tube using exact rational interval bounds. Only the original designated outer long rail segment is permitted to touch the tube. It must cross its full longitudinal extent and contain the chosen core strictly inside its endpoints; ALL other K arcs and segments must be disjoint from the closed tube.

The tubes around y=−6 and y=+6 have disjoint closed y ranges for every allowed r. They are also disjoint from **every closed removed Conway cube**, since all these cubes satisfy |y|≤3 whereas the tubes satisfy |y|≥23/4. This is independent of whether particular half-twist sectors have 0, positive or negative weights. The chosen source-candidate ribbon may pass through these knot neighborhoods *at the actual K boundary*; we make NO claim that the complete ribbon is currently a conforming triangle-face subcomplex of the 3D tubes.

## Adversarial proof controls

`tests/test_k_conforming_tetrahedral_tubes.py` verifies both atom signs, signed/zero sector vectors of n=2–4, exact radii 1/16,1/8,1/4, K as a real 3D mesh edge with circular star, independent 3-ball vertex links and volume, clearance from all other K segments and Conway solids, exact segment-vs-box clipping (including boundary tangency), deterministic mesh vertices/indices, and negative controls for deleted/repeated tetrahedra, changed source arc, forged status, collapsed links and r>0/valence-one unsupported geometry.

The two tube volumes are positive and **their physical 3D solids are disjoint**. They remain two separate local 3-balls, not one globally connected ambient source-tangle manifold. Their inner K mesh edge endpoints are subsegment points, **not** a new source Conway boundary piercing: there are exactly zero new source-derived Conway spheres.

## Remaining issue #60 acceptance

This is an essential **first K-conforming meshing kernel**, not its completion. A valid full source ambient must:

1. Preserve **every** physical K segment and its connection vertices as embedded 1-simplices in **one connected oriented 3-manifold** with the correct n Conway boundary components, not two local proper subsegments.
2. Preserve the original framed global twisted band and each source-defined plumbing patch as compatible oriented 2D face subcomplexes/refinements, with exact nonlocal tetrahedron intersection checks.
3. Slit the **source-chosen** patches, using the full edge/vertex-aware cutting algorithm of PR #59, and independently reconstruct the resulting **NEW** four-K-punctured Conway sphere components from actual segment–triangle intersections (PR #57), not simply naming the original cubical boundaries.
4. Transport the full signed/pure Conway mapping class and source equatorial arrows, and provide an independent positive 3D inverse and relative-boundary ambient-isotopy witness before claiming source correctness or any R± classification.

`UNVERIFIED_GEOMETRY` and all earlier source-UNVERIFIED status fields remain mandatory. No incomplete piece or absent witness proves the source construction false.

Primary source: Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.1–12.5, source image and page provenance detailed in Issues #52/#58/#60.

Run `python -m unittest discover -s tests -v` with hosted Python 3.11 and 3.12 and inspect the exact final head and post-merge main.
