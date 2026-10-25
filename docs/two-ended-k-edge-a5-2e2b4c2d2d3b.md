# A5.2e2b4c2d2d3b — Actual two-ended K-edge 3D attachment

**Scope:** An exact physically welded PL 3-ball containing one **complete original source-candidate K segment between two actual two-valent knot vertices**, with the original knot represented by a continuous tetrahedral mesh 1-chain. This is a *local* source-unverified model; it does not conform to all K arcs, the ribbon/Conway triangulation or a global Bonahon–Siebenmann (M,K) ambient, and derives no NEW four-punctured Conway sphere.

## Necessity

PR #68 physically welded the outgoing cap of a K corner P to the **near end** of a proper subsegment tube in the original K edge P→Q, but the tube's *far* end remained unsupported. PR #70 supplied the dual **positive-area incoming cap** of the actual NEXT K corner Q, avoiding a zero-area apex/nonmanifold point identification.

A genuine two-ended construction now needs to weld that incoming diamond disk to the **far positive-area square section of the SAME tube**. It must verify the entire physical mesh rather than concatenating proof labels.

## Exact physical construction

Input is the already verified PR #68 `WeldedKCornerToEdge` (21 physical vertices and 40 tetrahedra), plus the PR #70 `IncomingTruncatedKCorner` at precisely the **next true vertex** on the same original K polyline.

The near side's physical far-tube square is the real 5-vertex fan cross-section at original source K-edge parameter 3/4; the new incoming cap diamond center is at original parameter `1 − r_Q/2`, which is greater than 3/4 for the supported exact rational radius `r_Q≤1/4`. They have different rational transverse frames and generally lie in nonparallel 3D planes.

For each of four source K-centered fan sectors, connect one triangle of the square section to the corresponding triangle of the incoming diamond using a physical triangular prism, tetrahedralized into three coherently oriented 3-simplices. Try the four cyclic section rotations and both direction choices (eight possible dihedral correspondences), guarded by full geometric collision checks. The 12 new transition tetrahedra join the stored 40 near-side tetrahedra to the 16 incoming-corner tetrahedra, **welding equal Q³ mesh vertices and triangles** rather than merely naming equivalent source endpoints.

The intended finite certificate therefore has:

```
  21 near physical mesh vertices + 11 new incoming vertices = 32
  near corner + first transition + tube:                  40 tetrahedra
  second incoming diamond-to-square transition:          12 tetrahedra
  next K corner incoming octahedral cap:                   16 tetrahedra
  TOTAL:                                                    68 tetrahedra
  physical K 1-chain:          P→first cap→tube start→tube end→incoming cap→Q
  independent tetrahedra checked:                  C(68,2)=2278 pairs
  global exterior boundary:       1 connected PL 2-sphere, Euler χ=2
```

The physical knot chain must consist of **five actual tetrahedral 1-simplex edges**, all lying in strict increasing parameter order along the **same original source-candidate PL K edge** P→Q. The old source K vertices P and Q are not moved or substituted with diagrammatic labels.

## Independent geometric verification

The source-independent `verify_connected_tetrahedral_ball` checker recomputes all signed tetrahedral volumes, internal face incidence/orientation, dual connectedness, boundary Euler topology and every vertex link.

The **separate** exact 3D collision checker `audit_embedded_tetrahedra` enumerates all 2278 unordered pairs of physical closed simplices, deriving the bounded convex intersection by solving exact rational supporting-plane triples. Two cells may meet only on the actual shared mesh simplex (empty, shared point, complete shared edge or whole common triangular face), never on an invented contact, T-junction, nonlocal crossing or volume overlap.

Every face of the physically attached square and diamond fan disks must be a **true internal shared mesh face** of the combined 3-ball. The new transition tetrahedra must avoid every original K segment other than the claimed joined P→Q edge, checked by exact closed segment–tetrahedron clipping, and remain strictly outside every closed Conway solid ball.

At most eight exact-rational radius halvings with eight cyclic/dihedral alignments each are tried. If no compatible 3D geometry passes, report `UNVERIFIED_GEOMETRY`, **never** conclude that Bonahon–Siebenmann's source knot has no embedding or no isotopy.

## Tests and remaining source gate

`tests/test_two_ended_k_edge.py` checks a physical bent-to-bent source K edge and a bent-to-straight segment, including all geometry, K-chain and 2278 independent tetrahedral pairs. It rejects a fabricated source certification, wrong attachment sections/rotations, missing tetrahedra and impossible source corners.

Even a verified *single* two-ended K edge is not the complete 2n-arc K regular neighborhood or an ambient conforming to the complete ribbon and Conway sphere surfaces. The one-valent source puncture terminal stars, geometric assembly of all separate K arcs into an actual connected ambient, source Figure 12.10 plumbing patches and Figure 12.5-derived NEW four-punctured Conway spheres are outstanding. Therefore it is not source framed isotopy or ring R± certification. Keep `TWO_CORNER_K_EDGE_NOT_WHOLE_SOURCE_AMBIENT` and A5.1 `UNVERIFIED_GEOMETRY` explicit.

Hosted Python 3.11/3.12 exact PR-head CI and post-merge main checks required before claiming merge.
