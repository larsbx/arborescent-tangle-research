# A5.2d — Exact four-marked Conway spheres and complementary half-twist corridors

**Status: genuine partial 3D embedding; OPEN_CHAIN_NOT_ATOMIC.**

## Visual primary-source audit

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.2 printed p.142 **Fig.12.8** and §12.5 printed p.158 **Fig.12.29**, inspected as PDF page images (PDF indices 153 and 169 respectively):

https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The diagrams show **n Conway holes**, an outer **returning twisted band** with two long tracks and **n labelled half-twist regions** `a_1 ... a_n`. Fig.12.29 also shows **r closed rings near one end**. In particular, holes are not merely abstract port names, and the rings are not free-standing additional half-twist words. These structures must be embedded independently.

A5.2c built the two-strand twist regions in Q³ but had no geometric Conway boundaries. A5.2d supplies an explicit orientable four-marked PL 2-sphere for each port, then embeds all `a_i` local two-braid corridors between such disjoint boundary surfaces. We explicitly introduce **one auxiliary left anchor cube** in place of the not-yet-constructed source return band. This changes the geometry from the source: the resulting open chain **is not a realization of Fig.12.8 or Fig.12.29**, even for r=0.

## PL construction and exact audit

For centre `c ∈ Z`, `MarkedConwayCube(c)` is the boundary of the 3D cube:

```
c-2 ≤ x ≤ c+2,    -3 ≤ y,z ≤ 3.
```

Its four marked punctures, in the chosen cyclic equatorial order, are:

```
NW = (c-2,-2,0)   NE = (c+2,-2,0)
SW = (c-2,+2,0)   SE = (c+2,+2,0).
```

The boundary has 8 vertices, 18 edges and 12 outward oriented triangles. `audit` independently checks triangle orientation using exact cross products, cancellation of directed edges, each edge having exactly two incident faces, and Euler characteristic `8-18+12=2`. Its exact face-interior punctures are distinct, away from mesh vertices and edges. Here a **puncture** is an individually marked point on a real embedded sphere (the surface with that point removed is the four-punctured sphere); it is not yet a source-compatible parametrization of a Conway *seam*.

For an A5.1 atom with cyclic sector weights `(a₀, ..., aₙ₋₁)` and **r=0**, let `c₀=0` and

```
cᵢ₊₁ = cᵢ + 5 + 2|aᵢ|.
```

There are **n+1 cubes**: the first is an artificial left *anchor*, and the next n are each labelled by one genuine A5.1 `Port(atom.address, i)`. Thus every existing signed sector is preserved and the source's n free boundary holes are present as actual 3D PL surfaces. The anchor is an acknowledged topological incompleteness, NOT an extra source Conway port.

Each intervening signed sector `aᵢ` is realized by the genuine spatial two-braid `SpatialTwoBraid((sgn(aᵢ),)*|aᵢ|)`. Transform its exact local coordinates `(x,y,z)` by:

```
(x,y,z) -> (cᵢ+2+y, x, z).
```

The image connects the two punctures on the *east* face of cube i to the two punctures on the *west* face of cube i+1, using the signed crossings' actual 3D height separation. The intervening longitudinal gap has length `2|aᵢ|+1`. Even for `aᵢ=0`, its two disjoint straight segments have a positive-length open gap.

**Local embedding proof:** In any corridor, longitudinal coordinates strictly increase along both arcs. Internal points lie strictly between the adjacent cubes' facing faces; no other cube occupies that open gap. All curve points have `|y|≤2`, `|z|≤1`, and the previous A5.2c audit checks each pair's 3D crossing disjointness. Distinct corridors have disjoint open x-intervals. The cubical surface triangulations are oriented closed surfaces, each lying in a disjoint x-slab. Hence the embedded spatial arcs and the marked boundary surfaces are mutually disjoint except at their declared boundary endpoints. The two unused west marks on the anchor and the two unused east marks on the terminal source cube remain **OPEN**.

The whole chain is precisely serializable using rational-coordinate strings, independently replayable, and rejects source promotion through `require_source_atomic_realization()`.

## Validation

`tests/test_marked_conway_geometry.py` verifies sphere mesh topology/orientation and exact punctures, positive/negative/nonzero sectors, preservation of all original source ports, disjoint corridor gap ordering, open terminations, full deterministic coordinates, and adversarial deformation, port mutation, corridor reindexing, and fraudulent geometry status. Exhaustive small sector vectors of lengths 2 through 4 with entries from `{-2,-1,0,1,2}` test reproducibility, crossing counts and absence of overlap. Combined with A5.2c, run on hosted Python 3.11 and 3.12:

```sh
python -m unittest discover -s tests -v
```

## Crucial remaining A5.2 / Issue #33 gates

1. Replace the artificial anchor and the two open terminal ends by the **actual external long return band**, checking its 3D embedding and the global strand components. In Fig.12.8/12.29 the returning band is an essential part of the topology, not optional graphical decoration.
2. Replace/relate these chosen cubical puncture positions to the source **Conway-sphere parametrizations**; verify cyclic orientation and geometry-preserving equivalence under each source reparametrization.
3. Build the actual r=0 closed atomic model (no auxiliary Conway sphere) and independently certify the embedded curve collection and twisted band, including cyclic bond rotation and signs.
4. Embed the **r=1,2,3 ring loops** with exact disjointness/linking controls, mirror cases and sectors, in the source-prescribed location, not as plain braid additions.
5. Derive orientation-**reversing** gluing with the **infinite pure mapping-class data** proven necessary by A5.2b; track framing and V4 effects.
6. Only then build and check **both** R± transformations with relative-boundary ambient-isotopy certificates and source exceptional cases, and approach Theorem 12.21.

A verified local open-chain embedding is meaningful **positive geometric progress**, not proof that the source atomic model or general arborescent isotopy classification is correct. A5.1's `UNVERIFIED_GEOMETRY` gate stays unchanged.
