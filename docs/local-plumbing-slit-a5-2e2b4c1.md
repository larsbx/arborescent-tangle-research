# A5.2e2b4c1 — Actual local ambient-3D PL slit, with inverse reglue

**Positive scope:** a *local* three-dimensional triangulated PL product-ball can be cut along a physically embedded, doubly-arrowed two-triangle patch, and its two side-labelled boundary copies can be exactly rejoined. This is a real 3D cell construction, not a tree rewrite or a surface-edge deletion.

**Status:** `LOCAL_PRODUCT_3BALL_ONLY`, `LOCAL_SLIT_REGLUE_ONLY`, and the existing `UNVERIFIED_GEOMETRY`. **Nothing in this patch certifies the Bonahon–Siebenmann full Conway manifold or a source isotopy.**

## Primary source image and mathematical distinction

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.1, printed **p.140**, Figures **12.4–12.5**, inspected as PDF page images; §12.1 p.141 Figures 12.6–12.7 supply the subsequently required equatorial Conway arrow and plumbing conditions. Figure 12.5(a) marks **two transverse arrows** on a rectangular surface-plumbing patch, and Figure 12.5(b) illustrates the corresponding Conway-sphere parametrization. Section 12.1 explains that the converse construction involves **slitting open the ambient 3-manifold along the interiors of the plumbing patches** and parametrizing the resulting tangle boundaries.

Primary PDF: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The preceding A5.2e2b4a / A5.2e2b4b work added one equatorial arrow to the four-punctured Conway frame and separately identified the properly embedded **K ⊂ ∂B** within the completed ribbon-band boundary. Those results did **not** implement a cut of the ambient 3-manifold. This stage implements the **local product-ball** operation and explicitly keeps it separate from a global S³ slit and Conway 4-puncture reconstruction.

## Concrete oriented surface patch in the existing geometric candidate

A5.2e2b3b's actual Q³ exterior collar contains a planar, positive-area two-triangle roof rectangle, with exact ordered corners:

```
A=(w−¼,-2,4), B=(e+¼,-2,4),
C=(e+¼,0,4), D=(w−¼,0,4).
```

`collar_plumbing_chart(collar)` **verifies both real triangles (A,B,C) and (A,C,D) appear in the stored physical band mesh**. It places two oriented arrows at the common rectangle center:

```
O=((A.x+C.x)/2, -1, 4)
u: O → O+((B.x−A.x)/4,0,0)
v: O → O+(0,1/2,0).
```

The exact vectors satisfy `u·v=0` and `(u×v).z>0`. This is a **chosen local analogue** of Figure 12.5(a)'s two transverse surface arrows. It has not been identified as one of the source's actual framed plumbing patches. Reversed arrows, nonadjacent/degenerate corners or fabricated source-verification status fail.

## A physical PL 3-ball and its oriented midplane disk

For rational `h>0`, construct the embedded 3D product

```
P×[-h,h] = [A.x,B.x] × [A.y,D.y] × [4−h,4+h].
```

There are three layers, each with the same four rectangular vertices at `z=4−h`, `z=4`, and `z=4+h`: **12 exact Q³ vertices**. The rectangle has a diagonal A→C, so each layer-pair consists of two triangular prisms. Each triangular prism admits the canonical **three-tetrahedron staircase subdivision**, giving **12 tetrahedra** for the entire local product.

`TetrahedralProductNeighborhood.audit` independently checks every stored vertex, every tetrahedron's exact positive oriented volume, all internal face multiplicities and reversed orientations, the two midplane triangular faces occurring exactly twice (one tetrahedron on each side), and the closed orientable boundary having Euler characteristic 2. It additionally recomputes the **total exact tetrahedral volume** and compares it to the geometric rectangle area times `2h`:

```
Σ |det(v₁−v₀,v₂−v₀,v₃−v₀)|/6
    = (B.x−A.x)(D.y−A.y) · 2h.
```

A valid *local* three-dimensional triangulation therefore exists independently of any planar diagram or formal plumbing instruction. After cutting, an **additional independently reconstructed check** examines the actual tetrahedral faces of EACH of the two halves: exactly six tetrahedra per side, twelve boundary triangles, eight boundary vertices, every boundary edge with two incident boundary triangles, and boundary Euler characteristic two. It checks the two real cut triangles appear on the respective 3-ball boundaries; a missing tetrahedron or forged cut face fails this separate check.

## Slitting and exact reglue

`slit_product` assigns tetrahedra to the upper/lower spatial halves by their actual layer memberships. Four midplane vertices at `z=4` acquire a distinct combinatorial label on the LOWER side and on the UPPER side. The two **physically coincident** cut disk copies each have two triangles, but their **outward boundary orientations are opposite**: lower's cut normal points +z; upper's cut normal points −z.

`LocalSlitWitness.verify` independently reconstructs the canonical side-label assignments and both exact cut copies from the original stored 3D cell table, refuses crossed or forged vertex labels, compares corresponding Q³ coordinates, and removes the side distinctions only as a **replayable, explicitly specified inverse regluing operation**. The reconstructed tetrahedral incidence list must exactly equal the original.

This is a *local positive combinatorial PL slit/glue certificate on a geometrically embedded tetrahedral 3-ball*. It **does not** mean that the current surface's global embedding has been slit, that the entire 3-manifold complement was triangulated, that new Conway boundary spheres have four knot punctures, or that a relative-boundary isotopy to the source is known.

## Tests and adversarial cases

`tests/test_local_plumbing_slit.py` exercises source-labelled real collar fixtures from signed/zero sector candidates, patch-on-mesh membership, orthogonal arrows, positive rational neighborhood half-thicknesses `1/8,1/16,1/32`, exact 3D volume and Euler checks, distinct two-sided cut faces and inverse reconstruction, wrong arrows, reversed corner order, mutated tetrahedral coordinates and incidences, forged cut-side labels and `SOURCE_VERIFIED` rejection.

Run:

```sh
python -m unittest discover -s tests -v
```

Require GitHub-hosted Python 3.11/3.12 passing on **the exact PR head** and post-merge main.

## Remaining full issue #52 gates

This milestone **does not** discharge the complete geometric slitting or inverse theorem of Bonahon–Siebenmann §12.1. The next implementation must place all source-identified oriented plumbing disks in one **actual triangulated ambient 3-manifold**, cut their interiors, track new boundary 2-spheres and **exactly four transverse K punctures on each**, derive source-compatible equatorial arrows and full signed/pure Conway parametrizations, and replay inverse **geometric** plumbing. A small product 3-ball split into two halves is not already a 4-punctured Conway boundary construction.

Further tasks #49/#46/#33/#16 require source Fig12.10 band-preserving cyclic isotopy certificates, r>0 ring embeddings, and both R± relative-boundary witnesses. **Absent certificates remain UNVERIFIED, never nonisotopy.**
