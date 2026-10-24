# A5.2e2b2 — Physical Conway-boundary collars and a connected PL surface

**Positive scope:** an exactly specified, physically connected, 2-dimensional PL ribbon surface in a **chosen cubical Q³ embedding**, obtained from A5.2e2a's n separated ribbon disks. **Status: COLLARED_BAND_SOURCE_UNVERIFIED** (code literal `COLLARED_BAND_SOURCE_UNVERIFIED`).

**What this does not prove:** a source-canonical Bonahon–Siebenmann band, a correct boundary-marked Conway plumbing homeomorphism, ambient isotopy to Figure 12.10, either R± move, or any r>0 ring embedding.

## Primary-source distinction

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, printed p.142 Figure 12.8 (band), Figure 12.10 (plumbing patches), p.143 (band-preserving cyclic permutation), and p.158 Figure 12.29 (ringed model), inspected as primary PDF page images. The source defines a globally connected twisted band with n labelled sector half twists and n actual 4-punctured Conway boundaries. Nothing in the paper says our convenient **cubical roof-collar framing** is the specific source's boundary-parametrized band.

Primary PDF: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

A5.2e2b1 showed that abstract attachments cannot be guessed from the strand endpoints: the same patch may yield an annulus or Möbius band after an endpoint-width identification. The new A5.2e2b2 **does more than choose an abstract identification**: it creates actual polygonal surface faces **in Q³** that share their cap edges and vertices *at identical physical coordinates* with the already independently audited A5.2e2a ribbon disks.

## Exact construction

For every source Conway boundary cube with west face x=w and east face x=e and |y|,|z|≤3, create a three-rectangle U-shaped collar. Its longitudinal profile in the x–z plane is

```
(w,0) → (w,3) → (e,3) → (e,0).
```

At each profile station, use three cross-sectional y positions (-2,0,+2). Every profile vertex is an exact rational Q³ point. Triangulate each of the six rectangular half-width cells into two oriented triangles, giving **12 vertices and 12 nondegenerate triangular faces per collar**.

Every collar triangle lies on **one physical face of its cube**, specifically west, roof, or east. Therefore it does not enter the *interior* of the removed Conway ball. The west and east end caps are exactly

```
west: [(w,-2,0),(w,0,0),(w,+2,0)]
east: [(e,-2,0),(e,0,0),(e,+2,0)].
```

The A5.2e2a strips terminate at the same **three actual coordinates**; each of the two cap segments is shared by exactly one original strip face and exactly one new collar face.

The input disks are ordered:
- Patch 0 is the outer return strip (last cube east to first cube west).
- Patch i (i≥1) connects cube i−1 east to cube i west.

For cube i, physically attach the incoming patch's **end** cap at west and the outgoing patch's **start** cap at east. For cube 0 the incoming patch is patch 0; for the last cube the outgoing patch is patch 0.

In contrast to the earlier A5.2e2b1 abstract cap quotient, all these cap edges are actual shared geometric segments in Q³. There is no artificial second vertex-copy identification, nor an unimplemented cap-to-cap jump.

## Independent global mesh audit

`verify_spatial_surface(triangles)` is a standalone **coordinate-only** validator. Given nothing but an arbitrary collection of explicitly stored exact Point3 triangles, it:

1. Rejects non-typed, zero-area and repeated triangles; welds vertices by exact **3D coordinate equality**, not by source port names.
2. Reconstructs *all* triangle edges and requires incident-face multiplicity one (true global boundary) or two (internal seam).
3. Builds every triangle-corner link independently and requires **one** interval or cycle at each actual coordinate vertex. This rejects pinched or branched global surfaces.
4. Traces the degree-two boundary graph and counts actual boundary components. Computes connected components by face adjacency, Euler characteristic directly from welded vertices/edges/faces, and orientability by propagating signed edge-direction constraints.

The fixed construction separately checks geometric clearance: all collars lie **on disjoint x-slabs on the sphere boundary**, cannot penetrate their solid Conway cubes, and have no interior shared with each other. A5.2e2a's independently checked original patches touch cube boundaries only on their explicit cap chords. These two audited geometries therefore meet only on the declared matched end caps. This is a restricted construction proof; **the generic mesh verifier does not decide arbitrary triangle/triangle geometric intersection**.

The assembled physical Q³ surface is connected and has Euler characteristic 0. The checker determines, rather than assumes, whether it is an orientable annulus (b=2) or a nonorientable Möbius band (b=1). Both remain *candidate surfaces* subject to the source framing test.

## What requires future source verification

The roof collar **lies on the chosen cubical Conway boundary surface over an extended 2D region**, not merely at four transverse strand intersections. It provides an explicit and auditable *attachment scaffold* but **has not been isotoped off the Conway sphere** into a source-canonical collar or verified to preserve the source's Conway parametrization. That boundary transversality and pure mapping-class choice are essential for Figure 12.10.

The next A5.2e2b3 must: derive the actual Fig.12.10 *framed plumbing-patch attachments*; construct exact **off-boundary** ambient collars intersecting each Conway sphere in only the intended controlled arcs/patches; independently verify all PL triangle/triangle disjointness and crossing sign conventions; and exhibit replayable **band-preserving cyclic permutation** witnesses. Reject unsupported source interpretations. Then proceed to the r>0 rings of Figure 12.29 and independently verified R±.

The A5.1 `UNVERIFIED_GEOMETRY` gate stays mandatory and issues #43/#41/#39/#33/#16 remain open.

Run `python -m unittest discover -s tests -v` under hosted Python 3.11 and 3.12; inspect exact-head checks and post-merge main.
