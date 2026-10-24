# A5.2e2a — exact 2D ribbon patches on the r=0 3D return strands

**Progress toward #41; NOT its complete source-realization acceptance.**

This stage produces n genuine triangulated two-dimensional PL strip patches, one for each paired set of outer-return/interior strands constructed and globally collision-checked in A5.2e1 (#40). It is an improvement over the previous one-dimensional centerline candidate. Its scope is **RIBBON_PATCHES_SOURCE_UNVERIFIED**: it is not yet the Bonahon–Siebenmann source's global, correctly framed twisted band, and it certifies no source isotopy, Conway gluing, R± move or ringed atom.

Primary source context (previous page-image review): Bonahon–Siebenmann §12.2, printed p.142 Fig.12.8 and §12.5, printed p.158 Fig.12.29. PDF: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf . These figures motivated the existing return-strand candidate but *do not prove* that its arbitrary cubic source frames or global surface match the isotopy class of the diagram.

## Explicit mesh construction

For an independently embedded pair of source-candidate polygonal rails `A[0..m−1]` and `B[0..m−1]` (exact Q³ coordinates), build middle vertices `C[j]=(A[j]+B[j])/2` and replace two `C` vertices at the long outer arch's two crest stations (indices 2 and 3) by the *same x,y coordinates* and z=8. There are precisely **3m stored vertices** and **4(m−1) stored oriented triangles**; oriented triangles form the two halves of an actual width-carrying PL strip, rather than an abstract strand pairing. At the two outer U-turns, triangle diagonals are chosen individually to avoid genuine collinear/zero-area faces.

For each interior corridor patch, all sections have strictly increasing longitudinal x and distinct nonvanishing cross-section width vectors. Since adjacent intervals have disjoint open x slabs, and at each longitudinal x the cross-section is an embedded line segment, the entire sector patch is geometrically embedded.

For the outer return patch, the two rails extend beyond the leftmost and rightmost Conway cubes, while the long arch runs along y=±6. Its crest ridge has y=0,z=8. Each arch triangle is constrained to one signed-y half and its exact plane satisfies

`z = 8 − (4/3)|y|`.

Consequently, wherever `|y|≤2`, the arch has `z≥16/3`. Every interior braid patch lies in `|y|≤2, |z|≤1`, so this arch cannot intersect the interior patches. The remaining exterior triangles lie strictly before the first interior braid corridor or after the last. Interior corridors occupy disjoint open x slabs. These geometric separation inequalities are independently recomputed from the **stored mesh vertices**, rather than merely presumed from an A5.1 tree or port metadata.

## Independent surface checks

`TriangulatedRibbon.audit` checks every supplied triangle index, exact nonzero cross-product/area, every directed interior-edge cancellation, every interior edge appearing twice, connectedness, Euler characteristic `V−E+F=1`, and the actual one-cycle mesh boundary equal to **both already verified A5.2e1 physical strand polylines, plus the two end cap chords**. It rejects fake holes, face inversion, duplicate vertices, missing faces and collapsed cross sections.

`_clip_cube` independently uses exact rational Sutherland–Hodgman clipping of every triangle against all **six solid ball-box inequalities**. No triangle may penetrate *any* Conway cube: if a clipped intersection is nonempty, every resulting point must lie on that ribbon's explicitly checked end-cap edges, which attach to punctures on sphere face interiors. This is more rigorous than mesh-edge-only intersection checks. Source spheres remain chosen cubical parametrizations, **not yet proven equivalent to the source's Conway markings**.

`ZeroRingRibbonPatchSystem.audit` separately verifies alignment with independently embedded return strands, all signed sector metadata and ring_count=0, checks both ribbon/sphere disjointness and the prescribed exact arch clearance. Geometry and cyclic sector data are serialized without floats.

## Theorem (limited positive claim)

For a supported A5.1 atomic instruction with `ring_count=0` and `valence=n≥2`, the constructor returns n orientable PL disk **patches** with the exact A5.2e1 physical boundary strands and prescribed cap chords. Each disk has Euler characteristic 1, an explicit oriented triangulation, and nondegenerate triangles. The separate geometry audit certifies absence of intersection with the interiors of all n closed cubical Conway balls, and pairwise disjointness of the different patch interiors in the chosen exact coordinate model.

It does **not** prove these are the particular global embedded twisted surface and framing of the source diagrams. A surface with the same strand endpoints and even identical outer rail geometries may represent a different framing if an extra surface twist is introduced. The current canonical ridge selection rejects such unsupported variations rather than declaring them isotopic.

## Tests

`tests/test_pl_ribbon_patches.py` includes cases for n=2,3,4, signed and zero sectors, exhaustive `{-2,-1,0,1,2}^n` for n=2,3, exact triangle-vs-solid-cube clipping, deterministic canonical meshes, and adversarial mutation of triangles, face orientation, ridge heights and source status. Tested with the existing A5.2e1 independent strand intersection controls:

```sh
python -m unittest discover -s tests -v
```

Run on both hosted Python 3.11 and 3.12 before any merge.

## Remaining acceptance for issue #41

The **global supporting twisted band**, source-verified signed twist/framing identification, cyclic rotations, mirrors, valence-one special model, and independent Fig.12.8 diagram/witness comparison remain unsolved. The n independently embedded ribbons may need attachment or identification along Conway sphere boundary components to yield the exact source band. Establishing the source prescribed gluing and its mapping-class twist is not a consequence of these patches.

The r=1,2,3 loops of Fig.12.29, source-compatible orientation-reversing Conway sewing, both modified R± isotopy certificates, and Theorem 12.21 classification remain open under #33 and #16. Preserve A5.1 `UNVERIFIED_GEOMETRY` and do not close #41 on the basis of this PR.
