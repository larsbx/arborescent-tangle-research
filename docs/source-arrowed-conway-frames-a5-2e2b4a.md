# A5.2e2b4a — Figure 12.6/12.7 arrowed Conway frames and signed seam constraints

**Status: ARROW_FRAMING_SOURCE_UNVERIFIED.** This milestone audits *source-required extra marking data* that was absent from every previous PR; it is not a Bonahon–Siebenmann source isotopy proof.

## Primary source — newly checked at the actual page images

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.1, **printed p.141**, Figures **12.6(a),(b)** and **12.7**, and §12.2 **printed pp.142–143**, Figures 12.8 and 12.10 (PDF pages zero-based 152–154).

Primary PDF: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The key source detail was missing from A5.2a–b and A5.2e2b3b:

- An **arrowed equatorial arc**, not just four puncture labels, selects the *slope-zero Conway boundary parametrization*. The arrow from the paper's reference marked sphere corresponds to the standard arc between **(0,0,0) and (1,0,0)** (Figure 12.6). In general the source labels boundary arrows on whichever side is the tangle under consideration.
- For plumbing, Figure 12.7 requires **arrows along distinct adjacent equatorial arcs, emanating from the same puncture** (after allowing rotation/reflection of the picture). A four-point matching alone cannot test that specific condition.
- The source's plumbing joins **oriented** 3D boundary spheres. An orientation-preserving pure boundary mapping class is **not by itself** the orientation-reversing gluing between their induced boundary orientations.
- The band-preserving cyclic permutation statement on printed p.143 is a **tangle isomorphism respecting the global twisted band**, not just equality of rooted cyclic lists. This milestone does not claim that isomorphism.

## Exact typed layer

`ArrowedConwayBoundary(port, arrow_start, arrow_end, cyclic)` records a **directed adjacent arc**, not an unordered 4-point matching, with four standard labels `(P00,P10,P11,P01)`. The canonical zero-slope arrow is `P00→P10`, following Figure 12.6's standard oriented slope-zero arc.

`chosen_cubical_arrow(cube)` constructs a **real polygonal arrow on the existing Q³ cubical Conway sphere**. It maps the existing chosen cube labels `NW,NE,SE,SW` to `P00,P10,P11,P01` and uses the exact curve

```
(w,-2,0) → (w,-3,0) → (e,-3,0) → (e,-2,0).
```

Each open segment is on one *boundary face* (west, y=−3, or east). A chord running straight between the two physical endpoints would go *through the solid cube* and is rejected, as are reversed arrows, wrong endpoints, nonadjacent punctures and counterfeit SOURCE_VERIFIED metadata. This is an explicit **CHOSEN** cubical representative of the *combinatorial* arrow; the geometry has **not been demonstrated to realize the source's full Conway boundary parametrization**, nor are the source 3D equator arcs identified with a torus-cover coordinate automatically.

`SignedConwayBoundaryClass(a,b,c,d,tx,ty)` represents an exact signed integral 2×2 unimodular mapping class lift, with determinant +1 or **−1** and a binary affine 2-torsion shift. It independently derives the four branch-point images from the **full matrix modulo 2** (rather than using only V4 double-transpositions), supports exact group composition/inverse, and acts on primitive unoriented rational slopes. We require determinant −1 in a *proposed* orientation-reversing seam, whereas A5.2b only represented determinant +1 pure sphere self-maps.

`ArrowedConwaySeamObligation(incoming,outgoing,mapping)` checks the **source Figure 12.7 necessary combinatorics**: the image of the incoming arrow must start at the same outgoing marked puncture but head toward the *other* adjacent neighbor, as in two distinct arcs emanating from one point. It stores the full signed matrix and both generator slope images. It cannot claim spatial 3D Conway gluing or ambient isotopy. All statuses explicitly remain `SOURCE_UNVERIFIED`.

## New infinite negative control — arrow + four points still do not fix framing

For every integer n the signed matrices

```
     [0  1 ]
G_n= [1  2n]
```

have determinant −1 and **the same four-puncture map modulo 2**. All satisfy the same Figure 12.7 oriented-adjacent-arrow test, but send the primitive slope (0,1) to (1,2n), so the surface mapping classes are distinguishable. Therefore even adding the source *arrow* and orientation-reversing degree to the four punctures does not select the **pure mapping-class twist** required for a unique plumbing map. Source Fig 12.7 supplies additional geometric/homeomorphism conventions that must be extracted and actually realized; the constructor does not guess them.

A second family `T_n=[[1,2n],[0,1]]` fixes the horizontal primitive (1,0) and all punctures while twisting vertical slopes. It shows that **a slope-zero arrow alone does not eliminate every pure mapping-class ambiguity**. This is not a claim that the original monograph leaves the map ambiguous: the source model may fully determine it; our code does not yet provide that correspondence.

## Tests / next acceptance

`tests/test_source_arrowed_conway_frames.py` audits the chosen arrow on signed-sector 3D band candidate fixtures, rejects arrows passing through the Conway ball or reversing their source orientation, exhaustively tests small exact unimodular determinant ±1 matrices and all affine 2-torsion shifts for inverse/compose and branch-point actions, checks valid/invalid Fig 12.7 arrow compatibility, and exhibits 17 different pure-twist classes with the same endpoint/arrow conditions.

This PR **does not** discharge source Figure 12.10 signed framing or cyclic band-preserving isotopy; it makes these missing geometric choices explicit and refuses to self-certify. The next A5.2e2b4b must derive a **source-specific signed map for every seam** and construct its real 3D representative with verified transverse Conway marking. Only a replayable relative-boundary source isotopy (or equally rigorous independent diagrammatic witness) can promote a separately signed, specific positive certificate.

A5.1 `UNVERIFIED_GEOMETRY` and `OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED` remain unchanged. Issues #49, #46, #43, #41, #39, #33 and #16 remain open.

Run `python -m unittest discover -s tests -v` with hosted Python 3.11 and 3.12 and post-merge main.
