# A5.2e2b4b — Distinguish the proper tangle K from the entire band boundary

**Status: MARKED_KNOT_SOURCE_UNVERIFIED.** A positive geometric consistency test for the *chosen* candidate, not a source-isotopy certificate.

## New source interpretation that prevents a serious certification error

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.2 **printed p.142**, between Figures 12.9–12.10, expressly calls Figure 12.8's shaded surface **the surface with boundary containing the knot** (the *twisted band*). The authors do NOT say **the entire boundary of the twisted band equals the knot**.

Their arborescent tangle (M,K) is a knot/tangle pair, with every Conway sphere a boundary component of M and **four transverse marked punctures of K** on that sphere (see also §12.1 printed p.141, Figures 12.6–12.7, and the Introduction's Figure 0.2). Thus for an atomic n-free-bond tangle in its exterior 3-manifold with n boundary Conway spheres, the appropriate **proper knot marking** has 4n points at ∂M, grouped into **2n embedded knot arcs** with those endpoints. That is not the same object as the *entire* boundary of a completed auxiliary band surface.

Primary source: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

### Why this matters for previous PRs

A5.2e1 (#40) already constructed **2n physically embedded knot arcs** with 4n labelled Conway endpoints, independently checking all arc/solid-ball intersections.

A5.2e2a/b2/b3b then constructed physically completed **annulus/Möbius band surfaces** by adding collar faces. Those entire band surfaces have boundary consisting of **closed loops** (2 for orientable annulus or 1 for nonorientable Möbius band). It would be a **category error** to interpret these closed band-boundary loops themselves as the proper source-tangle knot K: their degree-two boundary graphs have *zero* one-valent ends.

The source statement is `K ⊂ ∂B` **as a separately distinguished, properly embedded 1-manifold**, not the assertion `K = ∂B`. In the chosen model, the additional collar-rim segments are part of `∂B`, but **are not marked K strands**. The new code makes this marking explicit and auditable.

## Independent exact PL marked-knot checker

`verify_distinguished_knot(triangles,cubes,knot_arcs)` reads only *actual rational 3D triangle coordinates*, *actual cubical Conway boundaries* and *candidate knot arc coordinates*. It does not trust source sector weights or a caller's declared topology.

1. Independently reconstruct every mesh edge and which edges occur in exactly one face, i.e. the **complete physical band boundary ∂B**. Verify ∂B is a union of 2-regular circuits, and compute its number of boundary loops.
2. Independently reuse A5.2e1's all-segment 3D arc disjointness and solid-ball intersection checker to confirm **exactly 2n knot arcs**, each with two distinct marked sphere endpoints, **all 4n physical punctures used once**.
3. Verify every single physical K arc edge is an **exposed real boundary mesh edge**, not a made-up chord between punctures, not a hidden interior edge of the band, and not repeated in any other K arc. Thus the assertion `K⊂∂B` is checked on actual geometry.
4. Independently rebuild the graph **K only**: its degree-1 vertex set must be *exactly* the 4n marked sphere punctures, with every other vertex having degree 2. Its connected component count must be **2n proper arcs**, not the band's whole closed boundary circuits.
5. At every physical Conway sphere puncture, require the K incident segment to approach the west/east cube face with a **nonzero outward longitudinal component**. This is a local PL boundary-transversality test, not a guessed puncture permutation.
6. Reject any attempt to assert that all band-boundary edges belong to K, or to promote the returned audit to SOURCE_GEOMETRY_VERIFIED. The extra collar/rim band-boundary edges are retained and explicitly counted as `complementary_band_edges`.

For supported r=0 configurations with `n≥2`, the independent checker derives

```
  |∂K| = 4n        proper Conway sphere punctures
  components(K) = 2n   physically embedded arcs
  components(∂B) = 2 or 1, depending on candidate twist parity
  K ⊊ ∂B          (strict subgraph, not equal)
```

This does not change the previously stored `OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED` band status or certify that the current band framing agrees with the source Figure 12.10 model.

## Integration and adversarial tests

`audit_source_marking_on_candidate(candidate)` integrates this independent checker with the existing `ExteriorCollaredBandCandidate` after validating that chosen embedded surface.

`tests/test_distinguished_knot_on_band.py` checks signs and positive/negative/zero sectors, knot-vs-band component counts, exact 4n puncture/2n arc properness, nonzero longitudinal transverse endpoints, fake chords through solid Conway balls, duplicated mesh faces, erroneous fabricated status, and **stable marked K data under auxiliary collar-surface refactoring** (comparing A5.2e2b2 roof surfaces and A5.2e2b3b exterior surfaces).

This correctly prevents an otherwise plausible source misclassification even if an auxiliary surface has χ=0 and passes all triangle/triangle geometry checks.

## Remaining #49 / #33 work

The displayed band *is a geometrically embedded, source-unverified candidate*, and the properness of its marked K arcs is an additional **necessary** consistency check. Still missing:

- The actual Figure 12.10 source plumbing-patch homeomorphisms, orientation and arc framing.
- The source Figure 12.6/12.7 arrow and pure Conway boundary mapping class realization in Q³ (A5.2e2b4a adds the necessary algebraic marking checks but not realization).
- A constructive **band-preserving cyclic sector rotation** isotopy relative to the source 4-punctured Conway boundary markings.
- Special n=1 source cases, r>0 rings, both actual R± moves and full arborescent isotopy classification.

Until then, **NO** `VERIFIED_TOPOLOGY`, source ambient-isotopy promotion, or real isotopy-class asymptotics. Issues #49/#46/#43/#41/#39/#33/#16 stay open.

Validation: `python -m unittest discover -s tests -v` on hosted Python 3.11/3.12 and post-merge main.
