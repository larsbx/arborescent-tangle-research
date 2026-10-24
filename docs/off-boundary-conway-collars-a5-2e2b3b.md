# A5.2e2b3b — Off-boundary exact PL Conway collars and full surface collision certification

**Status: OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED.** This is a *positive geometric construction in an explicitly chosen coordinate frame*, not a Bonahon–Siebenmann Figure 12.10 source-isotopy proof. In particular, no source 4-punctured Conway parametrization or pure mapping-class framing has yet been verified.

## Why the previous roof collar was insufficient

A5.2e2b2 (#45) built a physically connected, coordinate-welded PL ribbon surface but attached its patches via strips lying across **two-dimensional regions of the Conway cube surface**. That provided useful mesh construction evidence but not an embedded band in the *interior of the knot complement* meeting the Conway spheres only in controlled attachment arcs.

Bonahon–Siebenmann §12.2, **Figs 12.8 and 12.10, printed p.142**, describes the actual twisted band and its plumbing patches. §12.5, **Fig.12.29, printed p.158**, adds r ring components. The source also states (p.143) that cyclic permutations of sector weights are realized by band-preserving tangle isomorphisms. **Neither these source conclusions nor their specific maps have been reproduced by this chosen-coordinate construction.**

Primary PDF: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

## New physical collar geometry in exact Q³

For a Conway cube `C=[w,e]×[-3,3]×[-3,3]` with east/west faces `x=e,w`, construct a four-station ribbon profile in x–z:

```
(w,0) → (w−1,4) → (e+1,4) → (e,0).
```

At every station create width vertices with `y∈{−2,0,2}`. The exact 12 vertices and 12 triangles form three continuous 2-rectangle portions, with 2 triangles per half-width rectangle. The west/east end caps contain exactly the three coordinates of the existing inbound/outbound ribbon cap edges; those are genuinely **shared Q³ mesh edges**, not names identifying distant abstract vertices.

For any point in the **relative interior** of the western slope, `x<w` and `z=4(w−x)>0`; therefore it lies outside closed C. On the eastern slope, `x>e` and `z=4(x−e)>0`; on the central roof, `z=4>3`. The collar thus has no contact with its own solid cube except the two specified end-cap chords `x=w,e; z=0; -2≤y≤2`.

Other Conway cubes lie in disjoint longitudinal slabs. The nearest neighbor's face may be only **one unit away** for zero-weight sectors, but the collar reaches that neighboring slab only at `z=4`, strictly above the cube's `z≤3` roof. The exact `_clip_cube` algorithm is independently applied to **every triangle against every solid Conway cube**, and permits only points on the designated source-candidate end-cap chords.

## Independent *two-part* embedding verification

- `verify_spatial_surface(triangles)` (A5.2e2b2) independently welds exact coordinate vertices, checks triangle incidence, one path/cycle per vertex link, connectedness, Euler characteristic, orientability, and boundary-circle count. In the chosen model, the joined surface is **one connected χ=0 PL band**.
- `audit_embedded_triangle_mesh(triangles)` (A5.2e2b3a PR #47) independently checks **every pair of actual triangular faces** for geometric intersection, using exact rational coplanar polygon clipping or noncoplanar triangle-plane segment overlap. Only an identical common **declared endpoint** or entire common **declared edge** may meet; transverse interior crossings, same-plane overlapping triangles and unintended tangencies fail closed.

Together with the exact solid-cube clipping, the checks establish the *chosen model* is a physically connected, actually embedded PL surface in Q³ with all its open interiors outside the removed Conway balls. These are strictly **finite geometrical facts**, not source-isotopy statements.

The construction preserves the exact signed sector vector, genuine A5.1 port addresses and all 4n marked strand endpoints from previous stages. Tests exercise r=0 signed-sector fixtures, even/odd twist parity, physical end-cap equality, complete surface intersection audits and adversarial mesh or status mutation. The status remains `OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED`.

## Why the source-isotopy obligation is still open

Even a correct PL embedding with the expected abstract annulus/Möbius parity does **not** identify the precise Bonahon–Siebenmann twisted band and its **framed plumbing patches**, Conway boundary markings and *pure mapping-class twists*. Those source conventions are needed to prove the cyclic band-preserving rotation assertion and both modified R± ring-moves as relative-boundary isotopies. The present construction is **only a chosen geometric candidate**, with positive validation limited to its own coordinates.

Outstanding #46 items:
1. Figure-derived parametrizations including signed framing and all Conway pure twists; positive geometric comparison to the exact source.
2. Relative-boundary PL-isotopy/certified movie between this candidate and the source's global twisted band, including cyclic rotation and mirrors.
3. r=1,2,3 source ring embeddings with independent linking and then both R±.
4. Source exception regimes (including valence one) and complete §12.7 canonicality.

Do **not** close #46/#43/#41/#39/#33/#16 or promote A5.1 beyond `UNVERIFIED_GEOMETRY` from this exact-model accomplishment.

## Tests

```sh
python -m unittest discover -s tests -v
```

Require hosted Python 3.11 and 3.12 passing on the exact PR head and post-merge main.
