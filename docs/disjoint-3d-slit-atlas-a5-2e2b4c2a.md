# A5.2e2b4c2a — Simultaneous local 3D slits on every candidate Conway site

**Verified intermediate object:** a physically disjoint atlas of n local PL product 3-balls, each containing one *chosen* two-arrow oriented band patch, with an independently replayed local slit and inverse reglue. This is explicitly **NOT** the global ambient S³ triangulation or the Bonahon–Siebenmann Conway-tangle reconstruction promised in Issue #54.

## Source motivation and necessary distinctions

Bonahon–Siebenmann §12.1 printed p.140, Figures 12.4–12.5, discusses recovering boundary-marked Conway tangles by **slitting the ambient manifold along plumbing patch interiors**. The figure supplies two transverse oriented patch arrows. Printed p.141 Figures 12.6–12.7 then fixes an oriented equatorial arrow and describes the plumbing requirement. Our previous PR #53 verified a **single genuine local tetrahedral product-ball cut** on a two-face roof rectangle of the physically embedded r=0 auxiliary band. PR #53 made no claim to a globally connected ambient manifold or to producing any source boundary sphere.

Source: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

Here we implement the next concrete *geometric staging* requirement: put all supported local patch charts on the **same existing global candidate band**, verify that their entire **three-dimensional product neighborhoods are pairwise disjoint**, and independently replay each cut and reglue without reusing any geometric boundary port. Merely collecting cuts is not the source global slit: that requires a connected ambient triangulation and its actual 4-punctured boundary reconstruction.

## Exact construction and positive finite theorem

Given a supported A5.1 `AtomicInstruction` with `r=0` and `n≥2`, first use the existing A5.2e2b3b embedding to obtain n real source-port-labelled cubical Conway boundaries, n exterior-collared ribbon patches, and the distinguished 2n proper knot arcs. Let the chosen collar roof for boundary i be

```text
  x in [wᵢ−¼, eᵢ+¼],    y in [−2,0],    z=4.
```

The two-arrow source-inspired face chart on that rectangle is independently checked to be **two actual oriented faces** of the stored band mesh. For an **exact positive rational** `h∈(0,¼]`, let its local product 3-ball be

```text
  Pᵢ × [4−h,4+h].
```

Each local product neighborhood has **12 vertices, 12 positive-volume tetrahedra** and a verified closed Euler-2 boundary. Its exact volume is

```
  vol(Pᵢ × [4−h,4+h]) = (9/2) · 2 · (2h) = 18h.
```

The atlas's independent audit recomputes all actual 3D bounding boxes from the **stored tetrahedral vertices**, not from port IDs. Adjacent cubes' longitudinal gaps equal `1+2|aᵢ|` for the respective intervening source sector. Because their roof overhangs are only ¼ each, their product neighborhoods satisfy

```
  gap(Pᵢ,Pᵢ₊₁) = ½ + 2|aᵢ₊₁| ≥ ½ > 0.
```

Thus no two local 3D product neighborhoods meet, even in a zero-weight corridor. All have `z≥4−h≥15/4>3`, so the entire local product neighborhoods lie **strictly above every removed Conway cube** (whose z≤3). The distinguished knot K has physical vertices in `−1≤z≤1`, independently read from its actual Q³ arc table. Since every straight knot edge lies in its endpoint z-range, **K avoids every local product 3-ball and all their cut disks**.

For n candidate sites, the exact results are:

```
  local pre-slit 3-ball components: n
  local tetrahedra: 12n
  local post-slit 3-ball components: 2n
  tetrahedra after slitting: 12n
  two oriented cut-disk copies per site, with independent inverse
  exact total local product volume: 18nh
  punctures of K on these local cut disks: 0
  produced source 4-punctured Conway boundaries: 0
```

The atlas binds **every** source candidate `Port(atom.address,i)` to the corresponding real physical product geometry and typed `LocalSlitWitness`. Every local cut is checked by PR #53's independent face/volume and **both-side Euler-2 boundary** verifiers. Tampering with any port, chart, half-thickness, tetrahedron, cut side, or site order fails closed. `DisjointLocalSlitAtlas.canonical()` serializes complete exact physical coordinates and per-side cut incidence for replay.

### Research consequence: a concrete global-reconstruction obstruction

**The two locally cut product 3-ball boundary spheres are *not* the source Conway spheres.** Their cut-disk boundaries contain **zero** knot piercings in this chosen model, whereas source Conway boundaries must carry **four** transverse knot punctures apiece. Also, the atlas is a disjoint union of n 3D neighborhoods, not one ambient S³ or proper complement manifold.

Therefore a procedure that declares source-compatible 4-punctured Conway boundaries merely by applying n local product slits and renaming the 2n resulting 3-ball boundaries is **mathematically invalid**. The code enforces this distinction: `genuine_conway_boundary_spheres=0`, `knot_piercings_of_cut_disks=0`, and `DISJOINT_LOCAL_ATLAS_NOT_GLOBAL_MANIFOLD`.

This is **not** a proof of failure of the original source theorem; it is a precise failure of the *current incomplete reconstruction strategy*. The actual source construction requires the **whole ambient manifold and explicitly positioned plumbing cuts**, including correct K transversality, rather than a roof patch chosen for coding convenience.

## Acceptance evidence

`arborescent/disjoint_local_slit_atlas.py` implements the typed atlas and independent geometric checks. `tests/test_disjoint_local_slit_atlas.py` covers signed/mixed and zero sector vectors; n=2 and n=3; exact h=1/16,1/8,1/4; correct minimum gap and total volumes; exact candidate port-to-chart correspondence and deterministic canonical replay; adversarial site reordering, missing site, mismatched thickness, forged source status, incorrectly retagged cuts and r>0/valence-one rejection. A separate fixture confirms the physical proper-K z-range and **zero** source piercing for both cut copies.

Run the full test matrix on GitHub-hosted Python 3.11/3.12, and post-merge main.

## Next actual Issue #54 gate

The global source equivalence remains open. The next genuine geometric object must be a **single, oriented, tetrahedrally triangulated ambient manifold** containing both the actual candidate band and distinguished K as properly positioned embedded subcomplexes, with **source-identified Figure 12.10 plumbing patches**. A checked global cut must derive **new Conway sphere components**, each independently intersected in **four actual punctures** by K, and reconstruct the source-specific oriented equatorial arrow and full signed/pure mapping class. Inverse reglue needs a positive ambient cell/isotopy witness, not a local relabelling or a count.

Neither this atlas, the A5.1 plumbing plans nor the A5.2e2b3b band candidate may be promoted to source isotopy; `UNVERIFIED_GEOMETRY` remains the mandatory public gate. No ring relocation R± or Theorem 12.21 classification follows.
