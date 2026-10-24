# A5.2c — exact spatial embeddings for signed sector half-twist boxes

**Scope: verified local two-strand geometry, not the Figure 12.29 atomic tangle.**

## Image-checked source anchors

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.2 printed p.142, **Fig.12.8** draws atomic tangles as a band with signed integral half-twist regions `a_1,...,a_n`, including both crossing signs. §12.5 printed p.158, **Figs 12.28–12.29** repeats those twist boxes with a distinct left-hand chain of `r` rings. That figure also contains two long return strands, punctured gluing surfaces, and a surrounding band. It is **not** a concatenated two-strand braid, and its rings must not be simulated by adding extra braid crossings.

Primary source: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

## Geometric object constructed

`spatial_two_braid(Braid2Diagram(word))` outputs two explicitly labelled, disjoint polygonal arcs in exact `Q^3`, with fixed NW, NE, SW and SE endpoints. Coordinates: x horizontal, y increasing south, z height. For signed generator `s=±1` in cell `i`, the **left** incoming strand runs

`(-2,2i,0) -> (0,2i+1,s) -> (2,2i+2,0)`,

and the **right** strand runs

`(2,2i,0) -> (0,2i+1,-s) -> (-2,2i+2,0)`.

The projection has exactly one crossing per generator at `(0,2i+1)`; heights `+1` and `-1` separate the strands there in three dimensions. Between successive crossing levels, their projected x coordinates differ by linearity, so there are no other projected crossings. The terminal straight sections preserve the framed endpoints. The empty word yields two vertical disjoint arcs.

`SpatialTwoBraid.audit` independently checks every PL coordinate, over/under height, track swap and endpoint pairing. It accepts only this restricted geometric normal form, not arbitrary 3D arcs. `crossing_records` reads the height order from the spatial embedding, and `canonical_coordinates` serializes its **exact rational** vertices without floats.

`realize_unringed_sector_blocks(AtomicInstruction)` embeds the independent `a_i` half-twist boxes in the atom's preserved cyclic sector order. It is a component of a future source realization; it does **not** embed the whole atom, even when `r=0`. It rejects `ring_count>0` and valence-one closed atoms explicitly rather than silently deleting geometric information.

`verify_spatial_rii_certificate` independently rebuilds and audits the actual 3D arcs at **every step** of an A5.2a Reidemeister-II edit chain, checking labelled endpoints, signed writhe and local crossing metadata. The standard RII move gives a relative-boundary two-braid isotopy, but this is **not** an algorithm for general Conway tangles.

The `VERIFIED_LOCAL_PL_BRAID` status applies strictly to those local blocks. Any request to promote them to a source §12.5 atomic realization is rejected with `UNVERIFIED_GEOMETRY`, and the existing A5.1 `AtomicPlumbingPlan` still cannot certify anything topological.

## Evidence and next obligations

`tests/test_spatial_twist_blocks.py` exhaustively checks all two-braid words of length at most six, signs, endpoint matchings, crossing heights, zero cases, positive/negative angular sectors, valid RII replays and invalid edits. Adversarial fixtures corrupt both crossings and external endpoints and request fake source approval; all must fail closed. Run `python -m unittest discover -s tests -v` in Python 3.11 and 3.12.

The actual A5.2 (#33) still requires: (1) closing the outer bands and explicit four-punctured **3D Conway boundary spheres** in the correct Fig.12.29 cyclic order; (2) embedding **r=0,1,2,3 rings** including signed/mirror conventions and checkable linking; (3) orientation-reversing gluing including the infinite pure mapping-class information missing from V4, as proven in A5.2b; (4) independently checked ambient isotopy witnesses for **both** R± source orientations, with inverse, nonzero-sector, multiple- and adjacent-ring exception cases; and (5) complete source W/S/P/N canonicality before any Theorem 12.21 classification claim.

**This patch is genuine small 3D geometry, not a full Figure 12.29 realization and not a source ring-equivalence proof.** See issues #33 and #16.
