# A5.2e2b1 — Exact ribbon-surface quotient and the missing source attachment

**Status: ABSTRACT_SURFACE_SOURCE_UNVERIFIED.** This is an independent finite 2-complex **topology checker** on previously verified A5.2e2a PL disks, **not a 3D embedding of their identifications**, a source-complete twisted band, or a source R± isotopy proof.

## Source interpretation and gap

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.2 printed p.142 Figures 12.8 and 12.10 and p.143, and §12.5 printed p.158 Figure 12.29; source page images inspected.

Primary monograph: https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

Figure 12.8 distinguishes the globally connected **twisted band containing the knot** from just its boundary strands. Figure 12.10 marks **plumbing patches on that surface**, and p.143 states that cyclic permutations of a one-vertex weighted atomic presentation correspond to **band-preserving** tangle isomorphisms. The r=0 case of Figure 12.29 has the same source atomic construction as Figure 12.8. The source's framing/attachment choices are *geometric data*, not consequences of having the same set of labelled four-puncture endpoints.

A5.2e2a produced **n independently embedded PL disk patches**. Their union has n connected components, Euler characteristic n, and n circular surface boundaries. Those statements are about the actual meshes, **not** an already constructed whole band. To establish the source global band, we must specify and justify how the patches are attached around/along the Conway boundary surfaces, and then realize those attachments **geometrically in three dimensions**. A5.2e2b1 handles the first, purely **topological bookkeeping and obstruction** part.

## Typed input and independently computed topology

`AbstractRibbonSurface(patches,attachments)` takes exact independently audited `TriangulatedRibbon` meshes and explicit `CapAttachment(RibbonCap(...),RibbonCap(...),reverse_width)` entries. Each cap is a three-vertex end chord, with two oriented boundary segments. Each cap may be used at most once. Only original boundary edges may be paired. The `reverse_width` bit explicitly selects whether the two cap orientations match or reverse their width order.

The verifier treats every input triangle as a distinct 2-cell and constructs an **abstract quotient** of the input meshes by exact endpoint and edge identifications. It independently computes:

- Vertex and edge equivalence classes by separate disjoint-set structures. Every quotient edge must have exactly one or two incident triangles. No original nonboundary edge may be used as a cap.
- Each vertex's **combinatorial link** from triangle corners and incidences of paired edges. It must be one connected interval (surface boundary) or circle (interior). This rejects pinched vertices that Euler characteristic and boundary counts alone might miss.
- A degree-two boundary multigraph, including loops, traced into the exact number of boundary components.
- Face adjacency and orientation constraints, considering the *actual* edge direction under each proposed seam. Solves the ±1 constraints to decide if the resulting abstract surface is orientable.
- Connected components and Euler characteristic `χ = #vertex classes − #edge classes + #faces`.

No result is accepted on trust from caller-supplied Euler/Betti/genus metadata. Invalid/reused caps, inconsistent edge identification, nonmanifold vertex links, missing mesh faces, and fake status must fail closed.

### Exact obstruction: same patch, annulus versus Möbius band

Take **one and the same** already embedded rectangular triangulated PL disk from A5.2e2a. Identify its two *abstract* three-point end caps. If the first/second ends are paired in the same across-width order, the quotient has

`components=1, χ=0, boundary_components=2, orientable=True`,

i.e. abstract annulus topology.

If instead one cap is reversed in width, the quotient has

`components=1, χ=0, boundary_components=1, orientable=False`,

i.e. abstract Möbius-band topology.

Both quotients start with **exactly identical original physical strand coordinates and labelled endpoints** and differ solely by their explicit *attachment*. They therefore exhibit the obstruction:

> An independently verified collection of PL ribbon patches with correct boundary endpoints **does not determine the global supporting band topology** without further gluing/framing data.

This is a finite surface-topology counterpart to A5.2b's independent infinite-pure-Conway-mapping-class obstruction. Neither obstruction says the actual source construction is ambiguous: Bonahon–Siebenmann's source **does** prescribe a model; our code has not yet extracted and implemented its precise 3D band attachment rule.

### Second positive topology check

Abstractly glue the end cap of one ribbon disk to the start cap of a second, using either width identification. The independent validator computes a single **disk** (χ=1, one boundary circuit, orientable). A chain of n disk patches with n−1 abstract cap attachments has the same single-disk topology. This proves sound basic manifold accounting, *not* that the resulting disk agrees with the desired multiply punctured global source band.

## Tests and specific nonclaims

`tests/test_ribbon_surface_quotient.py` covers n=2,3,4 disconnected physical patches, one-cap gluing, annulus vs Möbius quotients from *the identical initial physical patch*, an exhaustive small sector/attachment sweep, illegal caps, reused seams, malformed/status-forged inputs, and mandatory source-certification refusal.

Run both GitHub-hosted Python 3.11 and 3.12 on the exact PR head, followed by post-merge `main`.

**No assertion of realized cap identification in Q³:** The original cap curves are generally *at distinct spatial positions*. The surface quotient identifies them *combinatorially only*; an ambient embedding requires explicitly constructed connecting/attachment geometry and full ribbon framing. No claimed "positive isotopy" follows.

**No source-derived gluing yet:** The source's Fig.12.10 twisted band still requires actual embedded cap-to-Conway-sphere attachment, figure-supported pure boundary mapping-class choices, global component/boundary topology, cyclic-rotation band-preserving isotopy, signed half-twist/framing controls and independent triangle/triangle intersection verification. Issues #43/#41/#39/#33/#16 remain open.

**No R± certificates and no isotopy-class counting:** The existing A5.1 `UNVERIFIED_GEOMETRY` boundary is unchanged.
