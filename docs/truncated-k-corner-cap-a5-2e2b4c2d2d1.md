# A5.2e2b4c2d2d1 — Positive-area PL interface on a true K corner

**Increment toward Issue #65; NOT yet a corner-to-edge 3D gluing.** Previous PR #64 constructs a genuine octahedral tetrahedral 3-ball around a source-candidate two-valent K vertex P. Its two K rays leave through isolated octahedron vertices (apices). Merely identifying such an apex with the center of PR #63's positive-area K edge-tube end would make a NONMANIFOLD point link.

This PR establishes one concrete, independent geometric prerequisite: **truncate the outgoing corner apex and retriangulate the remaining convex 3-ball while preserving the actual knot subcomplex and exposing a five-vertex / four-triangle disk of positive area.**

## Exact rational construction

Take the verified corner-star generators (u,v,w), with positive determinant

`Δ=det(u,v,w)>0`.

The original octahedron is the convex solid `P+xu+yv+zw` for `|x|+|y|+|z|≤1`. The outgoing knot ray is either +v (a genuine bend) or −u (collinear straight continuation). Cut the octahedron halfway to its outgoing apex, in the corresponding affine coordinate (e.g. `y≤1/2` for +v). The removed tip is a **true 3D pyramid with nonzero volume**, not a zero-area cap.

Keep the central knot vertex P and the five original outer octahedron vertices other than the removed outgoing apex. Add the new outgoing K vertex

`C=(P+outgoing_apex)/2`

and four diamond corners, each the midpoint of the old apex with one of its four neighboring octahedron vertices. These **11 rational vertices** triangulate the surviving convex solid:

- Four untouched original boundary triangles are coned to P.
- The other four original triangles are clipped into quadrilaterals; each is subdivided into TWO triangles coned to P.
- The newly exposed diamond is divided into FOUR triangles meeting at its true K center C, and these are also coned to P.

This yields **16 positively oriented 3D tetrahedra** comprising one convex PL 3-ball. The exposed physical cap has **five distinct vertices, eight edges, four triangles, Euler characteristic 1, and one boundary four-cycle**. Crucially each cap triangle is **actually an exterior face of the physical tetrahedral mesh**.

The outgoing K segment `P→C` is now a genuine internal tetrahedral edge with **four-cycle link**, while the incoming K ray is also a genuine edge. For a bend its adjacent original facet refinements create a SIX-cycle incoming edge link; for a collinear straight continuation the outgoing apex is opposite the incoming apex and its incoming link remains a FOUR-cycle. Both are independently checked rather than forcing a hardcoded original corner link after subdivision.

### Exact volume theorem

One full octahedron has volume `(4/3)Δ`. The removed similar tip pyramid above the halfway truncation plane has volume `Δ/12`. Hence

```
V_truncated = (4/3 − 1/12)Δ = (5/4)Δ
            = (15/16)·V_original.
```

The independent physical tetrahedral verifier recomputes oriented simplex volumes, all face incidences, 3D vertex-link disk/sphere topology, connected Euler-2 boundary and sum of 16 tetrahedral volumes. Source corner geometry is revalidated against **every physical K segment and Conway closed solid** via the existing PR #64 full convex-octahedron checker: truncating a subset cannot introduce new exterior intersection.

## Source limits and next positive interface-glue obligation

This milestone creates the **actual nondegenerate half-interface** needed before attempting any 3D attachment, but it does NOT claim to have joined the truncated corner cap to the differently framed positive-area square end of PR #63's K-edge tube. A future transition needs a proper, independently audited rational PL prism connecting both *actual physical section triangulations*, together with a common mesh refiner; matching only their K center coordinates or making the surfaces touch at one point is strictly forbidden.

The positive result is a local **K-conforming tetrahedral 3-ball with a positive-area outgoing disk**, not an embedded whole-K regular neighborhood. The knot's original corner and incident physical edges remain source-candidate geometries; Bonahon–Siebenmann Figure 12.5/12.10 source plumbing-patch identity, global S³ conformance, four NEW Conway puncture reconstruction, signed/pure boundary framing, R± and complete isotopy classification remain open.

`TRUNCATED_K_CORNER_CAP_NOT_GLUED` and the original A5.1 `UNVERIFIED_GEOMETRY` gates remain mandatory. Issue #65 remains OPEN, as do #62/#60/#58/#33/#16.

## Validation

`tests/test_truncated_k_corner_caps.py` exercises actual source-candidate signed bends and straight return continuation, 11-vertex/16-tetrahedron exact geometry, exposed five-vertex cap disk, independent incoming/outgoing link cycles, exact 15/16 volume, and adversarial cap-face, K-center, tetrahedron, source-status and unsupported ring/valence tampering. Both GitHub-hosted Python 3.11/3.12 must pass on final head and post-merge main.
