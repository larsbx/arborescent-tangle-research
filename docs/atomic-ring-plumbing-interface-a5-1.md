# A5.1 — Source-backed atomic-tangle plumbing *input contract*

**Status:** implementable, lossless **combinatorial instructions** for the
source's atomic tangles with rings. **NOT** a 3D Conway tangle realization,
a soundness proof of R±, or a completed tangle isotopy classifier.

## Exact primary source and diagrams audited

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots and the
Classification and Symmetries of Arborescent Knots* (2016), primary PDF:

https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The following printed pages were inspected as **images**, not inferred
from noisy text extraction:

| Printed page | PDF zero-based page | Source statement | Engineering consequence |
|---|---:|---|---|
| §12.4, **155** | 166 | Modified R± **relocates** the ±2—±1—±2 cap, retaining its weighted vertices and edges. | Do not equate a disappearing cap with a topology-preserving R± move. |
| §12.5, **157** | 168 | An **abbreviated tree** erases eligible ring subtrees *as a representation change* and assigns a ring number to each surviving vertex. | Preserve **per-vertex** multiplicity, not a global total or a boolean. |
| §12.5, **158**, Figs **12.28–12.29** | 169 | A weighted vertex with ring number r specifies a **model atomic tangle with r rings**; the weighted tree specifies how those atomic pieces are plumbed together. | Record cyclic boundary ports, sector integers, ring count and explicit instructions for seams. |
| §12.5, **159–160** | 170–171 | Proposition **12.16** characterizes abbreviated **positively canonical arborescent *pairs*** by conditions W,F,S,P,A, with negative analogue. | A *pair* screen must not be silently substituted for *tangle* canonicity. |
| §12.7.2, **166** | 177 | Tangle canonicality uses **W,S,P or N**; pair-only free-bond condition **F is absent**. | Do not require pair-F when assessing a tangle. |
| §12.7.3, **167**, Thm **12.21** | 178 | Canonical V4-labelled trees classify arborescent tangles under F1,F2,F3′ and modified ring moves, with source hypotheses. | Raw formal equality is **not** an isotopy certificate without these hypotheses and the actual tangle realization. |

### An important distinction

Section 12.5's Fig 12.29 **does provide a topological model** for an
atomic tangle with rings. Implementing only the sequence of sector
weights and the number of circles labelled r **does not yet draw or
certify that model**: one must also encode the arcs/crossing operations,
boundary 4-punctured spheres and their parametrizations, and 3D gluing
maps. The source explicitly states that plumbing the abbreviated
canonical tree produces the same tangle/pair as plumbing its original
canonical weighted tree. That statement is a theorem under the
source's conditions, not an excuse to skip constructing the actual
topological model in software.

## Typed *abstract* sewing plan

A5.1's new
`arborescent/atomic_plumbing_contract.py` converts an existing
`Abbreviated` presentation to a rooted instruction graph:

- **AtomicInstruction** at each retained vertex: canonical *rooted
  address*, its **full ordered sector-weight vector**, a
  nonnegative **ring multiplicity**, and the uniform source ring
  sign. It is an instruction for the still-unimplemented
  geometric atomic-tangle renderer.
- **Port** at every incident bond: an intrinsic rooted vertex
  address and an incident slot **after rotating the entire cyclic
  bond/sector list to its incoming parent bond**. The original
  vertex numbers and arbitrary bond names never enter the plan.
- **Seam** for every internal bond, represented by the two ports it
  joins. The plan deliberately supplies **NO** strand-pairing or
  3D orientation/gluing map; those are the missing geometric stage.
- **FreeBoundary** for each free bond, carrying its explicit
  `V4={0,ξ,η,ζ}` element and a single distinguished marked root.

All supplied pieces are immutable and type-checked. Every port must
appear **exactly once**, either at a seam or as a free boundary;
the seams must form the original rooted tree; exactly one free root
must occur at address `()` and port index zero. Unsupported
or malformed inputs raise a typed error instead of silently
discarding information.

### Theorem A5.1.1 — lossless round trip of **presentation data**

Let A be any syntactically valid, rooted, V4-labelled
`Abbreviated` instance. Then the typed sewing construction P(A)
is defined; its independently reconstructed marked rooted
instruction code is **exactly** the original
`Abbreviated.rooted_code(root)`:

```
P(A).rooted_instruction_code() == A.rooted_code(root).
```

Furthermore, the instruction plan and its code are invariant under
renaming interior edge identifiers, permuting stored vertex indices,
and simultaneously rotating the stored cyclic bond **and** sector
lists at any vertex.

**Proof.** A weighted plane tree is a connected acyclic graph in
which every interior bond appears at two vertices and every free
bond once. Root it at the marked free bond. There is a unique
incoming bond at every vertex, so cyclic rotation to that bond
defines its ordered sector vector independently of storage.
Assign each vertex the sequence of child incident slots along
its unique rooted path. Each internal bond then corresponds to
one unique parent-to-child seam, and each free bond to one
V4-labelled external port.

Conversely, walk the instruction tree from the marked root,
reading the sector vector and ring number at each atomic
vertex. At every other incident port, recurse through the
unique seam or read the labelled free port. This yields the
same recursive `(ring_number,sectors,ordered_children)`
definition as `Abbreviated.rooted_code(root)`. Every port is
accounted for exactly once by the audit, so neither a sector
weight nor a ring count or V4 free label can be omitted.
Storage renumbering affects none of the rooted addresses,
cyclically aligned sectors, or recursive code. □

This is a **representation isomorphism theorem**, not an
independent realization or isotopy soundness theorem.

## Proof obligation ledger: what is NOT implemented

The `AtomicPlumbingPlan.topology_status` field is hard-fixed to
`UNVERIFIED_GEOMETRY`; attempts to set
`VERIFIED_TOPOLOGY` are rejected by the constructor.
`require_verified_isotopy()` raises `NotImplementedError` rather
than promoting the abstract data to a false equivalence.

An actual *topological* A5.2 requires at minimum:

1. For each atomic instruction, construct the **explicit model
   atomic tangle with rings** of Fig 12.29: embedded strands,
   twists/sectors and r rings, including signs and crossings.
2. For each seam, specify and verify both 4-punctured Conway
   boundary parametrizations, the actual pairing of punctures,
   orientation conventions, and induced V4 action. A named
   internal edge by itself is NOT such a map.
3. Verify that the resulting manifold/tangle is the intended
   plumbing construction of §12.5, independently of the
   recursive rooted presentation encoding.
4. In that topological model, exhibit the source's two actual
   R± moves as 3D isotopy or explicit Reidemeister/Conway
   equivalences, with an independent checker.
5. Validate the source §12.7.2 **tangle** W/S/P/N and small/adjacent
   exceptions, and control flip interactions *before* claiming
   Theorem 12.21 identifies all canonical-tree isotopy classes.

Until these gates are met, the A4.x formal classes and the A3.x
square-root/length enumerations are **only about combinatorial
presentations**, not ambient isotopy of actual tangles.

## Reproducible acceptance

Run:

```sh
python -m unittest discover -s tests -v
```

`tests/test_atomic_plumbing_contract.py` independently compares the
new typed port/seam reconstruction to old rooted
`Abbreviated` codes for both signed corolla and degree-four
ring-bearing paths, with ring multiplicities 0–3, nonroot V4
labels, weighted/angular sectors and A3 fork-decoration input.
It also enumerates multiple unreduced ring placements collapsing
to one **identical atomic instruction input**, confirming the
*notational* §12.5 contraction and nothing more.

Independent negative controls mutate a boundary V4 label, sector
location and distribution of ring counts, and ensure the
reconstructed codes differ; malformed root/seam/port layouts fail
closed. A cyclic-storage/vertex/edge renumbering control must
produce *identical* plans. Every valid plan explicitly refuses
topological certification.
