# Arborescent tangle research

Research question: **which Catalan-type enumeration laws survive reduction and isotopy in explicitly specified arborescent tangle/link families?**

This repo begins with a deterministic **negative control**, not a proof of a general enumeration theorem.

## Three different counts

For an alphabet of m rational-tangle slopes and k ordered leaves:

1. Full ordered binary presentation trees: **C_(k-1) m^k**.
2. Presentations modulo formal associativity of horizontal sum: **m^k** ordered words. Catalan parenthesizations disappear.
3. Words modulo cyclic rotation and reversal: **B_m(k)** (bracelets), modeling cyclic **Montesinos link-closure presentations**, NOT isotopy classes of boundary-marked tangles.

The third counting series is logarithmic and non-algebraic for m>=2; Catalan syntax by itself does not guarantee a Catalan or algebraic series after quotienting.

## Run

Standard-library-only Python 3.11+:

```sh
python -m unittest discover -s tests -v
python -m arborescent.prototype --slopes 1/2 1/3 2/5 --max-k 8
```

The slopes must be distinct rational numbers in (0,1). Their values distinguish the atomic labels; the resulting presentation-orbit counts depend on the alphabet size, not the numerical values.

See [research program](docs/research-program.md), [exact prototype counts](docs/prototype.md), and [literature](docs/references.md). The program does not decide tangle isotopy, test non-rationality, calculate crossing numbers, or implement weighted arborescent-tree normal forms.

Next steps: establish the precise link-classification applicability and exceptional cases; specify and implement genuinely branching reduced weighted planar trees, including sector weights and the complete Bonahon–Siebenmann equivalence moves; test normal-form completeness; only then conjecture and attempt to prove a tree-type asymptotic for a specific nondegenerate grammar.

## A1 implemented: genuine branching weighted plane-tree syntax

[The A1 note](docs/weighted-tree-a1.md) specifies an explicit sector-weighted
planar rotation system with free bonds, checks **strict, nonexceptional**
Bonahon–Siebenmann (W), (F), (S) conditions, and proves exact counts for a
finite binary/ternary weighted-fork syntax grammar:

```text
F_m(x) = m + x*(F_m(x)^2 + F_m(x)^3)
a_n(m) = sum(j=0..n) (2n+j)! / ((n-j)! j! (n+j+1)!) * m^(n+j+1)
```

For m=2 and n=1..4, the syntax counts are 12, 192, 4080, 99264.
Exhaustive small-size enumeration, closed formulas, and W/F/S negative
controls have tests. This is still **not** the full Bonahon–Siebenmann
calculus, a topological canonical form, or an isotopy-class enumeration.

Run all tests: `python -m unittest discover -s tests -v`.

## A2.1: sourced elementary plumbing flips and exact restricted quotient

The [A2.1 flip-subgroup proof and experiment](docs/flip-subgroup-a2-1.md)
implements the Bonahon–Siebenmann **F1/F2** moves on weighted planar-tree
rotation systems, respecting angular sectors. On a connected tree these
moves generate a subgroup of size **at most four**: independent reflection
of the two vertex-bipartition classes. The quotient of the A1 finite
2/3-fork grammar is **exactly algebraic**, with fixed-point series
`E=m+x F(x²)(1+O)`, `O=m+x(E²+E³)`, `G=m+x F(x²)(1+G)`,
and `Q=(F+E+O+G)/4`. For m=2, the first positive-size F1/F2
presentation orbit counts are **9, 81, 1251, 26688, 669420**.
A square-root singularity proof shows the quotient's leading
asymptotic is one quarter of the raw syntax count.

These are **pair-isomorphism** moves from §12.3, not general
boundary-fixed-tangle isotopies. **A2 is not complete:** F3'/F3'', ring,
boundary and arithmetic moves, P/N canonicality and exceptional trees
are still outstanding. Reproduce with:

```sh
python -m arborescent.flips --max-n 6 --alphabet-size 2
python -m unittest discover -s tests -v
```

## A2.2: Klein-four boundary labels and a strict F3′ move

[The A2.2 tangle-calculus note](docs/v4-tangle-moves-a2-2.md) records a
critical distinction from A2.1: Bonahon–Siebenmann §12.7 requires **V4
labels on every free bond**. Its tangle F1 and F2 moves alter those labels
(ζ globally; ξ at even distance and η at odd distance). A restricted
W-preserving **F3′** weight transport is now implemented, including the
odd/even parity effect on the attached subtree and its boundary labels.

**New exact correction:** adding all four root-boundary labels makes
the F1/F2 action free on the A1 one-root-port grammar. Its correctly
V4-labelled presentation quotient has coefficients **a_n**, not the
unlabelled pair-presentation quotient **Q_n**. For m=2,n=2: raw=192,
unlabelled pair quotient=81, labelled tangle presentation quotient=192.
Neither is a full arborescent isotopy enumeration.

The §12.4 Free-Bond canonicity condition **F** pertains to the *pair*
calculus; §12.7.2 omits it for *tangles*. Tangle W/S/P/N and ring
equivalence remain open; the existing `wfs_failures()` must not be
interpreted as a tangle-canonicality certificate.

## A2.3: bounded replayable move search and new obstruction

[The A2.3 exact audit](docs/partial-move-closure-a2-3.md) fixes a
previously merged A2.2 test-fixture bug that caused hosted CI to fail.
It adds bounded exhaustive search under only implemented F1/F2/strict
F3′ moves on fixed indexed graphs, with replayable shortest move traces
and explicit incomplete status when the budget is reached.

**New exact obstruction:** The A1 grammar places nonzero weights only
at valence-one leaves. The strict F3′ weight-transport move requires two
distinct adjacent sectors at a weighted vertex. Therefore **F3′ is
inapplicable on every A1 syntax tree, even after F1/F2**. Its earlier
implementation must not be used to infer quotient behavior on A1.
Future nontrivial enumeration requires a different, predeclared grammar
with weighted internal forks and topology-valid tangle W/S/P/N conditions.

## A2.4: non-vacuous strict F3′ root-pivot quotient

[The exact root-pivot grammar and proof](docs/root-pivot-a2-4.md)
adds **one weighted internal root fork** to the A1 zero-weight-fork
language. Its strict root-F3′ transport is genuinely applicable,
unlike the old leaf-weight-only A1 grammar. With `m` distinct leaf
weights and `p` nonzero root weights, we have the *presentation*
generating functions:

```text
F = m + x(F² + F³)
H = p*x*(3F² + 4F³)        # all positions of the root weight
J = p*x*(F² + F³) = p*(F-m)  # quotient by root F3′ only
```

The quotient has an exact sector-zero normal form because the allowed
root-sector moves form a path. It preserves the same square-root
syntax asymptotic but does **not** establish a complete arborescent
isotopy-class enumeration or Bonahon–Siebenmann canonicality.
Reproduce with `python -m arborescent.root_pivot --max-n 6`.

## A2.5 — ring-signature obstruction

[The source-scoped necessary-pattern detector](docs/ring-signature-a2-5.md)
locates ±1-weight valence-three vertices with adjacent matching ±2
terminal weights. It does **not** implement R±. A1 and default A2.4
grammars have no ±1 vertex; since F1/F2/strict-F3′ preserve individual
vertex totals, they cannot acquire this necessary modified-ring motif.
The new unit-weight/±2-tip test alphabet offers a concrete next input
for a figure-exact ring rewrite with V4 boundary bookkeeping.

## A2.6 source correction: ring number is *not* a ring move

[The A2.6 ring-abbreviation audit](docs/correct-ring-abbreviation-a2-6.md)
corrects closed, **unmerged** PR #11. The source's modified R± moves
**relocate** an attached ±2–±1–±2 ring subtree in the plane; they do
**not delete its three vertices**. The source *separately* introduces
ring-numbered abbreviated trees (§12.5), erasing suitable **disjoint**
ring subtrees and storing their multiplicity at attachment vertices.

`ring_abbreviation.py` implements only a conservative, lossful
**notation change**, preserving weighted sectors, surviving V4 labels,
and explicit ring multiplicities. The exceptional case of one ring
attached to another is rejected, not silently normalized. No full R±
move, tangle isotopy classifier, or general Catalan quotient is claimed.

## A3: exact ring-numbered abbreviated tree grammar

[The A3 bivariate enumeration theorem](docs/ring-number-grammar-a3.md)
counts *syntactic* abbreviated presentations with an independent
ring number 0..R at each internal 2/3-ary fork, retaining explicit
V4-marked weighted-plane-tree codes. With `m` finite leaf weights:

```text
F_(m,R)(x,y) = m+x*(1+y+...+y^R)*(F_(m,R)^2+F_(m,R)^3)
[x^n y^q]F_(m,R) = a_n(m)*[y^q](1+y+...+y^R)^n
sum_q [x^n y^q]F_(m,R) = (R+1)^n*a_n(m)
```

Summing over ring counts retains the tree square-root exponent
`n^(-3/2)` while multiplying exponential growth by `R+1`.
Holding total ring number q fixed instead yields a
`n^(q-3/2)` polynomial prefactor. These are **combinatorial
presentation results**, not classification of canonical tangles or
a proven complete R± isotopy quotient.

## A3.1: finite-local ring compatibility alters growth

[The exact constrained grammar theorem](docs/local-ring-constraints-a3-1.md)
adds a **synthetic** (not source-canonical) local exclusion rule:
adjacent fork vertices may not both have positive ring numbers.
Its finite two-state algebraic equations are
`Z=m+x((Z+P)^2+(Z+P)^3)`,
`P=x(y+...+y^R)(Z^2+Z^3)`.
For `m=2,R=1,n=2` the guarded syntax count is **576**, versus
**768** under independent ring numbers. At *fixed total ring grade* q,
the leading factor remains `n^q/q!` (up to lower order);
when all grades are counted, a matching argument proves a **strict
exponential deficit** relative to the unrestricted grammar.
This is a falsifiable demonstration of finite-local sequence
classification, not a Bonahon–Siebenmann canonicality theorem.

## A3.2: exact exponential *presentation* fibers of ring abbreviation

[The explicit ring-location fiber theorem](docs/ring-placement-fibers-a3-2.md)
constructs a marked zero-sector base with `d` V4-labelled free legs
and `r` attached disjoint ±2–±1–±2 ring caps. There are exactly
`binom(d+r-1,r)` different rooted cyclic placement codes, but all
have **one** identical ring-numbered abbreviated code. For `d=r=n`,
this fiber is `binom(2n-1,n)~4^n/(2sqrt(pi n))`—exponential.
This proves **no polynomial-fiber bound for unreduced placements at
unbounded base valence**. It says nothing by itself about fibers of
already-reduced bounded-valence grammars or topological isotopy.

**Stronger A3.2 negative control:** an exact path-of-bases family has
**maximum bond valence four** and one ring cap at each of `k` bases.
Its two local cyclic cap positions per base produce **2^k** different
marked presentation codes, yet §12.5 ring-number abbreviation gives
**one** identical code with one ring at each base. Thus bounded
valence **alone** cannot give polynomial unreduced-placement fibers;
ring-location normalization is essential. This is not an isotopy count.

## A3.3: square-root theorem for the constrained ring grammar

[The A3.3 analytic and exact-algebraic proof](docs/guarded-critical-growth-a3-3.md)
shows the A3.1 **synthetic finite-local** ring exclusion retains
`K*mu^n*n^(-3/2)` total presentation growth for all fixed finite
leaf alphabet sizes `m>=1` and positive ring alphabet sizes `R>=1`.
Eliminating the two grammar states gives the scalar smooth implicit
function `G=H(x,G)`, with a positive nondegenerate critical pair.
For `m=2,R=1`, the exact critical discriminant contains a
specified degree-11 polynomial whose positive characteristic root
has `mu≈58.8432259437452`. Exact integer recurrence and rational
root-bracketing tests are included. This does **not** establish a
canonical tangle growth law or an isotopy quotient.

## A4.1 — source-motivated formal ring *relocation* (no cap deletion)

[The A4.1 proof and tested finite certificates](docs/ring-relocation-a4-1.md)
add a conservative, typed **adjacent ring-position slide** on
zero-sector weighted plane-tree presentations. It swaps a full
±2–±1–±2 ring subtree's stem with an adjacent ordinary bond **at
the same base**; it never deletes vertices or changes V4 labels,
and preserves the §12.5 abbreviated code. Overlapping rings, mixed
signs, nonzero base sectors, root-anchor crossings and ring-ring
swaps are rejected.

For the previously established A3.2 families, the finite formal
slide graph is connected on all binomially many corolla
placements and all `2^k` bounded-valence path placements, with
replayable shortest move witnesses. These are **presentation graph
moves**, not independently certified 3D R± tangle isotopies.
Complete Bonahon–Siebenmann Theorem 12.21 compatibility remains
open in issue #16.

## A4.2 — stable normalization of the strict formal ring slides

[The A4.2 proof and replay certificates](docs/ring-normal-form-a4-2.md)
orient the A4.1 formal adjacent ring/ordinary-bond slide as
`RO→OR`. The inversion count `#(R before O)` decreases exactly by
one on each step. Thus the finite slide system has a **unique
stable normal form** placing every ring after all ordinary
nonanchor bonds at its base, preserving their respective relative
orders. An implementation replays every typed move and checks
termination. For a base with `o` ordinary bonds and `r` identical
ring caps the rooted-presentation orbit contributes exactly
`binom(o+r,r)` positions; factors multiply over independent bases.

This is a mathematically exact **formal syntax normalizer**, not an
independently certified Bonahon–Siebenmann R± topological normal form
or knot/tangle isotopy classifier.

## A4.3 — rooted cyclic-storage invariant ring normal keys

[The A4.3 correction](docs/root-aligned-ring-normal-form-a4-3.md)
fixes a subtle presentation-isomorphism issue in A4.2:
`Vertex.slots[0]` is **not intrinsic** for a cyclic incident-bond
order. A cyclically rotated representation could be wrongly rejected
or normalized differently simply because a ring stem moved into
storage slot zero. A4.3 aligns every local cyclic order on the
unique **incoming bond toward the marked free root**, rotating its
sector weights simultaneously, before running A4.2's explicitly
anchored stable ring sort.

The root-alignment certificate distinguishes the **presentation
isomorphism** from the typed **formal ring slides**; tests verify
invariance under vertex/edge renumbering, cyclic starting slots,
both ring signs, arbitrary path/corolla placements and V4 labels.
Neither layer certifies actual tangle isotopy or a complete R±
classification.

## A4.4 — exact F1/F2-equivariance of restricted ring abbreviation

[The A4.4 commuting-square theorem](docs/abbreviated-flip-equivariance-a4-4.md)
implements the **source V4-labelled F1 and F2 flips** on the
ring-numbered abbreviated representation. For disjoint strict
±2–±1–±2 caps with zero sectors adjacent to each stem and no
adjacent-ring exceptional case, performing each source flip before
the §12.5 representation change yields **exactly** the same
abbreviated object as flipping afterward. For F2 the pivot is a
**surviving marked free bond**, so parity does not depend on erased
vertex indices. Ring numbers remain fixed at surviving vertices,
and all surviving V4 labels receive the correct ξ, η, ζ actions.
The four-action subgroup admits an identity-root V4 gauge and
rooted presentation key.

This is an exact formal **equivariance theorem**, **not** an
independent 3D isotopy proof, F3′ completion, or full canonical
tangle classification.

## A4.5 — guarded F3′/abbreviation commuting squares

[The A4.5 restricted naturality theorem](docs/abbreviated-f3prime-naturality-a4-5.md)
extends the A4.4 V4-labelled F1/F2 commuting-square checks to the
previously implemented strict **F3′** flip on a **surviving** internal
edge at a pivot identified by a **surviving free bond**. For disjoint
±2–±1–±2 caps on zero-sector attachment bases that remain strictly
recognized after F3′, abbreviation before or after the flip yields
**exactly** the same ring-numbered weighted planar-tree data.
Tests exercise both signs, odd/even pivot weights, parity-dependent
η changes to nonroot V4 labels, inverse moves, and fail-closed
guards for erased stems/unrecognized cap configurations.

This proves a restricted **formal equivariance**, not a complete
F3′/R± topology calculus or any claim of arborescent-tangle isotopy.

## A4.6 — joint F1/F2 V4 gauge and strict ring-position normalization

[The A4.6 formal quotient theorem and replay certificates]
(docs/joint-flip-ring-quotient-a4-6.md) combine the **source F1/F2**
four-action V4-labelled tangle flip subgroup with A4.1's
**strict formal** signed ring-position slides. First choose the
unique flip-group variant with identity root V4 label; then
root-align cyclic slots and apply A4.2's terminating `RO→OR`
ring-position normalizer.

For a fixed residual rooted skeleton and identical signed disjoint
ring caps on zero-sector bases, the resulting marked presentation
key is invariant under **both** generators and agrees, as an
equivalence criterion in the declared restricted domain, with the
independently computed §12.5 abbreviated F1/F2-orbit key.
Exact combined *indexed presentation* orbit counts are
`4·binom(d+r-1,r)` for a corolla and `16·3^(k−2)`
for a k-base bounded-valence path (k≥2).
Typed joint-move BFS supplies shortest replayable witnesses with
explicit cutoff status. **None of these counts are isotopy classes.**

## A4.7 — guarded mixed F1/F2/F3′/formal ring move graph

[The A4.7 replayable partial-calculus note](docs/mixed-partial-calculus-a4-7.md)
extends A4.6's strictly formal F1/F2/adjacent-ring actions with
the **source-scoped strict surviving-edge F3′** move, provided its
A4.5 ring-abbreviation commuting square passes. Every step rechecks
complete signed ring-cap recognition, per-base ring counts, V4
boundary data, and marked-root cyclic presentation alignment.
A bounded BFS supplies executable **positive transformation
certificates**, marking `complete=False` on search cutoff.
Unlike A4.6, **no complete normal form is claimed for the
F3′-extended graph**. A real seven-vertex witness verifies a
specific F3′/formal-ring-slide commuting diamond (both ring signs,
positive/negative odd/even root weights); zero-weight corollas
recover the earlier finite formal flip+slide orbit. No 3D
isotopy or source-complete R±/F3′ claim is made.

## A4.8 — necessary invariants and honest four-valued formal reachability

[The A4.8 proof and negative controls](docs/mixed-necessary-invariants-a4-8.md)
add an exact **rooted unordered weighted-ring skeleton fingerprint**:
retain each surviving vertex's *total* integer weight, per-base
ring number, incidence tree and count of incident unmarked free
legs, while forgetting cyclic slot positions and individual V4
labels that the declared F1/F2/F3′ moves may change. Every A4.7
**partial formal** move preserves this code, so unequal codes
supply a mathematically valid *necessary-invariant exclusion*
without search. Equality is neither sufficient nor evidence of
3D tangle isotopy.

A dedicated four-way API distinguishes `reachable_formally`
(with replayable typed trace), `excluded_by_necessary_invariant`,
`not_reached_in_complete_indexed_graph`, and
`undetermined_cutoff`. An independent negative control
proves why even a **complete indexed BFS** is not a presentation
isomorphism quotient: renaming internal edges and reordering
vertex storage can yield the *same* marked plane presentation
outside that indexed orbit. No isotopy-class nonreachability
or complete F3′/R± classification is claimed.

## A4.9 — partial formal search modulo rooted presentation isomorphism

[The A4.9 congruence proof and class-BFS audit](docs/rooted-class-quotient-a4-9.md)
repairs the remaining A4.8 distinction between a **complete
indexed** finite search and a complete search over **rooted
weighted plane-tree presentation-isomorphism classes**. A4.9
uses the existing faithful V4-labelled rooted plane code as
a hash key. Each class's stored representative is expanded
under all the *declared* A4.7 partial formal moves; tested
generator congruence under vertex/edge renaming and cyclic
slot rotation proves this class-level search is complete
when its own class budget is not exhausted.

Positive evidence consists of a replayable typed formal
move path **plus an independently checked final marked
presentation isomorphism**. A complete class graph can now
exclude the target *within the declared partial formal
calculus*; such exclusion is emphatically **not** evidence
of different arborescent-tangle isotopy classes. Unknown
budget cutoffs remain explicit.

## A4.10 — Gaussian distance ranks and shortest formal ring-slide witnesses

[The exact A4.10 metric theorem](docs/formal-ring-geodesics-a4-10.md)
upgrades A4.1/A4.2's **restricted formal** ring-position
connectivity to full distance distributions and explicit
shortest paths. For o ordinary bonds and r identical ring
caps, the distance-to-stable-form polynomial is the
**Gaussian binomial** `[o+r choose r]_q`, enumerating
ring-before-ordinary inversions. On the bounded-valence
k-base path the exact polynomial factors as
`(1+q)^2 (1+q+q^2)^(k-2)`, refining the earlier
`4·3^(k-2)` orbit cardinality. For same-index
presentations differing only by cap/ordinary interleavings,
the shortest slide distance is the total L1 displacement of
ordinary-bond counts before each ring. Code constructs
and replays a typed geodesic with that minimal length,
verified independently against exact BFS.

These are **formal presentation metrics**, not independent
Bonahon–Siebenmann R± or knot/tangle isotopy distances.

## A5.1 — source-backed atomic tangle *plumbing-input* contract

[The A5.1 source-image audit and lossless sewing proof](docs/atomic-ring-plumbing-interface-a5-1.md)
addresses the missing geometric bridge rather than inventing
additional formal quotient moves. Bonahon–Siebenmann §12.5
printed p.158, Figures 12.28–12.29, defines **model atomic
tangles with rings** at weighted vertices and plumbs them
according to internal edges. New typed code extracts a
storage-independent **abstract instruction graph**:
per-vertex cyclic integer sector weights and ring numbers,
one exact port for every incident bond, explicit seam
endpoints, and every free-bond V4 label. An independent
recursive reconstruction is proven/tested to equal the
original marked `Abbreviated.rooted_code`, with guards
against lost ports, multiplicities, sectors or labels.

**No 3D Conway strand diagram or boundary gluing is constructed
yet.** `topology_status=UNVERIFIED_GEOMETRY` is enforced
and any purported isotopy certification fails closed. This
is the correct input contract for an independently verified
Figure-12.29 realization, not an isotopy classifier.

## A5.2a — actual 4-ended braid diagrams, typed puncture maps and RII certificates

[The A5.2a exact limited-topology theorem and source comparison]
(docs/four-puncture-braid-kernel-a5-2a.md) introduce a **real
geometric diagram kernel** rather than another abstract tree
rewrite: explicit over/under crossing events for monotone
two-strand braids with four individually marked endpoints,
complete projected strand paths and deterministic SVG inspection
drawings. A typed **Reidemeister-II cancellation/insertion
certificate** is independently replayable and proves positive
diagrammatic isotopy relative to those endpoints when the two
braid exponents agree. The exact statement is scoped to
**two-strand braids B₂≅Z**, not all tangles.

The four cyclically named punctures `(NW,NE,SE,SW)`
admit an explicit faithful V4 permutation convention.
A proposed A5.1 Conway seam map is checked against
existing parent/child ports, but remains `ABSTRACT_ONLY`;
it does **not** implement orientation-correct 3D Conway
gluing. This kernel does **not yet realize** the source
§12.5 **Fig.12.29 model atomic tangle with rings** or
certify the modified R± moves. The A5.1 atomic sewing
plans retain `UNVERIFIED_GEOMETRY`.

## A5.2b — four puncture/V4 labels do NOT specify Conway seam topology

[The A5.2b pure mapping-class obstruction](docs/conway-pure-mapping-class-a5-2b.md)
makes precise why A5.2a's explicit puncture bijection cannot
yet count as geometric plumbing. The four-punctured sphere
has an **infinite pure mapping-class group**, represented
via the torus double cover by the congruence subgroup
`Gamma(2)/{±I}`. For every integer n the shear
`T_n=[[1,2n],[0,1]]` fixes **each** of the four
punctures, hence has the **same** V4/point permutation;
yet it sends primitive slope `(0,1)` to `(2n,1)`,
so all n give **different** marked-sphere mapping
classes. An exact typed witness can attach these
two distinct pure mapping classes to the **same
existing A5.1 parent-child seam** while holding its
four-puncture V4 matching fixed. The program
checks all branch points and primitive slopes.

This is a mathematically decisive **interface
incompleteness obstruction**: the actual
orientation-reversing Conway gluing and source's
sector/ring-dependent pure twist conventions
must be specified and independently validated.
The code still refuses all 3D isotopy certification;
actual §12.5 Fig.12.29 atoms and R± remain unimplemented.


## A5.2c/d — actual *partial* spatial geometry, not the source atomic tangle

[Exact local 3D sector-braid embeddings](docs/spatial-sector-geometry-a5-2c.md)
supply rational PL strands with true height-separated signed crossings and
independent Reidemeister-II step replay. The follow-up
[four-marked PL Conway boundary patches](docs/marked-conway-geometry-a5-2d.md)
adds genuine oriented cubical 2-spheres with four distinct punctures per
A5.1 port and disjoint exact two-strand braid corridors in their complement.
All signed sector weights are preserved in order and their spatial curves are
checked against the physically embedded boundary surfaces.

**NOT YET FIG.12.29:** the linear model has an *artificial anchor sphere*
and two open ends instead of the source's long returning band. It has no
actual ring loops, source-certified puncture parametrizations, complete
atomic embedding, orientation-reversing Conway seam, or verified modified
R± isotopy. These deliberate incompletenesses are encoded in
`OPEN_CHAIN_NOT_ATOMIC` and `UNVERIFIED_GEOMETRY`, not hidden behind a
formal equality. Issues #33 and #16 remain open.


## A5.2e1 — complete external *return strands*, incomplete source twisted ribbon

The [exact returning-strand model](docs/global-return-strands-a5-2e1.md)
removes A5.2d's artificial left sphere. It constructs precisely n
genuinely marked Conway 2-spheres, all n signed half-twist sectors,
and two disjoint long outer return arcs. Every one of the 4n
punctures is an endpoint of exactly one of 2n embedded PL strands.
An **independent exact-rational intersection checker** verifies every
curve against every solid Conway cube and other curve segment,
including self-intersections and false over/under collisions.

This is an important **partial** spatial realization of Fig.12.8/12.29:
the long *strands* are embedded, but the complete supporting
**twisted band/ribbon surface**, source framing, exact source
boundary-parametrization conventions and r>0 ring loops remain
unconstructed. The status is `RETURN_STRANDS_SOURCE_UNVERIFIED`.
The A5.1 `UNVERIFIED_GEOMETRY` gate and issues #39/#33/#16 stay open.


## A5.2e2a — actual 2D PL ribbon patches; source band still unverified

The [explicit oriented triangulated ribbon patch construction]
(docs/pl-ribbon-patches-a5-2e2a.md) thickens the A5.2e1
global return-strand candidate into **n genuine embedded
2D PL disk patches** with exact rational triangle tables.
Independent mesh validation checks oriented internal edges,
surface Euler characteristic, nonzero triangle area, the
complete two-rail-plus-cap boundary, exact clipping against
all n solid Conway balls, and exterior-arch/interior-sector
disjointness. Unsupported ridge/frame mutations are rejected.

**These are not yet the source's complete framed twisted
band or its Fig.12.8/12.29 isotopy class.** The source-specific
global band attachment, crossing/framing signs and boundary
parametrizations are not established; no ring embeddings,
Conway 3D sewing, R± equivalences or classification follows.
`RIBBON_PATCHES_SOURCE_UNVERIFIED` and A5.1
`UNVERIFIED_GEOMETRY` remain mandatory. Issue #41 stays open.


## A5.2e2b1 — source band attachment is missing mathematical data

[The independent abstract PL ribbon quotient validator]
(docs/abstract-ribbon-gluing-a5-2e2b1.md) checks proposed
three-vertex cap attachments between previously constructed
ribbon disks, deriving surface connectedness, Euler characteristic,
oriented-face consistency, individual vertex links and exact
boundary circuits from the supplied face tables and glue maps.

Crucially, gluing the *same* original patch ends with different
width orientation can produce an **annulus** or a **Möbius band**.
Thus correct signed braid strands, four marked points and real
individual ribbon disks do **not** alone determine the desired
global Fig.12.10 twisted-band topology. This is a computable
attachment-data obstruction; it does not say the source is ambiguous.
The quotient is **abstract only**: the glued caps remain physically
disjoint in Q³ until a source-checked 3D attaching construction is
provided. The result is `ABSTRACT_SURFACE_SOURCE_UNVERIFIED`,
never a Conway/plumbing or ambient-isotopy certificate.


## A5.2e2b2 — physically joined cubical Conway collar candidate (not source verified)

The [exact Q³ Conway roof-collar construction](docs/spatial-conway-collars-a5-2e2b2.md)
adds **real 2D triangle meshes** joining all n earlier ribbon
disks along the actual, exactly coincident end-cap edges on
their Conway cubes. Unlike A5.2e2b1's abstract cap quotient,
the new surface has physical shared 3D vertices and mesh edges.
A separate coordinate-only checker welds those points and
verifies manifold vertex links, connected components, Euler
characteristic, boundary circles and orientability. Its output
is a connected χ=0 **chosen coordinate-model band candidate**.

This does **not** establish the source's Fig.12.10 framing or
transverse Conway plumbing: the roof collars currently lie on
the cubical Conway boundary over **2D regions**. Their
source-canonical off-boundary realization, local triangle
intersection certification and band-preserving cyclic
rotation witnesses are open. `COLLARED_BAND_SOURCE_UNVERIFIED`
and A5.1 `UNVERIFIED_GEOMETRY` are unchanged by this progress.


## A5.2e2b3a — independent exact triangle-intersection validator

[The rational Q³ triangle-contact kernel](docs/exact-triangle-intersections-a5-2e2b3a.md)
independently classifies skew and coplanar triangle intersections as
disjoint, point, segment or overlapping positive area. Its all-pairs
mesh checker accepts only the *declared common* mesh vertex or whole
common edge, rejecting physically intersecting triangle interiors
that a valid combinatorial vertex-link/Euler audit might miss.
No floating tolerances are used; small rational separations are exact.

This adds a necessary test for the chosen A5.2e2b2 collared-band
candidate. It does **not** fix source Fig.12.10's Conway marking,
off-boundary plumbing collar geometry or ambient-isotopy obligation.
Results remain source-unverified; issue #46 and #33 stay open.


## A5.2e2b3b — off-boundary embedded Conway collar candidate

The [off-boundary exact-PL collar construction](docs/off-boundary-conway-collars-a5-2e2b3b.md)
replaces A5.2e2b2's roof strips lying on cubical Conway faces.
It physically connects each original ribbon to its neighbor using
an exact outer U-strip with roof height z=4 (> Conway-ball z=3),
meeting each Conway sphere **only at its two end-cap chords**.
An exact triangle-vs-solid-cube audit checks all collars against
every removed ball. The new independent triangle/triangle kernel
checks **all** physical mesh-face pairs in addition to the
vertex-link/Euler/orientability topology checker.

This is a *verified embedding of a chosen-coordinate PL band*,
**not a source-certified Fig.12.10 twisted band or plumbing
parametrization**. Source-relative isotopy, cyclic band-preserving
moves, r-ring embeddings and R± still require separate proofs.
The model stays `OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED`; A5.1
`UNVERIFIED_GEOMETRY` remains mandatory.


## A5.2e2b4a — source Fig12.6–12.7 oriented Conway arrows and signed boundary maps

[The primary-source arrow/framing audit](docs/source-arrowed-conway-frames-a5-2e2b4a.md)
adds a previously omitted boundary datum: Bonahon–Siebenmann §12.1's
**oriented equatorial slope-zero arrow**. Fig12.7 demands that both
plumbing arrows emanate from a common puncture along *distinct*
adjacent arcs. Code checks a real (CHOSEN) boundary arrow on each
Conway cube and an exact determinant −1 signed GL₂(Z) boundary
mapping-class proposal with full puncture/slope action. Pure-twist
negative controls demonstrate that the arrow, four points and gluing
orientation still do not fix the source's required 3D seam/framing.
No physical source homeomorphism or isotopy proof is inferred.
Statuses remain `ARROW_FRAMING_SOURCE_UNVERIFIED` and
`UNVERIFIED_GEOMETRY`. Issue #49 remains open.


## A5.2e2b4b — source K is a distinguished proper subgraph of band boundary

[The source knot-on-band properness audit](docs/distinguished-proper-knot-on-band-a5-2e2b4b.md)
enforces the precise wording of Bonahon–Siebenmann §12.2:
the twisted band has **boundary containing K**, not
necessarily boundary **equal** to K. The entire χ=0
annulus/Möbius auxiliary band boundary consists of closed
circuits, but the source (M,K) tangle has **2n proper
K-strand arcs terminating at all 4n marked Conway
punctures**. A separate exact physical checker now proves
that every K edge is an *exposed* triangle-mesh boundary
edge, the graph K alone has degree one **only** at those
punctures, and every puncture has a nonzero boundary-normal
PL approach. Extra collar-rim band edges stay unmarked.
This avoids incorrectly promoting a completed surface
boundary loop to a proper boundary-marked tangle.

The positive result concerns the *chosen-coordinate candidate*,
not the source's full framed Fig12.10 isomorphism or isotopy.
The `MARKED_KNOT_SOURCE_UNVERIFIED` and A5.1
`UNVERIFIED_GEOMETRY` guards remain mandatory.


## A5.2e2b4c1 — actual local ambient tetrahedral plumbing slit

[Source Figure 12.5 local 3D product-ball slit](docs/local-plumbing-slit-a5-2e2b4c1.md)
uses the paper's **two transverse directed arrows** to select
an oriented rectangular patch that is *actually two faces of
the existing PL ribbon mesh*. It triangulates an exact-rational
local product 3-ball around that patch into 12 **positive-volume
tetrahedra**, then independently replays a cut along the middle
two-triangle disk, duplicating its four vertices into distinct
upper/lower 3D boundary sides. A typed inverse reglue checks
reconstruction of the exact original tetrahedral incidence,
oriented cut faces and physical coordinates.

This is the **first actual local 3D ambient-cell slitting
operation** in the program, but **NOT** the source's complete
S³ slit into actual four-punctured Conway tangle components,
nor a source-identified Figure 12.5 patch or global isotopy.
The proof status is `LOCAL_SLIT_REGLUE_ONLY`; original
`UNVERIFIED_GEOMETRY` remains unchanged under issue #52.


## A5.2e2b4c2a — exact disjoint atlas of local 3D slits, NOT global S³

[The n-site tetrahedral slit-atlas construction](docs/disjoint-3d-slit-atlas-a5-2e2b4c2a.md)
now places **all n real local 3D product charts** on the SAME
verified chosen-coordinate band and independently checks their
physical separation, 3D volumes, candidate-port correspondence,
both local cut-boundary spheres, and exact inverse reglue.
Even for zero-weight sectors, chart gap is **at least 1/2** in
exact rational coordinates; with `h≤1/4` all charts lie above
both the Conway balls and the physically marked knot K.

This deliberately exposes the missing *global* theorem: these
are **n disconnected product 3-balls**, which slit into **2n
unpunctured local 3-balls**. Neither the local cut disks nor
their boundary spheres contain a K puncture. They therefore
are **not** the source's *four-punctured* Conway boundaries,
and cannot be treated as a recovered global tangle. The typed
status `DISJOINT_LOCAL_ATLAS_NOT_GLOBAL_MANIFOLD` and the
original `UNVERIFIED_GEOMETRY` remain enforced. Issue #54
still requires a globally connected ambient manifold, real
source plumbing patch positions and four actual knot piercings
per source Conway boundary.


## A5.2e2b4c2b — connected exact 3D tetrahedral roof ambient (not S³)

[The connected physical roof-slab construction](docs/connected-roof-ambient-slab-a5-2e2b4c2b.md)
bridges A5.2e2b4c2a's n separated 3D slit charts with **real
face-compatible oriented tetrahedral material**. One exact PL 3-ball
now contains all n source-candidate oriented roof patches as
two-sided **internal mesh faces**, with a separate coordinate-only
3D checker for signed tetrahedral volumes, interior face orientation,
Euler-2 closed boundary, connectivity and each vertex's *2-disk/
2-sphere link*. Its 2n−1 longitudinal slabs have 12(2n−1)
tetrahedra and exact volume `4h·(x_last−x_first)`.

However, this is only a **connected roof slab** far above the
physically embedded proper knot K, NOT an ambient triangulation
of the full source S³ or tangle manifold. Detaching a pair
of tetrahedral **dual faces** is not yet a complete edge-split
PL manifold slit; the chosen midplane disks have **zero K
piercings**, so no source four-punctured Conway boundaries
have been reconstructed. The typed flags
`CONNECTED_ROOF_SLAB_NOT_SOURCE_AMBIENT` and
`DUAL_FACE_BARRIERS_NOT_PL_MANIFOLD_SLIT` preserve the
source's `UNVERIFIED_GEOMETRY` requirement.


## A5.2e2b4c2c — independent geometric K punctures on existing Conway surfaces

[The exact rational K–triangle intersection audit](docs/physical-conway-piercings-a5-2e2b4c2c.md)
derives every *physical* transverse segment/triangle hit using
source-candidate Conway 2-spheres' actual oriented mesh faces,
including exact clipping to reject tangential or coplanar
contacts. **Only after** computing these spatial intersections
does it compare the points with the 4n named punctures on the
n already-existing cubes. It verifies every hit is a genuine
proper K terminal with nonzero outward face-normal and that
all 2n source-candidate knot arcs avoid the solid cubes.

This is a necessary future source-plumbing checker, but the
audited physical cubical Conway spheres **predate the global
slit** and were *not generated by* the Figure 12.5 operation.
The connected roof 3-ball of A5.2e2b4c2b is disjoint from K
and produces **zero** source 4-punctured boundaries. Never
conflate these objects: `EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED`
and `UNVERIFIED_GEOMETRY` remain mandatory.


## A5.2e2b4c2d1 — a FULL edge/vertex-aware slit of the connected roof 3-ball

[The complete rational PL midplane slit](docs/full-roof-midplane-pl-slit-a5-2e2b4c2d1.md)
replaces the earlier incomplete **tetrahedral dual-face barrier** with
a true two-sided cut along the entire separating z=4 disk in the
connected roof slab. Every cut-plane **vertex, edge and triangle**
has distinct LOWER/UPPER mesh identities and opposite outward
normals. An independent 3D PL checker reconstructs two genuine
3-balls, each with Euler-2 boundary and valid per-vertex link,
and verifies disk Euler-1, exact volume conservation and inverse
tetrahedral regluing.

**Only the slab is cut.** This full disk includes auxiliary
connector strips, is NOT the source's Figure 12.10 plumbing-patch
set, and it has **zero K piercings**: K lies outside the slab.
No source four-punctured Conway sphere, ambient framed isotopy or
R± certification follows. The flag
`FULL_ROOF_PL_SLIT_NOT_SOURCE_CONWAY` and existing
`UNVERIFIED_GEOMETRY` restrictions remain mandatory.


## A5.2e2b4c2d2a — actual K as an edge of an oriented tetrahedral 3-ball

[The exact knot-conforming 3D tube constructor](docs/k-conforming-tetrahedral-tubes-a5-2e2b4c2d2a.md)
now places **real physical subsegments of the two outer return
K rails** as bona fide 1-simplices in a 12-tetrahedron PL
3-ball apiece. Their center-to-center mesh edge has an independently
checked 4-cycle link. The separate manifold validator checks signed
volumes, oriented face incidence, a closed Euler-2 boundary and every
3D vertex link. Exact rational segment/solid clipping verifies each
tube is disjoint from **every other physical K segment** and every
removed Conway cube; the two tubes are mutually disjoint.

This is a deliberately **local** knot-conforming 3D meshing
primitive. It does NOT mesh the entire K, ribbon or source-sphere
surfaces into one connected ambient 3-manifold, and cannot certify
new four-punctured Conway spheres, source-isotopy or R±.
`LOCAL_K_EDGE_CONFORMING_NOT_GLOBAL_AMBIENT` and
`UNVERIFIED_GEOMETRY` remain in force under Issue #60.


## A5.2e2b4c2d2b — exact arbitrary-direction tetrahedral K subedge stars

[The rational K-segment frame and 3D tube kernel](docs/rational-arbitrary-k-edge-stars-a5-2e2b4c2d2b.md)
generalizes the two horizontal return-rail tetrahedral tubes in
PR #61 to **any selected existing physical source-candidate K
segment**, including transverse vertical rails and signed braid
crossings with nonzero z heights. Rational cross-product frames
construct a real oblique tetrahedral 3-ball (10 physical vertices,
12 positively oriented tetrahedra) with an exact K interior mesh
edge and audited four-cycle link. A separate exact 3D segment–
parallelepiped clip checks EVERY other K segment; a rational
separating-axis checker proves closed-solid clearance from EVERY
Conway cube. Radius shrinks deterministically until all geometry
passes or the bounded search reports UNVERIFIED.

These stars cover only proper **interiors** of original K edges.
They do not yet include vertex/corner junctions or form a connected
3D ambient whose mesh conforms to all K edges, ribbon and Conway
surfaces. They are source-unverified; Issue #62 and A5.1
`UNVERIFIED_GEOMETRY` remain open.


## A5.2e2b4c2d2c — true K vertices with TWO conforming 3D tetrahedral edges

[The exact knot-corner octahedral 3-ball](docs/rational-k-corner-stars-a5-2e2b4c2d2c.md)
adds a real convex 7-vertex, 8-tetrahedron **PL 3-ball
centered on an actual 2-valent knot vertex**. Both incoming
and outgoing physical K subsegments are independently checked
as tetrahedral mesh edges with their own circular 4-cycle
links. It works for bent and straight continuation vertices
with exact rational cross-product frames, positive tetrahedral
volumes, all vertex links, and exact 3D convex segment clipping
against EVERY other K edge and removed Conway cube.

These are still **independent local corner stars**, not
geometrically attached to PR #63's interior edge tubes or
assembled into one complete knot-containing 3-manifold.
No source-generated Conway four-punctured boundary, framed
isotopy, or R± theorem follows. `LOCAL_CORNER_STAR_NOT_WHOLE_K_AMBIENT`
and `UNVERIFIED_GEOMETRY` remain mandatory.


## A5.2e2b4c2d2d1 — genuine positive-area outgoing caps on K corners

[Truncated exact rational K corner 3-balls](docs/truncated-k-corner-cap-a5-2e2b4c2d2d1.md)
resolve the specific **zero-area apex** obstruction identified
in Issue #65: the PR #64 octahedral corner star is truncated
halfway along its outgoing actual K edge and independently
retriangulated into an embedded **11-vertex / 16-tetrahedron
PL 3-ball**. The newly exposed physical attachment cap is a
**five-vertex, four-triangle, positive-area disk** (χ=1)
with a real outgoing K mesh edge ending at the cap center.
Independent 3D manifold, oriented volume, exact cap-edge
incidence and incoming/outgoing knot edge-link audits validate
the construction. Its exact volume is 15/16 of the original
convex corner octahedron.

**A positive-area cap is not yet a 3D tube attachment.** An
explicit interface transition to the independently meshed
PR #63 knot tube remains open; no global knot ambient or
source Figure 12.10 Conway 4-puncture cut has been certified.
The `TRUNCATED_K_CORNER_CAP_NOT_GLUED` and
`UNVERIFIED_GEOMETRY` statuses remain fail-closed.


## A5.2e2b4c2d2d2 — exact 3D tetrahedron collision proof and conforming one-corner K weld

[The rational tetrahedral collision/transition implementation](docs/welded-k-corner-edge-a5-2e2b4c2d2d2.md)
independently computes the **physical convex intersection of every
pair of actual tetrahedra**, rejecting geometric crossings and
overlaps not explained by a literal shared mesh vertex, edge or
whole face. The verifier uses exact rational supporting-plane
triple intersections, not floating-point approximations or
combinatorial topology alone.

The construction joins PR #66's **exposed five-vertex diamond
corner cap** to PR #63's **five-vertex square edge-tube cap** by
four coherently triangulated three-dimensional sector prisms,
adding 12 transition tetrahedra. It welds the original 16 corner
and 12 tube tetrahedra by **actual physical coordinates**,
then independently checks all 3D manifold vertex links,
complete nonlocal tetrahedral intersections, source K 1-skeleton
continuity, exterior physical knot-segment contact and Conway
solid clearance. This is one knot-corner/outgoing-edge attachment,
**not** the entire K-conforming source manifold.

The current objects remain `ONE_K_CORNER_EDGE_GLUE_NOT_SOURCE_AMBIENT`;
full source Figure 12.10 ambient slitting, four NEW Conway
punctures, signed/pure framed isotopy and R± remain open.


## A5.2e2b4c2d2d3a — positive-area INCOMING K-corner attachment disk

[The incoming octahedral K cap construction](docs/incoming-k-corner-cap-a5-2e2b4c2d2d3a.md)
supplies the dual of PR #66's outgoing cap: a **11-vertex,
16-tetrahedron physically embedded rational PL 3-ball** clipped
at the *actual incoming knot ray*, with a **five-vertex,
four-triangle exposed positive-area disk** centered on the
same physical K segment. Both original K ray edges remain
genuine mesh 1-simplices with independently checked link
cycles. Exact volume is 15/16 of the original K-corner
octahedron, and all vertex links and cap face incidences
are independently audited.

This corrects the **zero-area next-corner apex** obstruction
preventing PR #68's outgoing K-edge tube from being attached
at its *far* end. **The incoming cap is not yet physically
welded** to that far square tube interface. The full two-ended
K edge and source Figure 12.10 ambient remain open in #69.
`INCOMING_K_CAP_NOT_TUBE_GLUED` and
`UNVERIFIED_GEOMETRY` remain mandatory.


## A5.2e2b4c2d2d3b — actual TWO-ENDED knot edge in one PL 3-ball

[The verified two-ended K-edge construction](docs/two-ended-k-edge-a5-2e2b4c2d2d3b.md)
joins the existing physical PR #68 corner/outgoing-edge 40-tet
3-ball to the actual *next K vertex's* PR #70 positive-area
incoming diamond cap. A four-sector 12-tetrahedron transition
welds the far tube square to that cap, producing one intended
**32-vertex / 68-tetrahedron PL 3-ball** with an actual
five-edge K 1-skeleton chain between two original K vertices.
Both attachment fans must become genuine internal faces.
An independent exact convex 3D collision checker audits all
2,278 tetrahedron pairs; all vertex links, signed volumes,
source K edge parameter order and closed Conway clearance
are separately checked.

This remains one source-candidate K **edge between two
corners**, NOT a mesh of every K arc or a connected
Bonahon–Siebenmann tangle (M,K). No full ribbon/Conway
face-conforming ambient, source Figure 12.5 slit producing
four NEW punctures, framed isotopy or R± proof follows.
`TWO_CORNER_K_EDGE_NOT_WHOLE_SOURCE_AMBIENT` and
`UNVERIFIED_GEOMETRY` remain mandatory.


## A5.2e2b4c2e1 — ONE true K corner core with TWO positive-area caps

[The doubly truncated shared K-corner core](docs/doubly-truncated-shared-k-corner-a5-2e2b4c2e1.md)
replaces the invalid strategy of overlapping an incoming-only and
outgoing-only 3D corner star at the same actual knot vertex.
One exact-rational convex octahedral K corner is instead clipped
at BOTH K rays, using a guarded cut fraction 3/4. The resulting
candidate contains **15 genuine Q³ vertices / 24 tetrahedra**
and **two disjoint, physically exposed 5-vertex/4-face cap disks**,
with distinct true K center-to-cap 1-simplex edges. Its exact volume
is **63/64 of the original octahedron**, and all 276 tetrahedral
pairs undergo independent rational convex intersection tests.

This is a **shared corner with two actual positive-area interfaces**,
NOT YET a two-edge conforming K neighborhood. True welding of both
adjacent K tubes, whole-K source ambient and source Figure 12.10
four-NEW-puncture reconstruction remain OPEN under #72.
`DUAL_CAP_K_CORE_NOT_EDGE_GLUED` / `UNVERIFIED_GEOMETRY`
remain fail-closed.


## A5.2e2b4c2e2 — two actual K edge tubes through ONE shared core

[The two-direction physically welded shared-corner construction](docs/shared-corner-two-k-rays-a5-2e2b4c2e2.md)
joins incoming and outgoing proper source-K segment tubes to the TWO
disjoint, positive-area cap disks of the **same single shared K corner**
from PR #73. A four-prism, 12-tetrahedron rational transition is
placed on each side. The intended connected mesh has **35 exact Q³
vertices / 72 oriented tetrahedra**, with all eight fan attachment
triangles welded as genuinely internal physical faces; the original
knot centerline is a conforming, strictly ordered 1-skeleton chain
across the actual original shared K vertex. Independent PL manifold
vertex-link checks, **2,556** exact physical tetrahedron-pair
intersections, source-K strand exclusion and Conway solid clearance
are mandatory.

Unlike a naive union of incoming-only/outgoing-only corner
octahedra, this construction never superposes two overlapping
3D solids at the same true K vertex. It is only a LOCAL two-ray
neighborhood, **not** complete source K, global ribbon/Conway
ambient or source Fig12.5 new four-punctured slitting/isotopy.
`TWO_K_RAYS_LOCAL_NOT_SOURCE_AMBIENT` and
`UNVERIFIED_GEOMETRY` stay fail-closed.


## A5.2e2b4c2e3 — ONE K tube between TWO doubly truncated true vertices

[The reusable shared-tube/two-corner PL 3-ball](docs/one-shared-k-tube-dual-corners-a5-2e2b4c2e3.md)
joins the outgoing four-triangle diamond cap of one genuine
source-candidate K vertex to a **SINGLE real original K tube**,
then joins its far square cap to the next K vertex's incoming
diamond. The 24+12+12+12+24 actual tetrahedra form an
independently checked **40-vertex / 84-tetrahedron PL 3-ball**.
Every one of **3,486 physical tetrahedron pairs** is checked
for nonconforming collision, and the source K segment is a
five-edge properly ordered tetrahedral **1-skeleton chain**.
The two unused opposite K-corner caps remain exposed for
future larger-chain construction.

This avoids the subtle mistake of unioning two individually
certified PR #75 neighborhoods at successive K vertices:
both would contain the SAME original K tube material, so
their tetrahedra would **overlap in positive volume**.
No whole-K or source-ribbon/Conway ambient result follows,
and no Fig. 12.5 newly slit four-K-punctured Conway sphere,
framed isotopy or ring R± has been established. The
`ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT` and
`UNVERIFIED_GEOMETRY` gates remain strict.
