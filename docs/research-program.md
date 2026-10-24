# Catalan syntax versus arborescent isotopy

**Status:** precise counting theorems for presentations, a negative control for raw Catalan factors, and separately labeled conjectures. No general arborescent-link isotopy classifier has been implemented or proved.

## Objects

A *boundary-fixed two-string tangle* has four individually marked endpoints; isotopies fix the boundary. A *cyclic Montesinos link closure* is a different object: its symmetries may include cyclic permutation and reversal of rational pieces. Never identify their equivalence relations.

Fix a finite alphabet L of m>=2 pairwise distinct rational slopes in (0,1). Let Syn(k,L) be the full rooted ordered binary trees with k ordered leaves in L, where internal nodes denote the chosen horizontal gluing operation. Syntactic size is k-1 internal nodes, not crossing number. A *reduced weighted arborescent tree* must later be defined independently of parenthesization and with actual sector and cyclic-order data.

## Proven presentation-level statements

**Theorem 1 (raw syntax).** |Syn(k,L)| = C_(k-1) m^k, with C_n = binom(2n,n)/(n+1). This follows from the bijection with plane full binary trees and independent leaf labels.

**Theorem 2 (associativity counterexample).** Quotient the raw syntax by formal associativity (X+Y)+Z = X+(Y+Z). The quotient is the set of m^k ordered leaf words. For every fixed word, its fiber of parenthesizations has size C_(k-1), exponential in k, not polynomial. This is a formal quotient and does not assert that distinct words yield distinct isotopy classes. The claim that raw Catalan presentations have polynomial fibers is therefore false.

**Theorem 3 (dihedral presentation quotient).** Further quotient ordered words by cyclic rotation and reversal **as a model for cyclic link-closure presentations**. The orbit count B_m(k) is given exactly by Burnside's lemma:

```
B_m(k) = [R_m(k) + F_m(k)] / (2k)
R_m(k) = sum_(d|k) phi(d)*m^(k/d)
F_m(k) = k*m^((k+1)/2)                 if k is odd
       = (k/2)*(m^(k/2+1) + m^(k/2))  if k is even
```

For fixed m>=2, only the identity rotation contributes m^k; every other symmetry fixes at most O(m^((k+2)/2)) words. Hence B_m(k) = m^k/(2k) + O(poly(k)*m^((k+2)/2)). The generating series G_m(x)=sum_(k>=3) B_m(k)x^k has a singular term -(1/2)log(1-mx); the remainder extends analytically beyond |x|<1/m. It is **not algebraic over Q(x)** because algebraic functions have local Puiseux expansions without logarithms.

These are theorems about syntactic equivalence classes, **not** an unconditional identification with Montesinos-link isotopy classes. Using Bonahon–Siebenmann's link classification needs an exact match of closure conventions and exceptional cases.

The second quotient has fibers of size at most 2k over ordered words, but by then the Catalan factor is gone. A polynomial-fiber theorem stated at the wrong stage does not restore Catalan growth.

## Conjectural positive program

Choose in advance a **specific finite** admissible family of reduced arborescent weighted planar trees with truly branching vertices, local weights, cyclic orders and precise topological equivalence moves.

* **A — canonicalization:** a stated reduction system terminates and is confluent on that subclass, or another effectively checkable canonical certificate agrees with the relevant isotopy theorem. Exceptional cases and marked symmetries must be handled, not waved away.
* **B — restricted tree-type growth:** for an explicit irreducible, aperiodic, genuinely branching finite-local grammar whose generating equations satisfy a square-root singularity theorem, the reduced-form counts have a_n ~ K*mu^n*n^(-3/2), with explicit K,mu>0. A finite grammar alone does **not** imply this exponent: path-like languages, periodicity, and group quotients may yield other laws.
* **C — topological quotient:** independently prove a uniform or polynomial fiber bound for a precise map from reduced **normal forms** to isotopy classes, if that map is not already one-to-one by a valid classification theorem.
* **D — geometric comparison:** compare reduced-tree size with crossings, bridge number or other geometric complexities only after the syntactic classification is established.

**Boundary question:** characterize the transition from arborescent to polyhedral Conway decomposition pieces. Non-arborescent does not automatically mean the global decomposition graph is nonplanar or nonbinary. No faster-than-Catalan dichotomy is assumed.

## Any unbounded sequence?

Drop the unconstrained universality claim. Fix the grammar and a specific state-sum/skein/TQFT-style representation **before** asking which invariant or counting sequences are realizable. The dihedral quotient already shows that finite local syntax can have a non-algebraic quotient series.

## Validation ladder

S0: verify Catalan raw syntax. S1: verify reassociation to words. S2: verify bracelet Burnside against independent orbit enumeration (this PR). S3: verify the applicable Montesinos-link theorem and exceptions before treating counts as isotopy classes. S4: implement actual reduced weighted-tree calculus and test soundness/completeness. S5: perform singularity analysis and independent geometric verification. **The present PR stops at S2.**

## Stage A1 — completed exact combinatorial platform

The [A1 weighted-plane-tree note](weighted-tree-a1.md) gives a literal
rotation/sector/free-bond representation, a marked-plane-tree *presentation*
isomorphism code, strict nonexceptional screens for (W),(F),(S), and an
independently checked finite binary/ternary branching grammar. The exact
syntax generating series is F=m+x(F²+F³), with a cycle-lemma coefficient
formula by number of ternary forks. **Neither passing W/F/S nor matching
presentation codes establishes an arborescent isotopy equivalence.**
A2 (full (P)/(N) conditions, exceptional cases, flips/ring/boundary
plumbing moves and certified quotient computation) remains open.

## A2.1 — rigorous partial plumbing-move quotient (implemented)

[Proof and exact Burnside results](flip-subgroup-a2-1.md): F1/F2
reflect cyclic orders on weighted planar trees and, by the bipartition
of any connected tree, generate a group with at most four elements.
This subsystem preserves W/F/S, admits exact bounded-orbit enumeration,
and has **algebraic** quotient generating series with the same
`μ^n n^(-3/2)` growth as the syntactic grammar, up to a factor `1/4`
for fixed leaf alphabet size at least two. Distinguish this theorem
from the still-open growth problem for the **full** Bonahon–Siebenmann
equivalence and from pointwise boundary-fixing of tangles.

Issue #4 tracks the remaining sourced, topology-certified moves,
positive/negative canonicity and exceptions. The A2 completion bar is
**not** met by this partial group action.

## Stage A2.2 — boundary-label repair and strict F3′ (implemented)

[Full definitions, source audit, exact theorem and tests](v4-tangle-moves-a2-2.md).
The §12.7 *tangle* calculus requires a V4 element on each free bond,
and F1/F2 **change** those elements. The ordinary unlabelled A2.1
quotient Q is a **pair-presentation** result and not a
boundary-parametrized tangle quotient. With exactly one marked free
root bond and all four labels, F1/F2 acts freely; the finite labelled
quotient count is exactly the A1 raw syntax count a_n.

A source-restricted W-preserving F3′ moves a single nonzero sector
weight across an internal edge, and when odd reflects selected
subtree vertices and updates free-bond V4 labels by distance parity.
The complete F3′/ring equivalence is **not** implemented.

Importantly, pair-canonicity condition F is **absent** from tangle
canonicity in §12.7.2. Do not apply the A1 pair W/F/S screen to certify
a tangle. This stage corrects the mathematical interpretation of A2.1
without retracting its valid syntactic unlabelled enumeration result.

## A2.3 — evidence correction and partial-move exploration

[A2.3 note](partial-move-closure-a2-3.md) records that A2.2's hosted
Python tests initially **failed after merge**, from an incorrectly
constructed `Fork` fixture. The fix is part of A2.3.

The exact structural result A2.3.1 demonstrates that **strict F3′ is
vacuous on the original A1 leaf-weight-only grammar**, including its
F1/F2 orbit. This blocks treating A2.2's implemented F3′ as new
quotient evidence for the A1 growth law. The source-scoped finite
move explorer supplies shortest, replayable certificates on fixed
indexed graphs, including an actual weighted-internal-vertex example.
The results are proofs about formal presentation moves; full tangle
canonicality, modified ring moves, and isotopy quotient remain open.

## Stage A2.4 — non-vacuous strict-F3′ presentation quotient

[The A2.4 proof and exact results](root-pivot-a2-4.md) addresses
A2.3's obstruction by fixing a finite grammar with precisely one
nonzero sector at a root fork and zero-weight descendant forks.
This gives a nontrivial, reversible strict-F3′ move on a genuinely
weighted internal vertex. The exact syntax series are
`H=p*x*(3F²+4F³)` and `J=p*x*(F²+F³)=p*(F-m)` for the quotient
under root-sector strict F3′ moves **alone**. Independent exhaustive
normalization checks small cases. The same square-root exponent
survives this *restricted presentation* quotient.

This is **not** the full Theorem 12.21 quotient: independent checks of
positive/negative tangle canonicality, exceptional forms, R± and
boundary-parametrized isotopy remain outstanding.

## A2.5 — necessary ring signature, not an R± rewrite

[The A2.5 note](ring-signature-a2-5.md) audits the modified R±
diagrams and derives an exact *unit vertex weight* obstruction for
all partial F1/F2/F3′-reachable grammars omitting ±1. A permissive
detector returns necessary ±1 center/±2-tip graph motifs. Source
diagram sector locations and direct R± replacement remain future
verification obligations; the present detector is not sufficient
to declare a topologically allowed ring move.

## Corrected A2.6 — ring-number abstraction (not R± deletion)

[Primary-source image audit and restricted abbreviation](correct-ring-abbreviation-a2-6.md)
retracts the false claim in closed, unmerged PR #11 that R± deletes a
three-vertex cap. The §12.4 diagram relocates the ring; §12.5 creates
a **different ring-numbered abbreviated tree notation**, where
non-overlapping caps may be erased *as an abstraction* and counted at
the base. The new software provides explicit RingCount data and
refuses source-described adjacent-ring exceptions. Full R± diagram
moves, topology soundness/completeness, and true quotient growth remain
open in issue #4.

## A3 — exact bivariate counting laws for ring-numbered syntax

[Proof and reproducible exact data](ring-number-grammar-a3.md):
assign bounded independent nonnegative ring numbers to vertices of
the A1 2/3-ary fork grammar, as a restricted family of **abbreviated
presentations** in the sense of §12.5. Its bivariate generating series
is `F=m+x(1+...+y^R)(F²+F³)` and exact coefficient factorization
is `a_n(m)[y^q](1+...+y^R)^n`. This is a predeclared structural
multiplicative count, not the earlier vacuous arbitrary sequence
factor. For fixed q, the polynomial prefactor changes to
`n^(q-3/2)`; summing all q keeps `n^(-3/2)` but changes
exponential growth to `(R+1)*mu_m`. None of these claims yet
identifies isotopy classes under the full modified ring calculus.

## A3.1 — local grammar compatibility, beyond independent labels

[Exact coupled equations, independently checked counts and proofs]
(local-ring-constraints-a3-1.md): impose the predeclared
**synthetic** restriction that adjacent forks cannot both carry a
positive ring number. Two-type recursions produce an algebraic
bivariate series and invalidate the naive independent
`(R+1)^n a_n` factorization already at two forks. Despite this
change, for any *fixed* ring grade q the leading coefficient still
behaves as `a_n n^q/q!` by an edge-exclusion union bound. Summing all
ring grades instead yields a strict exponential upper-rate gap
from the independent model, using a matching of size ≥(n−1)/7.
No claim is made that this local exclusion reproduces §12.5's
actual exceptional adjacent-ring theorem or defines isotopy classes.

## A3.2 — raw ring-location ambiguity can be exponential

[The A3.2 proof and exact exhaustive witnesses](ring-placement-fibers-a3-2.md)
construct a corolla with `d` identical V4 free legs (one marked root)
and `r` identical disjoint ring caps. Exactly
`binom(d+r-1,r)` marked plane presentation codes collapse to a
**single ring-numbered abbreviated presentation code**. At `d=r=n`,
this is `~4^n/(2sqrt(pi n))`. This is a **negative control for
unreduced presentation-to-abbreviation fibers** and demonstrates
why ring-location normalization must precede
any polynomial-fiber conjecture. It is NOT a calculation of tangle
isotopy classes and does not refute a bounded-valence reduced
normal-form hypothesis.

**A3.2 stronger obstruction (bounded valence):** a path of `k` base
vertices carrying one disjoint ring cap and one free bond at every
base, with all degrees ≤4 and all weights in `{0,sign,2*sign}` for either fixed sign,
has a selected binary family of exactly `2^k=2^{N/4}` distinct marked
plane-tree presentation codes but one ring-numbered abbreviated code
(the unrestricted fiber may be larger). This strictly strengthens
the single high-valence corolla counterexample and proves **bounded
valence alone does not imply polynomial raw-to-abbreviation fibers**.
Only *after* a genuine ring-location normalization could a polynomial
fiber theorem be plausible for a canonical grammar. Full topology
classification remains separately open.

## A3.3 — square-root universality for a nonindependent local grammar

[The A3.3 proof](guarded-critical-growth-a3-3.md) applies the
Flajolet–Sedgewick/Meir–Moon smooth implicit-function schema to the
A3.1 synthetic positive-ring adjacency exclusion after setting y=1.
The two-state polynomial specification reduces to
`G=H(x,G)` where `A=m+x(G²+G³)` and
`H=A+R*x(A²+A³)`. For fixed `m>=1,R>=1`, the coefficientwise
solution is positive and aperiodic, with a unique positive
nondegenerate characteristic singularity. Thus total guarded
**syntax** counts satisfy `K*mu^n*n^(-3/2)`, with algebraic
`mu` specified by the critical system. An exact degree-11
polynomial and rational root bracket are included for m=2,R=1.
Nothing here transfers the exponent to genuine arborescent isotopy.

## A4.1 — finite formal ring-position connectivity, not isotopy

[Proof, guards and replay certificates](ring-relocation-a4-1.md):
A2.6's lost ring locations can be explored using a **conservative
formal adjacent interchange** of a strict ring cap's full attachment
stem with a neighboring ordinary bond at a zero-sector base.
Unlike the retracted PR #11, the full cap remains present.
The restricted §12.5 abbreviated code is verified unchanged.
On the A3.2 corolla the presentation-move graph connects all
`binom(d+r-1,r)` placements; on the bounded-valence k-base
path it connects all `2^k` placements. These results are exact
**syntactic** connectivity theorems and deliberately do **not**
certify that each interchange is the published topological R±
move for an arbitrary noncanonical presentation.
The remaining direct plumbing and boundary proof is tracked in #16.

## A4.2 — terminating and confluent strict *formal* ring normal form

[The inversion-measure proof](ring-normal-form-a4-2.md) orients the
A4.1 signed zero-sector, nonoverlapping-ring *presentation*
interchange `RO↔OR` as `RO→OR`. At every base the total
inversion number decreases by exactly one, and each anchored
incident list reaches the **unique** stable form
`O_1…O_o R_1…R_r`. Deterministic typed replay certificates
verify every step. For structurally identical cap subtrees and
a fixed rooted residual skeleton, the resulting marked
normal-form code corresponds exactly to the restricted §12.5
ring-numbered abbreviation code, up to the declared
presentation-isomorphism relation. These are **syntactic**
equivalence classes, not verified topological tangle isotopy.

## A4.3 — corrected rooted cyclic isomorphism before formal ring sorting

[Rooted cyclic-normalization proof and tests](root-aligned-ring-normal-form-a4-3.md):
A4.2's anchored ring bubble-sort was dependent on the chosen first
**storage** slot at each cyclic vertex; that first slot is not a
mathematical plane-tree invariant. A4.3 first cyclically rotates
the bond **and sector** lists at every vertex to the unique incoming
edge toward a marked root free bond, as a **presentation
isomorphism**, then replays the typed strict ring slides.
The resulting marked rooted normal-form code is invariant under
cyclic slot starts, bond renaming and vertex indexing, retaining
all V4 labels. Still not an isotopy theorem or full R± calculus.

## A4.4 — commutation of §12.5 abbreviation and F1/F2 V4 tangle flips

[Primary-source proof and exhaustive small-instance tests]
(abbreviated-flip-equivariance-a4-4.md): Bonahon–Siebenmann §12.7.1
requires that F1/F2 flips also transform every free-bond V4 label.
The restricted §12.5 abbreviation **commutes exactly with F1 and
F2** on disjoint, uniformly signed ring caps with zero-sector
attachment bases and no adjacent-ring exception. F2 pivot
correspondence uses a named *surviving free bond*, not a stale
integer vertex index. The induced four-element action on abbreviated
ring-numbered data has an exact identity-root V4 gauge, preserving
ring multiplicities and boundary labels. The proof and regression
fixtures are for the source-defined **formal presentation moves**;
neither genuine 3D isotopy nor complete F3′/R± canonicality follows.

## A4.5 — guarded F3′ naturality under ring-number abbreviation

[The A4.5 commuting-square proof](abbreviated-f3prime-naturality-a4-5.md)
handles a strict source-modified F3′ move across a **surviving**
weighted-tree edge. Parity-sensitive branch reversal and ξ/η
V4 transformations agree before/after source-style §12.5
abbreviation because pendant ring erasure does not change graph
distances between surviving vertices. It requires unchanged strict
recognition of every ring center, zero-sector bases, an explicitly
surviving free-bond pivot, and unchanged per-base ring counts;
unsupported cases fail closed. Nontrivial even/odd pivot witnesses
verify the result as **indexed presentation data**. Full F3′,
source canonicality exceptions and independent tangle isotopy
are still open.

## A4.6 — joint formal F1/F2 plus ring-placement quotient

[The A4.6 two-layer proof](joint-flip-ring-quotient-a4-6.md)
integrates the source V4-labelled F1/F2 four-group with the
rooted, terminating formal ring-location normalizer. A unique
identity-root-V4 flip group variant is followed by cyclic parent
alignment and stable `RO→OR` ring placement sorting. The resulting
marked presentation key is invariant under both flip generators
and all admitted strict formal ring slides. On a fixed residual
skeleton with indistinguishable cap subtrees it agrees with the
F1/F2 root-gauged §12.5 ring-number abbreviation criterion.
Small explicit joint orbits have `4*binom(d+r-1,r)` corolla
states and `16*3^(k-2)` bounded-valence path states. These are
rigorous **formal indexed presentation** orbits, not
Bonahon–Siebenmann topology-certified isotopy classes.

## A4.7 — partial mixed flip/F3′/ring-slide graph with replay witnesses

[The source-scoped A4.7 note](mixed-partial-calculus-a4-7.md)
constructs a typed finite-index **formal** transition graph for
V4-labelled F1/F2, strict F3′ on surviving pivot edges, and
the restricted A4.1 ring-attachment interchange. Each accepted
source-defined flip must pass the independently implemented
§12.5-abbreviation commuting square; every intermediate
presentation is root-cyclic-aligned and must preserve the
per-base signed ring-count vector. The graph supports replayable
shortest positive witnesses and honest `complete=False` cutoffs.
A signed seven-vertex example verifies one actual strict-F3′ /
**formal** ring-slide diamond; it does not establish general
confluence or genuine tangle isotopy. The F3′-extended calculus
has **no asserted complete normal form**; A4.6's stable formal
key remains limited to its earlier subsystem.

## A4.8 — necessary rooted weighted-ring invariant and evidence types

[The A4.8 exact proof](mixed-necessary-invariants-a4-8.md)
adds a rooted, **unordered** residual weighted-ring skeleton
fingerprint to the partial A4.7 formal graph. Every F1/F2,
strict surviving-edge F3′ and *formal* ring-position step
preserves the rooted tree incidence, the total integer
weight at each vertex and its §12.5 ring multiplicity.
The fingerprint deliberately ignores cyclic sector
locations and mutable V4 free labels. Different fingerprints
thus give a valid *restricted formal nonreachability* witness,
but equal fingerprints do not imply any equivalence.
The search API now distinguishes **replayable reachability,
necessary-invariant exclusion, complete *indexed* target
nonreachability, and unresolved cutoff**. Even fully completed
indexed BFS is not a topological-isotopy decision procedure,
or even a quotient under arbitrary vertex/bond renaming.

## A4.9 — exact partial *formal* quotient modulo rooted isomorphism

[The A4.9 proof and independent regressions](rooted-class-quotient-a4-9.md)
replace fixed-index A4.7 BFS's storage-sensitive equality
with the faithful rooted V4-labelled plane-tree presentation
code. All coded F1/F2, strict surviving-edge F3′ and
**formal** ring-slide successors respect marked presentation
isomorphism: moving source vertex/edge indices and cyclic
starting slots commutes with enumerating available moves.
One representative per code therefore suffices for exact
bounded BFS on *presentation-isomorphism classes*, with a
qualified `complete=True` closure result for the declared
partial move rules. Positive witnesses combine an executable
typed move trace and a checked final rooted isomorphism;
a complete formal-class exclusion still does **not** establish
different 3D tangle isotopy classes or source-complete
canonical equivalence.
