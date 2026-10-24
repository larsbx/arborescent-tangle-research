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

## A3.2 — raw ring-location ambiguity can be exponential

[The A3.2 proof and exact exhaustive witnesses](ring-placement-fibers-a3-2.md)
construct a corolla with `d` identical V4 free legs (one marked root)
and `r` identical disjoint ring caps. Exactly
`binom(d+r-1,r)` marked plane presentation codes collapse to a
**single ring-numbered abbreviated presentation code**. At `d=r=n`,
this is `~4^n/(2sqrt(pi n))`. This is a **negative control for
unreduced presentation-to-abbreviation fibers** and demonstrates
why bounded-valence and reduction must be named before making
any polynomial-fiber conjecture. It is NOT a calculation of tangle
isotopy classes and does not refute a bounded-valence reduced
normal-form hypothesis.
