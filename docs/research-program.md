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
