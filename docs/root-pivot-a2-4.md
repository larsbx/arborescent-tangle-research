# A2.4 — A genuinely non-vacuous strict-F3′ root-pivot grammar

**Result type:** exact formal presentation theorem, not a knot or tangle
isotopy-class enumeration. This stage repairs the *mathematical*
vacuity identified in [A2.3](partial-move-closure-a2-3.md), without
silently asserting positive/negative topological canonicity.

## Source and scope

The published Bonahon–Siebenmann calculus allows (F3′) to relocate
one nonzero angular sector weight across an adjacent internal bond
while reversing the affected branch at odd distances **if the weight is
odd**, and modifying its free-bond V4 labels. This is the precise
nonexceptional restricted implementation already present in
`arborescent/tangle_labels.py`.

Reference: Bonahon–Siebenmann (2016), §12.7.3 (printed p. 167,
PDF index p. 178), Theorem 12.21 on the same page, and §12.7.2 for the
extra *tangle* canonicality hypotheses.
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

This experiment deliberately uses **only strict F3′ at the root**,
not F1/F2 or ring moves. The root has one free bond `f0` labelled
`ι`, and that free-bond label is unchanged by root-F3′ because it is
outside every moved child branch. There are no other free bonds.

## Fixed, bounded-local grammar

Fix independent finite alphabets:

- `L`: m distinct terminal leaf weights, each integer with |w| ≥ 2;
- `P`: p distinct nonzero root sector weights, each integer with |a| ≥ 2.

Let `F=F_m(x)` be the A1 ordered 2/3-ary zero-internal-weight fork
grammar with weighted terminal leaves:

```
F = m + x(F²+F³).
```

The new grammar selects exactly one weighted **root fork**, with
either two or three ordered child subtrees from F. Its root has
respectively three or four angular sectors (including the one
adjacent to free bond f0). Exactly one sector carries a root weight
`a in P`; every other sector at the root and at descendant forks
has weight zero. All terminal vertices carry their original A1
weights.

Unlike A1, a generated root can carry nonzero weight at valence ≥ 3,
so the strict F3′ transport is **genuinely applicable**.

The grammar satisfies the sector sparsity condition **W**. We have
not verified the complete *tangle* (S) and (P)/(N) canonicity
conditions or exceptions. The strict pair-only F condition cannot
be used to exclude this grammar as a tangle.

## Theorem A2.4.1 — exact weighted-root syntax generating function

If `h_n(m,p)` counts all presentations, with n internal fork nodes
and one distinguished weighted root fork, then

```
H_(m,p)(x) = p*x*(3F(x)² + 4F(x)³).
```

**Proof.** A binary root has three angular sectors, so for each
root weight it admits three placements. A ternary root has four
sectors, giving four placements. The ordered child subtrees are
independent; the root contributes one internal node.

## Theorem A2.4.2 — restricted F3′ root-sector normal form

For a presentation with root child bonds `e0,...,e(r-2)`, cyclic
order `(f0,e0,...,e(r-2))` and root valence `r = 3 or 4`, the
angular sectors `0,1,...,r-1` form a **path** under allowed strict
root-F3′ moves:

```
sector 0 --e0-- sector 1 --e1-- ... --e(r-2)-- sector r-1.
```

Across a child bond, the weight `a` is moved one sector left or
right; if odd, the source-prescribed parity reflections are applied
to *that* child subtree. Each move is reversible by the opposite
transport. No other internal vertex in this grammar carries a
nonzero sector (other than terminal leaves, which have only one
sector), so **no additional strict F3′ moves exist** on descendants.

**Consequently** each root-sector orbit contains precisely **r**
different **indexed** presentations, one for each sector, and
there is a unique `sector=0` representative. Distinct normalized
root presentations have different rooted plane-presentation codes
unless related by ordinary presentation isomorphism. An odd-weight
subtree reflection does not leave the grammar: it permutes the
children of a zero-weight descendant fork.

This gives the exact restricted quotient generating series:

```
J_(m,p)(x) = p*x*(F(x)² + F(x)³)
           = p*(F(x)-m).
```

**Important:** this quotient is by **root-sector strict F3′ only**.
Adding F1/F2, R±, or broader canonical-tree equivalences requires a
*new* proof. The quotient theorem is not an isotopy classification
statement.

## Exact results

For `m=p=2`, including one positive odd and one positive even
root weight (P={3,4}):

| Forks n | All weighted-root presentations | Strict root-F3′ orbits |
|--:|--:|--:|
| 1 | 88 | 24 |
| 2 | 1,440 | 384 |
| 3 | 30,816 | 8,160 |
| 4 | 752,256 | 198,528 |
| 5 | 19,853,568 | 5,229,312 |
| 6 | 551,725,056 | 145,135,104 |

The counting code evaluates `H` and `J` with exact integer
coefficients from the A1 recursion. An independent small-n test
generates every decorated weighted-root presentation, executes
source-guarded strict F3′ transports, normalizes to sector zero,
and confirms the same numbers.

## Asymptotic consequence for *presentations*

For fixed `m≥2,p≥1`, the A1 series has
`[x^n]F ~ K_m μ_m^n n^(-3/2)`, with `K_m, μ_m>0`.
Since `J=p(F-m)` exactly,

```
[x^n]J ~ p K_m μ_m^n n^(-3/2).
```

Thus the **restricted non-vacuous F3′ quotient** retains precisely
the A1 square-root exponent and exponential growth constant.

For `τ=F(ρ)`, the ratio of raw to reduced leading coefficients
is

```
lim(n→∞) h_n / j_n = (6τ+12τ²)/(2τ+3τ²),
```

which lies strictly between 3 and 4. It need not equal a
constant finite orbit size because binary and ternary roots have
different numbers of sector positions.

## Reproduce

```sh
python -m unittest discover -s tests -v
python -m arborescent.root_pivot --max-n 6 --leaf-labels 2 --root-labels 2
```

## Remaining blockers to an isotopy theorem

A2.4 only counts weighted-tree **syntax modulo a specified strict
F3′ subsystem**. To use Bonahon–Siebenmann Theorem 12.21 for true
isotopy classes, the allowed tree subclass would also need verified
*positive or negative tangle canonicality* including exceptional
cases, the full F1/F2/F3′ and modified R± equivalences, and the
theorem's boundary and irreducibility hypotheses. The actual
diagrammatic plumbing map has not been independently implemented.
These remain open in issue #4.
