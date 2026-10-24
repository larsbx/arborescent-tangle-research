# A3.2 — Exact exponential fiber of ring-numbered abbreviation

**Status:** explicit negative control for **combinatorial presentation
fibers**, not a theorem about marked-tangle isotopy.

The previous Catalan-skeleton discussion proposed polynomial fibers
from reduced decorated trees to topological classes. That question is
not decided here. This note establishes that a *different map*, from
unreduced locations of ring subtrees to source-style ring-numbered
abbreviated trees, can have **exponentially large fibers** when the
base's degree increases with presentation size. It therefore warns
against transferring any polynomial bound between these distinct
representation levels.

## Primary source

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots and
the Classification and Symmetries of Arborescent Knots*,
§12.4 printed p.155 (modified ring relocation) and §12.5 printed
pp.157–158 (abbreviation retaining ring multiplicities), 2016:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The source says the **location** of an attached ring subtree is
irrelevant *under its ring-move and canonicity hypotheses*. The
abbreviated representation erases suitable ring subtrees while
counting their multiplicities at the surviving attachment base.
Its exceptional case of adjacent rings is not part of this example:
all our ring subtrees are pairwise disjoint and have the same base.

This experiment uses only the *safe, highly restricted* abbreviation
implementation already merged in PR #12. In particular it does
**not** implement the topology-preserving R± location move or prove
the input presentations are canonical arborescent tangles.

## A fixed explicit family

Let `d>=3` and `r>=0` be integers. Construct a rooted weighted
plane-tree with one **zero-sector base** of bond valence `d+r`:

- `d` free boundary bonds `f0,...,f(d-1)`, all labelled `ι=0`
  in V4. The bond `f0` is the marked root.
- `r` identical positive ring caps attached to the base by separate
  internal bonds. Each ring cap has exactly three vertices: a central
  trivalent vertex with sectors `(0,+1,0)` in the order
  `(stem,tip1,tip2)`, and two terminal vertices weighted `+2`.
  For the negative family replace every cap sign by minus.
- The base's incident cyclic order is arbitrary subject to the
  marked `f0` and a choice of interleaving the `r` caps among the
  `d-1` remaining indistinguishable free bonds. All base angular
  sectors are zero.

All caps are internally attached, disjoint and satisfy the
restrictive §12.5 abbreviation recognizer. The base has bond
valence at least three, so it survives any abbreviation. Cyclic
starting-point rotations are accounted for by using `f0` as the
anchor, and edge names are **not** treated as invariants.

## Theorem A3.2.1 — exact marked-presentation fiber

Write `P(d,r)` for the number of distinct **marked weighted-plane
presentation-isomorphism codes** in the family, and `B(d,r)`
for the number of distinct source-style abbreviated presentation
codes. Then

```
P(d,r) = binom(d+r-1,r)
B(d,r) = 1.
```

**Proof.** Anchor the cyclic incident list at the marked `f0`.
What remains is a length-`d-1+r` word on two different symbols:
`F` (an ordinary free leg with V4 label zero), appearing `d-1`
times, and `R` (a three-vertex ring cap of the chosen sign),
appearing `r` times. Because the cap has nonzero weighted vertices
and the free leg does not, distinct such words have different
rooted plane-tree codes. Conversely all permutations within each
symbol type are forgotten by graph/bond renaming. The number of
two-symbol words is `binom(d+r-1,r)`.

In the §12.5 abbreviation every cap is erased **as a change of
representation**, and the base is assigned ring number `r`.
After erasure the remaining base has exactly the same `d` free
legs and zero sectors in every source presentation, and its ring
number is the same. Thus all abbreviations have one identical
rooted ring-numbered code. □

This proof is about presentation encodings only, and does not assume
the source's complete tangle-isotopy theorem.

### Explicit examples

| Free legs d | Ring caps r | Distinct marked locations | Abbreviated codes |
|--:|--:|--:|--:|
| 3 | 1 | 3 | 1 |
| 3 | 2 | 6 | 1 |
| 4 | 2 | 10 | 1 |
| 5 | 3 | 35 | 1 |
| 8 | 8 | 6435 | 1 |

The implementation `placement_fibers.py` constructs all exact
weighted-plane presentations with V4 boundary labels and **executes
the source-restricted abbreviation**. Tests compare closed binomial
counts against independent exhaustive rooted codes for
`3<=d<=6`, `0<=r<=3` and both signs.

## Corollary A3.2.2 — exponentially large raw syntax fibers

Choose `d=r=n>=3`. Then

```
P(n,n) = binom(2n-1,n)
       = (1/2) binom(2n,n)
       ~ 4^n/(2 sqrt(pi*n)).
```

This is exponential in the number of caps and in the source
presentation's vertex count `3n+1`. **Therefore there cannot be
a polynomial uniform fiber bound for the map from these raw
marked placements to this abbreviated syntax language when base
valence grows without bound.**

**Important qualifications:**

- This does **not** refute a polynomial-fiber conjecture for an
  *already reduced, fixed bounded-valence* grammar. Our base's
  valence `d+r=2n` grows with `n`; that violates a bounded-valence
  hypothesis.
- This is not a count of different topological tangles mapping to
  the same tangle; it is precisely a count of **presentation codes**
  lost by changing representation. The source's own equivalence
  theorem would need to be invoked, with all hypotheses checked,
  before a topological equivalence statement can be made.
- Ring numbers only at specific vertices and free-port V4 labels
  are indispensable in the abbreviated code; ignoring them would
  create further artificial identification.

## Next: isolate bounded-valence canonical subclasses

This negative control narrows the research objective: predeclare
a canonical weighted planar-tree subclass with bounded valence,
precise tangle W/S/P/N and exceptional hypotheses, and an explicit
R± relocation calculus. Count representations **after** §12.5
abbreviation, not their discarded raw ring placements; test the
full topological quotient separately. The exact A3/A3.1 generating
functions remain statements about syntax, independent of these
topology obligations.
