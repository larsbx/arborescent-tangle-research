# A3.1 — finite-local compatibility, exact counts, and sharp growth bounds

**Status:** exact, independently testable theorem on a **synthetic finite-state
grammar of ring-numbered presentations**. It is deliberately *not* an
assertion that Bonahon–Siebenmann’s adjacent-ring exception equals the
simple graph-adjacency prohibition below. The actual source §12.5 has
additional topological and exceptional hypotheses; this guarded grammar
is a classification experiment, not an isotopy-class count.

## Fixed grammar and algebraicity

Start with the A3 bounded ring-number grammar: each rooted, ordered
2/3-ary fork carries an integer `r_v∈{0,...,R}`, leaves carry no
ring number, and each leaf is labelled from an alphabet of size `m`.
Impose one additional **local** rule, fixed *before* counting:

> Along every edge joining two fork vertices, at least one of the two
> ring numbers must equal zero.

Equivalently, the positive-ring fork vertices form an independent set
in the induced tree on internal vertices. This resembles a hard-core
constraint, but is **not** claimed to be the source's canonicality
predicate for actual arborescent tangles.

Let `Z(x,y)` enumerate a leaf or a tree whose root fork has no rings;
let `P(x,y)` enumerate a tree with a positive-ring root fork; let
`G=Z+P` and `D_+(y)=y+...+y^R`. Then

```
Z = m + x (G² + G³)
P = x D_+(y) (Z² + Z³)
G = Z + P.
```

**Theorem A3.1.1.** The exact bivariate counting series `G` is
algebraic over `Q(x,y)`.

**Proof.** A ring-free root allows either type of child, while a
positive-ring root requires every child root to be ring-free.
Binary/ternary child orders are distinguished. The displayed system
is a finite polynomial system with a unique coefficientwise formal
power-series solution at `x=0`. Eliminating the auxiliary
series yields an algebraic equation for `G`. □

The implemented integer convolution recurrences match the equations,
and independent enumeration of explicit `Abbreviated` presentations
plus an edge-by-edge ring-label filter checks the first sizes.

## Exact counts and a falsifiable break from independent decoration

For `m=2,R=1`, the count vector `[x^n y^q]G` at fixed `n` is:

| forks n | q=0,1,...,n | total guarded | independent A3 total |
|--:|---|--:|--:|
| 0 | 2 | 2 | 2 |
| 1 | 12, 12 | 24 | 24 |
| 2 | 192, 384, 0 | 576 | 768 |
| 3 | 4080, 12240, 4080, 0 | 20,400 | 32,640 |
| 4 | 99264, 397056, 297792, 17856, 0 | 811,968 | 1,588,224 |

**Immediate counterexample to independent-label factorization:**
two forks are adjacent in a two-fork tree. They cannot both carry a
positive ring number, so the `n=2,q=2` coefficient is **zero** rather
than the unrestricted value **192**. Thus the exact factorization
`a_n(m)(R+1)^n` is **not universal** for bounded-local grammars.

For `m=2,R=2,n=2`, the guarded vector is
`(192,384,384,0,0)`, total **960**; the independent A3 grammar has
`(192,384,576,384,192)`, total **1,728**.

## Theorem A3.1.2 — fixed ring grade has universal leading power

Let `c_(n,q)^(R)(m)=[x^n y^q]G` for fixed `m≥2,R≥1,q≥0`.
For *each fixed q*,

```
c_(n,q)^(R)(m)
= a_n(m) * (n^q/q! + O(n^(q-1))),
```

uniformly over all A1 tree shapes. In particular, using the proven A1
square-root asymptotic `a_n~K_m μ_m^n n^(-3/2)`,

```
c_(n,q)^(R)(m) ~ (K_m/q!) μ_m^n n^(q-3/2).
```

**Proof.** Each such tree has `n` internal fork vertices and
exactly `n−1` edges joining internal forks. For assignments giving
exactly one ring to each of `q` distinct forks, there are
`binom(n,q)` subsets. At most
`(n−1)binom(n−2,q−2)=O(n^(q−1))` include an adjacent pair, by
the union bound over the `n−1` internal edges (for `q<2` the
subtraction is zero). This gives the uniform lower bound
`n^q/q!+O(n^(q−1))`.

For an upper bound, ignore adjacency and the per-fork cap `R`:
the `q` rings can be distributed among `n` distinguishable forks
in at most `binom(n+q−1,q)=n^q/q!+O(n^(q−1))` ways. Sandwich
the admissible count between these two uniform bounds, and multiply
by `a_n(m)`. □

**Implication:** bounded local exclusion **does not change** the
leading polynomial prefactor at any *fixed* total ring grade.
This is subtler than the summation-over-grades result.

## Theorem A3.1.3 — exponential loss when summing all grades

Let `C_n(m,R)=sum_q c_(n,q)^(R)(m)` and fix `R≥1`.
Define

```
gamma_R = (1+2R)/(R+1)²,   0 < gamma_R < 1.
```

Then

```
a_n(m) <= C_n(m,R)
       <= a_n(m) (R+1)^n gamma_R^((n-1)/7).
```

**Proof.** The induced fork graph is a tree on `n` vertices, with
maximum degree at most 4: a fork has at most 3 children and 1 parent.
Greedily selecting disjoint edges of this graph produces a matching
of size `h≥(n−1)/7`, since each selected edge eliminates at most
`2·4−1=7` edges. Each matched pair would have `(R+1)²`
independent labels without the rule, but has only `1+2R` legal
choices because two positive ring numbers are forbidden. Every other
fork has at most `R+1` choices. Thus any fixed tree has at most
`(R+1)^n gamma_R^h` admissible assignments. Sum over A1 tree
shapes. The lower bound assigns ring number zero everywhere. □

Consequently, for the exponential growth upper rate,

```
limsup_(n→∞) C_n(m,R)^(1/n)
   <= μ_m (R+1) gamma_R^(1/7)
   <  μ_m (R+1).
```

This is a rigorous **strict exponential gap** from independent
decoration, without claiming a full square-root asymptotic for
this coupled grammar or a theorem about link isotopy.

## Interpretation, sources and obligations

The exact model answers the earlier universality criticism:
bounded-local tree labels produce constrained, algebraic growth
sequences, but the naive independent correction factor is destroyed
by even one local compatibility constraint.

The ring-number vocabulary is motivated by Bonahon–Siebenmann,
*New Geometric Splittings of Classical Knots*, §12.5, printed
pp.157–158:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The **adjacency guard is a synthetic research restriction**,
not a transcription of the section's exceptional two-ring
configuration. Before any tangle-topology interpretation, the
program still needs diagram-verified R± moves, the source's
full exceptional canonicality, ring-numbered representation
completeness, and an independently validated plumbing map.

Run:

```sh
python -m unittest discover -s tests -v
python -m arborescent.ring_constraints --max-n 5 --m 2 --max-rings 1
```

The regression checks compare the polynomial dynamic program against
a completely separate exhaustive tree-generation-and-filter route,
and verify strict loss from independent decoration for `n>=2`.
