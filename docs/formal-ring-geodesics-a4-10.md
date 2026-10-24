# A4.10 — Gaussian-binomial rank and exact geodesics in the formal ring-slide graph

**Status:** an exact enumerative and algorithmic theorem on the
**restricted A4.1 formal ring-position calculus**, not an isotopy theorem,
a general R± implementation, or a quotient of the full tangle category.

## Motivation and scope

The corrected Bonahon–Siebenmann source (§12.4, printed p.155)
**relocates** a signed ±2—±1—±2 ring subtree; §12.5, printed
pp.157–158, separately introduces the **lossy** ring-numbered
abbreviated representation. Our A4.1 model conservatively formalizes
only adjacent swaps of a strict ring cap stem R and an ordinary
nonanchor bond O at the same zero-sector base.

Primary source:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

A4.2/A4.3 established a terminating normal form for this **formal**
rewrite, and A4.6 combined it with the source F1/F2 V4 action.
A4.7/A4.9 extended **partial** formal search to strict F3′
without a complete normal form. A4.10 returns to the precise A4.1
subsystem and derives the **entire shortest-path and distance-rank
distribution**, rather than only its orbit cardinality.

## Theorem A4.10.1 — Gaussian-binomial inversion enumerator

Fix a rooted zero-sector base with `o≥0` ordered ordinary
attachments O and `r≥0` indistinguishable strictly recognized
ring cap stems R (both sequences internally ordered as stored).
The rooted cyclic anchor is fixed, and only adjacent
`RO ↔ OR` interchanges are allowed; `RR` and `OO` are not moves.

For each interleaving word `w`, define

```text
I(w) = # {(i,j) : i<j, w_i=R, w_j=O}.
```

This is exactly its A4.2 shortest distance **to the sorted stable
normal form `O^o R^r`**. The exact rank-generating polynomial is

```text
G_(o,r)(q) = sum_interleavings q^I(w)
           = [o+r choose r]_q.
```

Here the right-hand side is the **Gaussian binomial polynomial**,
not the ordinary binomial count. It obeys the exact recurrence

```text
G_(o,0)=G_(0,r)=1
G_(o,r)=G_(o,r-1) + q^r G_(o-1,r)
```

for positive o and r.

**Proof.** Every legal adjacent `RO→OR` step removes exactly one
R-before-O inversion, so the A4.2 bubble sort realizes a path of
length `I(w)`. Any path to the sorted form must change the
inversion count by `I(w)), with each step changing it by at most
one; hence the distance is minimal.

To enumerate the words, partition them according to their last
symbol. A word ending R contributes the same inversion count as
its prefix of o O's and r−1 R's, giving `G(o,r−1)`. A word
ending O has exactly r inversions contributed by that final
ordinary bond after all r rings; its prefix contains o−1 O's and
r R's, giving `q^r G(o−1,r)`. This proves the recurrence. □

At q=1 the polynomial specializes to the ordinary orbit
cardinality `binom(o+r,r)`. Its degree is or, and its coefficients
are symmetric: reversing a binary word while exchanging the
initial/final sorted positions sends `I → or−I`. In particular
the maximum normal-form distance is or.

## Theorem A4.10.2 — bounded-valence path rank distribution

Consider A3.2/A4.1's k-base path with k≥2. The rooted cyclic
attachments at either endpoint admit one ordinary nonanchor
bond and one ring cap, giving `G_(1,1)=1+q`.
Each interior degree-four base has two ordinary nonanchor bonds
and one cap, giving `G_(2,1)=1+q+q²`.

The attachment bases are independent under the strict *formal*
slides (no cap moves between bases). Thus the full distance
distribution to the common formal stable-normal representative is

```text
H_k(q) = (1+q)^2 (1+q+q²)^(k-2).
```

**Proof.** The global inversion statistic is the sum of each
base's local inversion counts, and every permitted move acts at
exactly one base. Therefore the generating polynomial factors
as the product of local `G(o,r)` polynomials. □

At q=1,

```text
H_k(1) = 4*3^(k-2),
```

recovering (and refining) the A4.1 full formal-orbit count.
The polynomial is palindromic of degree `2k-2` and gives
**the number of formal presentations at every exact distance**
to the common stable normal form. For example,

```text
H_2(q) = 1+2q+q²
H_3(q) = 1+3q+4q²+3q³+q⁴
```

These are formal presentation counts, *not isotopy strata*.

## Theorem A4.10.3 — shortest path between arbitrary supported states

Let U and V be two **same-index, same-bond, same-V4-label**
rooted presentations with identical cap stems and their order,
identical ordinary bond subsequences at every affected base,
identical entire nonbase vertices, and identical zero angular
sectors at every affected base. They may differ only in the
**interleaving** of R and O incident bonds. Both are cyclically
aligned toward the same marked free root before comparison.

For ring stem `R_j` in its fixed order at a base, let
`b_U(j)` and `b_V(j)` denote the number of ordinary bonds
preceding it in U and V respectively. Then

```text
d_formal(U,V) = sum_bases sum_j |b_U(j)-b_V(j)|.
```

**Proof (lower bound).** An allowed adjacent interchange moves
one indexed ring stem past exactly one ordinary bond. It changes
exactly one integer `b(j)` by ±1. Hence at least the displayed
L¹ distance many steps are required.

**Proof (construction).** Scan each target anchored local list
from left to right. If the next stored bond differs from the
target bond, find that target bond farther right and successively
move it left through the intervening bonds. Preservation of
both relative subsequences ensures every crossed bond is of
the opposite R/O type, so every interchange is legal.
After each settled prefix, continue to the next target slot.
All bases can be treated independently because the moves
preserve the underlying incidence tree. No ring/ordinary pair
crosses more than necessary, so the constructed path performs
exactly the displayed L¹ distance. □

Implementation `formal_slide_geodesic` returns a
`FormalSlideGeodesic` certificate containing the source and
target rooted cyclic alignments (presentation-isomorphism steps,
**not** slides), the exact integer minimal distance, and a
replayable sequence of typed A4.1 `RingSlide` moves. Its verifier
checks both alignments, replays every formal step, and recomputes
the independent L¹ lower bound.

Inputs that differ by vertex indices, internal bond names,
nonbase weights, boundary V4 labels, ring sign/multiplicities
or ordinary branch order **fail closed**, since this theorem
concerns the explicit fixed-index strict ring-slide graph.
A4.9 separately treats **presentation-isomorphism classes**
and does not silently conflate them with typed ring slides.

## Independent exhaustive verification

The standard-library implementation
`arborescent/ring_geodesics.py` computes Gaussian-binomial
coefficients by the displayed exact integer recurrence.

`tests/test_ring_geodesics.py` compares:

1. The Gaussian polynomial coefficients to exhaustive concrete
   A3.2 corolla interleavings, computing the **actual**
   A4.2 inversion count at each presentation;
2. The complete A4.1 BFS orbit on bounded-valence k-paths
   (k=2..6, both signs) to the entire coefficient vector
   of `(1+q)^2(1+q+q²)^(k−2)`;
3. The returned geodesic's length to the **independent shortest
   distance computed by A4.1 BFS** for multiple corolla and
   bounded-valence path families, and all ordered pairs of
   4-free-leg/two-cap presentations;
4. Rooted cyclic slot rotations at source/target (zero-cost
   presentation isomorphism, not a ring move); and
5. Strict invalid-input guards, including changed ring sign,
   cap count, bond naming, sector data and integer arguments.

Reproduce:

```sh
python -m unittest discover -s tests -v
```

## Topological boundary

A4.10 is a precise metric theorem for a **restricted formal
presentation calculus**. It does not independently establish
that its adjacent ring slides are actual source-valid R± tangle
isotopies. It adds no general F3′, exceptional adjacent-ring
normalization, true W/S/P/N canonicality or independent
three-dimensional Conway/plumbing proof. Therefore neither its
Gaussian-binomial coefficients nor geodesics enumerate or
certify ambient-isotopy classes of tangles. The full
topological acceptance remains open in issue #16.
