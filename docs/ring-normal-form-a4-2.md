# A4.2 — Stable normal form for the restricted **formal** ring-slide calculus

**Status:** constructive terminating and confluent presentation rewrite
for the strictly guarded A4.1 subsystem. **NOT** a canonical arborescent
tangle, knot-isotopy classifier, or complete Bonahon–Siebenmann R± move.

## Source and scope

The source Bonahon–Siebenmann monograph §12.4 (printed p.155)
draws modified R± *ring relocations* that retain a signed
`±2—±1—±2` ring cap. Section 12.5 (printed p.157) separately defines
ring-numbered **abbreviated trees** which erase suitable ring caps in
a change of representation and record the count at the base. Actual
topological equivalence applies within the source's positively or
negatively canonical classes, including exceptional cases.

Primary source:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

A4.1 [formal ring-position slides](ring-relocation-a4-1.md) allow only
a strictly recognized signed ring cap and swap its attachment stem
with an adjacent ordinary incident bond at the **same** base when:

* all sectors at the base are zero;
* a specified noncap anchor occupies base slot zero, which no slide
  may cross;
* the central cap has the exact `(0,±1,0)` wedge and matching ±2 tips;
* the chosen sign is uniform and every candidate ring signature is
  recognized; no overlapping/adjacent ring centers or unknown
  arithmetic rearrangements are allowed; and
* free-bond V4 labels and all vertex weight totals remain unchanged.

The present stage **adds no topology claims**. It analyzes the
well-defined A4.1 abstract rewrite, which is an explicitly restricted
subset of possible presentation rearrangements.

## Theorem A4.2.1 — termination and a unique stable normal form

At a base, call any strictly recognized signed ring cap stem `R`,
and every other nonanchor bond `O`. Ignore particular bond names only
for this word analysis (the software retains them in certificates).

A4.1 allows `R O ↔ O R` for adjacent nonanchor slots; it explicitly
does **not** allow `R R ↔ R R` or changing the order of two ordinary
bonds. Define the inversion number

```text
I(w) = # {(i,j) : i<j, w_i=R, w_j=O},
I(tree) = sum_base I(base word).
```

**Theorem.** The right-moving orientation `RO → OR` terminates
after exactly `I(tree)` moves, and the result is the unique stable
normal form in the indexed-letter class:

```text
[anchor, O_1,...,O_o, R_1,...,R_r]
```

where the relative orders among ordinary bonds and among ring
stems agree with the input. Distinct reductions of the same valid
input have the same result, and every finite A4.1 undirected slide
orbit reaches that stable normal form.

**Proof.** Each `RO→OR` swap decreases `I` by exactly one; it
cannot change any other `R/O` pair's order, so termination is
immediate. A word admits no decreasing move iff no `R` appears
before an `O`, i.e. iff every `O` precedes every `R`. All swaps
preserve the subsequence of `O` bonds and the subsequence of `R`
stems, so only one such indexed sequence exists. Moves at different
bases commute (each changes only that base's incident rotation);
therefore the global normal form is unique. Every undirected
elementary slide preserves those two subsequences, so applying the
normalization to either endpoint yields the same result. □

The deterministic implementation scans incident slots from right
to left and bubble-sorts rings rightward; its certificate includes
every typed `RingSlide` action. `RingNormalForm.verify()` replays
each move and independently checks that `I` decreases by exactly
one, and finally reaches zero. Insufficient explicit move budgets
raise `ValueError`, never silently produce a partial normal form.

## Corollary A4.2.2 — exact finite presentation orbits

At a fixed base with `o` ordered ordinary nonanchor bonds and
`r` indistinguishable structurally identical rings, there are
exactly

```text
binom(o+r,r)
```

different rooted plane presentation codes in its formal A4.1 orbit.
For several independent bases these factors multiply.

**Proof.** Allowed moves are adjacent swaps between `R` and `O`.
They generate every interleaving of the two preserved sequences.
Choose the `r` positions for the `R` symbols; all choices are
reachable from the stable form, and the rooted planar order records
each distinct interleaving. Distinct indices of otherwise
indistinguishable ring caps disappear only at the
presentation-*isomorphism code*, not from the replay trace. □

This includes the previously verified families:

* **Corolla:** `o=d-1` and `r` caps at one base, giving
  `binom(d+r-1,r)`.
* **Bounded-degree path:** k bases with one cap each, two options
  at either end and three at each interior base, giving the exact
  full orbit size `4·3^(k-2)`. The earlier `2^k` family is a proper
  subfamily when `k≥3`.

This is an *abstract* orbit-counting theorem. It does not state that
the orbit consists of tangles belonging to a particular canonical
topological category.

## Theorem A4.2.3 — abbreviated key agrees with stable position

Within this same strict, anchor-fixed, all-zero-sector,
indistinguishable-cap subclass, the rooted marked normal-form
presentation-isomorphism key is constant on every formal slide orbit.
The §12.5-style restricted abbreviated code is also constant,
because all caps remain attached to the same base in every slide,
their counts are unchanged, and the ordinary branch order, all
residual sector weights and all V4 boundary labels survive.

For presentations with **the same rooted residual skeleton and
fixed-cap type**, equality of restricted abbreviated codes is
equivalent to equality of stable-normal-form rooted codes.

**Proof.** Abbreviation discards only the positions of these strictly
matching caps, recording their number at each base. The stable
normalizer relocates every such cap after the ordered ordinary
bonds. Hence the residual skeleton and ring counts uniquely
determine the stabilized rooted code, and conversely one can erase
the stable-position caps to recover the abbreviated code. No
other information is forgotten in the *declared subclass*. □

This equivalence is **not** an independent theorem equating such
codes to topological isotopy classes. Bonahon–Siebenmann's complete
W/S/P/N canonicality, F1/F2/F3′, exceptional ring cases, V4 boundary
transport and genuine R± diagram soundness are still required.

## Reproducible computational obligations

Implementation: `arborescent/ring_normalizer.py` exports
`normalize_ring_positions`,
`ring_position_inversions`,
`strict_ring_normal_key`, and the replayable
`RingNormalForm`.

Tests independently verify the fixed-endpoint `I` certificate,
all small corolla interleavings (both signs), full A4.1 BFS
orbits, bounded-valence path orbits, sign/wedge/anchor/sector
fail-closed guards, and preservation of every marked free-bond
V4 label. They also check that no false normal form is reported
when the explicit move budget is insufficient.

```sh
python -m unittest discover -s tests -v
```

## Topological acceptance still open

A4.2 supplies a candidate **normalization interface** for a narrow
formal model. It does *not* prove the constructed presentation swap
is a source-valid R± move for arbitrary noncanonical plumbing data.
Issue #16 remains open for independent 3D Conway/plumbing
construction and figure-level isotopy certificates; exceptional
ring-center handling and tangle W/S/P/N canonicality; boundary V4
transport; and complete topology-level quotient computations.

The A3.x algebraic `n^(-3/2)` asymptotics are about generated
**presentations** only and are not transferred by this theorem
to genuine tangle-isotopy classes.
