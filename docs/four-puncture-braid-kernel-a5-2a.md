# A5.2a — Four-ended braid diagrams and independently replayable RII isotopies

**Status:** first genuine *diagrammatic* (not merely weighted-tree)
isotopy certificates in the repo, for an intentionally small class:
**monotone two-strand braids with four individually marked endpoints**.
This stage does **not** implement the Bonahon–Siebenmann §12.5
Figure 12.29 *atomic tangles with rings* or its actual R± plumbing
moves. The source-specific topological classification remains open.

## The precisely defined geometric category

Use a rectangle `[0,2]×[0,n+1]` with marked endpoints
`NW=(0,0)`, `NE=(2,0)`, `SW=(0,n+1)`,
`SE=(2,n+1)`. A two-braid diagram is a sequence
`w=(ε_1,...,ε_n)`, `ε_i∈{+1,−1}`.

At level i its two strands cross in projection at
`(1,i+1/2)`. We record **which physical strand passes over**:
- `ε_i=+1`: the strand on the left immediately above
  level i passes **over** the right strand;
- `ε_i=−1`: the strand on the right immediately above
  level i passes **over** the left strand.

Both strands exchange left/right positions at each crossing.
Their original identities remain anchored at the individually
marked NW and NE endpoints. The north/south matching therefore
depends only on parity: for even n, NW–SW and NE–SE; for odd n,
NW–SE and NE–SW.

`Braid2Diagram.crossing_events()` computes every over/under
pair; `projection_paths()` gives both complete polygons from
north to south, meeting at the correctly indexed crossing
projection coordinates. `render_svg()` produces a small
deterministic inspection drawing with a gap in the underpass
and the overpass drawn last. These are bona fide planar
**two-braid projections**, unlike A5.1's abstract sewing graph.

**Important category restriction:** neither the SVG renderer
nor the exact braid-word checker can interpret arbitrary
Conway/tangle diagrams, winding ring components, nonbraided
rational tangles or 3D plumbing seams. The monotone two-braid
category is a proper subclass of the topology required by
Fig.12.29.

## Theorem A5.2a.1 — constructive relative-endpoint RII certificates

Define `σ` as the positive elementary crossing and
`σ⁻¹` as the negative one. A local adjacent inverse pair

```
    ... σ σ⁻¹ ...  <---->  ...
    ... σ⁻¹ σ ...  <---->  ...
```

can be removed/inserted by a **Reidemeister-II** move inside
the rectangle, fixing the four endpoints. Crucially, the
signed over/under data is part of the crossing event; this is
not simply erasing a pair from a symbol string.

Let `w(D)=Σ ε_i` be the braid exponent.

**Proposition:** For two such diagrams `D,D′` with
`w(D)=w(D′)`, `prove_two_braid_equivalence(D,D′)`
constructs a **replayable sequence of relative-boundary
Reidemeister-II diagram moves**, taking one diagram to the other.

**Proof:** Read the word left to right, removing adjacent
inverse crossings until none remain. Every removal is a legal
Reidemeister-II move and decreases diagram crossing count by
exactly two without changing `w` or the endpoint matching.
A word with no adjacent inverse symbols in the two-letter
alphabet `{σ,σ⁻¹}` consists entirely of `σ` or entirely
of `σ⁻¹`; its length and sign are determined by `w`.
Thus both input diagrams reduce to exactly the same canonical
word. Reverse the second reduction as RII insertions.
The resulting sequence takes the first diagram to the
second, and `Braid2IsotopyCertificate.verify()` replays
every edit, enforcing the pair signs, crossing positions,
matching and invariant at every stage. □

The familiar presentation of the **two-strand braid group**
`B₂≅Z` says that `w` is complete for *two-strand braid
isotopy*. We deliberately limit negative claims to that
category. A failed or mismatching exponent does **not**
constitute proof that arbitrary 3D tangles are nonisotopic
under different boundary conventions. The checker always
raises instead of issuing a negative topology certificate.

## Four labelled punctures and explicit, still abstract seams

Conway boundary components are 4-punctured spheres. Before
attempting a geometric seam, A5.2a fixes a typed **clockwise
puncture enumeration**

```
    (NW, NE, SE, SW).
```

The Klein four-group `V4` is represented explicitly by
permutations of these four indices:

```
    id   = (NW)(NE)(SE)(SW)
    ξ    = (NW NE)(SW SE)
    η    = (NW SW)(NE SE)
    ζ    = (NW SE)(NE SW).
```

They compose as bitwise XOR, and each element is its own
inverse. This is a **declared coordinate convention**
consistent with `ξ η = ζ`, **not independently verified**
to be the particular Fig.12.29 and §12.7 Conway boundary
parametrization adopted by Bonahon–Siebenmann.

`FourPunctureFrame` attaches this enumeration to a typed
A5.1 `Port`. `proposed_seam_puncture_matching`
checks that the named parent and child ports are endpoints
of an actual A5.1 `Seam`, then returns a four-point
bijection for a declared V4 element. The resulting
`AbstractPunctureMatching` is explicitly marked
`ABSTRACT_ONLY`: it does **not** specify an
orientation-reversing 4-punctured-sphere homeomorphism, framing,
strand seaming or the source's actual plumbing gluing.
`certify_seam_isotopy()` raises `NotImplementedError`.

Thus the repository now has **two distinct evidence layers**:
a genuinely checked RII-isotopy certificate for a narrowly
defined four-ended braid **diagram**, and a merely abstract
four-punctured boundary matching proposal for an A5.1
atomic sewing port. They must not be conflated.

## Independently checkable test coverage

`tests/test_conway_braid_diagrams.py` checks:

1. V4 acts *faithfully* on four labelled punctures;
   all 16 compositions match XOR and all maps are
   involutions with four distinct images.
2. One exactly labelled four-puncture frame per A5.1
   port; proposed seam maps are accepted **only** for real
   directed A5.1 parent/child seam endpoints. Unsupported
   or falsely certified gluing is rejected.
3. Full over/under strand identities, projection crossing
   coordinates, top/bottom endpoint matching and SVG drawing
   metadata for explicit positive, negative and zero-crossing
   examples.
4. **All** two-braid crossing words of length ≤6 (126
   nonempty words, plus the empty word) reduce with independently
   replayed RII deletions, correct parity/writhe and inverse
   edit reconstruction.
5. Pairwise positive equivalence witnesses built from two
   independent RII normalizing paths for words with the same
   exponent; words with unequal exponent produce **no
   overclaimed tangle nonisotopy certificate**.
6. Strict rejection of malformed generators, crossings,
   RII edits, boundary maps and fake 3D seam approvals.
7. A5.1 ring-numbered sewing plans remain marked
   `UNVERIFIED_GEOMETRY` even when an unrelated two-strand
   braid certificate succeeds.

Run:

```sh
python -m unittest discover -s tests -v
```

## Source comparison and the remaining geometry bridge

Primary source checked as **page images**:

Bonahon–Siebenmann,
*New Geometric Splittings of Classical Knots and the
Classification and Symmetries of Arborescent Knots*,
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

- **§12.4 printed p.155**, modified R± ring diagrams:
  a signed three-vertex ring subtree is **relocated**, not
  removed.
- **§12.5 printed p.157**: ring-number abbreviation erases
  certain caps as a **representation change**.
- **§12.5 printed p.158, Fig.12.29**: the *actual model
  atomic tangle with r rings*, drawn as a winding,
  multiply punctured piece with explicit twists `a_i`.
  **This is not a two-strand braid tile.**

What remains for the later stages of
[Issue #33](https://github.com/larsbx/arborescent-tangle-research/issues/33):
an independent exact **Fig.12.29 atomic tile** with both
through-strands and ring components, positive/negative
sector twists, a fully oriented 4-punctured Conway
seam homeomorphism, and explicit Reidemeister or
3D-isotopy witnesses for the actual R± diagrams
with V4 boundary parametrizations and source-canonical
W/S/P-or-N exceptions.

**No present A5.2a certificate may be attached to a
Bonahon–Siebenmann ring move or used as an isotopy-class
growth theorem.**
