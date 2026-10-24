# A5.2b — A four-puncture/V4 gluing map is **not a Conway seam**

**Status:** exact algebraic obstruction to an unsound geometric
shortcut. It is independent of our previous weighted-tree normal
forms and uses the topology of the **four-punctured sphere**.

**Result:** Even specifying an explicit bijection of all four
individually labelled punctures, including its V4 permutation,
**does not determine a Conway sewing homeomorphism up to isotopy**.
There are infinitely many nonisotopic puncture-fixing sphere
mapping classes that induce the **same identity permutation**
on all punctures. We therefore cannot turn the A5.2a finite
V4 permutation into a geometry certificate.

## Why this is the next source-critical obligation

Bonahon–Siebenmann §12.5, printed **p.158**, Fig.12.29,
provides a genuine model atomic tangle with `r` rings at
every weighted vertex. The associated tree's internal bonds
indicate which **four-punctured Conway sphere boundaries**
must be sewn. Section 12.7 carries explicit V4 free-boundary
parametrizations.

Our merged A5.1 `AtomicPlumbingPlan` retains all weighted-tree
sector and ring-number data and declares its missing geometry
`UNVERIFIED_GEOMETRY`. Merged A5.2a built real elementary
two-braid diagram RII witnesses, plus a typed *proposed* finite
V4 bijection on the four punctures of each A5.1 seam. However,
a four-point bijection is only a boundary *marking* datum,
not a complete 4-punctured-sphere **mapping class**.
This note supplies the precise obstruction.

Primary source visually checked:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf
(§12.5 printed p.158, PDF zero-based p.169;
§12.7 pp.166–167, PDF pp.177–178).

## The double branched cover of a four-punctured sphere

The two-fold branched covering
`T² → S²` with four branch points gives the
well-known pillowcase model of `S²_{0,4}`.
The four branch points correspond to the **2-torsion
points** of the torus:

```
   (0,0), (1,0), (0,1), (1,1)  in (Z/2Z)².
```

A linear torus homeomorphism represented by
`M ∈ SL₂(Z)` descends through its
deck involution `x ↦ −x` and acts on the
branch points by the matrix `M mod 2`.

**The pure congruence subgroup** is

```
Gamma(2) = { M ∈ SL₂(Z) : M ≡ I (mod 2) }.
```

Its elements fix all four labelled punctures individually.
On the quotient pillowcase the matrices `M` and `−M`
induce the same mapping class; the appropriate pure
projective group is `Gamma(2)/{±I}`.

The mapping class of the four-punctured sphere remembers
how an essential **slope class** is transported, not merely
the four marked puncture positions. Unoriented primitive
slopes are represented by coprime pairs `(p,q)`
modulo simultaneous sign, and the action is

```
   (p,q) ↦ projective_primitive( M · (p,q) ).
```

No 3-dimensional tangle isotopy claim enters this
construction: it is an exact model of **boundary surface
mapping classes** only.

## Theorem A5.2b.1 — infinitely many undetected pure seam twists

For each signed integer `n ∈ Z`, set

```
               [ 1  2n ]
    T_n    =   [ 0   1 ].
```

Then all of the following statements hold:

1. `T_n ∈ Gamma(2)`; it fixes **each** of the four
   punctures, not just their unordered set.
2. Its projective action on a primitive vertical
   slope is `T_n(0,1)=(2n,1)`, interpreted as an
   unoriented primitive pair.
3. If `m ≠ n`, `T_m` and `T_n` represent
   **different** mapping classes of the marked
   four-punctured sphere.
4. The **same** proposed four-puncture V4 matching,
   for **every** `v4∈{0,1,2,3}`, is compatible
   with all these distinct pure twists.

**Proof.** The determinant of `T_n` is one,
and it reduces to the identity matrix modulo two.
Thus it fixes the four branch points in the
double cover and all four marked sphere punctures.
Matrix multiplication maps `(0,1)` to `(2n,1)`.
For distinct n, these are distinct primitive
unoriented slopes, so the induced actions on
essential curve classes differ, hence their
sphere mapping classes differ. Yet their action on
puncture labels is identical, since each is pure.
Composing a fixed proposed puncture permutation
with any of these pure actions therefore leaves
that four-point map unchanged while varying its
underlying surface mapping class. □

A useful exact witness is the determinant between
the vertical slope and its twisted image:

```
  |det((0,1),(2n,1))| = 2 |n|.
```

This integer is a transparent algebraic way to
distinguish each nonzero twist from the identity.
The actual minimal intersection number of essential
unoriented curves on the four-punctured sphere is
related to twice the determinant, but our tests
intentionally certify only the **integer homology
determinant**, not a geometric curve renderer.

**Corollary A5.2b.2.** Any interface that records
only four labelled punctures and a V4 element is
**insufficient** for uniquely specifying even a
single 4-punctured Conway seam **mapping class**.
Adding all integer sectors and ring numbers of the
incident weighted vertices does not itself supply
the *missing choice of gluing homeomorphism* without
an explicit source construction connecting them.

The source's plumbing conventions may indeed select
such a canonical homeomorphism — but the program must
implement and verify that rule from the actual
Fig.12.29/§12.7 boundary parametrizations rather than
infer it from a four-point permutation.

## Exact implementation and proof-replay

`arborescent/conway_mapping_class.py` provides:

- `PureConwayMappingClass(a,b,c,d)`: exact
  `SL₂(Z)` determinant-one matrices required to
  reduce to identity modulo two. A projective key
  identifies `M` and `−M`.
- `branch_point_action()`: calculates the
  **actual** four 2-torsion images modulo two
  instead of merely returning a claimed property.
- `compose`, `inverse` and exact action on
  primitive unoriented integer slopes.
- `pure_twist(n)`: the explicit infinite
  sequence `T_n`.
- `seam_ambiguity_witness(plan,parent,child,v4,n)`:
  requires a **real A5.1 typed parent-child seam**,
  holds the chosen four-puncture matching fixed,
  and returns two distinct pure surface mapping
  classes with a primitive slope witnessing their
  difference.
- `require_verified_conway_gluing()` deliberately
  raises `NotImplementedError`: an algebraic
  pure mapping class alone does not orient, frame
  or realize a 3D tangle seam.

The tests independently check matrices for all
powers `−24…24`, signed composition and
inversion, projective ±I equivalence, slope
separation, determinant preservation, *all*
four 2-torsion branch points, every V4 coordinate
permutation and actual signed A5.1 k-base seam
fixtures. Invalid nonprimitive slopes, non-pure
matrices, nonexistent seams and false geometry
approval all fail closed.

Run:

```sh
python -m unittest discover -s tests -v
```

## A5.2c acceptance contract — still-open geometry

This stage does not render Fig.12.29, and
the A5.2a two-braid RII certificate remains
restricted to two-braid geometry. To construct a
**source-sound** Conway seam we still need to:

1. Determine the exact boundary *markings* and
   pure mapping-class twists selected by the
   source's Fig.12.29 atomic model, including
   the sector integers and rings.
2. Specify the actual **orientation-reversing**
   gluing homeomorphism between the two oriented
   four-punctured Conway spheres, including
   framing and the induced action on all strand
   endpoints. A `Gamma(2)` element by itself is
   orientation-preserving and **not** a seam gluing.
3. Validate the resulting assembled embedded
   tangle against independent planar diagrams
   and relative-boundary isotopy certificates.
4. Only then attempt both signed actual R±
   source diagrams with canonicality hypotheses
   and exceptional adjacent-ring rules.
5. Preserve A5.1's explicit
   `UNVERIFIED_GEOMETRY` status until an
   independently checked *specific geometric
   proof object* exists.

This is a **new rigorous obstruction** identifying
information missing from the interface, not a
negative result about the original conjecture's
isotopy-class growth. Issues #33 and #16 remain
open.
