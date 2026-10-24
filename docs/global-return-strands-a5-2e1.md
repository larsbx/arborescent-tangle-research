# A5.2e1 — Global returning strands with n genuine Conway boundaries

**Topological scope:** exact, independently audited PL *curve complement* candidate for the **r=0** Bonahon–Siebenmann atomic model. This realizes the two long returning **strands**, with no artificial anchor boundary and no unpaired strand endpoint. The underlying **twisted band / ribbon surface, its framing, and source-equivalence have NOT been constructed**. The code remains `RETURN_STRANDS_SOURCE_UNVERIFIED`.

## Primary source: figures inspected as page images

Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots*, §12.2 **Figure 12.8**, printed p.142, PDF index 153, and §12.5 **Figure 12.29**, printed p.158, PDF index 169:

https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

Figure 12.8 shows *n* actual Conway boundaries aligned in the body of an atomic tangle, *n* signed half-twist regions `a₁,…,aₙ`, and two very long exterior return paths. Figure 12.29 adds *r* closed rings at the left end and is identical for `r=0` at the depicted level. A5.2d modeled a **different** arrangement of *n+1* 3D cubes with an artificial, unmarked left cube and four exposed extreme endpoints. That was an intermediate linear construction, **not** Figure 12.8's strand geometry. A5.2e1 eliminates the artificial sphere and supplies return curves.

This construction is a **chosen spatial convention**, not yet a proof that it is isotopic relative to the source's boundary parametrizations. In particular, source images alone do not certify the placement and framing of a missing complete ribbon surface. The named `Atom.sign` ring polarity is carried without being used to infer unproved signed crossing conventions.

## Exact coordinate model

For source cyclic integer sector vector `(a₀,…,aₙ₋₁)`, `n≥2`, `ring_count=0`:

- Initial unmarked left entrance plane `x=0`, with two *non-boundary*, intermediate joining vertices `(0,−2,0)` and `(0,+2,0)`. These are *not* Conway sphere punctures.
- Actual n Conway cubes, one for each genuine A5.1 `Port(atom.address,i)`. Cube 0 has centre `c₀=3+2|a₀|`. For i≥1, `cᵢ=cᵢ₋₁+5+2|aᵢ|`. Each cube has x half-width 2 and y/z half-width 3, with marked punctures on the left/right x faces at y=±2, z=0.
- Signed half-twist region `a₀` connects the entrance plane to the *west face* of cube 0. Its longitudinal gap is `2|a₀|+1`. For i≥1 signed `aᵢ` connects the *east face* of cube i−1 to the *west face* of cube i, again with a `2|aᵢ|+1` positive-length corridor. They use the exact A5.2c signed rational spatial two-strand geometry.
- The two *long return rails* independently join the two points of the entrance plane to the **east-face punctures of cube n−1** by routing out to y=−6 (upper) and y=+6 (lower), respectively. They are piecewise rectilinear and stay outside *all* solid cubes because the cubes have |y|≤3. At the left end they detour through x=−3, and at the right end through `x=last.east+3`.
- Splice each return rail to its corresponding `a₀` entrance strand. There are **two complete long outer arcs**, from the last cube's east face to the first cube's west face, and **two short arcs per interior gap**. The global embedded 1-manifold therefore consists of exactly `2n` disjoint polygonal arcs, with **each of the 4n marked punctures appearing exactly once as a terminal point**. No source-boundary endpoint is left unmatched.

The first outer pair has parity-controlled endpoint connections: even `|a₀|` sends final east NE to first west NW and final east SE to first west SW; odd `|a₀|` swaps the west destinations. Each interior pair independently has the corresponding parity-controlled matching. These statements follow from the *constructed PL paths*, not from the A5.1 rooted code.

## Independent collision validator

`independent_pl_arc_audit(cubes, arcs)` accepts only physical PL triangle-boundaries and polygonal arcs. It is **not** a call to a constructor's layout method. It rederives the entire boundary matching by reading the **actual global terminal coordinates**, checks that the endpoints are exactly the specified 4n cube punctures, and validates:

1. Every segment has positive length, every curve is typed rational 3D geometry, and every source Conway puncture is a terminal point of exactly one arc.
2. Every segment against every *closed solid cube*, not only mesh edges, using **exact interval clipping**. A whole interval through/on the cube fails; an isolated contact is allowed only if it coincides with a *terminal marked puncture*.
3. Every unordered pair of segments, including segments of the *same* arc, with exact rational segment intersection. Skew pairs pass; point contacts pass **only** for consecutive segments sharing their one intended vertex; collinear positive-length overlap, unanticipated tangency and all other intersections fail.
4. Source sector indexing/weights, spatial corridor words, all genuine port addresses, cube placement, and canonical return rail layout are also audited in the typed `ZeroRingReturnStrandCandidate`.

The independent checker uses no floating tolerances, and distinguishes source-boundary contact from false intersections in the 2D planar projection (which are disjoint at over/under z-heights). Its positive result is **a concrete valid embedded curve collection with 4n boundary marks**, not a classification theorem.

## Validation and falsification

`tests/test_atomic_return_strands.py` includes exact segment intersection negative cases (skew, parallel, collinear overlap, tangencies, rational crossings), exact closed-solid-cube clipping, both atom signs, even/odd crossing parities, signed/nonzero sector weights, complete 4n endpoint coverage, deterministic coordinates, and an exhaustive n=2,3 sweep over weights in `{-2,-1,0,1,2}`. Adversarial mutation tries to force a strand through a cube, reuse a boundary endpoint, reindex corridors, exchange the two return rails, forge a status, or request unsupported `r>0` and valence-one construction. All such requests must fail closed.

```sh
python -m unittest discover -s tests -v
```

Require hosted Python 3.11 and 3.12 on the exact pull-request head and post-merge main.

## Pending A5.2e2, A5.2f and #33 obligations

**A5.2e2:** Construct the **complete oriented, embedded supporting twisted band/ribbon** that carries this candidate curve model. Verify it with an independent PL surface-mesh audit (edge incidences, strip orientation and two boundary tracks). Derive the precise mapping from the *chosen spatial* port frames and sector signs to the *source Figure 12.8's* convention. Test cyclic rotations, mirrors, and n=1 separately. Only then can any specific r=0 model be called source verified.

**A5.2f:** Model r=1,2,3 closed ring loops from Figure 12.29 and independently check exact embedded separation and **linking**. Do not reinterpret closed rings as ordinary extra braid generators.

**Beyond:** Source-correct orientation-reversing 4-punctured Conway gluing including the nontrivial infinite pure mapping-class ambiguity of A5.2b; both source-modified R± relative-boundary isotopy certificates with adjacent-ring exceptions; only thereafter Theorem 12.21 and isotopy-class enumerations.

**No public API may turn this object or any A5.1 plumbing plan into `VERIFIED_TOPOLOGY`.** This patch establishes precise positive 3D embedding statements only for the constructed curve candidate.
