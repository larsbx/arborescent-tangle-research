# Independent primary-diagram audit of A2.6 (2026-10-23)

Source: Bonahon–Siebenmann, *New Geometric Splittings of Classical Knots and the Classification and Symmetries of Arborescent Knots*, 2016 edition:
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

This audit is of the **printed-page diagrams**, not an inference from
the previous PR prose.

| Location | Direct observation | Consequence for code |
| --- | --- | --- |
| §12.4, printed p.155, PDF index 166 | Both R+ and R− replacements show the signed three-vertex ±2—±1—±2 subtree on **both sides** of the figure, moved from one side of its attachment to the other. | Vertex **deletion cannot be represented as an R± rewrite**. The retracted PR #11 contraction must not be reintroduced as such. |
| §12.5, printed p.157, PDF index 168 | The authors explicitly erase ring subtrees **and the corresponding attachment bond** at vertices of valence at least 3, recording the number erased as the **ring number** at the surviving vertex. Three placements of the +1 central weight are depicted as eligible *ring subtrees* up to the earlier flip; the sample includes multiplicities exceeding 1. | The current matcher covers only the strictly oriented between-tips case and requires zero sectors on both sides of the stem at the base. This is intentionally **incomplete** with respect to all source ring subtrees. The abbreviation is a distinct, information-losing representation. |
| §12.5, printed p.157 bottom and p.158 top, PDF indices 168–169 | When two ring centers are adjacent, erasing one destroys the other; the source gives a special tiny-tree representative instead of assigning independent ring numbers. | Current code correctly **rejects** the exception. Implementing the special representative remains open. |
| §12.5, printed p.158, PDF index 169 | The ring-numbered tree is given a plumbing interpretation, and the text asserts agreement of the underlying arborescent pair for its applicable canonical trees. | This is **not** a verification of the entire plumbing realization or its boundary-parameterized variant by this implementation. |
| §12.5, printed p.159, PDF index 170 | Lemma 12.14 and Theorem 12.15 describe a classification of suitably canonical trees via abbreviated trees and flips, not equality of arbitrary labelled presentation codes. | Do not promote `Abbreviated.rooted_code()` equality into an isotopy theorem. |

## Adversarial distinction

The added `tests/test_ring_abbreviation_adversarial.py` checks two
genuinely different **marked plane-tree presentations** with a changed
cyclic placement of rings relative to a nonring weight-3 branch. The
abbreviation codes agree, demonstrating information loss. It is
deliberately *not* asserted that code equality establishes a valid
R± path or labelled-tangle isotopy in an unrestricted category.

The tests also check a changed **non-root** V4 free-bond label cannot
disappear under abbreviation, both signs of the adjacent-ring
exception raise, and the non-between-tips signed angular variants are
left intact rather than silently normalized. The p.155 diagrams do
not encode the full V4 boundary bookkeeping: retaining labels in the
software is a syntactic preservation check, **not** a new source
proof of boundary-fixed isotopy.

## Acceptance and open scope

PR #12's exact proposed head
`f975ce81993877902f06ee0b382912f2487bf945` passed the
GitHub-hosted `python -m unittest discover -s tests -v` matrix:
67 tests on each of Python 3.11 and 3.12, workflow run 40034626922.
PR #12 merged as `343fa2d3c6b8cab9ae50b1c3d6c5097992882c37`
and PR #10 was closed as superseded. These checks establish software
acceptance at that head; they do **not** prove topological soundness.

Open under A2 issue #4: genuine source-exact R± relocation and inverses,
full F3 variants, exceptional canonicality, independent plumbing
soundness, marked boundary equivalence, and classification of full
isotopy classes or Catalan growth.
