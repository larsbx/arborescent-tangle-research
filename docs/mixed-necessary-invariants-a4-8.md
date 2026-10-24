# A4.8 — Necessary invariants and four-valued reachability evidence

**Scope:** the explicitly implemented *partial formal presentation*
move graph A4.7, with source-defined F1/F2 and surviving-edge strict
F3′, plus A4.1's **formal** ring-stem slides. There is **no claim of
ambient tangle isotopy**, a complete Bonahon–Siebenmann calculus, or
negative isotopy decisions.

## The problem

A4.7 introduced replayable mixed-move paths and honest cutoffs, but
the mathematical interpretation of a failed search needs care:

1. A positive, replayed path proves **formal reachability** under
   the declared moves.
2. An incomplete search proves no negative result whatsoever.
3. A complete search proves nonreachability only of an **exact
   indexed presentation target** in its explored move graph.
4. Different vertex-storage indices and arbitrary internal-bond
   names are not moves in that graph, although the same marked
   weighted plane-tree presentation may have several such encodings.

The fourth point is an independent obstruction to reading
`complete=True` as inequivalence of unindexed presentations,
let alone inequivalence of actual tangles.

We therefore add a genuinely mathematical **necessary invariant**
that can exclude pairs *without any move search* and distinguish
those exclusions from inconclusive enumerations.

## Definition: rooted unordered weighted-ring skeleton

For a full V4-labelled input T, require complete, uniform-sign
strict ring-cap recognition under A2.6/A4.7, and let `A_s(T)`
be the restricted §12.5-style **abbreviated representation**.

Fix an existing surviving marked free root bond `f`. Root the
residual tree at the vertex incident to f. For each vertex v,
define recursively

```text
C(v,parent) = (
    W(v),
    R(v),
    sort( { ("free",) for each other free incident bond,
            ("branch", C(child,v)) for each child bond } )
)
I_f(T) = (s, C(root,f)).
```

Here `W(v)` is the integer *total* of its angular-sector weights
and `R(v)` is the nonnegative integer ring number stored at the
surviving vertex. The multiset-sort deliberately discards:

- cyclic child order and arbitrary starting slots;
- individual sector *locations* (but **not** their integer totals);
- V4 boundary labels (which F1/F2/F3′ can update);
- arbitrary vertex and internal-bond identifiers;
- unmarked free-leg names (which do not themselves define intrinsic
  topological identifiers).

The root remains individually marked. This choice makes the code
weak enough to be preserved by each of the declared generators,
but strong enough to separate many weighted/incidence structures
that have the same number of vertices or the same weight multiset.

## Proposition A4.8.1 — preservation under every accepted move

Let T be in the strict A4.7 source-pattern domain and let
`T→T′` be one accepted F1, F2, strict surviving-edge F3′, or
formal ring-position slide, with the same root f and ring sign s.
Then

```text
I_f(T) = I_f(T′).
```

**Proof.** F1 and F2 reflect the cyclic lists at designated
vertices; sectors are permuted, so their *sums* are fixed.
They do not change graph incidence or free-bond ownership.
Their ξ/η/ζ transformations may alter individual V4 labels,
but the invariant explicitly omits those labels.

A strict F3′ moves one integer angular-sector weight between
two sectors **of the same pivot vertex**; any odd-weight reflection
again permutes sectors in its selected branch. Thus every surviving
vertex's total integer weight and the rooted incidence graph
remain fixed. The allowed F3′ is required to preserve strict
ring recognition and each surviving base's ring number by its
A4.5 commuting-square guard.

An A4.1 formal ring slide swaps the cyclic positions of two
attachments at one base. It does not change the underlying
incidence tree, any total vertex weight, any V4 label, or any
strict §12.5 ring-number attachment. A4.3 rooted cyclic
re-alignment is a presentation isomorphism, and also preserves
all these data. Recursion from the same marked free root,
sorting the unordered child subcodes at every vertex, therefore
returns the same tuple in each case. □

**Corollary A4.8.2 — rigorous restricted negative certificate.**
If `I_f(T)≠I_f(U)`, no sequence of the *implemented partial
formal moves* can connect T and U, regardless of search depth.

**Proof.** Every accepted step preserves the invariant.
An imagined connecting path would preserve it by induction,
contradicting unequal endpoint codes. □

This corollary does **not** prove that the associated 3D tangles
are nonisotopic. Missing source moves could relate them.
Conversely, `I_f(T)=I_f(U)` does not establish reachability under
even the partial formal rules.

## Four distinct outcomes: API contract

`decide_partial_formal_reachability(start,target,sign,root,max_states)`
returns one of four statuses:

| Status | Precise meaning |
|---|---|
| `reachable_formally` | A replayable typed A4.7 move path reaches the **exact root-aligned indexed target**. |
| `excluded_by_necessary_invariant` | The exact rooted unordered weighted-ring codes differ; no A4.7 *formal* path exists, without enumerating any states. |
| `not_reached_in_complete_indexed_graph` | Full finite BFS exhausted its **indexed graph** without that exact stored target; renumbered or topologically equivalent presentations may still exist. |
| `undetermined_cutoff` | State budget stopped search before exhaustion; no negative statement follows. |

The response also reports the pair of fingerprints, states
actually examined, completeness and (only for reached states)
a replayable `MixedTrace`. The boolean
`FormalReachabilityDecision.verify_positive()` verifies only
the positive path; it does not turn any negative status into an
isotopy oracle.

## Independent negative controls and exact examples

Regressions in `tests/test_mixed_invariants.py` include:

1. Every one-step generator emitted on signed corollas, k=2
   bounded-valence ring paths, and nontrivial signed odd/even
   strict F3′ witnesses preserves the independently computed
   rooted unordered fingerprint.
2. Two trees with the **same total-weight multiset** but
   different rooted branching shapes have different fingerprints.
   This demonstrates a stronger obstruction than a raw
   vertex-count or total-weight check.
3. Changing one weighted leaf yields a provable necessary-invariant
   exclusion without any BFS.
4. Adjacent ring locations produce a positive shortest typed
   formal path with exact replay.
5. A one-state cutoff reports `undetermined_cutoff` rather than
   inequivalence.
6. A graph-isomorphic presentation produced by **reversing the
   stored vertex array and renaming every internal edge** has the
   same fingerprint and `labelled_rooted_code` as its source,
   but is not reached by a complete *fixed-index* move search.
   This independently verifies why completeness of indexed BFS
   is weaker than completeness modulo presentation isomorphism.
7. Changed individual V4 boundary labels, F1 reflection and
   cyclic storage renaming do not change the necessary fingerprint.
   This is deliberate: the fingerprint is a **necessary** rather
   than a complete invariant.
8. Unsupported roots/signs and invalid budgets fail closed.

Reproduce:

```sh
python -m unittest discover -s tests -v
```

## What remains open

[Issue #16](https://github.com/larsbx/arborescent-tangle-research/issues/16)
still requires an independent Conway/plumbing model, source-complete
positively/negatively canonical hypotheses and exceptions, actual
modified R± moves, full F3′, V4 boundary-parametrized isotopy
and a proof of a complete **topological** equivalence quotient.

This stage intentionally strengthens the **epistemic status
discipline** of the code. It neither claims that a source-defined
partial flip is a verified 3D isotopy nor that failure of any
formal search implies different isotopy classes.
