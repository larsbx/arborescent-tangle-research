# A2.4: bounded-local weighted internal forks

## Source audit and scope

Primary source: Bonahon–Siebenmann, *New Geometric Splittings of Classical
Knots*, January 26, 2016,
https://dornsife.usc.edu/francis-bonahon/wp-content/uploads/sites/205/2023/06/BonSieb-compressed.pdf

The actual PDF text and diagrams were inspected: printed p.153 (PDF page
165), p.166 (PDF page 178), and §12.7.3, printed p.167. Page indices here
are one-based. The four P diagrams on p.153 are **singleton sticks**,
not a sequence of adjacent negative weights: isolated -1, isolated -2,
and one- or two-attached singleton -2. The analogous N diagrams reverse
sign. Section 12.7.2 exempts sticks carrying a free bond. No pair F
condition belongs in the tangle screen.

`internal_forks.tangle_screen` implements W, the tangle S condition,
and P/N on trees with at least one valence-at-least-three vertex.
Sticks are components after deleting such vertices, counting free bonds
in valence. Zero weights are allowed only at free end vertices; nonzero
weights alternate after zeros are removed; end ±1 requires a free bond.
P/N rejects singleton -2/+2 without free bonds. Isolated ±1 whole trees
and exceptional zero whole trees are outside this implementation:
**all unbranched whole trees return unresolved**, not canonical.
Ring moves, canonical reduction, and independent plumbing/isotopy
soundness remain open under issue #4. Passing this screen verifies the
implemented structural source conditions only.

## Witness and minimal production set

The PR #7 four-vertex witness has a weight-3 valence-three root fork,
a zero-weight valence-three second fork, two weighted tips 2 and 3,
and two free ports. It passes positive W/S/P. It fails N at its weight-2
singleton tip stick. It fails the inappropriate pair F screen at the
weighted root free port. This is expected, not a tangle defect.

Fork vertices are removed before S/P/N are evaluated, so their integer
weights are unrestricted by those conditions. W allows at most one
nonzero sector at each. For the binary grammar below, every remaining
stick is a singleton tip: tips |a|>=2 satisfy S, positive canonicity
excludes -2, negative canonicity excludes +2. Internal weights can be
positive, negative or zero. The default experiment uses {0,3}.

Use productions `F -> D(T,T)`, `T -> Tip(a) | Port | F`, with two
ordered children and incoming bond at each F. D selects zero (one
encoding), or an integer weight and one of three angular sectors.
The root incoming bond is f0. Port terminals create additional free
bonds; Tip creates a valence-one weighted vertex. All free bonds get
explicit identity V4 labels. Binary forks alone suffice for the witness;
ternary forks are not needed for this first extension. This is a minimal
production vocabulary for this witness, not a uniqueness/minimality theorem
among all conceivable grammars.

Generation fixes n forks and p additional ports. With m tip labels and
D decorations, its exact count is C_n binomial(n+1,p) m^(n+1-p) D^n.
Compilation is preorder with deterministic edge/port names. The witness
is generated at n=2,p=1; its graph agrees with PR #7 up to vertex/edge
renaming. The regression tests exercise the same branch transport.

## Exact bounded rewrite experiment

Counts below are **intersections of seed sets with indexed fixed-graph
orbits**, not counts modulo graph isomorphism and not isotopy counts.
All seeds have identity V4 labels. Searches can leave the grammar slice.
Every full F1/F2/strict-F3′ search must complete; cutoff raises instead of
returning a quotient count. Both F2 pivot parity classes are used in the
baseline; these indexed actions must not be confused with the four-action
root gauge slice in A2.2.

| seed family | forks | seeds | F1/F2 classes | partial-move classes |
|---|---:|---:|---:|---:|
| A1 binary restriction, tips {2,3}, internal 0, one root port | 1 | 4 | 4 | 4 |
| same | 2 | 16 | 16 | 16 |
| new binary grammar, tips {3}, internal {0,3}, two ports | 1 | 8 | 8 | 6 |
| same | 2 | 96 | 96 | 50 |

The A1 control is an explicitly matched **binary subgrammar**, not the
full A1 binary/ternary enumeration. The new rows use a different tip/port
slice and cannot support a growth-rate comparison. All new seeds pass
positive W/S/P. The two-seed witness/transport-target control has 2 F1/F2
classes and 1 partial class; a replayable path proves an identification
absent from F1/F2 even on the exact indexed representation. No claim is
made here about novelty modulo arbitrary graph isomorphism or topological
isotopy. Every recorded new seed identification replays through the existing
move implementation.

## Validation and acceptance

Reproduce with `python -m unittest discover -s tests -v`. New tests cover
independent combinatorial counts, duplicates, witness structure, all reachable
witness screens, replay of new identifications, cutoff rejection, A1 negative
control, invalid grammar inputs, P/N sign differences, tangle free endpoints,
S zero/sign/end failures, and unresolved whole-tree exceptions.

PR #7 is now merged at 197a69a7b7a5d9959b03bd81e94e0ecd66e7da9e.
Its hosted acceptance debt is resolved by Tests run 39993723571:
both `unit (3.11)` and `unit (3.12)` succeeded. This independently verifies
the hosted fixture repair; it does not transfer acceptance to this PR.
https://github.com/larsbx/arborescent-tangle-research/actions/runs/39993723571

This PR retains the hosted 3.11/3.12 gates unchanged. Full exceptional
canonicity, modified rings, marked-boundary plumbing validation, and an
independent topological justification remain separate obligations. In
particular, Theorem 12.21 is not being promoted to a new implemented
isotopy classifier or a Catalan isotopy counting theorem.
