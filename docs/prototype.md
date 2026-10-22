# Explicit negative control for Catalan survival

These are **presentation-level** counts, not a computation of topological isotopy classes. The alphabet consists of distinct rational slopes in (0,1). Parenthesized syntax -> ordered words -> cyclic dihedral presentation orbits.

## Exact output: two labels (1/2, 1/3)

| k | C(k-1) | raw binary trees | ordered words | dihedral presentations |
|--:|--:|--:|--:|--:|
| 3 | 2 | 16 | 8 | 4 |
| 4 | 5 | 80 | 16 | 6 |
| 5 | 14 | 448 | 32 | 8 |
| 6 | 42 | 2688 | 64 | 13 |
| 7 | 132 | 16896 | 128 | 18 |
| 8 | 429 | 109824 | 256 | 30 |

## Exact output: three labels (1/2, 1/3, 2/5)

| k | C(k-1) | raw binary trees | ordered words | dihedral presentations |
|--:|--:|--:|--:|--:|
| 3 | 2 | 54 | 27 | 10 |
| 4 | 5 | 405 | 81 | 21 |
| 5 | 14 | 3402 | 243 | 39 |
| 6 | 42 | 30618 | 729 | 92 |
| 7 | 132 | 288684 | 2187 | 198 |
| 8 | 429 | 2814669 | 6561 | 498 |

The tests check Catalan values, rational slope normalization, fixed-word reparenthesization, rotation/reversal invariance and independent exhaustive orbit enumeration against Burnside for m=1,2,3 and k=1..8.

Even though the map **words -> dihedral orbits** has fibers of size at most 2k, **binary syntax -> ordered words** has fibers C(k-1), so the Catalan factor is lost. The quotient sequence B_m(k) is asymptotic to m^k/(2k), with a logarithmic generating-function singularity.

## Topology limitation

The Montesinos-link classification reported in Futer–Kalfagianni–Purcell §8.1 refers to cyclic link **closures** and depends on explicit hypotheses and conventions; its converse is not programmed here. In particular, do not read `dihedral-orbits` as `boundary-fixed tangle isotopy`, or even as `link isotopy` without verifying those hypotheses. No genuinely branching reduced arborescent weighted tree is enumerated yet.

Run: `python -m unittest discover -s tests -v` and `python -m arborescent.prototype --slopes 1/2 1/3 2/5 --max-k 8`.
