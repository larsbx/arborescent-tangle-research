# Prototype M0: exact counts

Reproduce with:

```sh
python -m arborescent_tangles.counts --max-k 10 --alphabet-size 3
python -m unittest discover -s tests -v
```

Alphabet: `L={1/5,2/5,1/7}` with three distinct proper rational slopes; integral twist `e=0`. The dihedral-orbit column represents Montesinos **link** isotopy classes only in the restricted `k>=3` setting under the classification hypotheses explained in [M0](montesinos-prototype.md). `k=1,2` are **not** topology claims.

| k | Raw binary syntax | Ordered words | Dihedral orbits |
|--:|--:|--:|--:|
| 1 | 3 | 3 | 3 |
| 2 | 9 | 9 | 6 |
| 3 | 54 | 27 | 10 |
| 4 | 405 | 81 | 21 |
| 5 | 3402 | 243 | 39 |
| 6 | 30618 | 729 | 92 |
| 7 | 288684 | 2187 | 198 |
| 8 | 2814669 | 6561 | 498 |
| 9 | 28146690 | 19683 | 1219 |
| 10 | 287096238 | 59049 | 3210 |

**Exact:** raw / ordered = `C_(k-1)`; binary parenthesizations disappear under associativity. Burnside yields the bracelet count, and the small cases have been cross-checked against direct D_k orbit enumeration.

**Asymptotic:** `B_k(3) ~ 3^k/(2k)`. This explicit family does **not** retain the naive Catalan factor.

**Not tested:** an independent diagrammatic ambient-isotopy implementation, canonical weighted plumbing-tree moves, crossing numbers or a universal `n^(-3/2)` law. The observed integer sequence is not by itself a proof of any topological classification; that implication relies on the published Montesinos theorem.
