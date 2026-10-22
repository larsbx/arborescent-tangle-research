# Prototype M0: restricted Montesinos link enumeration

## What is being counted

Fix the alphabet of distinct proper rational slopes

```text
L = {1/5, 2/5, 1/7}
```

and fix the integer twist parameter `e = 0`. For `k >= 3` take standard cyclic Montesinos link presentations `M(0; q_1,...,q_k)` with `q_i in L`. We study **unoriented links in oriented S³ under orientation-preserving ambient isotopy**, not boundary-fixed rational tangles. The numerical program computes parameter-word orbits; only the published Montesinos classification permits the topological interpretation.

An additional raw binary-sum *syntax* object includes a full planar binary tree parenthesizing the `k` ordered labels. `n = k-1` is the internal-node count of this **syntax tree**, not the crossing number or reduced arborescent plumbing-tree size.

Four stages must remain separate:

1. **Raw binary syntax**: `C_(k-1) m^k` with alphabet size `m`.
2. **Associativity of horizontal tangle sum**: erase parentheses, leaving `m^k` ordered words. Different words need not be different boundary-fixed tangles.
3. **Cyclic link presentation**: quotient by rotations and reversal of the rational-piece order, obtaining `D_k` orbits (bracelets).
4. **Link isotopy**: in the standard non-exceptional Montesinos classification the slope fractional parts up to dihedral order and total slope/twist classify these links. Our alphabet has each slope in `(0,1)`, integral term fixed at zero, and each denominator at least 5. Thus `sum(1/denominator(q_i)) <= k/5 <= k-2` for `k >= 3`, meeting the conventional non-exceptional restriction. The tuple orbit then determines the total slope as well.

The precise theorem normalization, signs of the global twist, exceptional cases, mirrors, and boundary markings must be rechecked against [Bonahon–Siebenmann and subsequent expositions](references.md) before extending the alphabet. The dihedral quotient is **not** an algorithm for testing tangle isotopy. Rows `k=1,2` are combinatorial controls only.

## Exact combinatorics

Let `C_n = binom(2n,n)/(n+1)`. Define

```text
S_k(m) = C_(k-1) * m^k                raw binary syntax
W_k(m) = m^k                          ordered words
B_k(m) = (R_k(m) + F_k(m)) / (2k)     dihedral word orbits
R_k(m) = sum(j=0..k-1) m^gcd(j,k)
F_k(m) = k*m^((k+1)/2)                             if k odd
       = (k/2)*(m^(k/2+1)+m^(k/2))                 if k even
```

Burnside's lemma proves the formula: rotation by `j` has `gcd(j,k)` letter cycles. Every odd-length reflection has `(k+1)/2` cycles; for even `k` half have `k/2+1` and half `k/2` cycles.

For fixed `m >= 2`, `B_k(m) ~ m^k/(2k)`. The calculation uses exact integers; a separate brute-force enumeration of canonical dihedral word keys cross-checks small cases.

## A concrete negative result

For a fixed ordered length-`k` word, all `C_(k-1)` binary parenthesizations produce the same horizontal tangle sum by associativity. Consequently a **polynomial bound on fibers from unreduced binary syntax to isotopy is false**, even before closure. This is an explicit obstruction, not a numerical guess.

After parenthesization reduction, `W_k(m)/B_k(m) ~ 2k`: the remaining word-to-bracelet ambiguity is polynomial, but **the Catalan growth factor has vanished**. This demonstrates why a polynomial-fiber result *alone* cannot prove Catalan-led isotopy growth.

## Next target

The bona fide arborescent object is the weighted planar **plumbing tree** of Bonahon–Siebenmann, with vertex/sector weights, local cyclic orders, and published moves (Chapter 12, especially ``12.2–12.7). Build a separate grammar and certified reduction for **these** trees before inferring any Catalan enumeration for arborescent isotopy classes.
