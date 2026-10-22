# Arborescent tangle research

A reproducible research program on **Catalan syntax, arborescent plumbing trees, and the effect of isotopy on enumeration**. The initial repository establishes an exact, restricted Montesinos **link** experiment rather than claiming a Catalan classification of arbitrary non-rational tangles.

## First result: binary Catalan shapes disappear in a Montesinos test

For `k` ordered rational-tangle labels from a fixed alphabet of size `m`:

- Unreduced binary parenthesizations: `C_(k-1) * m^k`.
- After associativity of fixed planar sum: `m^k` ordered words.
- After cyclic rotation and reversal of a standard Montesinos **link** presentation: `B_k(m)`, the number of dihedral word orbits.

For fixed `m >= 2`, `B_k(m) ~ m^k/(2k)`. Thus the Catalan factor is eliminated by associativity in this test. A purported universal Catalan isotopy-growth law already fails here.

The dihedral number is an isotopy count **only under the Montesinos classification theorem's hypotheses and equivalence conventions**. The restricted family fixes `L = {1/5, 2/5, 1/7}`, integral twist `e = 0`, and `k >= 3`; see [the model and qualifications](docs/montesinos-prototype.md).

## Reproduce the experiment

Python >= 3.11; no third-party dependencies:

```sh
python -m unittest discover -s tests -v
python -m arborescent_tangles.counts --max-k 10 --alphabet-size 3
```

The Burnside dihedral orbit count is independently checked against brute-force orbit enumeration for small parameters.

## Project documents

- [Research conjectures, mathematical scope, falsification, roadmap](docs/research-program.md)
- [Montesinos prototype and theorem assumptions](docs/montesinos-prototype.md)
- [Exact results](docs/prototype-results.md)
- [Primary literature](docs/references.md)

**Limits:** this code does not construct knot diagrams or test ambient isotopy, enumerate Bonahon–Siebenmann reduced plumbing trees, compute crossing numbers, or classify boundary-fixed tangles. Those are explicitly separate phases.
