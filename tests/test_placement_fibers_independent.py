"""Independent incidence checks of the selected A3.2 binary family."""
import unittest
from itertools import product

from arborescent.placement_fibers import bounded_valence_chain_presentations
from arborescent.ring_abbreviation import Abbreviated, abbreviate_disjoint_rings
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def decode_choices(source, sign):
    """Recover choices from rooted incidence, without codes or cap recognizers.

    Zero-sector vertices are bases; the unique nonzero neighboring subtree
    is the cap. Follow the base path from the owner of the marked free leg.
    Verify the cap's complete local shape rather than trusting vertex indices.
    """
    vertices = source.tree.vertices
    owners = {}
    for i, vertex in enumerate(vertices):
        for bond in vertex.slots:
            owners.setdefault(bond, []).append(i)
    bases = {i for i, v in enumerate(vertices) if not any(v.sectors)}
    current, = owners['f0']
    incoming, visited, bits = 'f0', set(), []
    while True:
        assert current in bases and current not in visited
        visited.add(current)
        vertex = vertices[current]
        start = vertex.slots.index(incoming)
        slots = vertex.slots[start:] + vertex.slots[:start]
        kinds = []
        outgoing = None
        for bond in slots[1:]:
            if bond.startswith('f'):
                kinds.append('F')
                continue
            other, = (i for i in owners[bond] if i != current)
            if other in bases:
                kinds.append('P')
                outgoing = (other, bond)
                continue
            kinds.append('R')
            center = vertices[other]
            j = center.slots.index(bond)
            assert center.sectors[j:] + center.sectors[:j] == (0, sign, 0)
            assert len(center.slots) == 3
            for tip in center.slots:
                if tip != bond:
                    leaf, = (i for i in owners[tip] if i != other)
                    assert vertices[leaf] == Vertex((tip,), (2 * sign,))
        patterns = (('P', 'R'), ('R', 'P')) if incoming == 'f0' else (
            (('P', 'F', 'R'), ('P', 'R', 'F')) if outgoing else
            (('F', 'R'), ('R', 'F'))
        )
        bits.append(patterns.index(tuple(kinds)))
        if outgoing is None:
            break
        current, incoming = outgoing
    assert visited == bases
    return tuple(bits)


def expected_path(k, sign):
    """Construct the output directly, without deleting caps from a fixture."""
    vertices = []
    for i in range(k):
        slots = ((f'e{i-1}',) if i else ()) + (f'f{i}',)
        if i < k - 1:
            slots += (f'e{i}',)
        # Interior source order is incoming, outgoing, free, after erasure.
        if 0 < i < k - 1:
            slots = (f'e{i-1}', f'e{i}', f'f{i}')
        vertices.append(Vertex(slots, (0,) * len(slots)))
    return Abbreviated(V4LabeledTree(WeightedPlaneTree(tuple(vertices)),
                                  tuple((f'f{i}', 0) for i in range(k))),
                       (1,) * k, sign)


def reindex(source):
    """Change storage and names, preserving root and cyclic orientation."""
    names = {b for v in source.tree.vertices for b in v.slots}
    rename = {b: ('f0' if b == 'f0' else b[0] + str(i + 100))
              for i, b in enumerate(sorted(names))}
    vertices = []
    for v in reversed(source.tree.vertices):
        slots = tuple(rename[b] for b in v.slots)
        vertices.append(Vertex(slots[1:] + slots[:1],
                               v.sectors[1:] + v.sectors[:1]))
    labels = tuple((rename[b], value) for b, value in source.labels)
    return V4LabeledTree(WeightedPlaneTree(tuple(vertices)), labels)


class IndependentBoundedFiberTests(unittest.TestCase):
    def test_decode_all_choices_and_exact_outputs(self):
        for k in range(2, 8):
            for sign in (-1, 1):
                with self.subTest(k=k, sign=sign):
                    expected = expected_path(k, sign)
                    sources = bounded_valence_chain_presentations(k, sign)
                    decoded, codes = set(), set()
                    for source in sources:
                        decoded.add(decode_choices(source, sign))
                        codes.add(source.labelled_rooted_code())
                        self.assertEqual(len(source.tree.vertices), 4 * k)
                        self.assertLessEqual(max(len(v.slots) for v in
                                                 source.tree.vertices), 4)
                        self.assertLessEqual({w for v in source.tree.vertices
                                              for w in v.sectors},
                                             {0, sign, 2 * sign})
                        ab = abbreviate_disjoint_rings(source, sign)
                        self.assertEqual(ab.rings, expected.rings)
                        self.assertEqual(ab.labeled.labels, expected.labeled.labels)
                        self.assertEqual(ab.rooted_code(), expected.rooted_code())
                        self.assertTrue(all(not any(v.sectors)
                                            for v in ab.labeled.tree.vertices))
                    self.assertEqual(decoded, set(product((0, 1), repeat=k)))
                    self.assertEqual(len(codes), 2**k)

    def test_codes_ignore_storage_names_and_cyclic_start(self):
        for sign in (-1, 1):
            for source in bounded_valence_chain_presentations(3, sign):
                changed = reindex(source)
                self.assertEqual(decode_choices(changed, sign),
                                 decode_choices(source, sign))
                self.assertEqual(source.labelled_rooted_code(),
                                 changed.labelled_rooted_code())
                self.assertEqual(abbreviate_disjoint_rings(source, sign).rooted_code(),
                                 abbreviate_disjoint_rings(changed, sign).rooted_code())

    def test_nonroot_boundary_and_per_base_counts_are_not_discarded(self):
        source = bounded_valence_chain_presentations(3)[0]
        baseline = abbreviate_disjoint_rings(source, 1)
        changed = V4LabeledTree(source.tree,
                               tuple((b, 2 if b == 'f2' else value)
                                     for b, value in source.labels))
        self.assertNotEqual(source.labelled_rooted_code(), changed.labelled_rooted_code())
        self.assertNotEqual(baseline.rooted_code(),
                            abbreviate_disjoint_rings(changed, 1).rooted_code())
        # Same total ring count, different distribution along the rooted path.
        moved = Abbreviated(baseline.labeled, (2, 1, 0), 1)
        self.assertNotEqual(baseline.rooted_code(), moved.rooted_code())

    def test_wrong_wedge_does_not_abbreviate_as_a_valid_cap(self):
        for sign in (-1, 1):
            source = bounded_valence_chain_presentations(3, sign)[0]
            vertices = list(source.tree.vertices)
            center = vertices[3]
            vertices[3] = Vertex(center.slots, (sign, 0, 0))
            changed = V4LabeledTree(WeightedPlaneTree(tuple(vertices)), source.labels)
            with self.assertRaises(AssertionError):
                decode_choices(changed, sign)
            ab = abbreviate_disjoint_rings(changed, sign)
            self.assertEqual(sum(ab.rings), 2)
            self.assertEqual(len(ab.labeled.tree.vertices), 6)
            self.assertNotEqual(ab.rooted_code(), expected_path(3, sign).rooted_code())


if __name__ == '__main__':
    unittest.main()
