"""Independent counterexamples: encoding equality is not move/isotopy equality."""
import unittest
from dataclasses import replace
from arborescent.ring_abbreviation import Abbreviated, abbreviate_disjoint_rings
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def fixture(sign=1, placement=1, count=1):
    # Three individually V4-labelled free ports; caps in any base sector.
    caps = tuple(f'e{3*i}' for i in range(count))
    slots = list(('f0', 'f1', 'f2'))
    slots[placement:placement] = caps
    vertices = [Vertex(tuple(slots), (0,) * len(slots))]
    for i in range(count):
        a, b, c = (f'e{3*i+j}' for j in range(3))
        vertices.extend((Vertex((a,b,c), (0,sign,0)),
                         Vertex((b,), (2*sign,)), Vertex((c,), (2*sign,))))
    return V4LabeledTree(WeightedPlaneTree(tuple(vertices)),
                        (('f0',0), ('f1',1), ('f2',2)))


class AdversarialRingTests(unittest.TestCase):
    def test_distinct_placements_collapse_only_in_notation(self):
        for sign in (-1,1):
            sources = [fixture(sign, p) for p in range(4)]
            self.assertNotEqual(sources[0].labelled_rooted_code(),
                                sources[1].labelled_rooted_code())
            results = [abbreviate_disjoint_rings(s,sign) for s in sources]
            self.assertEqual(len({r.rooted_code() for r in results}),1)
            for s,r in zip(sources,results):
                self.assertEqual(len(s.tree.vertices),4)  # source never mutated
                self.assertEqual(len(r.labeled.tree.vertices),1)
                self.assertEqual(r.labeled.labels,s.labels)
                # Bare cap deletion fails to retain the ring invariant.
                self.assertNotEqual(r.rooted_code(),
                                    Abbreviated(r.labeled,(0,),sign).rooted_code())

    def test_multiplicity_and_every_boundary_label_are_detectable(self):
        results = [abbreviate_disjoint_rings(fixture(count=n),1)
                   for n in (0,1,2,3)]
        self.assertEqual(len({r.rooted_code() for r in results}),4)
        for n,r in enumerate(results):
            self.assertEqual(r.rings,(n,))
        r = results[1]
        for port in ('f0','f1','f2'):
            labels = tuple((p,v ^ 3 if p == port else v)
                           for p,v in r.labeled.labels)
            changed = replace(r,labeled=V4LabeledTree(r.labeled.tree,labels))
            self.assertNotEqual(changed.rooted_code(),r.rooted_code())

    def test_distinct_bases_keep_separate_counts_and_nonzero_sectors(self):
        source = V4LabeledTree(WeightedPlaneTree((
            Vertex(('f0','e0','e1','f1'),(0,0,7,0)),
            Vertex(('e1','e4','f2','f3'),(0,0,9,0)),
            Vertex(('e0','e2','e3'),(0,1,0)),
            Vertex(('e2',),(2,)),Vertex(('e3',),(2,)),
            Vertex(('e4','e5','e6'),(0,1,0)),
            Vertex(('e5',),(2,)),Vertex(('e6',),(2,)),
        )),tuple((f'f{i}',i) for i in range(4)))
        result = abbreviate_disjoint_rings(source,1)
        self.assertEqual(result.rings,(1,1))
        self.assertEqual(tuple(v.weight for v in result.labeled.tree.vertices),(7,9))
        self.assertEqual(result.labeled.labels,source.labels)
        # Equal total multiplicity does not identify different base assignments.
        self.assertNotEqual(result.rooted_code(),replace(result,rings=(2,0)).rooted_code())

    def test_opposite_sign_adjacent_centers_must_not_be_partially_erased(self):
        source = V4LabeledTree(WeightedPlaneTree((
            Vertex(('e0','e1','e2'),(0,1,0)),
            Vertex(('e0','e3','e4'),(0,-1,0)),
            Vertex(('e1',),(2,)),Vertex(('e2',),(2,)),
            Vertex(('e3',),(-2,)),Vertex(('e4',),(-2,)),
        )),())
        for sign in (-1,1):
            with self.subTest(sign=sign):
                with self.assertRaisesRegex(ValueError,'rings attached to rings'):
                    abbreviate_disjoint_rings(source,sign)

    def test_nonzero_base_sector_is_preserved_without_abbreviation(self):
        source = fixture()
        base = source.tree.vertices[0]
        for j in (0,1):
            sectors = tuple(5 if k == j else 0 for k in range(4))
            changed = replace(source,tree=WeightedPlaneTree((replace(base,sectors=sectors),
                                                          *source.tree.vertices[1:])))
            r = abbreviate_disjoint_rings(changed,1)
            self.assertEqual(r.labeled,changed)
            self.assertEqual(sum(r.rings),0)

if __name__ == '__main__':
    unittest.main()
