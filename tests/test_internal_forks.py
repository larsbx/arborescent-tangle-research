import unittest
from arborescent.internal_forks import (
    Tip,Port,InternalFork,compile_fork,generate,tangle_screen,
    compare_partial_orbits,elementary_states,
)
from arborescent.partial_moves import Move,apply,explore
from arborescent.weighted_trees import Vertex,WeightedPlaneTree,enumerate_specs,to_weighted_tree
from arborescent.tangle_labels import label_free_ports

class InternalForkTests(unittest.TestCase):
    def witness(self):
        return compile_fork(InternalFork((InternalFork((Tip(3),Port())),Tip(2)),3))

    def test_exact_old_witness_and_new_equivalence(self):
        s=self.witness()
        self.assertEqual(s.tree,WeightedPlaneTree((
            Vertex(('f0','e0','e2'),(3,0,0)),
            Vertex(('e0','e1','f1'),(0,0,0)),
            Vertex(('e1',),(3,)),Vertex(('e2',),(2,)),
        )))
        t=apply(s,Move('F3prime',0,'e0','before'))
        self.assertNotIn(t,elementary_states(s))
        a,b,traces=compare_partial_orbits((s,t))
        self.assertEqual((a,b),(2,1))
        self.assertTrue(traces)
        self.assertTrue(all(trace.verify() for trace in traces))
        for state in explore(s).predecessors:
            self.assertTrue(tangle_screen(state).passed)
        self.assertTrue(tangle_screen(s).passed)
        self.assertIn(('N',3),tangle_screen(s,'-').failures)
        self.assertTrue(s.tree.wfs_failures()) # F is inappropriate here
        self.assertIn(s,generate(2))
        renamed=label_free_ports(WeightedPlaneTree(tuple(
            Vertex(tuple('e'+str(int(x[1:])+100) if x.startswith('e') else x
                         for x in v.slots[1:]+v.slots[:1]),
                   v.sectors[1:]+v.sectors[:1])
            for v in s.tree.vertices[::-1])))
        self.assertEqual(s.labelled_rooted_code(),renamed.labelled_rooted_code())
        self.assertTrue(tangle_screen(renamed).passed)

    def test_generation_counts_independent_formula_and_replays(self):
        # Catalan(n)*choose(n+1,p)*m**(n+1-p)*D**n.
        from math import comb
        from arborescent.prototype import catalan
        for n,expected in ((1,(8,6)),(2,(96,50))):
            seeds=generate(n,tips=(3,),weights=(0,3),extra_ports=1)
            self.assertEqual(len(seeds),catalan(n)*comb(n+1,1)*4**n)
            self.assertEqual(len(set(seeds)),len(seeds))
            self.assertTrue(all(tangle_screen(s).passed for s in seeds))
            a,b,traces=compare_partial_orbits(seeds)
            self.assertEqual((a,b),expected)
            self.assertTrue(all(t.verify() for t in traces))
        with self.assertRaises(ValueError):
            compare_partial_orbits((self.witness(),),max_states=1)

    def test_a1_matched_binary_control(self):
        for n in (1,2):
            seeds=generate(n,tips=(2,3),weights=(0,),extra_ports=0)
            old=[label_free_ports(to_weighted_tree(s)) for s in enumerate_specs(n) ]
            old_codes={s.labelled_rooted_code() for s in old}
            self.assertTrue(all(s.labelled_rooted_code() in old_codes for s in seeds))
            a,b,traces=compare_partial_orbits(seeds)
            self.assertEqual(a,b)
            self.assertFalse(traces)

    def test_guards_and_exception_status(self):
        for kwargs in ({'forks':True},{'forks':1,'extra_ports':3},
                       {'forks':1,'tips':(True,)},{'forks':1,'weights':(0,0)}):
            with self.assertRaises(ValueError): generate(**kwargs)
        for spec in (Tip(3),InternalFork((Tip(3),)),
                     InternalFork((Tip(3),Tip(3)),0,1),
                     InternalFork((Tip(True),Tip(3)))):
            with self.assertRaises(ValueError): compile_fork(spec)
        single=label_free_ports(WeightedPlaneTree((Vertex(('f0',),(0,)),)))
        self.assertFalse(tangle_screen(single).passed)
        self.assertTrue(tangle_screen(single).unresolved)
        with self.assertRaises(ValueError): tangle_screen(single,'bad')

    def test_tangle_free_end_sticks_and_sign_guards(self):
        def stick(weights):
            vertices=[Vertex(('f0','e0','e1'),(0,0,0)),Vertex(('e1',),(3,))]
            for j,w in enumerate(weights):
                slots=(f'e{0 if j==0 else j+1}',f'e{j+2}' if j+1<len(weights) else 'f1')
                vertices.append(Vertex(slots,(w,0)))
            return label_free_ports(WeightedPlaneTree(tuple(vertices)))
        for ws in ((3,0),(3,-1),(3,-3,0)):
            self.assertTrue(tangle_screen(stick(ws)).passed)
        for ws,kind in (((0,3),'S-zero'),((3,3),'S-sign'),((1,-3),'S-end')):
            self.assertIn(kind,{k for k,_ in tangle_screen(stick(ws)).failures})
        negative=compile_fork(InternalFork((Tip(-2),Tip(3)),3))
        self.assertIn('P',{k for k,_ in tangle_screen(negative).failures})
        self.assertTrue(tangle_screen(negative,'-').passed)
