"""A5.2e1: global r=0 return arcs, exact collision checker, source-claim gate."""
import unittest
from dataclasses import replace
from fractions import Fraction
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.spatial_twist_blocks import Point3
from arborescent.atomic_return_strands import (
    _segment_contact, _segment_cube_contact, independent_pl_arc_audit,
    realize_zero_ring_return_strands,
)
from arborescent.marked_conway_geometry import MarkedConwayCube


class ExactGlobalSegmentAuditTests(unittest.TestCase):
    def test_independent_exact_3d_segment_intersections(self):
        p=lambda *t: Point3(*t)
        cases=(
            ((p(0,0,0),p(2,2,0),p(0,2,0),p(2,0,0)), "point",p(1,1,0)),
            ((p(0,0,0),p(2,2,0),p(0,2,1),p(2,0,1)), "none",None),
            ((p(0,0,0),p(2,0,0),p(1,0,0),p(3,0,0)), "overlap",None),
            ((p(0,0,0),p(1,0,0),p(1,0,0),p(2,0,0)), "point",p(1,0,0)),
            ((p(0,0,0),p(1,0,0),p(2,0,0),p(3,0,0)), "none",None),
            ((p(0,0,0),p(1,0,0),p(2,1,0),p(3,1,0)), "none",None),
            ((p(0,0,0),p(0,0,2),p(0,0,1),p(0,0,3)), "overlap",None),
            ((p(0,0,0),p(0,0,1),p(0,0,2),p(0,0,3)), "none",None),
            ((p(0,0,0),p(2,2,2),p(0,2,2),p(2,0,0)), "point",p(1,1,1)),
            ((p(0,0,0),p(2,2,2),p(0,2,3),p(2,0,1)), "none",None),
        )
        for (a,b,c,d),kind,point in cases:
            with self.subTest(case=(a,b,c,d)):
                self.assertEqual(_segment_contact(a,b,c,d),(kind,point))
                self.assertEqual(_segment_contact(c,d,a,b),(kind,point))
        h=Fraction(1,2)
        self.assertEqual(_segment_contact(p(0,0,0),p(1,0,0),
                                          p(h,-1,0),p(h,1,0)),
                         ("point",p(h,0,0)))
        with self.assertRaises(ValueError):
            _segment_contact(p(0,0,0),p(0,0,0),p(1,0,0),p(2,0,0))

    def test_solid_cube_clipping_not_just_mesh_edges(self):
        c=MarkedConwayCube(10, Port((),0))
        p=lambda *a: Point3(*a)
        self.assertEqual(_segment_cube_contact(p(0,0,0),p(8,-2,0),c),
                         ("point",p(8,-2,0)))
        self.assertEqual(_segment_cube_contact(p(9,0,0),p(13,0,0),c)[0],
                         "overlap")
        self.assertEqual(_segment_cube_contact(p(0,4,0),p(25,4,0),c)[0],
                         "none")
        self.assertEqual(_segment_cube_contact(p(0,0,0),p(7,0,0),c)[0],
                         "none")


class ReturnStrandCandidateTests(unittest.TestCase):
    def test_two_n_boundary_arcs_not_an_auxiliary_anchor(self):
        fixtures=((0,0),(1,0),(-1,2),(3,-2,0),(0,-1,2,-2),
                  (2,-2,1,0,-1),(0,)*6)
        for sign in (-1,1):
            for weights in fixtures:
                atom=AtomicInstruction((1,2),weights,0,sign)
                obj=realize_zero_ring_return_strands(atom)
                self.assertEqual(obj.status,"RETURN_STRANDS_SOURCE_UNVERIFIED")
                self.assertEqual(len(obj.cubes),len(weights))
                self.assertEqual(len(obj.embedded_arcs()),2*len(weights))
                self.assertEqual(tuple(c.port for c in obj.cubes),
                                 tuple(Port((1,2),i) for i in range(len(weights))))
                self.assertEqual(len(obj.return_rails),2)
                self.assertEqual(len(obj.internal),len(weights)-1)
                self.assertEqual(obj.entry.writhe,weights[0])
                self.assertEqual(tuple(c.local.writhe for c in obj.internal),
                                 weights[1:])
                self.assertEqual(obj.boundary_pairings(),
                                 independent_pl_arc_audit(obj.cubes,obj.embedded_arcs()))
                self.assertEqual(obj.canonical_coordinates(),
                                 realize_zero_ring_return_strands(atom).canonical_coordinates())
                seen=[p for pair in obj.boundary_pairings() for p in pair]
                self.assertEqual(len(seen),4*len(weights))
                self.assertEqual(len(set(seen)),4*len(weights))
                self.assertEqual({p[0] for p in seen},
                                 {Port((1,2),i) for i in range(len(weights))})
                self.assertEqual(sum(len(c.local.word) for c in obj.internal)
                                 +len(obj.entry.word),
                                 sum(abs(w) for w in weights))
                with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                    obj.require_source_atomic_realization()

    def test_endpoint_connectivity_and_parity_are_computed_from_geometry(self):
        for w0 in (-4,-3,-2,-1,0,1,2,3,4):
            for w1 in (-2,0,2):
                obj=realize_zero_ring_return_strands(
                    AtomicInstruction((),(w0,w1),0,1))
                pairs=obj.boundary_pairings()
                self.assertEqual(pairs[0][0],(Port((),1),"NE"))
                self.assertEqual(pairs[1][0],(Port((),1),"SE"))
                if abs(w0)%2:
                    self.assertEqual(pairs[0][1],(Port((),0),"SW"))
                    self.assertEqual(pairs[1][1],(Port((),0),"NW"))
                else:
                    self.assertEqual(pairs[0][1],(Port((),0),"NW"))
                    self.assertEqual(pairs[1][1],(Port((),0),"SW"))
                self.assertEqual({p[0] for p in pairs[2:]},
                                 {(Port((),0),"NE"),(Port((),0),"SE")})
                self.assertEqual({p[1] for p in pairs[2:]},
                                 {(Port((),1),"NW"),(Port((),1),"SW")})

    def test_exhaustive_small_signed_sectors_and_disjointness(self):
        for n in (2,3):
            for weights in product((-2,-1,0,1,2),repeat=n):
                obj=realize_zero_ring_return_strands(
                    AtomicInstruction((),weights,0,1))
                self.assertEqual(len(obj.boundary_pairings()),2*n)
                self.assertEqual(len(obj.cubes),n)
                self.assertEqual(len(obj.internal),n-1)

    def test_malformed_return_paths_and_independent_collision_detection(self):
        obj=realize_zero_ring_return_strands(
            AtomicInstruction((),(2,-1,0),0,1))
        arcs=list(obj.embedded_arcs())
        # Force the first strand to enter the last closed Conway ball.
        first=list(arcs[0])
        first[1]=Point3(obj.cubes[-1].center_x,-2,0)
        arcs[0]=tuple(first)
        with self.assertRaisesRegex(ValueError,"sphere|intersect"):
            independent_pl_arc_audit(obj.cubes,tuple(arcs))
        arcs=list(obj.embedded_arcs())
        arcs[1]=arcs[0]
        with self.assertRaises(ValueError):
            independent_pl_arc_audit(obj.cubes,tuple(arcs))
        with self.assertRaises(ValueError):
            replace(obj,return_rails=(obj.return_rails[1],obj.return_rails[0]))
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,atom=AtomicInstruction((),(2,1,0),0,1))
        with self.assertRaises(ValueError):
            replace(obj,cubes=obj.cubes[:-1])
        with self.assertRaises(ValueError):
            replace(obj,internal=obj.internal[::-1])
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            realize_zero_ring_return_strands(AtomicInstruction((),(0,0),1,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            realize_zero_ring_return_strands(AtomicInstruction((),(0,),0,1))
        with self.assertRaises(ValueError):
            realize_zero_ring_return_strands("not an atom")


if __name__=="__main__":
    unittest.main()
