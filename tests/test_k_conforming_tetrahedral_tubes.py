"""Local K-conforming tetrahedron stars, not a source ambient tangle."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.k_conforming_tetrahedral_tubes import (
    _clip_closed_box, _core_edge_link, _fan_tetrahedra,
    build_k_core_tubes, KEdgeConformingTube,
)
from arborescent.spatial_twist_blocks import Point3


class KConformingTubesTests(unittest.TestCase):
    def test_two_real_k_edges_are_3d_mesh_simplices(self):
        for sign in (-1,1):
            for weights in ((0,0),(1,-1),(0,2,-2),(-1,0,1,0)):
                with self.subTest(sign=sign,weights=weights):
                    atlas=build_k_core_tubes(
                        AtomicInstruction((2,),weights,0,sign))
                    self.assertEqual(atlas.status,
                        "TWO_DISJOINT_K_TUBES_NOT_SOURCE_AMBIENT")
                    self.assertEqual(atlas.first.source_arc,0)
                    self.assertEqual(atlas.second.source_arc,1)
                    self.assertEqual(atlas.first.core()[0].y,-6)
                    self.assertEqual(atlas.second.core()[0].y,6)
                    for tube in (atlas.first,atlas.second):
                        result=tube.audit()
                        self.assertEqual((result.vertices,result.tetrahedra),
                                         (10,12))
                        self.assertEqual((result.boundary_euler,
                                          result.boundary_components),(2,1))
                        self.assertEqual(result.volume,
                            (tube.core()[1].x-tube.core()[0].x)*
                            (2*tube.radius)**2)
                        self.assertEqual(_core_edge_link(tube.tetrahedra),(4,4))
                        self.assertTrue(any(
                            0 in tet and 5 in tet
                            for tet in tube.tetrahedra))
                        with self.assertRaisesRegex(NotImplementedError,
                                                     "UNVERIFIED_GEOMETRY"):
                            tube.require_source_ambient()
                    self.assertLess(atlas.first.bounds()[1][1],
                                    atlas.second.bounds()[1][0])
                    with self.assertRaisesRegex(NotImplementedError,
                                                 "UNVERIFIED_GEOMETRY"):
                        atlas.require_source_global_cut()

    def test_exact_rational_radii_and_reproducibility(self):
        atom=AtomicInstruction((),(0,1),0,1)
        for h in (Fraction(1,16),Fraction(1,8),Fraction(1,4)):
            obj=build_k_core_tubes(atom,h)
            second=build_k_core_tubes(atom,h)
            self.assertEqual(obj.first.vertices,second.first.vertices)
            self.assertEqual(obj.first.tetrahedra,second.first.tetrahedra)
            self.assertEqual(obj.second.vertices,second.second.vertices)
            self.assertEqual(obj.second.tetrahedra,second.second.tetrahedra)
            self.assertEqual(obj.first.audit().volume,
                             obj.second.audit().volume)
            self.assertEqual(obj.first.bounds()[2],(-h,h))
            self.assertGreater(obj.first.core()[1].x,obj.first.core()[0].x)

    def test_exact_segment_box_clipping_and_no_touching_shortcuts(self):
        p=lambda *q:Point3(*q)
        box=((Fraction(0),Fraction(2)),
             (Fraction(-1),Fraction(1)),
             (Fraction(-1),Fraction(1)))
        self.assertEqual(_clip_closed_box(p(-1,0,0),p(3,0,0),box),
                         (Fraction(1,4),Fraction(3,4)))
        self.assertEqual(_clip_closed_box(p(1,0,0),p(1,3,0),box),
                         (Fraction(0),Fraction(1,3)))
        self.assertIsNone(_clip_closed_box(p(-1,2,0),p(3,2,0),box))
        # A point touching the closed boundary is geometrically a hit.
        self.assertEqual(_clip_closed_box(p(-1,1,0),p(3,1,0),box),
                         (Fraction(1,4),Fraction(3,4)))
        with self.assertRaises(ValueError):
            _clip_closed_box(p(0,0,0),p(0,0,0),box)
        with self.assertRaises(ValueError):
            _clip_closed_box(p(0,0,0),p(1,0,0),
                             ((Fraction(0),Fraction(0)),)+box[1:])

    def test_mutated_tetrahedra_false_k_edges_and_source_status_rejected(self):
        obj=build_k_core_tubes(AtomicInstruction((),(1,0),0,1))
        tube=obj.first
        with self.assertRaises(ValueError):
            replace(tube,status="SOURCE_AMBIENT_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_CONWAY_SLIT")
        with self.assertRaises(ValueError):
            replace(obj,first=obj.second)
        with self.assertRaises(ValueError):
            replace(tube,source_arc=1)
        with self.assertRaises(ValueError):
            replace(tube,radius=Fraction(1,2))
        with self.assertRaises(ValueError):
            replace(tube,vertices=tube.vertices[::-1])
        with self.assertRaises(ValueError):
            replace(tube,tetrahedra=tube.tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(tube,tetrahedra=(tube.tetrahedra[0],)*12)
        with self.assertRaises(ValueError):
            _core_edge_link(((0,1,2,5),))
        with self.assertRaises(ValueError):
            build_k_core_tubes("not source atom")
        with self.assertRaises(ValueError):
            build_k_core_tubes(AtomicInstruction((),(0,0),0,1),
                               Fraction(0))
        with self.assertRaises(ValueError):
            build_k_core_tubes(AtomicInstruction((),(0,0),0,1),
                               Fraction(1,2))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_k_core_tubes(AtomicInstruction((),(0,0),1,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_k_core_tubes(AtomicInstruction((),(0,),0,1))


if __name__=="__main__":
    unittest.main()
