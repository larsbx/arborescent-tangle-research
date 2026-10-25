"""Test a truly two-ended K edge welded between two real corner 3-balls."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.two_ended_k_edge import (
    build_two_ended_k_edge, _assemble, TwoEndedKVertexEdge,
)
from arborescent.exact_tetrahedron_contacts import audit_embedded_tetrahedra


class TwoEndedKConformingEdgeTests(unittest.TestCase):
    def test_real_bent_to_bent_original_k_edge(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_two_ended_k_edge(atom,0,1)
        topology=obj.audit()
        self.assertEqual(obj.status,"TWO_CORNER_K_EDGE_NOT_WHOLE_SOURCE_AMBIENT")
        self.assertEqual((topology.vertices,topology.tetrahedra),(32,68))
        self.assertEqual((topology.boundary_euler,
                          topology.boundary_components),(2,1))
        self.assertEqual(len(obj.bridge_tetrahedra),12)
        self.assertEqual(len(obj.near_faces),4)
        self.assertEqual(len(obj.far_faces),4)
        self.assertEqual(obj.far.original.corner,2)
        self.assertEqual(obj.near.tube.segment,1)
        self.assertEqual(audit_embedded_tetrahedra(
            obj.vertices,obj.tetrahedra).pair_count,2278)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            obj.require_source_plumbing_certificate()

    def test_real_bent_to_straight_k_edge(self):
        atom=AtomicInstruction((),(0,0),0,-1)
        obj=build_two_ended_k_edge(atom,0,4)
        topology=obj.audit()
        self.assertEqual((topology.vertices,topology.tetrahedra),(32,68))
        self.assertEqual((topology.boundary_euler,
                          topology.boundary_components),(2,1))
        self.assertEqual(obj.far.original.corner,5)

    def test_adversarial_nonconforming_3d_faces_and_wrong_source_edge(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_two_ended_k_edge(atom,0,1)
        for attrs in (
            {"status":"SOURCE_PLUMBING_VERIFIED"},
            {"vertices":obj.vertices[::-1]},
            {"tetrahedra":obj.tetrahedra[:-1]},
            {"bridge_tetrahedra":obj.bridge_tetrahedra[:-1]},
            {"near_faces":obj.far_faces},
            {"far_faces":obj.near_faces},
            {"rotation":4},
        ):
            with self.subTest(attrs=attrs):
                with self.assertRaises(ValueError):
                    replace(obj,**attrs)
        with self.assertRaises(ValueError):
            _assemble(obj.near,obj.far,4,False)
        with self.assertRaises(ValueError):
            _assemble(obj.near,obj.far,0,1)
        with self.assertRaises(ValueError):
            build_two_ended_k_edge(atom,False,1)
        with self.assertRaises(ValueError):
            build_two_ended_k_edge(atom,0,0)
        with self.assertRaises(ValueError):
            build_two_ended_k_edge(atom,0,1,Fraction(0))


if __name__=="__main__":
    unittest.main()
