"""Exact physical pair-of-K-edges joined through ONE shared 3D corner."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.exact_tetrahedron_contacts import audit_embedded_tetrahedra
from arborescent.shared_corner_two_k_rays import (
    _assemble, build_shared_corner_two_k_directions
)


class SharedCornerTwoKDirectionsTests(unittest.TestCase):
    def test_one_actual_bent_k_corner_has_two_genuinely_welded_edges(self):
        model=build_shared_corner_two_k_directions(
            AtomicInstruction((),(0,0),0,1),0,1)
        top=model.audit()
        self.assertEqual(model.status,"TWO_K_RAYS_LOCAL_NOT_SOURCE_AMBIENT")
        self.assertEqual((top.vertices,top.tetrahedra),(35,72))
        self.assertEqual((top.boundary_components,top.boundary_euler),(1,2))
        self.assertEqual((len(model.incoming_bridge),
                          len(model.outgoing_bridge)),(12,12))
        self.assertEqual(audit_embedded_tetrahedra(
            model.vertices,model.tetrahedra).pair_count,2556)
        self.assertEqual(model.corner.original.corner,1)
        self.assertEqual(model.incoming.segment,0)
        self.assertEqual(model.outgoing.segment,1)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            model.require_source_full_k()

    def test_true_straight_k_continuation_one_nonoverlapping_corner(self):
        obj=build_shared_corner_two_k_directions(
            AtomicInstruction((),(0,0),0,-1),0,5)
        top=obj.audit()
        self.assertEqual((top.vertices,top.tetrahedra),(35,72))
        self.assertEqual((top.boundary_euler,top.boundary_components),(2,1))
        self.assertEqual(obj.incoming.segment,4)
        self.assertEqual(obj.outgoing.segment,5)
        self.assertEqual(obj.corner.original.corner,5)

    def test_mesh_status_faces_and_source_labels_cannot_fake_welding(self):
        obj=build_shared_corner_two_k_directions(
            AtomicInstruction((),(0,0),0,1),0,1)
        for mutation in (
            {"status":"SOURCE_PLUMBING_CERTIFIED"},
            {"vertices":obj.vertices[::-1]},
            {"tetrahedra":obj.tetrahedra[:-1]},
            {"incoming_bridge":obj.incoming_bridge[:-1]},
            {"outgoing_bridge":obj.outgoing_bridge[:-1]},
            {"incoming_square":obj.incoming_diamond},
            {"outgoing_square":obj.outgoing_diamond},
            {"incoming_rotation":4},
        ):
            with self.subTest(mutation=mutation):
                with self.assertRaises(ValueError):
                    replace(obj,**mutation)
        with self.assertRaises(ValueError):
            _assemble(obj.corner,obj.incoming,obj.outgoing,4,0,False,False)
        with self.assertRaises(ValueError):
            build_shared_corner_two_k_directions(
                AtomicInstruction((),(0,0),0,1),True,1)
        with self.assertRaises(ValueError):
            build_shared_corner_two_k_directions(
                AtomicInstruction((),(0,0),0,1),0,0)
        with self.assertRaises(ValueError):
            build_shared_corner_two_k_directions(
                AtomicInstruction((),(0,0),0,1),0,1,Fraction(0))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_shared_corner_two_k_directions(
                AtomicInstruction((),(0,0),1,1),0,1)


if __name__=="__main__":
    unittest.main()
