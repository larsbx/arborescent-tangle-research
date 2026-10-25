"""A5.2e2b4c2e3 — one genuine physical K tube joins TWO shared corners."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.shared_tube_between_k_vertices import (
    SharedTubeBetweenKVertices,build_one_shared_k_tube,_assemble,
)
from arborescent.exact_tetrahedron_contacts import audit_embedded_tetrahedra
from arborescent.connected_roof_ambient_slab import verify_connected_tetrahedral_ball


class OneSharedTubeKVertexTests(unittest.TestCase):
    def test_bent_to_bent_actual_k_edge_has_one_shared_mesh_tube(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_one_shared_k_tube(atom,0,1)
        # Construction already performs an independent 3,486-pair
        # exact physical collision audit. Recheck topology separately
        # instead of replaying the expensive collision proof twice.
        topology=verify_connected_tetrahedral_ball(obj.vertices,obj.tetrahedra)
        self.assertEqual(obj.status,"ONE_SHARED_K_TUBE_NOT_WHOLE_SOURCE_AMBIENT")
        self.assertEqual((topology.vertices,topology.tetrahedra),(40,84))
        self.assertEqual((topology.boundary_components,
                          topology.boundary_euler),(1,2))
        self.assertEqual((len(obj.near_bridge),len(obj.far_bridge)),(12,12))
        self.assertEqual(obj.left.original.corner,1)
        self.assertEqual(obj.tube.segment,1)
        self.assertEqual(obj.right.original.corner,2)
        self.assertEqual(audit_embedded_tetrahedra(
            obj.vertices,obj.tetrahedra).pair_count,3486)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            obj.require_source_ambient()

    def test_bent_to_straight_true_source_k_edge(self):
        obj=build_one_shared_k_tube(
            AtomicInstruction((),(0,0),0,-1),0,4)
        self.assertEqual(len(obj.vertices),40)
        self.assertEqual(len(obj.tetrahedra),84)
        self.assertEqual(obj.left.original.corner,4)
        self.assertEqual(obj.right.original.corner,5)
        topology=verify_connected_tetrahedral_ball(
            obj.vertices,obj.tetrahedra)
        self.assertEqual((topology.boundary_components,
                          topology.boundary_euler),(1,2))

    def test_corrupt_shared_tube_faces_source_status_and_geometry_fail(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_one_shared_k_tube(atom,0,1)
        for attrs in (
            {"status":"SOURCE_FIG1210_CERTIFIED"},
            {"vertices":obj.vertices[::-1]},
            {"tetrahedra":obj.tetrahedra[:-1]},
            {"near_bridge":obj.near_bridge[:-1]},
            {"far_bridge":obj.far_bridge[:-1]},
            {"near_square":obj.far_square},
            {"far_diamond":obj.near_diamond},
            {"left_rotation":4},
        ):
            with self.subTest(attrs=attrs):
                with self.assertRaises(ValueError):
                    replace(obj,**attrs)
        with self.assertRaises(ValueError):
            _assemble(obj.left,obj.tube,obj.right,4,0,False,False)
        with self.assertRaises(ValueError):
            build_one_shared_k_tube(atom,False,1)
        with self.assertRaises(ValueError):
            build_one_shared_k_tube(atom,0,0)
        with self.assertRaises(ValueError):
            build_one_shared_k_tube(atom,0,1,Fraction(0))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_one_shared_k_tube(
                AtomicInstruction((),(0,0),1,1),0,1)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_one_shared_k_tube(
                AtomicInstruction((),(0,),0,1),0,1)


if __name__=="__main__":
    unittest.main()
