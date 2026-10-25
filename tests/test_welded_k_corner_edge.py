"""Exact geometric tetrahedron contact and welded K corner-to-edge tests."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.exact_tetrahedron_contacts import (
    tetrahedron_intersection_vertices,audit_embedded_tetrahedra,
)
from arborescent.welded_k_corner_edge import (
    build_welded_k_corner_to_edge,_clip_segment_tet,_as_meshes,
)
from arborescent.spatial_twist_blocks import Point3


class ExactTetrahedralContactTests(unittest.TestCase):
    def test_closed_contact_at_exact_faces_edges_vertices(self):
        p=lambda *a:Point3(*a)
        a=(p(0,0,0),p(2,0,0),p(0,2,0),p(0,0,2))
        face=(p(0,0,0),p(0,2,0),p(2,0,0),p(0,0,-2))
        self.assertEqual(set(tetrahedron_intersection_vertices(a,face)),
                         {a[0],a[1],a[2]})
        result=audit_embedded_tetrahedra(
            (p(0,0,0),p(2,0,0),p(0,2,0),p(0,0,2),p(0,0,-2)),
            ((0,1,2,3),(0,2,1,4)))
        self.assertEqual(result.pair_count,1)
        self.assertEqual(result.legitimate_face_pairs,1)
        disjoint=(p(4,0,0),p(6,0,0),p(4,2,0),p(4,0,2))
        self.assertFalse(tetrahedron_intersection_vertices(a,disjoint))
        near=Fraction(1,10**12)
        self.assertFalse(tetrahedron_intersection_vertices(
            a,(p(2+near,0,0),p(3,0,0),p(3,1,0),p(3,0,1))))
        vertex=(p(0,0,0),p(-2,0,0),p(0,-2,0),p(0,0,-2))
        self.assertEqual(tetrahedron_intersection_vertices(a,vertex),
                         (p(0,0,0),))

    def test_disjoint_mesh_metadata_does_not_hide_3d_overlap(self):
        p=lambda *a:Point3(*a)
        v=(p(0,0,0),p(2,0,0),p(0,2,0),p(0,0,2),
           p(Fraction(1,4),Fraction(1,4),Fraction(1,4)),
           p(1,Fraction(1,4),Fraction(1,4)),
           p(Fraction(1,4),1,Fraction(1,4)),
           p(Fraction(1,4),Fraction(1,4),1))
        with self.assertRaisesRegex(ValueError,"cross or overlap"):
            audit_embedded_tetrahedra(v,((0,1,2,3),(4,5,6,7)))
        with self.assertRaises(ValueError):
            audit_embedded_tetrahedra(v,((0,1,2,3),)*2)
        with self.assertRaises(ValueError):
            audit_embedded_tetrahedra((),())
        with self.assertRaises(ValueError):
            tetrahedron_intersection_vertices(
                (p(0,0,0),)*4,tuple(v[:4]))
        with self.assertRaises(ValueError):
            audit_embedded_tetrahedra(v,((0,1,2,100),))

    def test_exact_seg_to_closed_tet_clipping(self):
        p=lambda *a:Point3(*a)
        tetra=(p(0,0,0),p(2,0,0),p(0,2,0),p(0,0,2))
        self.assertEqual(_clip_segment_tet(
            p(-1,0,0),p(3,0,0),tetra),(Fraction(1,4),Fraction(3,4)))
        self.assertIsNone(_clip_segment_tet(
            p(3,3,3),p(4,3,3),tetra))
        with self.assertRaises(ValueError):
            _clip_segment_tet(p(0,0,0),p(0,0,0),tetra)


class WeldedKCornerToTubeTests(unittest.TestCase):
    def test_actual_bent_corner_conforming_bridge_on_existing_k(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_welded_k_corner_to_edge(atom,0,1)
        self.assertEqual(obj.status,"ONE_K_CORNER_EDGE_GLUE_NOT_SOURCE_AMBIENT")
        result=obj.audit()
        self.assertEqual(result.vertices,21)
        self.assertEqual(result.tetrahedra,40)
        self.assertEqual((result.boundary_components,result.boundary_euler),(1,2))
        self.assertEqual(len(obj.bridge_tetrahedra),12)
        self.assertEqual(len(obj.near_faces),4)
        self.assertEqual(len(obj.far_faces),4)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            obj.require_source_ambient()

    def test_actual_straight_continuation_on_existing_k(self):
        atom=AtomicInstruction((),(0,0),0,-1)
        obj=build_welded_k_corner_to_edge(atom,0,5)
        result=obj.audit()
        self.assertEqual((result.vertices,result.tetrahedra),(21,40))
        self.assertEqual(len(obj.near_faces),len(obj.far_faces))

    def test_wrong_seam_and_status_fail_closed(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_welded_k_corner_to_edge(atom,0,1)
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,bridge_tetrahedra=obj.bridge_tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(obj,near_faces=obj.far_faces)
        with self.assertRaises(ValueError):
            replace(obj,vertices=obj.vertices[::-1])
        with self.assertRaises(ValueError):
            _as_meshes(obj.corner,obj.tube,4)
        with self.assertRaises(ValueError):
            build_welded_k_corner_to_edge(atom,False,1)
        with self.assertRaises(ValueError):
            build_welded_k_corner_to_edge(atom,0,1,Fraction(0))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_welded_k_corner_to_edge(
                AtomicInstruction((),(0,0),1,1),0,1)


if __name__=="__main__":
    unittest.main()
