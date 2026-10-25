"""Positive-area INCOMING K-corner disk is a real 3D tetrahedral boundary."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.incoming_k_corner_caps import (
    build_incoming_k_corner_cap,truncate_incoming_k_corner,
    _incoming_geometry,
)
from arborescent.truncated_k_corner_caps import _cap_disk, _link_cycle


class IncomingKCornerCapsTests(unittest.TestCase):
    def test_real_bent_and_signed_source_k_corners(self):
        fixtures=(((0,0),0,1),((1,-1),1,2),((-1,0,1),4,1))
        for sign in (-1,1):
            for weights,arc,corner in fixtures:
                with self.subTest(sign=sign,weights=weights,arc=arc):
                    item=build_incoming_k_corner_cap(
                        AtomicInstruction((),weights,0,sign),arc,corner)
                    result=item.audit()
                    self.assertEqual(item.status,"INCOMING_K_CAP_NOT_TUBE_GLUED")
                    self.assertEqual((result.vertices,result.tetrahedra),(11,16))
                    self.assertEqual((result.boundary_components,
                                      result.boundary_euler),(1,2))
                    self.assertEqual(_cap_disk(
                        item.vertices,item.cap_faces,item.cap_center),(5,8,4))
                    self.assertEqual(_link_cycle(
                        item.tetrahedra,item.k_center,item.cap_center),(4,4))
                    self.assertEqual(_link_cycle(
                        item.tetrahedra,item.k_center,item.outgoing_vertex),(6,6))
                    self.assertEqual(result.volume,
                                     Fraction(15,16)*item.original.audit().volume)
                    self.assertEqual(
                        (item.vertices,item.tetrahedra,item.cap_faces,
                         item.outgoing_vertex,item.cap_center,item.k_center),
                        _incoming_geometry(item.original))
                    with self.assertRaisesRegex(
                        NotImplementedError,"UNVERIFIED_GEOMETRY"):
                        item.require_full_source_plumbing()

    def test_straight_k_vertex_opposite_ray_retains_4_cycle(self):
        item=build_incoming_k_corner_cap(
            AtomicInstruction((),(0,0),0,1),0,5)
        self.assertEqual((item.audit().vertices,item.audit().tetrahedra),(11,16))
        self.assertEqual(_link_cycle(
            item.tetrahedra,item.k_center,item.outgoing_vertex),(4,4))
        self.assertEqual(_link_cycle(
            item.tetrahedra,item.k_center,item.cap_center),(4,4))
        cap_points={item.vertices[j] for f in item.cap_faces for j in f}
        self.assertEqual(len(cap_points),5)

    def test_forged_geometry_source_and_interfaces_fail_closed(self):
        item=build_incoming_k_corner_cap(
            AtomicInstruction((),(0,0),0,1),0,1)
        for kwargs in (
            {"status":"VERIFIED_SOURCE_PLUMBING"},
            {"vertices":item.vertices[::-1]},
            {"tetrahedra":item.tetrahedra[:-1]},
            {"cap_faces":item.cap_faces[:-1]},
            {"cap_center":item.outgoing_vertex},
            {"k_center":item.cap_center},
        ):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    replace(item,**kwargs)
        with self.assertRaises(ValueError):
            truncate_incoming_k_corner("not an existing physical K vertex")
        with self.assertRaises(ValueError):
            build_incoming_k_corner_cap(AtomicInstruction((),(0,0),0,1),0,0)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_incoming_k_corner_cap(AtomicInstruction((),(0,0),1,1),0,1)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_incoming_k_corner_cap(AtomicInstruction((),(0,),0,1),0,1)


if __name__=="__main__":
    unittest.main()
