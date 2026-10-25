"""Exact positive-area truncated octahedral K corner interface tests."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.rational_k_corner_stars import build_rational_k_corner_star
from arborescent.truncated_k_corner_caps import (
    _canonical_geometry, _cap_disk, _link_cycle,
    truncate_source_k_corner, build_truncated_k_corner
)


class TruncatedKCornerTests(unittest.TestCase):
    def test_bent_physical_k_corners_with_real_positive_area_attachment(self):
        for sign in (-1,1):
            for weights,arc,corner in (
                ((0,0),0,1),
                ((1,-1),1,2),
                ((-1,0,1),4,1),
            ):
                with self.subTest(sign=sign,weights=weights):
                    atom=AtomicInstruction((),weights,0,sign)
                    obj=build_truncated_k_corner(atom,arc,corner)
                    before=obj.original.audit()
                    after=obj.audit()
                    self.assertEqual(obj.status,
                        "TRUNCATED_K_CORNER_CAP_NOT_GLUED")
                    self.assertEqual((len(obj.vertices),len(obj.tetrahedra)),
                                     (11,16))
                    self.assertEqual((after.vertices,after.tetrahedra),
                                     (11,16))
                    self.assertEqual((after.boundary_euler,
                                      after.boundary_components),(2,1))
                    self.assertEqual(_cap_disk(
                        obj.vertices,obj.cap_faces,obj.cap_center),(5,8,4))
                    self.assertGreater(after.volume,0)
                    self.assertEqual(after.volume,
                        Fraction(15,16)*before.volume)
                    self.assertEqual(_link_cycle(
                        obj.tetrahedra,obj.k_center,obj.cap_center),(4,4))
                    self.assertEqual(_link_cycle(
                        obj.tetrahedra,obj.k_center,obj.incoming_vertex),(6,6))
                    self.assertEqual(obj.vertices[obj.k_center],
                        obj.original.candidate.patches.strands.embedded_arcs()[
                            arc][corner])
                    self.assertEqual(
                        (obj.vertices,obj.tetrahedra,obj.cap_faces,
                         obj.incoming_vertex,obj.cap_center,obj.k_center),
                        _canonical_geometry(obj.original))
                    with self.assertRaisesRegex(
                        NotImplementedError,"UNVERIFIED_GEOMETRY"):
                        obj.require_source_tube_attachment()

    def test_straight_return_continuation_has_noncollapsed_outgoing_disk(self):
        atom=AtomicInstruction((),(0,0),0,1)
        source=build_rational_k_corner_star(atom,0,5)
        obj=truncate_source_k_corner(source)
        result=obj.audit()
        self.assertEqual((result.vertices,result.tetrahedra),(11,16))
        self.assertEqual(_link_cycle(
            obj.tetrahedra,obj.k_center,obj.cap_center),(4,4))
        self.assertEqual(_link_cycle(
            obj.tetrahedra,obj.k_center,obj.incoming_vertex),(4,4))
        # The diamond's 4 boundary vertices are distinct and do not
        # collapse onto the actual outgoing K-cap center.
        ring={v for face in obj.cap_faces for v in face}
        self.assertEqual(len(ring),5)
        self.assertEqual(len({obj.vertices[i] for i in ring}),5)

    def test_mutated_physical_mesh_and_fake_gluing_refused(self):
        obj=build_truncated_k_corner(
            AtomicInstruction((),(0,0),0,1),0,1)
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_CONWAY_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,vertices=obj.vertices[::-1])
        with self.assertRaises(ValueError):
            replace(obj,tetrahedra=obj.tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(obj,cap_faces=obj.cap_faces[:-1])
        with self.assertRaises(ValueError):
            replace(obj,cap_faces=tuple(tuple(reversed(t))
                                             for t in obj.cap_faces))
        with self.assertRaises(ValueError):
            replace(obj,cap_center=obj.incoming_vertex)
        with self.assertRaises(ValueError):
            replace(obj,k_center=1)
        with self.assertRaises(ValueError):
            truncate_source_k_corner("not a K corner star")
        with self.assertRaises(ValueError):
            _link_cycle(obj.tetrahedra,0,0)
        with self.assertRaises(ValueError):
            _cap_disk(obj.vertices,(),obj.cap_center)
        with self.assertRaisesRegex(NotImplementedError,
                                      "UNVERIFIED_GEOMETRY"):
            build_truncated_k_corner(
                AtomicInstruction((),(0,0),1,1),0,1)
        with self.assertRaisesRegex(NotImplementedError,
                                      "UNVERIFIED_GEOMETRY"):
            build_truncated_k_corner(
                AtomicInstruction((),(0,),0,1),0,1)


if __name__=="__main__":
    unittest.main()
