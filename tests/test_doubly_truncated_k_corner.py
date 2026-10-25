"""A5.2e2b4c2e1: two disjoint real cap disks on ONE shared K 3-ball."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.doubly_truncated_k_corner import (
    C, build_doubly_truncated_k_corner, truncate_both_k_rays,
    _geometry, _clip_clean,
)
from arborescent.incoming_k_corner_caps import build_incoming_k_corner_cap
from arborescent.truncated_k_corner_caps import (
    _cap_disk, _link_cycle, build_truncated_k_corner,
)
from arborescent.exact_tetrahedron_contacts import audit_embedded_tetrahedra


class DoublyTruncatedKCornerTests(unittest.TestCase):
    def test_signed_bent_and_straight_source_k_corners(self):
        fixtures=(
            ((0,0),0,1),
            ((0,0),0,5),
            ((1,-1),1,2),
            ((-1,0,1),4,1),
        )
        for sign in (-1,1):
            for weights,arc,corner in fixtures:
                with self.subTest(sign=sign,weights=weights,arc=arc,corner=corner):
                    atom=AtomicInstruction((3,),weights,0,sign)
                    obj=build_doubly_truncated_k_corner(atom,arc,corner)
                    result=obj.audit()
                    self.assertEqual(obj.status,"DUAL_CAP_K_CORE_NOT_EDGE_GLUED")
                    self.assertEqual((result.vertices,result.tetrahedra),(15,24))
                    self.assertEqual((result.boundary_components,
                                      result.boundary_euler),(1,2))
                    self.assertEqual(_cap_disk(
                        obj.vertices,obj.incoming_faces,obj.incoming_center),(5,8,4))
                    self.assertEqual(_cap_disk(
                        obj.vertices,obj.outgoing_faces,obj.outgoing_center),(5,8,4))
                    self.assertEqual(_link_cycle(
                        obj.tetrahedra,0,obj.incoming_center),(4,4))
                    self.assertEqual(_link_cycle(
                        obj.tetrahedra,0,obj.outgoing_center),(4,4))
                    self.assertFalse(
                        {i for t in obj.incoming_faces for i in t}
                        & {i for t in obj.outgoing_faces for i in t})
                    self.assertEqual(
                        result.volume,obj.original.audit().volume*Fraction(63,64))
                    self.assertEqual(
                        audit_embedded_tetrahedra(
                            obj.vertices,obj.tetrahedra).pair_count,276)
                    self.assertEqual(
                        (obj.vertices,obj.tetrahedra,obj.incoming_faces,
                         obj.outgoing_faces,obj.incoming_center,
                         obj.outgoing_center),
                        _geometry(obj.original))
                    self.assertEqual(C,Fraction(3,4))
                    with self.assertRaisesRegex(
                        NotImplementedError,"UNVERIFIED_GEOMETRY"):
                        obj.require_global_source_ambient()

    def test_naive_old_incoming_plus_outgoing_corner_solids_overlap_in_volume(self):
        atom=AtomicInstruction((),(0,0),0,1)
        outgoing=build_truncated_k_corner(atom,0,1)
        incoming=build_incoming_k_corner_cap(atom,0,1)
        shared=build_doubly_truncated_k_corner(atom,0,1)
        base=shared.original.audit().volume
        # Both previous single-capped convex solids were subsets
        # of the SAME full source octahedron. Inclusion-exclusion
        # gives a lower bound 7/8 of the original 3D volume of
        # OVERLAP. A naive union double-counts this physical region.
        overlap_lower=outgoing.audit().volume+incoming.audit().volume-base
        self.assertEqual(overlap_lower,Fraction(7,8)*base)
        self.assertGreater(overlap_lower,0)
        self.assertLess(shared.audit().volume,base)
        self.assertNotEqual(shared.audit().volume,
                            incoming.audit().volume+outgoing.audit().volume)

    def test_exact_two_plane_polygon_clipping_and_disjoint_caps(self):
        one=Fraction(1)
        zero=Fraction(0)
        facet=((one,zero,zero),(zero,one,zero),(zero,zero,one))
        x=_clip_clean(facet,0,1)
        xy=_clip_clean(x,1,1)
        self.assertEqual(len(x),4)
        self.assertEqual(len(xy),5)
        self.assertTrue(all(v[0]<=C and v[1]<=C for v in xy))
        self.assertEqual(_clip_clean(((one,zero,zero),),0,1),())
        for args in ((3,1),(0,0),(1,2)):
            with self.assertRaises(ValueError):
                _clip_clean(facet,*args)

    def test_corrupted_source_mesh_and_fabricated_geometry_fail_closed(self):
        atom=AtomicInstruction((),(0,0),0,1)
        obj=build_doubly_truncated_k_corner(atom,0,1)
        for changes in (
            {"status":"SOURCE_ISOTOPY_VERIFIED"},
            {"vertices":obj.vertices[::-1]},
            {"tetrahedra":obj.tetrahedra[:-1]},
            {"incoming_faces":obj.outgoing_faces},
            {"outgoing_center":obj.incoming_center},
            {"incoming_center":0},
        ):
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    replace(obj,**changes)
        with self.assertRaises(ValueError):
            truncate_both_k_rays("not a K corner")
        with self.assertRaises(ValueError):
            build_doubly_truncated_k_corner(atom,0,0)
        with self.assertRaisesRegex(
            NotImplementedError,"UNVERIFIED_GEOMETRY"
        ):
            build_doubly_truncated_k_corner(
                AtomicInstruction((),(0,0),1,1),0,1)
        with self.assertRaisesRegex(
            NotImplementedError,"UNVERIFIED_GEOMETRY"
        ):
            build_doubly_truncated_k_corner(
                AtomicInstruction((),(0,),0,1),0,1)


if __name__=="__main__":
    unittest.main()
