"""A5.2e2b4c1: genuinely tetrahedral local 3D slit/reglue, source-unverified."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.exterior_conway_collars import build_exterior_collared_band
from arborescent.local_plumbing_slit import (
    LocalPlumbingPatch, CutVertex, collar_plumbing_chart,
    product_neighborhood, slit_product, _vec, _cross,
)
from arborescent.spatial_twist_blocks import Point3


class ExactLocalPlumbingSlitTests(unittest.TestCase):
    def _fixture(self, weights=(0,0)):
        band=build_exterior_collared_band(
            AtomicInstruction((),weights,0,1))
        patch=collar_plumbing_chart(band.collars[0])
        neighborhood=product_neighborhood(patch)
        return band,patch,neighborhood

    def test_patch_is_actual_pair_of_faces_and_two_transverse_arrows(self):
        for weights in ((0,0),(1,0),(-1,2),(2,-2,0)):
            band,patch,solid=self._fixture(weights)
            physical={frozenset(face) for face in band.collars[0].triangles()}
            self.assertTrue(all(frozenset(face) in physical
                                for face in patch.triangles()))
            self.assertEqual(len(patch.corners),4)
            self.assertEqual(_cross(_vec(patch.origin,patch.u_tip),
                                    _vec(patch.origin,patch.v_tip))[:2],
                             (0,0))
            self.assertGreater(_cross(_vec(patch.origin,patch.u_tip),
                                     _vec(patch.origin,patch.v_tip))[2],0)
            with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                patch.require_source_plumbing_patch()
            with self.assertRaises(ValueError):
                replace(patch,u_tip=patch.v_tip)
            with self.assertRaises(ValueError):
                replace(patch,corners=patch.corners[::-1])
            with self.assertRaises(ValueError):
                replace(patch,status="SOURCE_VERIFIED")

    def test_twelve_real_tetrahedra_form_exact_local_3ball(self):
        _,patch,solid=self._fixture()
        self.assertEqual(len(solid.vertices),12)
        self.assertEqual(len(solid.tetrahedra),12)
        self.assertEqual(solid.half_thickness,Fraction(1,16))
        self.assertEqual(solid.midplane_faces(),((4,5,6),(4,6,7)))
        solid.audit()
        for h in (Fraction(1,8),Fraction(1,16),Fraction(1,32)):
            obj=product_neighborhood(patch,h)
            self.assertEqual(len(obj.tetrahedra),12)
            obj.audit()
            self.assertEqual(obj.vertices[4].z,patch.corners[0].z)
            self.assertEqual(obj.vertices[0].z,patch.corners[0].z-h)
            self.assertEqual(obj.vertices[8].z,patch.corners[0].z+h)
        with self.assertRaises(ValueError):
            product_neighborhood(patch,Fraction(0))
        with self.assertRaises(ValueError):
            product_neighborhood(patch,Fraction(-1,8))
        with self.assertRaises(ValueError):
            product_neighborhood(patch,1)
        with self.assertRaises(ValueError):
            replace(solid,vertices=solid.vertices[::-1])
        with self.assertRaises(ValueError):
            replace(solid,tetrahedra=solid.tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(solid,status="SOURCE_3D_CERTIFIED")

    def test_two_distinct_cut_sides_and_exact_reverse_reglue(self):
        for sector in ((0,0),(-1,1),(1,0,-2)):
            _,patch,solid=self._fixture(sector)
            cut=slit_product(solid)
            self.assertTrue(cut.verify())
            self.assertEqual(cut.status,"LOCAL_SLIT_REGLUE_ONLY")
            self.assertEqual(len(cut.cut_tetrahedra),12)
            self.assertEqual(len(cut.lower_cut),2)
            self.assertEqual(len(cut.upper_cut),2)
            lower={v for f in cut.lower_cut for v in f}
            upper={v for f in cut.upper_cut for v in f}
            self.assertEqual(len(lower),4)
            self.assertEqual(len(upper),4)
            self.assertFalse(lower & upper)
            self.assertEqual({v.original_index for v in lower},
                             {4,5,6,7})
            self.assertEqual({v.original_index for v in upper},
                             {4,5,6,7})
            # Cut face normal orientations are opposite in true 3D.
            low_normal=_cross(_vec(*patch.triangles()[0][:2]),
                              _vec(patch.triangles()[0][0],
                                   patch.triangles()[0][2]))
            self.assertGreater(low_normal[2],0)
            # Retagging only changes incidence, not the spatial
            # coordinates of any tetrahedron or plane face.
            self.assertEqual(tuple(tuple(v.original_index for v in t)
                                   for t in cut.cut_tetrahedra),
                             solid.tetrahedra)
            with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                cut.require_source_conway_slitting()

    def test_forged_slit_boundary_orientation_and_identity_rejected(self):
        _,patch,solid=self._fixture((0,0))
        cut=slit_product(solid)
        bad_tets=list(cut.cut_tetrahedra)
        first=list(bad_tets[0])
        first[0]=CutVertex(first[0].original_index,"UPPER")
        bad_tets[0]=tuple(first)
        with self.assertRaises(ValueError):
            replace(cut,cut_tetrahedra=tuple(bad_tets))
        wrong_upper=list(cut.upper_cut)
        wrong_upper[0]=cut.lower_cut[0]
        with self.assertRaises(ValueError):
            replace(cut,upper_cut=tuple(wrong_upper))
        with self.assertRaises(ValueError):
            replace(cut,status="SOURCE_ISOTOPY")
        with self.assertRaises(ValueError):
            CutVertex(4,"WRONG")
        with self.assertRaises(ValueError):
            CutVertex(True,"LOWER")
        with self.assertRaises(ValueError):
            slit_product("not a 3D manifold")
        with self.assertRaises(ValueError):
            collar_plumbing_chart("not a collar")


if __name__=="__main__":
    unittest.main()
