"""A5.2e2b4c2a: all physical local 3D slits remain a DISCONNECTED atlas."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.disjoint_local_slit_atlas import (
    DisjointLocalSlitAtlas, LocalSlitAtlasMetrics,
    _bounds, build_disjoint_local_slit_atlas,
)
from arborescent.local_plumbing_slit import slit_product


class DisjointLocalSlitAtlasTests(unittest.TestCase):
    def test_every_candidate_port_has_distinct_physical_product_3ball(self):
        for sign in (-1,1):
            for weights in ((0,0),(1,0),(-1,1,0),(0,-1,2)):
                with self.subTest(sign=sign,weights=weights):
                    atom=AtomicInstruction((2,),weights,0,sign)
                    atlas=build_disjoint_local_slit_atlas(atom)
                    result=atlas.audit()
                    n=len(weights)
                    self.assertEqual(result.local_product_components,n)
                    self.assertEqual(result.local_cut_3balls,2*n)
                    self.assertEqual(result.tetrahedra_before_slit,12*n)
                    self.assertEqual(result.tetrahedra_after_slit,12*n)
                    self.assertEqual(result.genuine_conway_boundary_spheres,0)
                    self.assertEqual(result.knot_piercings_of_cut_disks,0)
                    self.assertEqual(result.total_local_volume,
                                     Fraction(18*n,16))
                    self.assertEqual(result.minimum_longitudinal_gap,
                        min(Fraction(1,2)+2*abs(w) for w in weights[1:]))
                    self.assertEqual(tuple(s.port for s in atlas.sites),
                                     tuple(Port((2,),i) for i in range(n)))
                    self.assertEqual(atlas.status,
                        "DISJOINT_LOCAL_ATLAS_NOT_GLOBAL_MANIFOLD")
                    self.assertTrue(all(site.slit.verify() for site in atlas.sites))
                    with self.assertRaisesRegex(NotImplementedError,
                                                 "UNVERIFIED_GEOMETRY"):
                        atlas.require_global_source_slitting()

    def test_exact_finite_rational_thickness_and_canonical_replay(self):
        atom=AtomicInstruction((),(0,1),0,1)
        for h in (Fraction(1,16),Fraction(1,8),Fraction(1,4)):
            atlas=build_disjoint_local_slit_atlas(atom,h)
            self.assertEqual(atlas.audit().total_local_volume,18*len(
                atlas.sites)*h)
            self.assertEqual(atlas.canonical(),build_disjoint_local_slit_atlas(
                atom,h).canonical())
            for site in atlas.sites:
                self.assertEqual(site.solid.half_thickness,h)
                self.assertEqual(site.slit,slit_product(site.solid))
                self.assertEqual(site.spatial_bounds()[2],
                                 (Fraction(4)-h,Fraction(4)+h))

    def test_fail_closed_against_fake_global_sphere_or_cut(self):
        atlas=build_disjoint_local_slit_atlas(
            AtomicInstruction((),(0,0,1),0,-1))
        with self.assertRaises(ValueError):
            replace(atlas,status="SOURCE_ISOTOPY_VERIFIED")
        with self.assertRaises(ValueError):
            replace(atlas,sites=atlas.sites[::-1])
        with self.assertRaises(ValueError):
            replace(atlas,sites=atlas.sites[:-1])
        with self.assertRaises(ValueError):
            replace(atlas,thickness=Fraction(1,2))
        with self.assertRaises(ValueError):
            replace(atlas,thickness=Fraction(1,8))
        with self.assertRaises(ValueError):
            replace(atlas,sites=(
                replace(atlas.sites[0],port=Port((),1)),
                *atlas.sites[1:]))
        with self.assertRaises(ValueError):
            replace(atlas.sites[0],solid=atlas.sites[1].solid)
        with self.assertRaises(ValueError):
            replace(atlas.sites[0].slit,
                    cut_tetrahedra=atlas.sites[0].slit.cut_tetrahedra[:-1])
        with self.assertRaises(ValueError):
            LocalSlitAtlasMetrics(1,12,12,2,4,0,Fraction(0),Fraction(1))
        with self.assertRaises(ValueError):
            _bounds(())
        with self.assertRaises(ValueError):
            build_disjoint_local_slit_atlas("not an atomic instruction")
        with self.assertRaises(ValueError):
            build_disjoint_local_slit_atlas(
                AtomicInstruction((),(0,0),0,1),Fraction(0))
        with self.assertRaises(ValueError):
            build_disjoint_local_slit_atlas(
                AtomicInstruction((),(0,0),0,1),Fraction(1,2))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_disjoint_local_slit_atlas(
                AtomicInstruction((),(0,0),1,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_disjoint_local_slit_atlas(
                AtomicInstruction((),(0,),0,1))

    def test_two_cut_boundary_spheres_per_patch_have_zero_k_punctures(self):
        atlas=build_disjoint_local_slit_atlas(
            AtomicInstruction((),(1,-2,0),0,1))
        k_arcs=atlas.candidate.patches.strands.embedded_arcs()
        self.assertEqual(len(k_arcs),6)
        self.assertEqual(len(atlas.sites),3)
        all_k_z=[Fraction(p.z) for arc in k_arcs for p in arc]
        self.assertLessEqual(max(all_k_z),Fraction(1))
        for site in atlas.sites:
            minimum_z=site.spatial_bounds()[2][0]
            self.assertGreaterEqual(minimum_z,Fraction(15,4))
            self.assertGreater(minimum_z,max(all_k_z))
            self.assertEqual(len(site.slit.lower_cut),2)
            self.assertEqual(len(site.slit.upper_cut),2)
            lower={v for f in site.slit.lower_cut for v in f}
            upper={v for f in site.slit.upper_cut for v in f}
            self.assertFalse(lower & upper)
        # DO NOT claim the 2n local product boundary spheres are
        # the source's desired FOUR-punctured Conway boundaries.
        self.assertEqual(atlas.audit().genuine_conway_boundary_spheres,0)
        self.assertEqual(atlas.audit().knot_piercings_of_cut_disks,0)


if __name__=="__main__":
    unittest.main()
