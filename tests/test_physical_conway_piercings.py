"""A5.2e2b4c2c: exact K-versus-physical Conway triangles, not source cut."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.exterior_conway_collars import build_exterior_collared_band
from arborescent.connected_roof_ambient_slab import build_connected_roof_ambient_slab
from arborescent.physical_conway_piercings import (
    segment_triangle_hit, audit_existing_cube_piercings,
    audit_existing_band_conway_piercings,
)
from arborescent.spatial_twist_blocks import Point3


class PhysicalConwayPiercingsTests(unittest.TestCase):
    def test_exact_segment_triangle_contacts_and_nontransverse_cases(self):
        p=lambda *a:Point3(*a)
        face=(p(0,-2,-2),p(0,2,-2),p(0,0,2))
        proper=segment_triangle_hit(p(-1,0,0),p(1,0,0),face)
        self.assertEqual(proper.kind,"TRANSVERSE")
        self.assertEqual(proper.point,p(0,0,0))
        reverse=segment_triangle_hit(p(1,0,0),p(-1,0,0),face)
        self.assertEqual(reverse,proper)
        self.assertEqual(segment_triangle_hit(
            p(1,0,0),p(2,0,0),face).kind,"DISJOINT")
        self.assertEqual(segment_triangle_hit(
            p(-1,0,0),p(-2,0,0),face).kind,"DISJOINT")
        self.assertEqual(segment_triangle_hit(
            p(0,0,-1),p(0,0,1),face).kind,"COPLANAR_SEGMENT")
        self.assertEqual(segment_triangle_hit(
            p(0,-2,-2),p(0,-4,-4),face).kind,"COPLANAR_POINT")
        self.assertEqual(segment_triangle_hit(
            p(Fraction(-1,2),0,0),p(Fraction(1,2),0,0),face).point,p(0,0,0))
        gap=Fraction(1,10**10)
        self.assertEqual(segment_triangle_hit(
            p(gap,0,0),p(2,0,0),face).kind,"DISJOINT")
        with self.assertRaises(ValueError):
            segment_triangle_hit(p(0,0,0),p(0,0,0),face)
        with self.assertRaises(ValueError):
            segment_triangle_hit(p(-1,0,0),p(1,0,0),(
                p(0,0,0),p(0,1,0),p(0,2,0)))

    def test_all_source_candidate_cubes_have_four_derived_geometric_hits(self):
        for sign in (-1,1):
            for weights in ((0,0),(1,-1),(-2,0,1),(0,1,-1,2)):
                with self.subTest(sign=sign,weights=weights):
                    band=build_exterior_collared_band(
                        AtomicInstruction((1,),weights,0,sign))
                    audit=audit_existing_band_conway_piercings(band)
                    n=len(weights)
                    self.assertEqual(audit.existing_spheres_checked,n)
                    self.assertEqual(audit.source_cut_generated_spheres,0)
                    self.assertEqual(audit.status,
                                     "EXISTING_CUBES_NOT_SOURCE_CUT_DERIVED")
                    self.assertEqual(len(audit.piercings),4*n)
                    for i in range(n):
                        hits=[r for r in audit.piercings if r.port==Port((1,),i)]
                        self.assertEqual({r.puncture for r in hits},
                                         {"NW","NE","SE","SW"})
                        self.assertEqual(len({r.point for r in hits}),4)
                        self.assertEqual(len({r.knot_arc for r in hits}),4)
                    with self.assertRaisesRegex(NotImplementedError,
                                                 "UNVERIFIED_GEOMETRY"):
                        audit.require_source_conway_boundaries()

    def test_closed_roof_ambient_still_does_not_generate_source_spheres(self):
        atom=AtomicInstruction((),(0,1,-1),0,1)
        slab=build_connected_roof_ambient_slab(atom)
        actual=audit_existing_cube_piercings(
            slab.atlas.candidate.patches.strands.cubes,
            slab.atlas.candidate.patches.strands.embedded_arcs())
        self.assertEqual(len(actual.piercings),12)
        self.assertEqual(actual.source_cut_generated_spheres,0)
        self.assertEqual(slab.face_barriers().actual_knot_piercings,0)

    def test_mutations_and_fabricated_source_3d_claims_are_rejected(self):
        band=build_exterior_collared_band(AtomicInstruction((),(0,1),0,1))
        cubes=band.patches.strands.cubes
        arcs=band.patches.strands.embedded_arcs()
        positive=audit_existing_cube_piercings(cubes,arcs)
        with self.assertRaises(ValueError):
            replace(positive,status="SOURCE_FIGURE_CERTIFIED")
        with self.assertRaises(ValueError):
            replace(positive,source_cut_generated_spheres=2)
        with self.assertRaises(ValueError):
            replace(positive,existing_spheres_checked=0)
        with self.assertRaises(ValueError):
            audit_existing_cube_piercings(cubes,arcs[:-1])
        fake=list(arcs)
        first=list(arcs[0])
        first[1]=Point3(cubes[-1].center_x,0,0)
        fake[0]=tuple(first)
        with self.assertRaises(ValueError):
            audit_existing_cube_piercings(cubes,tuple(fake))
        fake=list(arcs)
        fake[0]=(arcs[0][0],Point3(cubes[-1].center_x,0,0),
                 arcs[0][-1])
        with self.assertRaises(ValueError):
            audit_existing_cube_piercings(cubes,tuple(fake))
        with self.assertRaises(ValueError):
            audit_existing_cube_piercings((),arcs)
        with self.assertRaises(ValueError):
            audit_existing_cube_piercings(cubes,())
        with self.assertRaises(ValueError):
            audit_existing_band_conway_piercings("claim")


if __name__=="__main__":
    unittest.main()
