"""A5.2c exact 3D local geometry, deliberately NOT atomic-ring isotopy."""
import unittest
from dataclasses import replace
from fractions import Fraction
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.conway_braid_diagrams import (
    Braid2Diagram, Braid2IsotopyCertificate, RIIEdit,
    prove_two_braid_equivalence,
)
from arborescent.spatial_twist_blocks import (
    Point3, spatial_two_braid, realize_unringed_sector_blocks,
    verify_spatial_rii_certificate,
)


class SpatialTwistTests(unittest.TestCase):
    def test_real_3d_strands_and_crossings(self):
        for n in range(7):
            for word in product((-1, 1), repeat=n):
                diagram = Braid2Diagram(tuple(word))
                shape = spatial_two_braid(diagram)
                shape.audit()
                self.assertEqual(shape.crossing_count, n)
                self.assertEqual(shape.writhe, sum(word))
                self.assertEqual(shape.endpoint_matching(), diagram.endpoint_matching())
                self.assertEqual(shape.crossing_records(), tuple(
                    (i, event.over, event.under)
                    for i, event in enumerate(diagram.crossing_events())
                ))
                self.assertEqual(len(shape.nw), 2*n+2)
                self.assertEqual(len(shape.ne), 2*n+2)
                self.assertEqual(shape.canonical_coordinates()[0], tuple(word))
                self.assertEqual(shape.canonical_coordinates()[1][0], ('-2', '0', '0'))
                with self.assertRaisesRegex(NotImplementedError, "UNVERIFIED_GEOMETRY"):
                    shape.require_source_atomic_realization()

    def test_exact_coordinates_and_tampering_are_rejected(self):
        geometry = spatial_two_braid(Braid2Diagram((1, -1, 1)))
        with self.assertRaisesRegex(ValueError, "exact rational"):
            Point3(1.0, 0, 0)
        with self.assertRaisesRegex(ValueError, "exact rational"):
            Point3(True, 0, 0)
        self.assertEqual(Point3(Fraction(1,2), 2, 0).serial(), ('1/2', '2', '0'))
        for field in ('nw', 'ne'):
            points = list(getattr(geometry, field))
            for idx, distorted in (
                (1, Point3(0, 1, 0)),             # Lose 3D overpass separation
                (1, Point3(Fraction(1,2), 1, 1)), # No transverse crossing
                (2, Point3(1, 2, 0)),             # Broken end of twist cell
                (-1, Point3(0, 7, 0)),            # Corrupted free endpoint
            ):
                modified = points.copy()
                modified[idx] = distorted
                with self.assertRaises(ValueError):
                    replace(geometry, **{field: tuple(modified)})
        with self.assertRaises(ValueError):
            replace(geometry, geometry_status='VERIFIED_SOURCE_ATOM')
        with self.assertRaises(ValueError):
            replace(geometry, word=(1, 1, 1))
        with self.assertRaises(ValueError):
            spatial_two_braid('not a braid')

    def test_signed_sector_boxes_not_full_ringed_atomic_tangles(self):
        for sign in (-1, 1):
            for weights in ((0, 0), (3, -2, 0, 1), (-1, -2, 2)):
                a = AtomicInstruction((), weights, 0, sign)
                blocks = realize_unringed_sector_blocks(a)
                self.assertEqual(tuple(x.signed_half_twists for x in blocks), weights)
                self.assertEqual(tuple(x.sector_index for x in blocks), tuple(range(len(weights))))
                self.assertTrue(all(x.source_scope == 'LOCAL_SECTOR_ONLY' for x in blocks))
                self.assertEqual(tuple(x.geometry.writhe for x in blocks), weights)
                self.assertTrue(all(x.atom_address == () for x in blocks))
                with self.assertRaises(ValueError):
                    replace(blocks[0], signed_half_twists=10)
            for rings in (1, 2, 3):
                with self.assertRaisesRegex(NotImplementedError, 'UNVERIFIED_GEOMETRY'):
                    realize_unringed_sector_blocks(AtomicInstruction((), (0,0,0), rings, sign))
        with self.assertRaisesRegex(NotImplementedError, 'UNVERIFIED_GEOMETRY'):
            realize_unringed_sector_blocks(AtomicInstruction((), (3,), 0, 1))

    def test_relative_boundary_rii_each_step_is_spatially_audited(self):
        fixtures = (
            ((1,-1,1,-1,1), (1,)),
            ((-1,1,-1,1), ()),
            ((1,-1,1,-1), (-1,1)),
            ((-1,-1,1,1,-1), (-1,)),
        )
        for source, target in fixtures:
            cert = prove_two_braid_equivalence(Braid2Diagram(source), Braid2Diagram(target))
            self.assertTrue(cert.verify())
            self.assertTrue(verify_spatial_rii_certificate(cert))
            reverse = prove_two_braid_equivalence(cert.target, cert.start)
            self.assertTrue(verify_spatial_rii_certificate(reverse))
        self.assertTrue(verify_spatial_rii_certificate(
            Braid2IsotopyCertificate(Braid2Diagram((1,)), Braid2Diagram((1,)), ())
        ))
        self.assertFalse(verify_spatial_rii_certificate('invented certificate'))
        self.assertFalse(verify_spatial_rii_certificate(Braid2IsotopyCertificate(
            Braid2Diagram((1,)), Braid2Diagram((-1,)), ()
        )))
        self.assertFalse(verify_spatial_rii_certificate(Braid2IsotopyCertificate(
            Braid2Diagram((1,1)), Braid2Diagram(()), (RIIEdit('delete',0,1),)
        )))
        self.assertFalse(verify_spatial_rii_certificate(Braid2IsotopyCertificate(
            Braid2Diagram((1,-1)), Braid2Diagram(()), (RIIEdit('delete',1,1),)
        )))


if __name__ == '__main__':
    unittest.main()
