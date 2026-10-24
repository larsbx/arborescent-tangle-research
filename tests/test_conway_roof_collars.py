"""A5.2e2b2 exact physical roof-collar surfaces in Q³ (not source isotopy)."""
import unittest
from dataclasses import replace
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.spatial_twist_blocks import Point3
from arborescent.conway_roof_collars import (
    build_spatial_collared_band, roof_collar, verify_spatial_surface,
    SpatialCollaredBandCandidate, ConwayRoofCollar,
)


class PhysicalConwayCollarTests(unittest.TestCase):
    def test_a_real_surface_not_an_abstract_quotient(self):
        for sign in (-1,1):
            for sector in ((0,0),(1,0),(-1,2),(2,-2,1),(0,1,-1,2)):
                obj=build_spatial_collared_band(
                    AtomicInstruction((),sector,0,sign))
                obj.audit()
                t=obj.topology()
                self.assertEqual(t.components,1)
                self.assertEqual(t.euler,0)
                self.assertIn(t.boundary_components,(1,2))
                self.assertEqual(t.orientable,t.boundary_components==2)
                self.assertEqual(len(obj.collars),len(sector))
                self.assertEqual(t.faces,len(obj.triangles()))
                self.assertEqual(t.euler,t.vertices-t.edges+t.faces)
                self.assertEqual(obj.status,"COLLARED_BAND_SOURCE_UNVERIFIED")
                self.assertEqual(obj.canonical_coordinates(),
                    build_spatial_collared_band(
                        AtomicInstruction((),sector,0,sign)).canonical_coordinates())
                with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                    obj.require_source_atomic_realization()

    def test_conway_face_geometry_and_glued_cap_equality(self):
        obj=build_spatial_collared_band(AtomicInstruction((),(2,-1,0),0,1))
        for i,c in enumerate(obj.collars):
            self.assertEqual(c.cube,obj.patches.strands.cubes[i])
            self.assertEqual(len(c.vertices),12)
            self.assertEqual(len(c.faces),12)
            self.assertEqual(c.west_cap(),
                             (Point3(c.cube.west,-2,0),
                              Point3(c.cube.west,0,0),
                              Point3(c.cube.west,2,0)))
            self.assertEqual(c.east_cap(),
                             (Point3(c.cube.east,-2,0),
                              Point3(c.cube.east,0,0),
                              Point3(c.cube.east,2,0)))
            self.assertEqual(c,roof_collar(c.cube))
            # Along the actual surface, no triangle penetrates the open ball.
            for triangle in c.triangles():
                self.assertTrue(
                    any(all(getattr(p,axis)==fixed for p in triangle)
                        for axis,fixed in (("x",c.cube.west),
                                          ("x",c.cube.east),("z",3)))
                )

    def test_signed_small_vectors_and_topological_parity(self):
        for n in (2,3):
            for sector in product((-1,0,1),repeat=n):
                obj=build_spatial_collared_band(
                    AtomicInstruction((1,),sector,0,1))
                t=obj.topology()
                self.assertEqual((t.components,t.euler),(1,0))
                self.assertIn(t.boundary_components,(1,2))

    def test_bad_surface_components_are_rejected_independently(self):
        obj=build_spatial_collared_band(
            AtomicInstruction((),(1,-2,0),0,-1))
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,collars=obj.collars[:-1])
        with self.assertRaises(ValueError):
            replace(obj,collars=obj.collars[::-1])
        first=obj.collars[0]
        vertices=list(first.vertices)
        vertices[4]=Point3(vertices[4].x,0,-3)
        with self.assertRaises(ValueError):
            replace(first,vertices=tuple(vertices))
        with self.assertRaises(ValueError):
            replace(first,faces=first.faces[:-1])
        with self.assertRaises(ValueError):
            replace(first,faces=(first.faces[0],)*len(first.faces))
        with self.assertRaises(ValueError):
            roof_collar("fake")
        triangles=obj.triangles()
        with self.assertRaisesRegex(ValueError,"duplicate"):
            verify_spatial_surface(triangles+(triangles[0],))
        with self.assertRaises(ValueError):
            verify_spatial_surface(triangles[:-1])
        with self.assertRaises(ValueError):
            verify_spatial_surface(())

    def test_unsupported_rings_and_valence_fail_closed(self):
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_spatial_collared_band(AtomicInstruction((),(0,0),1,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_spatial_collared_band(AtomicInstruction((),(0,),0,1))
        with self.assertRaises(ValueError):
            build_spatial_collared_band("invalid")


if __name__=="__main__":
    unittest.main()
