"""A5.2e2b3b: exact off-boundary collar coordinates and geometric contacts."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.spatial_twist_blocks import Point3
from arborescent.pl_ribbon_patches import _clip_cube, _in_segment
from arborescent.marked_conway_geometry import MarkedConwayCube
from arborescent.exterior_conway_collars import (
    ExteriorConwayCollar, exterior_collar, build_exterior_collared_band,
)
from arborescent.triangle_intersections import audit_embedded_triangle_mesh


class ExteriorBandTests(unittest.TestCase):
    def test_exact_closed_ball_avoidance_except_two_attachment_caps(self):
        for center in (-10,0,7):
            cube=MarkedConwayCube(center,Port((1,),0))
            collar=exterior_collar(cube)
            collar.audit()
            self.assertEqual(len(collar.vertices),12)
            self.assertEqual(len(collar.faces),12)
            self.assertEqual(collar.west_cap(),(
                Point3(cube.west,-2,0),Point3(cube.west,0,0),
                Point3(cube.west,2,0)))
            self.assertEqual(collar.east_cap(),(
                Point3(cube.east,-2,0),Point3(cube.east,0,0),
                Point3(cube.east,2,0)))
            for triangle in collar.triangles():
                points=_clip_cube(triangle,cube)
                for p in points:
                    self.assertTrue(any(
                        _in_segment(p,a,b)
                        for cap in (collar.west_cap(),collar.east_cap())
                        for a,b in zip(cap,cap[1:])
                    ))
            self.assertEqual(exterior_collar(cube),collar)

    def test_source_candidate_is_a_real_embedded_surface(self):
        for sign in (-1,1):
            for sectors in ((0,0),(1,0),(0,-1),(-1,1),(2,-1),
                            (0,0,0),(1,0,-1)):
                with self.subTest(sign=sign,sectors=sectors):
                    b=build_exterior_collared_band(
                        AtomicInstruction((),sectors,0,sign))
                    t=b.topology()
                    self.assertEqual(b.status,"OFF_BOUNDARY_BAND_SOURCE_UNVERIFIED")
                    self.assertEqual((t.components,t.euler),(1,0))
                    self.assertEqual(t.orientable,sum(sectors)%2==0)
                    self.assertEqual(t.boundary_components,
                                     2 if sum(sectors)%2==0 else 1)
                    self.assertEqual(t.faces,len(b.triangles()))
                    # Even for zero-weight sectors, roofs have a
                    # positive 1/2 longitudinal separation.
                    for left,right in zip(b.collars,b.collars[1:]):
                        self.assertGreaterEqual(
                            right.vertices[3].x-left.vertices[6].x,
                            Fraction(1,2))
                    audit=b.geometric_intersections()
                    self.assertEqual(audit.checked_pairs,
                                     t.faces*(t.faces-1)//2)
                    if sign==1 and sectors==(0,0):
                        self.assertEqual(b.canonical_coordinates(),
                            build_exterior_collared_band(
                                AtomicInstruction((),sectors,0,sign)
                            ).canonical_coordinates())
                    with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                        b.require_source_atomic_realization()

    def test_adversarial_geometry_and_status_fail_closed(self):
        obj=build_exterior_collared_band(
            AtomicInstruction((),(1,-1),0,1))
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,collars=obj.collars[::-1])
        with self.assertRaises(ValueError):
            replace(obj,collars=obj.collars[:-1])
        collar=obj.collars[0]
        verts=list(collar.vertices)
        verts[4]=Point3(verts[4].x,verts[4].y,0)
        with self.assertRaises(ValueError):
            replace(collar,vertices=tuple(verts))
        with self.assertRaises(ValueError):
            replace(collar,faces=collar.faces[:-1])
        with self.assertRaises(ValueError):
            replace(collar,status="SOURCE_REALIZED")
        with self.assertRaises(ValueError):
            exterior_collar("fake")
        # A triangle can have valid combinatorial mesh incidence yet
        # penetrate another: the independent geometric checker rejects.
        faces=obj.triangles()
        extra=(Point3(0,0,0),Point3(0,0,2),Point3(1,2,1))
        with self.assertRaises(ValueError):
            audit_embedded_triangle_mesh(faces+(faces[0],))
        # Source ring_count and valence-one remain unsupported.
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_exterior_collared_band(AtomicInstruction((),(0,0),2,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_exterior_collared_band(AtomicInstruction((),(0,),0,1))
        with self.assertRaises(ValueError):
            build_exterior_collared_band("not an atom")


if __name__=="__main__":
    unittest.main()
