"""A5.2e2b4c2b: connected 3D ambient roof-slab checks, not source S3."""
import unittest
from dataclasses import replace
from fractions import Fraction
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.connected_roof_ambient_slab import (
    build_connected_roof_ambient_slab, verify_connected_tetrahedral_ball,
)
from arborescent.spatial_twist_blocks import Point3


class ConnectedRoofSlabTests(unittest.TestCase):
    def test_global_3ball_across_every_chart_and_positive_connector_gap(self):
        for sign in (-1,1):
            for weights in ((0,0),(1,0),(-1,0,2),(0,1,-2,0)):
                with self.subTest(sign=sign,weights=weights):
                    obj=build_connected_roof_ambient_slab(
                        AtomicInstruction((1,),weights,0,sign))
                    n=len(weights)
                    stats=obj.topology()
                    self.assertEqual(obj.status,
                                     "CONNECTED_ROOF_SLAB_NOT_SOURCE_AMBIENT")
                    self.assertEqual(stats.vertices,12*n)
                    self.assertEqual(stats.tetrahedra,12*(2*n-1))
                    self.assertEqual(stats.boundary_components,1)
                    self.assertEqual(stats.boundary_euler,2)
                    self.assertEqual(stats.boundary_triangles,
                                     12*(2*n-1)+8)
                    self.assertEqual(stats.volume,(obj.x_slices[-1]-
                         obj.x_slices[0])*4*obj.atlas.thickness)
                    barriers=obj.face_barriers()
                    self.assertEqual(len(barriers.paired_faces),2*n)
                    self.assertEqual(barriers.actual_knot_piercings,0)
                    # The connector 3D material leaves the ambient
                    # dual cell graph CONNECTED even after face barriers.
                    self.assertEqual(barriers.dual_components_after_detachment,1)
                    self.assertEqual(tuple(f.port for f in barriers.paired_faces
                                          if f==barriers.paired_faces[2*f.site]),
                                     tuple(s.port for s in obj.atlas.sites))
                    self.assertEqual(obj.audit(),stats)
                    with self.assertRaisesRegex(NotImplementedError,
                                                 "UNVERIFIED_GEOMETRY"):
                        obj.require_source_global_slit()

    def test_spatial_connector_faces_are_actually_shared(self):
        obj=build_connected_roof_ambient_slab(
            AtomicInstruction((),(0,1,0),0,1))
        xs=obj.x_slices
        self.assertEqual(len(xs),6)
        self.assertTrue(all(a<b for a,b in zip(xs,xs[1:])))
        self.assertEqual(xs[2]-xs[1],Fraction(1,2)+2)
        self.assertEqual(xs[4]-xs[3],Fraction(1,2))
        for i in range(1,len(xs)-1):
            self.assertEqual(len(tuple(v for v in obj.vertices if v.x==xs[i])),6)
        self.assertEqual(len(obj.face_barriers().paired_faces),6)

    def test_exact_thickness_values_and_reproducible_3d_coordinates(self):
        atom=AtomicInstruction((),(1,-1),0,1)
        for h in (Fraction(1,16),Fraction(1,8),Fraction(1,4)):
            a=build_connected_roof_ambient_slab(atom,h)
            self.assertEqual(a.audit().volume,
                             (a.x_slices[-1]-a.x_slices[0])*4*h)
            self.assertEqual(a.vertices[0].z,4-h)
            self.assertEqual(a.vertices[2].z,4)
            self.assertEqual(a.vertices[4].z,4+h)
            b=build_connected_roof_ambient_slab(atom,h)
            self.assertEqual(a.vertices,b.vertices)
            self.assertEqual(a.tetrahedra,b.tetrahedra)
            self.assertEqual(a.face_barriers(),b.face_barriers())

    def test_bounded_signed_3site_sweep_and_framing_nonclaim(self):
        for weights in product((-1,0,1),repeat=3):
            a=build_connected_roof_ambient_slab(
                AtomicInstruction((),weights,0,1))
            self.assertEqual(a.audit().tetrahedra,60)
            self.assertEqual(a.face_barriers().actual_knot_piercings,0)

    def test_geometry_mutations_and_missing_global_cut_fail_closed(self):
        obj=build_connected_roof_ambient_slab(
            AtomicInstruction((),(0,0),0,1))
        with self.assertRaises(ValueError):
            replace(obj,status="SOURCE_GLOBAL_ISOTOPY_VERIFIED")
        with self.assertRaises(ValueError):
            replace(obj,x_slices=obj.x_slices[:-1])
        with self.assertRaises(ValueError):
            replace(obj,vertices=obj.vertices[:-1])
        with self.assertRaises(ValueError):
            replace(obj,tetrahedra=obj.tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(obj,tetrahedra=(obj.tetrahedra[0],)*len(obj.tetrahedra))
        v=list(obj.vertices)
        v[0]=Point3(v[0].x,v[0].y,v[0].z+Fraction(1,8))
        with self.assertRaises(ValueError):
            replace(obj,vertices=tuple(v))
        with self.assertRaises(ValueError):
            verify_connected_tetrahedral_ball(obj.vertices,())
        with self.assertRaises(ValueError):
            verify_connected_tetrahedral_ball(
                obj.vertices,obj.tetrahedra+(obj.tetrahedra[0],))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_connected_roof_ambient_slab(
                AtomicInstruction((),(0,0),1,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_connected_roof_ambient_slab(
                AtomicInstruction((),(0,),0,1))


if __name__=="__main__":
    unittest.main()
