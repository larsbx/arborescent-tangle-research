"""Exact-Q3 triangle interior collision tests (A5.2e2b3a)."""
import unittest
from fractions import Fraction

from arborescent.spatial_twist_blocks import Point3
from arborescent.triangle_intersections import (
    TriangleContact, triangle_contact, audit_embedded_triangle_mesh,
)


class ExactTriangleIntersectionTests(unittest.TestCase):
    def setUp(self):
        p=lambda *args:Point3(*args)
        self.p=p
        self.a=(p(0,0,0),p(4,0,0),p(0,4,0))

    def check_contact(self,a,b,kind,points=None):
        forward=triangle_contact(a,b)
        backward=triangle_contact(b,a)
        self.assertEqual(forward.kind,kind)
        self.assertEqual(backward.kind,kind)
        self.assertEqual(set(forward.points),set(backward.points))
        if points is not None:
            self.assertEqual(set(forward.points),set(points))

    def test_coplanar_disjoint_point_edge_and_overlap(self):
        p=self.p
        self.check_contact(self.a,(p(5,5,0),p(6,5,0),p(5,6,0)),"DISJOINT")
        self.check_contact(self.a,(p(0,0,0),p(-2,0,0),p(0,-2,0)),
                           "POINT",(p(0,0,0),))
        self.check_contact(self.a,(p(0,0,0),p(4,0,0),p(1,-2,0)),
                           "SEGMENT",(p(0,0,0),p(4,0,0)))
        self.check_contact(self.a,(p(1,1,0),p(2,1,0),p(1,2,0)),"AREA")
        self.check_contact(self.a,(p(0,1,0),p(4,1,0),p(2,-1,0)),"AREA")
        self.check_contact(self.a,(p(4,0,0),p(4,2,0),p(6,0,0)),
                           "POINT",(p(4,0,0),))

    def test_noncoplanar_true_crossing_skew_and_touch(self):
        p=self.p
        self.check_contact(self.a,(p(0,0,0),p(4,0,0),p(0,0,2)),
                           "SEGMENT",(p(0,0,0),p(4,0,0)))
        self.check_contact(self.a,(p(1,1,-1),p(1,1,1),p(3,1,0)),
                           "SEGMENT",(p(1,1,0),p(3,1,0)))
        self.check_contact(self.a,(p(1,1,0),p(1,0,2),p(0,1,2)),
                           "POINT",(p(1,1,0),))
        self.check_contact(self.a,(p(1,1,1),p(3,1,1),p(1,3,1)),"DISJOINT")
        self.check_contact(self.a,(p(0,0,2),p(0,2,1),p(3,0,1)),"DISJOINT")
        self.check_contact(self.a,(p(0,0,0),p(-1,-1,1),p(-1,1,1)),
                           "POINT",(p(0,0,0),))

    def test_exact_rational_subdivision_not_floating_tolerance(self):
        p=self.p
        h=Fraction(1,2)
        self.check_contact(self.a,(p(h,h,-1),p(h,h,1),p(h+1,h,0)),
                           "SEGMENT",(p(h,h,0),p(h+1,h,0)))
        gap=Fraction(1,10**12)
        self.check_contact(self.a,(p(1,1,gap),p(2,1,gap),p(1,2,gap)),
                           "DISJOINT")
        self.check_contact(self.a,(p(1,1,0),p(2,1,0),p(1,2,0)),"AREA")
        for bad in (
            (p(0,0,0),p(0,0,0),p(0,1,0)),
            (p(0,0,0),p(1,1,0),p(2,2,0)),
            (p(0,0,0),p(1,0,0)),
        ):
            with self.assertRaises(ValueError):
                triangle_contact(bad,self.a)

    def test_disjoint_mesh_vs_interior_crossing(self):
        p=self.p
        t1=(p(0,0,0),p(2,0,0),p(0,2,0))
        t2=(p(2,0,0),p(2,2,0),p(0,2,0))
        legal=audit_embedded_triangle_mesh((t1,t2))
        self.assertEqual(legal.checked_pairs,1)
        self.assertEqual(legal.allowed_shared_edges,1)
        with self.assertRaisesRegex(ValueError,"triangle interiors"):
            audit_embedded_triangle_mesh((t1,t2,
                (p(1,1,-1),p(1,1,1),p(3,1,0))))
        with self.assertRaisesRegex(ValueError,"triangle interiors"):
            audit_embedded_triangle_mesh((t1,t2,
                (p(1,1,0),p(3,1,0),p(2,3,0))))
        with self.assertRaisesRegex(ValueError,"triangle interiors"):
            audit_embedded_triangle_mesh((t1,t1))
        with self.assertRaises(ValueError):
            audit_embedded_triangle_mesh(())
        with self.assertRaises(ValueError):
            TriangleContact("DISJOINT",(p(0,0,0),))
        with self.assertRaises(ValueError):
            TriangleContact("SEGMENT",(p(0,0,0),p(0,0,0)))

    def test_noncoplanar_shared_vertex_legal_but_tangent_penetration_not(self):
        p=self.p
        a=self.a
        c=(p(0,0,0),p(-2,0,1),p(0,-2,1))
        good=audit_embedded_triangle_mesh((a,c))
        self.assertEqual(good.allowed_shared_vertices,1)
        with self.assertRaisesRegex(ValueError,"triangle interiors"):
            audit_embedded_triangle_mesh((
                a,(p(0,0,0),p(2,0,-1),p(0,2,1))
            ))


if __name__ == "__main__":
    unittest.main()
