"""A5.2e2b4c2d2b exact oblique rational K edge 3D manifold stars."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.exterior_conway_collars import build_exterior_collared_band
from arborescent.k_conforming_tetrahedral_tubes import _core_edge_link
from arborescent.rational_k_segment_stars import (
    _frame, _clip_k_segment, _axes_vertices, _box_disjoint_from_cube,
    build_rational_k_segment_star, RationalKSegmentStar,
)
from arborescent.spatial_twist_blocks import Point3


class RationalKSegmentStarsTests(unittest.TestCase):
    def test_transverse_frames_for_x_y_z_and_rational_skew(self):
        p=lambda *a:Point3(*a)
        for a,b in (
            (p(0,0,0),p(4,0,0)),
            (p(0,0,0),p(0,5,0)),
            (p(0,0,0),p(0,0,3)),
            (p(1,2,3),p(3,5,4)),
            (p(Fraction(1,3),0,0),p(1,Fraction(1,4),-2)),
        ):
            d,u,v=_frame(a,b)
            self.assertEqual(sum(x*y for x,y in zip(d,u)),0)
            self.assertEqual(sum(x*y for x,y in zip(d,v)),0)
            self.assertEqual(sum(x*y for x,y in zip(u,v)),0)
            self.assertEqual(max(abs(x) for x in u),1)
            self.assertEqual(max(abs(x) for x in v),1)
            self.assertGreater(sum(d[i]*(u[(i+1)%3]*v[(i+2)%3]-
                                  u[(i+2)%3]*v[(i+1)%3])
                                   for i in range(3)),0)
            vertices=_axes_vertices(a,b,Fraction(1,8))
            self.assertEqual(len(vertices),10)
            self.assertEqual(_clip_k_segment(a,b,a,b,Fraction(1,8)),
                             (Fraction(1,4),Fraction(3,4)))
        with self.assertRaises(ValueError):
            _frame(p(0,0,0),p(0,0,0))

    def test_all_directions_on_actual_signed_braid_candidate(self):
        atom=AtomicInstruction((),(1,-1),0,1)
        paths=build_exterior_collared_band(atom).patches.strands.embedded_arcs()
        # Choose one real source-K segment in each present direction;
        # diagonal sectors include nonzero projected and z changes.
        kinds={}
        for ai,arc in enumerate(paths):
            for si,(p,q) in enumerate(zip(arc,arc[1:])):
                changes=tuple(qv-pv for pv,qv in zip(
                    (p.x,p.y,p.z),(q.x,q.y,q.z)))
                if changes[2] != 0 and changes[0] != 0:
                    kind="braid-skew"
                elif changes[0]==0:
                    kind="vertical-return"
                elif changes[1]==0 and changes[2]==0:
                    kind="horizontal-return"
                else:
                    kind="other"
                kinds.setdefault(kind,(ai,si))
        self.assertIn("braid-skew",kinds)
        self.assertIn("vertical-return",kinds)
        self.assertIn("horizontal-return",kinds)
        for name,(ai,si) in kinds.items():
            with self.subTest(kind=name,arc=ai,segment=si):
                star=build_rational_k_segment_star(atom,ai,si)
                audited=star.audit()
                self.assertEqual(star.status,
                    "ARBITRARY_EDGE_STAR_NOT_GLOBAL_K_AMBIENT")
                self.assertEqual((audited.vertices,audited.tetrahedra),(10,12))
                self.assertEqual((audited.boundary_euler,
                                  audited.boundary_components),(2,1))
                self.assertEqual(_core_edge_link(star.tetrahedra),(4,4))
                self.assertLessEqual(star.radius,Fraction(1,4))
                self.assertGreater(star.radius,0)
                with self.assertRaisesRegex(NotImplementedError,
                                             "UNVERIFIED_GEOMETRY"):
                    star.require_source_ambient()

    def test_exact_segment_box_clearance_near_k_corner(self):
        p=lambda *a:Point3(*a)
        a,b=p(0,0,0),p(4,0,0)
        r=Fraction(1,8)
        self.assertEqual(_clip_k_segment(a,b,a,b,r),
                         (Fraction(1,4),Fraction(3,4)))
        self.assertIsNone(_clip_k_segment(
            p(0,0,0),p(0,1,0),a,b,r))
        self.assertIsNone(_clip_k_segment(
            p(4,0,0),p(4,1,0),a,b,r))
        self.assertEqual(_clip_k_segment(
            p(2,0,0),p(2,1,0),a,b,r),
            (Fraction(0),Fraction(1,8)))
        self.assertIsNone(_clip_k_segment(
            p(1,2,0),p(3,2,0),a,b,r))
        self.assertEqual(_clip_k_segment(
            p(2,0,r),p(2,1,r),a,b,r),
            (Fraction(0),Fraction(1,8)))
        gap=Fraction(1,10**11)
        self.assertIsNone(_clip_k_segment(
            p(2,0,r+gap),p(2,1,r+gap),a,b,r))

    def test_nonmutation_and_fabricated_status_fail_closed(self):
        atom=AtomicInstruction((),(0,0),0,1)
        star=build_rational_k_segment_star(atom,0,2)
        self.assertEqual(star.audit().tetrahedra,12)
        with self.assertRaises(ValueError):
            replace(star,status="SOURCE_FIGURE_ISOTOPY_VERIFIED")
        with self.assertRaises(ValueError):
            replace(star,radius=Fraction(1,2))
        with self.assertRaises(ValueError):
            replace(star,vertices=star.vertices[::-1])
        with self.assertRaises(ValueError):
            replace(star,tetrahedra=star.tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(star,segment=100)
        with self.assertRaises(ValueError):
            replace(star,arc=-1)
        with self.assertRaises(ValueError):
            build_rational_k_segment_star(atom,True,0)
        with self.assertRaises(ValueError):
            build_rational_k_segment_star(atom,0,100)
        with self.assertRaises(ValueError):
            build_rational_k_segment_star(atom,0,1,Fraction(1,2))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_rational_k_segment_star(
                AtomicInstruction((),(0,0),1,1),0,1)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_rational_k_segment_star(
                AtomicInstruction((),(0,),0,1),0,1)


if __name__=="__main__":
    unittest.main()
