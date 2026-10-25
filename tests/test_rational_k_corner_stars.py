"""A5.2e2b4c2d2c: actual K two-valent corner tetrahedral 3-ball stars."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.exterior_conway_collars import build_exterior_collared_band
from arborescent.rational_k_corner_stars import (
    _generators,_vertices,_tetrahedra,_edge_link,
    _convex_clip_segment,build_rational_k_corner_star,
)
from arborescent.spatial_twist_blocks import Point3


class RationalKCornerStarTests(unittest.TestCase):
    def test_bent_corner_and_straight_continuation_exact_geometry(self):
        p=lambda *args:Point3(*args)
        cases=(
            (p(-1,0,0),p(0,0,0),p(0,2,0),False),
            (p(0,0,0),p(1,2,3),p(2,3,5),False),
            (p(-2,0,0),p(0,0,0),p(3,0,0),True),
            (p(0,0,-2),p(0,0,0),p(0,0,3),True),
        )
        for before,corner,after,straight in cases:
            with self.subTest(before=before,after=after):
                u,v,w,first,second=_generators(
                    before,corner,after,Fraction(1,8))
                self.assertEqual((first,second),(1,2 if straight else 3))
                self.assertGreater(sum(u[j]*(
                    v[(j+1)%3]*w[(j+2)%3]-
                    v[(j+2)%3]*w[(j+1)%3]) for j in range(3)),0)
                vertices=_vertices(corner,u,v,w)
                tetrahedra=_tetrahedra(vertices)
                self.assertEqual((len(vertices),len(tetrahedra)),(7,8))
                self.assertEqual(vertices[0],corner)
                self.assertEqual(_edge_link(tetrahedra,first),(4,4))
                self.assertEqual(_edge_link(tetrahedra,second),(4,4))
                for source in (before,after):
                    clip=_convex_clip_segment(corner,source,corner,(u,v,w))
                    self.assertIsNotNone(clip)
                    self.assertEqual(clip[0],0)
                    self.assertTrue(0<clip[1]<1)
                self.assertIsNone(_convex_clip_segment(
                    p(10,10,10),p(20,10,10),corner,(u,v,w)))
        with self.assertRaises(ValueError):
            _generators(p(1,0,0),p(0,0,0),p(2,0,0),Fraction(1,8))
        with self.assertRaises(ValueError):
            _generators(p(0,0,0),p(0,0,0),p(0,1,0),Fraction(1,8))
        with self.assertRaises(ValueError):
            _generators(p(1,0,0),p(0,0,0),p(0,1,0),Fraction(0))

    def test_true_candidate_k_vertices_have_tetrahedral_incident_edges(self):
        for weights,chosen in (
            ((0,0),((0,1),(0,2),(1,3))),
            ((1,-1),((0,1),(1,2),(2,1))),
            # Sector a1=0 has two-vertex straight connector arcs
            # (2,3) with NO internal 2-valent vertices. Test true
            # corners of the signed a2=1 sector (arcs 4,5) instead.
            ((-1,0,1),((0,2),(4,1),(5,1))),
        ):
            atom=AtomicInstruction((),weights,0,1)
            for arc,corner in chosen:
                with self.subTest(weights=weights,arc=arc,corner=corner):
                    star=build_rational_k_corner_star(atom,arc,corner)
                    stats=star.audit()
                    self.assertEqual(
                        star.status,"LOCAL_CORNER_STAR_NOT_WHOLE_K_AMBIENT")
                    self.assertEqual((stats.vertices,stats.tetrahedra),(7,8))
                    self.assertEqual((stats.boundary_euler,
                                      stats.boundary_components),(2,1))
                    self.assertGreater(stats.volume,0)
                    self.assertGreater(star.radius,0)
                    self.assertLessEqual(star.radius,Fraction(1,4))
                    original=star.candidate.patches.strands.embedded_arcs()
                    self.assertEqual(star.vertices[0],original[arc][corner])
                    with self.assertRaisesRegex(
                        NotImplementedError,"UNVERIFIED_GEOMETRY"
                    ):
                        star.require_global_source_ambient()

    def test_exact_adversarial_geometry_source_restrictions(self):
        atom=AtomicInstruction((),(0,0),0,1)
        star=build_rational_k_corner_star(atom,0,1)
        with self.assertRaises(ValueError):
            replace(star,status="FIGURE_1210_SOURCE_CERTIFIED")
        with self.assertRaises(ValueError):
            replace(star,arc=True)
        with self.assertRaises(ValueError):
            replace(star,corner=0)
        with self.assertRaises(ValueError):
            replace(star,radius=Fraction(1,2))
        with self.assertRaises(ValueError):
            replace(star,vertices=star.vertices[::-1])
        with self.assertRaises(ValueError):
            replace(star,tetrahedra=star.tetrahedra[:-1])
        with self.assertRaises(ValueError):
            replace(star,tetrahedra=(star.tetrahedra[0],)*8)
        with self.assertRaises(ValueError):
            _edge_link(star.tetrahedra,0)
        with self.assertRaises(ValueError):
            _edge_link(((0,1,3,5),),1)
        with self.assertRaises(ValueError):
            build_rational_k_corner_star(atom,0,0)
        # An untwisted two-vertex corridor carries one straight K
        # segment, not a middle K corner. Reject fake intermediate
        # vertices instead of fabricating triangulation geometry.
        with self.assertRaises(ValueError):
            build_rational_k_corner_star(atom,2,1)
        with self.assertRaises(ValueError):
            build_rational_k_corner_star(atom,0,1,Fraction(0))
        with self.assertRaises(ValueError):
            build_rational_k_corner_star("not an atomic instruction",0,1)
        with self.assertRaisesRegex(
            NotImplementedError,"UNVERIFIED_GEOMETRY"
        ):
            build_rational_k_corner_star(
                AtomicInstruction((),(0,0),1,1),0,1)
        with self.assertRaisesRegex(
            NotImplementedError,"UNVERIFIED_GEOMETRY"
        ):
            build_rational_k_corner_star(
                AtomicInstruction((),(0,),0,1),0,1)


if __name__=="__main__":
    unittest.main()
