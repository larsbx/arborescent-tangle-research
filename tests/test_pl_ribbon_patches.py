"""A5.2e2a: PL surface meshes and independent Conway solid-ball clipping."""
import unittest
from dataclasses import replace
from fractions import Fraction
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.spatial_twist_blocks import Point3
from arborescent.marked_conway_geometry import MarkedConwayCube
from arborescent.pl_ribbon_patches import (
    TriangulatedRibbon, ZeroRingRibbonPatchSystem, _clip_cube, _ribbon,
    build_zero_ring_ribbon_patches,
)


class ExactRibbonMeshTests(unittest.TestCase):
    def test_ribbon_surface_euler_boundary_and_port_attachments(self):
        vectors=((0,0),(1,-1),(-2,3),(0,1,-2,2),
                 (1,0,-1,2,0),(-1,)*7)
        for sign in (-1,1):
            for weights in vectors:
                atom=AtomicInstruction((1,2),weights,0,sign)
                assembly=build_zero_ring_ribbon_patches(atom)
                self.assertEqual(assembly.status,"RIBBON_PATCHES_SOURCE_UNVERIFIED")
                self.assertEqual(len(assembly.ribbons),len(weights))
                self.assertEqual(sum(len(r.faces) for r in assembly.ribbons),
                                 sum(4*(len(r.first)-1) for r in assembly.ribbons))
                self.assertEqual(tuple(r.kind for r in assembly.ribbons),
                                 ("EXTERIOR_RETURN",)+
                                 ("SECTOR_CORRIDOR",)*(len(weights)-1))
                physical=assembly.strands.embedded_arcs()
                for i,r in enumerate(assembly.ribbons):
                    r.audit()
                    self.assertEqual((r.first,r.second),physical[2*i:2*i+2])
                    self.assertEqual(len(r.boundary_paths()),4)
                    self.assertEqual(len(r.vertices),3*len(r.first))
                    self.assertEqual(len(r.faces),4*(len(r.first)-1))
                    self.assertTrue(all(
                        isinstance(p,Point3) for p in r.vertices))
                self.assertEqual(assembly.canonical_meshes(),
                                 build_zero_ring_ribbon_patches(atom).canonical_meshes())
                with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                    assembly.require_source_atomic_realization()

    def test_closed_solid_conway_cube_triangle_clipping(self):
        cube=MarkedConwayCube(10,Port((),0))
        p=lambda *x: Point3(*x)
        self.assertEqual(_clip_cube((p(0,0,0),p(0,1,0),p(0,0,1)),cube),())
        on_face=_clip_cube((p(12,-2,0),p(12,2,0),p(15,0,0)),cube)
        self.assertEqual(set(on_face),{p(12,-2,0),p(12,2,0)})
        penetrate=_clip_cube((p(6,-2,0),p(10,0,0),p(14,2,0)),cube)
        self.assertNotEqual(penetrate,())
        self.assertTrue(all(8<=t.x<=12 and -3<=t.y<=3 and -3<=t.z<=3
                            for t in penetrate))
        rational=_clip_cube((p(Fraction(23,2),0,-5),
                             p(Fraction(23,2),0,5),p(14,0,0)),cube)
        self.assertNotEqual(rational,())
        self.assertTrue(all(t.x<=12 and abs(t.z)<=3 for t in rational))

    def test_exhaustive_short_signed_sectors(self):
        for n in (2,3):
            for weights in product((-2,-1,0,1,2),repeat=n):
                obj=build_zero_ring_ribbon_patches(AtomicInstruction((),weights,0,1))
                self.assertEqual(len(obj.ribbons),n)
                self.assertEqual(len(obj.strands.boundary_pairings()),2*n)
                self.assertTrue(all(len(r.faces)>0 for r in obj.ribbons))

    def test_surface_tampering_collapsed_and_reversed_triangulations(self):
        obj=build_zero_ring_ribbon_patches(AtomicInstruction((),(1,-2,0),0,-1))
        outer=obj.ribbons[0]
        with self.assertRaises(ValueError):
            replace(obj,status="VERIFIED_SOURCE_ATOM")
        with self.assertRaises(ValueError):
            replace(obj,ribbons=obj.ribbons[:-1])
        with self.assertRaises(ValueError):
            replace(obj,ribbons=obj.ribbons[::-1])
        with self.assertRaises(ValueError):
            replace(outer,faces=outer.faces[:-1])
        faces=list(outer.faces)
        faces[0]=(faces[0][1],faces[0][0],faces[0][2])
        with self.assertRaises(ValueError):
            replace(outer,faces=tuple(faces))
        with self.assertRaises(ValueError):
            replace(outer,faces=(outer.faces[0],)*len(outer.faces))
        with self.assertRaises(ValueError):
            replace(outer,ridge=outer.first)
        ridge=list(outer.ridge)
        ridge[2]=Point3(ridge[2].x,ridge[2].y,0)
        # Noncanonical ridge either collapses triangles or violates
        # the independently proved global arch-height separation.
        with self.assertRaises(ValueError):
            changed=replace(outer,ridge=tuple(ridge))
            replace(obj,ribbons=(changed,)+obj.ribbons[1:])
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_zero_ring_ribbon_patches(AtomicInstruction((),(1,0),2,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            build_zero_ring_ribbon_patches(AtomicInstruction((),(1,),0,1))
        with self.assertRaises(ValueError):
            build_zero_ring_ribbon_patches("invented")

    def test_independent_geometry_not_just_matching_ends(self):
        obj=build_zero_ring_ribbon_patches(AtomicInstruction((),(1,-1),0,1))
        a,b=obj.ribbons[1].first,obj.ribbons[1].second
        ridge=tuple(Point3(x.x,x.y,-4) for x in obj.ribbons[1].ridge)
        with self.assertRaises(ValueError):
            modified=_ribbon("SECTOR_CORRIDOR",a,b,ridge)
            replace(obj,ribbons=(obj.ribbons[0],modified))
        # A genuine triangulated ribbon includes a 2D face *interior*;
        # its verified boundary consists exactly of the original
        # polygonal physical strands plus two sphere-face caps.
        for r in obj.ribbons:
            cap0,cap1=r.boundary_paths()[2:]
            self.assertEqual(cap0[0],r.first[0])
            self.assertEqual(cap0[-1],r.second[0])
            self.assertEqual(cap1[0],r.first[-1])
            self.assertEqual(cap1[-1],r.second[-1])


if __name__=="__main__":
    unittest.main()
