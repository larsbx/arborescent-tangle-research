"""A5.2e2b4c2d1: true full midplane cut, not just detached dual faces."""
import unittest
from dataclasses import replace
from fractions import Fraction

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.full_roof_midplane_slit import (
    SlitVertex, _disk_topology, _outward_faces,
    slit_complete_roof_midplane,
)
from arborescent.connected_roof_ambient_slab import (
    verify_connected_tetrahedral_ball,
)


class FullRoofMidplaneSlitTests(unittest.TestCase):
    def test_both_3balls_actual_cut_vertex_stars_and_inverse(self):
        fixtures=((0,0),(1,-1),(-2,0,1),(1,0,-1,2))
        for sign in (-1,1):
            for weights in fixtures:
                atom=AtomicInstruction((1,),weights,0,sign)
                with self.subTest(sign=sign,weights=weights):
                    cut=slit_complete_roof_midplane(atom)
                    stats=cut.audit()
                    n=len(weights)
                    m=2*n-1
                    self.assertEqual(cut.status,"FULL_ROOF_PL_SLIT_NOT_SOURCE_CONWAY")
                    self.assertEqual(stats.lower_tetrahedra,6*m)
                    self.assertEqual(stats.upper_tetrahedra,6*m)
                    self.assertEqual(stats.duplicated_midplane_vertices,2*(m+1))
                    self.assertEqual(stats.cut_triangles_per_side,2*m)
                    self.assertEqual(stats.each_ball_euler,2)
                    self.assertEqual(stats.cut_disk_euler,1)
                    self.assertEqual(stats.knot_piercings,0)
                    self.assertEqual(stats.source_cut_four_punctured_spheres,0)
                    self.assertEqual(stats.volume_lower+stats.volume_upper,
                                     cut.slab.topology().volume)
                    self.assertEqual(stats.volume_lower,stats.volume_upper)
                    for a,b in zip(cut.lower_faces,cut.upper_faces):
                        self.assertEqual({v.original for v in a},
                                         {v.original for v in b})
                        self.assertFalse(set(a)&set(b))
                    self.assertEqual(_disk_topology(cut.lower_faces),
                                     (4*n,8*n-3,4*n-2))
                    with self.assertRaisesRegex(NotImplementedError,
                                                 "UNVERIFIED_GEOMETRY"):
                        cut.require_source_conway_certificate()

    def test_exact_3ball_surface_topology_and_rational_scaling(self):
        atom=AtomicInstruction((),(0,1),0,1)
        for h in (Fraction(1,16),Fraction(1,8),Fraction(1,4)):
            cut=slit_complete_roof_midplane(atom,h)
            stats=cut.audit()
            self.assertEqual(stats.volume_lower,cut.slab.topology().volume/2)
            self.assertEqual(stats.volume_upper,cut.slab.topology().volume/2)
            self.assertEqual(len(_outward_faces(cut.lower)),8*3+4)
            self.assertEqual(len(_outward_faces(cut.upper)),8*3+4)
            self.assertEqual(cut.lower_faces,slit_complete_roof_midplane(
                atom,h).lower_faces)

    def test_forged_cut_mapped_edges_and_cells_are_rejected(self):
        cut=slit_complete_roof_midplane(
            AtomicInstruction((),(1,0),0,1))
        with self.assertRaises(ValueError):
            replace(cut,status="SOURCE_VERIFIED")
        with self.assertRaises(ValueError):
            replace(cut,lower=cut.lower[:-1])
        with self.assertRaises(ValueError):
            replace(cut,lower_faces=cut.lower_faces[:-1])
        with self.assertRaises(ValueError):
            replace(cut,upper_faces=cut.lower_faces)
        bad=list(cut.lower)
        t=list(bad[0])
        first=t[0]
        t[0]=SlitVertex(first.original,"UPPER")
        bad[0]=tuple(t)
        with self.assertRaises(ValueError):
            replace(cut,lower=tuple(bad))
        with self.assertRaises(ValueError):
            SlitVertex(True,"LOWER")
        with self.assertRaises(ValueError):
            SlitVertex(-1,"EXTERIOR")
        with self.assertRaises(ValueError):
            SlitVertex(4,"SOURCE_ISOTOPY")
        with self.assertRaises(ValueError):
            _disk_topology(())
        with self.assertRaises(ValueError):
            _disk_topology((cut.lower_faces[0],cut.lower_faces[0]))
        # A full independent check must reject deleting a tetrahedron:
        # there is no special cap-name trust to bypass the 3D link audit.
        with self.assertRaises(ValueError):
            verify_connected_tetrahedral_ball(
                cut.slab.vertices,cut.slab.tetrahedra[:-1])

    def test_zero_k_and_unsupported_ring_or_valence_do_not_certify_source(self):
        cut=slit_complete_roof_midplane(
            AtomicInstruction((),(0,0,1),0,1))
        k=cut.slab.atlas.candidate.patches.strands.embedded_arcs()
        self.assertLess(max(p.z for arc in k for p in arc),
                        min(p.z for p in cut.slab.vertices))
        self.assertEqual(cut.audit().knot_piercings,0)
        self.assertEqual(cut.audit().source_cut_four_punctured_spheres,0)
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            slit_complete_roof_midplane(AtomicInstruction((),(0,0),1,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            slit_complete_roof_midplane(AtomicInstruction((),(0,),0,1))
        with self.assertRaises(ValueError):
            slit_complete_roof_midplane("not an atom")


if __name__ == "__main__":
    unittest.main()
