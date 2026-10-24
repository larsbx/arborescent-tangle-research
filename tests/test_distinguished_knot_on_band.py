"""A5.2e2b4b: distinguish proper source K arcs from all band boundary loops."""
import unittest
from dataclasses import replace

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.exterior_conway_collars import build_exterior_collared_band
from arborescent.spatial_twist_blocks import Point3
from arborescent.distinguished_knot_on_band import (
    verify_distinguished_knot, audit_source_marking_on_candidate
)


class DistinguishedKnotTests(unittest.TestCase):
    def test_actual_source_marked_k_on_candidate_boundary(self):
        for sign in (-1,1):
            for weights in ((0,0),(1,0),(1,-1),(-1,0,2),(2,-1,0)):
                with self.subTest(sign=sign,weights=weights):
                    candidate=build_exterior_collared_band(
                        AtomicInstruction((1,),weights,0,sign))
                    outcome=audit_source_marking_on_candidate(candidate)
                    self.assertEqual(outcome.status,"MARKED_KNOT_SOURCE_UNVERIFIED")
                    self.assertEqual(outcome.distinguished_knot_arcs,2*len(weights))
                    self.assertEqual(outcome.marked_conway_endpoints,4*len(weights))
                    self.assertEqual(outcome.band_boundary_loops,
                                     candidate.topology().boundary_components)
                    self.assertEqual(outcome.band_boundary_loops,
                                     1 if sum(weights)%2 else 2)
                    self.assertGreater(outcome.complementary_band_edges,0)
                    self.assertGreater(outcome.knot_edges,0)
                    self.assertEqual(outcome.band_boundary_edges,
                                     outcome.knot_edges+outcome.complementary_band_edges)
                    # Physically complete annulus/Mobius band ≠ closed knot!
                    self.assertNotEqual(outcome.band_boundary_loops,
                                        outcome.distinguished_knot_arcs)

    def test_chosen_cube_knot_is_proper_not_whole_band_boundary(self):
        obj=build_exterior_collared_band(AtomicInstruction((),(1,-1,0),0,1))
        cubes=obj.patches.strands.cubes
        arcs=obj.patches.strands.embedded_arcs()
        positive=verify_distinguished_knot(obj.triangles(),cubes,arcs)
        self.assertEqual(positive.distinguished_knot_arcs,6)
        self.assertEqual(positive.marked_conway_endpoints,12)
        self.assertEqual(positive.band_boundary_loops,2)
        self.assertNotEqual(positive.knot_edges,positive.band_boundary_edges)
        with self.assertRaises(ValueError):
            replace(positive,status="SOURCE_GEOMETRY_CERTIFIED")
        with self.assertRaises(ValueError):
            replace(positive,complementary_band_edges=0)

    def test_bad_terminal_transversality_and_nonboundary_knot_edges_rejected(self):
        obj=build_exterior_collared_band(AtomicInstruction((),(0,1),0,1))
        cubes=obj.patches.strands.cubes
        arcs=obj.patches.strands.embedded_arcs()
        triangles=obj.triangles()
        with self.assertRaises(ValueError):
            verify_distinguished_knot(triangles,cubes,arcs[:-1])
        with self.assertRaises(ValueError):
            verify_distinguished_knot(triangles,cubes,arcs[::-1]+(arcs[0],))
        # A direct chord through a removed Conway ball is not K.
        fake=list(arcs)
        fake[0]=(arcs[0][0],arcs[0][-1])
        with self.assertRaises(ValueError):
            verify_distinguished_knot(triangles,cubes,tuple(fake))
        # Edit the normal direction of an external puncture approach
        # to zero; the source K must meet the Conway face transversely.
        fake=list(arcs)
        strand=list(arcs[0])
        strand[1]=Point3(strand[0].x,strand[0].y+1,strand[0].z)
        fake[0]=tuple(strand)
        with self.assertRaises(ValueError):
            verify_distinguished_knot(triangles,cubes,tuple(fake))
        # A mesh face added over a legal K edge changes its status
        # from exposed boundary to interior or nonmanifold.
        extra=(triangles[0][0],triangles[0][1],triangles[0][2])
        with self.assertRaises(ValueError):
            verify_distinguished_knot(triangles+(extra,),cubes,arcs)
        with self.assertRaises(ValueError):
            verify_distinguished_knot(triangles,(),arcs)
        with self.assertRaises(ValueError):
            audit_source_marking_on_candidate("fabricated isotopy")

    def test_same_marked_k_arcs_survive_collar_refactoring(self):
        """The meaning of K is independent of how a band collar is thickened.

        Previous A5.2e2b2 roof surfaces and A5.2e2b3b exterior surfaces
        have different embedded 2D faces, yet retain the same K arcs.
        A source proof may replace an auxiliary band WITHOUT silently
        changing the boundary-marked knot/tangle.
        """
        from arborescent.conway_roof_collars import build_spatial_collared_band
        atom=AtomicInstruction((),(1,0,-1),0,1)
        old=build_spatial_collared_band(atom)
        newer=build_exterior_collared_band(atom)
        arcs=newer.patches.strands.embedded_arcs()
        self.assertEqual(old.patches.strands.embedded_arcs(),arcs)
        old_mark=verify_distinguished_knot(
            old.triangles(),old.patches.strands.cubes,arcs)
        new_mark=audit_source_marking_on_candidate(newer)
        self.assertEqual(old_mark.distinguished_knot_arcs,
                         new_mark.distinguished_knot_arcs)
        self.assertEqual(old_mark.knot_edges,new_mark.knot_edges)
        self.assertEqual(old_mark.marked_conway_endpoints,
                         new_mark.marked_conway_endpoints)


if __name__=="__main__":
    unittest.main()
