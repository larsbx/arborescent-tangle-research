"""Source Fig12.6-7 arrow frames and signed Conway algebra; no source-isotopy claim."""
import unittest
from dataclasses import replace
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.conway_mapping_class import pure_twist, BRANCH_COORDINATES
from arborescent.marked_conway_geometry import MarkedConwayCube
from arborescent.exterior_conway_collars import build_exterior_collared_band
from arborescent.spatial_twist_blocks import Point3
from arborescent.source_arrowed_conway_frames import (
    ArrowedConwayBoundary, CubicalArrowRepresentative,
    SignedConwayBoundaryClass, ArrowedConwaySeamObligation,
    chosen_cubical_arrow, proposed_arrowed_seam
)


class SourceArrowedConwayFrameTests(unittest.TestCase):
    def test_source_zero_slope_arrow_and_all_existing_physical_ports(self):
        for sign in (-1,1):
            for weights in ((0,0),(1,0),(-2,1,0),(1,-1,2,-2)):
                b=build_exterior_collared_band(
                    AtomicInstruction((2,),weights,0,sign))
                for cube in b.patches.strands.cubes:
                    arrow=chosen_cubical_arrow(cube)
                    self.assertEqual(arrow.status,"CUBICAL_CHOICE_SOURCE_UNVERIFIED")
                    self.assertEqual(arrow.frame.port,cube.port)
                    self.assertTrue(arrow.frame.source_zero_slope_arrow)
                    self.assertEqual(arrow.arc[0],Point3(cube.west,-2,0))
                    self.assertEqual(arrow.arc[-1],Point3(cube.east,-2,0))
                    self.assertEqual(arrow.arc[1:3],(
                        Point3(cube.west,-3,0),Point3(cube.east,-3,0)
                    ))
                    with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                        arrow.frame.require_source_parametrization()

    def test_bad_arrow_direction_marking_and_solid_cube_shortcuts_fail(self):
        cube=MarkedConwayCube(5,Port((),1))
        arrow=chosen_cubical_arrow(cube)
        with self.assertRaises(ValueError):
            ArrowedConwayBoundary(Port((),1),"P00","P11")
        with self.assertRaises(ValueError):
            ArrowedConwayBoundary(Port((),1),"P10","P10")
        with self.assertRaises(ValueError):
            ArrowedConwayBoundary(Port((),1),status="SOURCE_VERIFIED")
        with self.assertRaises(ValueError):
            replace(arrow,arc=tuple(reversed(arrow.arc)))
        with self.assertRaises(ValueError):
            replace(arrow,arc=(arrow.arc[0],arrow.arc[-1]))
        with self.assertRaises(ValueError):
            replace(arrow,arc=(arrow.arc[0],
                               Point3(cube.center_x,-2,0),
                               arrow.arc[-1]))
        with self.assertRaises(ValueError):
            replace(arrow,frame=ArrowedConwayBoundary(Port((),2)))
        with self.assertRaises(ValueError):
            replace(arrow,status="SOURCE_ISOTOPY")
        with self.assertRaises(ValueError):
            chosen_cubical_arrow(MarkedConwayCube(5,None))

    def test_signed_gl2_integer_group_inverse_and_four_branch_points(self):
        tested=0
        for a,b,c,d in product(range(-2,3),repeat=4):
            det=a*d-b*c
            if det not in (-1,1):
                continue
            for tx,ty in product((0,1),repeat=2):
                f=SignedConwayBoundaryClass(a,b,c,d,tx,ty)
                actual=f.puncture_permutation()
                self.assertEqual(set(actual),set(range(4)))
                for i,(x,y) in enumerate(BRANCH_COORDINATES):
                    image=((a*x+b*y+tx)%2,(c*x+d*y+ty)%2)
                    self.assertEqual(BRANCH_COORDINATES[actual[i]],image)
                inverse=f.inverse()
                self.assertEqual(f.compose(inverse).puncture_permutation(),
                                 (0,1,2,3))
                self.assertEqual(inverse.compose(f).puncture_permutation(),
                                 (0,1,2,3))
                self.assertEqual(f.compose(inverse).image_slope(1,0),(1,0))
                self.assertEqual(f.compose(inverse).image_slope(0,1),(0,1))
                self.assertEqual(f.compose(inverse).determinant,1)
                tested+=1
        self.assertGreater(tested,50)
        for matrix in (
            (1,0,0,2,0,0), # Non-unimodular
            (True,0,0,1,0,0),
            (1,0,0,1,2,0),
            (1,0,0,1,0,-1)
        ):
            with self.assertRaises(ValueError):
                SignedConwayBoundaryClass(*matrix)
        for slope in ((0,0),(2,2),(True,1)):
            with self.assertRaises(ValueError):
                SignedConwayBoundaryClass(0,1,1,0).image_slope(*slope)

    def test_figure_127_adjacent_emanating_arrows_reject_wrong_maps(self):
        incoming=ArrowedConwayBoundary(Port((),0))
        outgoing=ArrowedConwayBoundary(Port((),1))
        swap=SignedConwayBoundaryClass(0,1,1,0)
        positive=proposed_arrowed_seam(incoming,outgoing,swap)
        self.assertEqual(positive.mapping.determinant,-1)
        self.assertEqual(positive.mapping.image_puncture("P00"),"P00")
        self.assertEqual(positive.mapping.image_puncture("P10"),"P01")
        self.assertEqual(positive.source_alignment_record()[2],(0,3,2,1))
        with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
            positive.require_ambient_isotopy()
        with self.assertRaises(ValueError):
            proposed_arrowed_seam(incoming,outgoing,
                                 SignedConwayBoundaryClass(1,0,0,1))
        with self.assertRaises(ValueError):
            proposed_arrowed_seam(incoming,outgoing,
                                 SignedConwayBoundaryClass(1,0,0,-1))
        with self.assertRaises(ValueError):
            proposed_arrowed_seam(incoming,outgoing,
                                 SignedConwayBoundaryClass(0,1,1,0,1,0))
        with self.assertRaises(ValueError):
            proposed_arrowed_seam(incoming,incoming,swap)
        with self.assertRaises(ValueError):
            replace(positive,status="GEOMETRY_VERIFIED")

    def test_same_arrows_and_endpoint_matching_hide_infinite_pure_twists(self):
        incoming=ArrowedConwayBoundary(Port((),0))
        outgoing=ArrowedConwayBoundary(Port((),1))
        witnesses=set()
        for n in range(-8,9):
            # Slope matrices [[0,1],[1,2n]] have det=-1.
            # All reduce to the same F2 puncture action and
            # satisfy the SAME source Fig12.7 adjacent-arrow check.
            m=SignedConwayBoundaryClass(0,1,1,2*n)
            proof=proposed_arrowed_seam(incoming,outgoing,m)
            self.assertEqual(proof.source_alignment_record()[2],(0,3,2,1))
            self.assertEqual(m.image_slope(1,0),(0,1))
            witness=m.image_slope(0,1)
            self.assertEqual(witness,(1,2*n) if n>=0 else (1,2*n))
            witnesses.add(witness)
            self.assertTrue(proof.status.endswith("SOURCE_UNVERIFIED"))
        self.assertEqual(len(witnesses),17)
        self.assertEqual(pure_twist(1).image_slope(1,0),(1,0))
        # A source arrow along slope-zero does NOT kill all
        # pure mapping classes; a shear stabilizing its slope exists.
        self.assertEqual(pure_twist(6).image_slope(1,0),(1,0))
        self.assertNotEqual(pure_twist(6).image_slope(0,1),
                            pure_twist(0).image_slope(0,1))


if __name__=="__main__":
    unittest.main()
