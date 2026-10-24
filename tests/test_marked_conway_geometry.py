"""A5.2d: source-marked PL spheres, complementary braid corridors, failure gates."""
import unittest
from dataclasses import replace
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction, Port
from arborescent.spatial_twist_blocks import Point3
from arborescent.marked_conway_geometry import (
    MarkedConwayCube, TwoStrandCorridor, OpenAtomicConwayChain,
    build_open_unringed_atomic_chain,
)


class MarkedConwayGeometryTests(unittest.TestCase):
    def test_closed_outward_oriented_marked_sphere(self):
        for center in (-13, 0, 4, 100):
            sphere = MarkedConwayCube(center, Port((1,2), 0))
            sphere.audit()
            self.assertEqual(len(sphere.vertices()), 8)
            self.assertEqual(len(sphere.triangles()), 12)
            self.assertEqual(tuple(name for name,_ in sphere.punctures()),
                             ("NW", "NE", "SE", "SW"))
            self.assertEqual(sphere.punctures()[0][1], Point3(center-2,-2,0))
            self.assertEqual(sphere.punctures()[1][1], Point3(center+2,-2,0))
            self.assertEqual(sphere.punctures()[2][1], Point3(center+2,2,0))
            self.assertEqual(sphere.punctures()[3][1], Point3(center-2,2,0))
        with self.assertRaises(ValueError):
            MarkedConwayCube(True)
        with self.assertRaises(ValueError):
            MarkedConwayCube(0, "bad port")

    def test_sectors_preserved_with_real_disjoint_corresponding_arcs(self):
        fixtures = (
            (0,0), (1,-2), (-1,3,0), (0,1,-1,2),
            (4,-3,1,0,2), (0,)*6
        )
        for sign in (-1, 1):
            for weights in fixtures:
                atom = AtomicInstruction((1,2), weights, 0, sign)
                chain = build_open_unringed_atomic_chain(atom)
                chain.audit()
                self.assertEqual(chain.status, "OPEN_CHAIN_NOT_ATOMIC")
                self.assertEqual(len(chain.cubes), len(weights)+1)
                self.assertEqual(len(chain.corridors), len(weights))
                self.assertEqual(chain.cubes[0].port, None)
                self.assertEqual(tuple(c.port for c in chain.cubes[1:]),
                                 tuple(Port((1,2), i) for i in range(len(weights))))
                self.assertEqual(tuple(c.signed_weight for c in chain.corridors),weights)
                self.assertEqual(tuple(c.local.writhe for c in chain.corridors), weights)
                self.assertEqual(len(chain.unused_anchor_endpoints()), 2)
                self.assertEqual(len(chain.unused_terminal_endpoints()), 2)
                self.assertEqual(chain.canonical_coordinates()[0], atom.address)
                self.assertEqual(chain.canonical_coordinates()[2], weights)
                for i,corridor in enumerate(chain.corridors):
                    self.assertEqual(corridor.left, chain.cubes[i])
                    self.assertEqual(corridor.right, chain.cubes[i+1])
                    self.assertEqual(corridor.right.west-corridor.left.east,
                                     2*abs(weights[i])+1)
                    arcs = corridor.strands()
                    self.assertEqual(len(arcs),2)
                    self.assertEqual(arcs[0][0],Point3(corridor.left.east,-2,0))
                    self.assertEqual(arcs[1][0],Point3(corridor.left.east,2,0))
                    self.assertEqual({arc[-1] for arc in arcs},
                                     {Point3(corridor.right.west,-2,0),
                                      Point3(corridor.right.west,2,0)})
                with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                    chain.require_source_atomic_realization()

    def test_exhaustive_small_signed_weight_words(self):
        for length in range(2,5):
            for weights in product((-2,-1,0,1,2), repeat=length):
                atom = AtomicInstruction((),tuple(weights),0,1)
                chain = build_open_unringed_atomic_chain(atom)
                # Canonical exact coordinates are deterministic.
                self.assertEqual(chain.canonical_coordinates(),
                                 build_open_unringed_atomic_chain(atom).canonical_coordinates())
                # Crossing number is the sum of sector half-twist counts.
                self.assertEqual(sum(c.local.crossing_count for c in chain.corridors),
                                 sum(abs(w) for w in weights))
                self.assertEqual(len({cube.center_x for cube in chain.cubes}),
                                 length+1)

    def test_adversarial_distortion_and_false_promotion_rejected(self):
        atom = AtomicInstruction((),(2,-1,0),0,-1)
        chain = build_open_unringed_atomic_chain(atom)
        with self.assertRaises(ValueError):
            replace(chain,status="VERIFIED_ATOMIC_TANGLE")
        cubes=list(chain.cubes)
        cubes[1]=replace(cubes[1], center_x=cubes[1].center_x+1)
        with self.assertRaises(ValueError):
            replace(chain,cubes=tuple(cubes))
        cubes=list(chain.cubes)
        cubes[1]=replace(cubes[1], port=Port((),2))
        with self.assertRaises(ValueError):
            replace(chain,cubes=tuple(cubes))
        moved=list(chain.corridors)
        with self.assertRaises(ValueError):
            replace(moved[0], signed_weight=1)
        with self.assertRaises(ValueError):
            replace(moved[0], left=chain.cubes[1])
        moved[0]=replace(moved[0], sector_index=10)
        with self.assertRaises(ValueError):
            replace(chain,corridors=tuple(moved))
        with self.assertRaises(ValueError):
            replace(chain,atom=AtomicInstruction((),(2,0,0),0,-1))
        with self.assertRaisesRegex(NotImplementedError, "UNVERIFIED_GEOMETRY"):
            build_open_unringed_atomic_chain(AtomicInstruction((),(0,0),2,1))
        with self.assertRaisesRegex(NotImplementedError, "UNVERIFIED_GEOMETRY"):
            build_open_unringed_atomic_chain(AtomicInstruction((),(1,),0,1))
        with self.assertRaises(ValueError):
            build_open_unringed_atomic_chain("fabricated atom")
        with self.assertRaises(ValueError):
            TwoStrandCorridor(chain.cubes[0],chain.cubes[2],0,2,chain.corridors[0].local)


if __name__ == "__main__":
    unittest.main()
