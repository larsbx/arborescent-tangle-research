"""A5.2e2b1: finite PL surface quotient, no claim of source geometric gluing."""
import unittest
from dataclasses import replace
from itertools import product

from arborescent.atomic_plumbing_contract import AtomicInstruction
from arborescent.pl_ribbon_patches import build_zero_ring_ribbon_patches
from arborescent.ribbon_surface_quotient import (
    RibbonCap, CapAttachment, AbstractRibbonSurface,
)


def patches(*weights):
    obj=build_zero_ring_ribbon_patches(AtomicInstruction((),weights,0,1))
    return obj.ribbons


class RibbonQuotientTests(unittest.TestCase):
    def test_disconnected_source_candidate_patches_are_not_global_band(self):
        for weights in ((0,0),(1,-1),(2,-2,0),(0,1,2,-1)):
            r=patches(*weights)
            surface=AbstractRibbonSurface(r,())
            topo=surface.analyze()
            self.assertEqual((topo.components,topo.euler,topo.boundary_components),
                             (len(r),len(r),len(r)))
            self.assertTrue(topo.orientable)
            self.assertEqual(topo.face_count,sum(len(p.faces) for p in r))
            self.assertEqual(topo.vertex_classes-topo.edge_classes
                             +topo.face_count,topo.euler)
            with self.assertRaisesRegex(NotImplementedError,"UNVERIFIED_GEOMETRY"):
                surface.require_source_band_certificate()

    def test_connect_two_disks_along_one_abstract_cap(self):
        for weights in ((0,0),(1,-1),(2,0),(0,2)):
            r=patches(*weights)
            for flip in (False,True):
                attach=CapAttachment(RibbonCap(0,"end"),RibbonCap(1,"start"),flip)
                x=AbstractRibbonSurface(r,(attach,))
                t=x.analyze()
                self.assertEqual((t.components,t.euler,t.boundary_components),(1,1,1))
                self.assertTrue(t.orientable)

    def test_same_endpoint_count_allows_annulus_or_mobius(self):
        # Both surfaces start from *identical embedded disk mesh*,
        # with identical 4 corner/end points in their chosen frames.
        # Abstractly identifying the two cap chains with or without
        # width reversal is extra data, not the already stored V4.
        r=(patches(0,0)[0],)
        for reverse,expected_boundaries,orientable in (
                (False,2,True),(True,1,False)):
            x=AbstractRibbonSurface(
                r,(CapAttachment(RibbonCap(0,"start"),
                                  RibbonCap(0,"end"),reverse),))
            t=x.analyze()
            self.assertEqual((t.components,t.euler),(1,0))
            self.assertEqual(t.boundary_components,expected_boundaries)
            self.assertEqual(t.orientable,orientable)
            self.assertEqual(t.vertex_classes-t.edge_classes+t.face_count,t.euler)

    def test_unsupported_and_nonmanifold_gluings_fail_closed(self):
        r=patches(1,-1,0)
        c0=RibbonCap(0,"start")
        c1=RibbonCap(1,"start")
        c2=RibbonCap(2,"start")
        with self.assertRaises(ValueError):
            RibbonCap(-1,"start")
        with self.assertRaises(ValueError):
            RibbonCap(0,"middle")
        with self.assertRaises(ValueError):
            RibbonCap(True,"start")
        with self.assertRaises(ValueError):
            CapAttachment(c0,c0)
        with self.assertRaises(ValueError):
            CapAttachment(c0,c1,1)
        with self.assertRaises(ValueError):
            CapAttachment(c0,c1,False,"TOPOLOGY_VERIFIED")
        with self.assertRaises(ValueError):
            AbstractRibbonSurface(r,(CapAttachment(c0,c1),
                                     CapAttachment(c0,c2)))
        with self.assertRaises(ValueError):
            AbstractRibbonSurface(r,(CapAttachment(c0,RibbonCap(5,"end")),))
        with self.assertRaises(ValueError):
            AbstractRibbonSurface(r,(), "SOURCE_ISOTOPY")
        with self.assertRaises(ValueError):
            AbstractRibbonSurface((),())
        with self.assertRaises(ValueError):
            AbstractRibbonSurface((r[0],"fabricated mesh"),())

    def test_bounded_gluing_sweep_requires_independent_topology_check(self):
        for n in (2,3):
            for w in product((-1,0,1),repeat=n):
                r=patches(*w)
                for reverse in (False,True):
                    join=tuple(CapAttachment(
                        RibbonCap(i,"end"),RibbonCap(i+1,"start"),
                        reverse if i%2 else False)
                        for i in range(n-1))
                    t=AbstractRibbonSurface(r,join).analyze()
                    self.assertEqual(t.components,1)
                    self.assertEqual(t.euler,1)
                    self.assertEqual(t.boundary_components,1)
                    self.assertTrue(t.orientable)


if __name__=="__main__":
    unittest.main()
