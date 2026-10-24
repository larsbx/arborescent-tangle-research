"""A5.2a: genuine two-braid RII witnesses and typed four puncture frames.

This suite never claims the source's Fig.12.29 atomic rings were drawn.
"""
import unittest
from itertools import product

from arborescent.conway_braid_diagrams import (
    PUNCTURES, PERMUTATIONS, FourPunctureFrame,
    PuncturePermutation, AbstractPunctureMatching,
    Braid2Diagram, Crossing, RIIEdit, Braid2IsotopyCertificate,
    frames_for_plan, proposed_seam_puncture_matching,
    reduce_two_braid, prove_two_braid_equivalence,
)
from arborescent.atomic_plumbing_contract import (
    Port, atomic_plumbing_plan,
)
from arborescent.placement_fibers import (
    placement_presentations, bounded_valence_chain_presentations,
)
from arborescent.ring_abbreviation import abbreviate_disjoint_rings


class PunctureGeometryTests(unittest.TestCase):
    def test_v4_is_an_explicit_faithful_action_on_four_punctures(self):
        self.assertEqual(PUNCTURES, ("NW", "NE", "SE", "SW"))
        self.assertEqual(len(set(PERMUTATIONS)), 4)
        for a, b in product(range(4), repeat=2):
            lhs = PuncturePermutation(a).compose(PuncturePermutation(b))
            self.assertEqual(lhs.v4, a ^ b)
            self.assertEqual(lhs.inverse().compose(lhs).v4, 0)
            self.assertEqual(
                tuple(lhs.apply(p) for p in PUNCTURES),
                tuple(PUNCTURES[k] for k in PERMUTATIONS[a ^ b]),
            )
        for a in (1, 2, 3):
            self.assertEqual(len(set(PuncturePermutation(a).image_indices)), 4)
            self.assertNotEqual(PuncturePermutation(a).image_indices,
                                PERMUTATIONS[0])

    def test_atomic_port_frames_have_exact_four_declared_punctures(self):
        for sign in (-1, 1):
            for k in (2, 3, 5):
                source = bounded_valence_chain_presentations(k, sign)[0]
                plan = atomic_plumbing_plan(
                    abbreviate_disjoint_rings(source, sign)
                )
                frames = frames_for_plan(plan)
                self.assertEqual(
                    len(frames),
                    2 * len(plan.seams) + len(plan.free),
                )
                self.assertEqual(len({f.port for f in frames}), len(frames))
                self.assertEqual({f.cyclic_punctures for f in frames},
                                 {PUNCTURES})
                self.assertEqual(plan.topology_status, "UNVERIFIED_GEOMETRY")

    def test_abstract_four_puncture_matching_is_not_certified_gluing(self):
        matching = AbstractPunctureMatching(
            FourPunctureFrame(Port((), 1)),
            FourPunctureFrame(Port((1,), 0)),
            PuncturePermutation(3),
        )
        self.assertEqual(
            matching.pairs(),
            (("NW", "SE"), ("NE", "SW"),
             ("SE", "NW"), ("SW", "NE")),
        )
        self.assertEqual(matching.status, "ABSTRACT_ONLY")
        with self.assertRaisesRegex(NotImplementedError, "UNVERIFIED_GEOMETRY"):
            matching.certify_seam_isotopy()
        with self.assertRaises(ValueError):
            AbstractPunctureMatching(matching.parent, matching.parent,
                                    PuncturePermutation(0))
        with self.assertRaises(ValueError):
            AbstractPunctureMatching(matching.parent, matching.child,
                                    matching.mapping, "TOPOLOGY_PROVED")
        with self.assertRaises(ValueError):
            FourPunctureFrame(Port((), 0), ("NW", "SE", "NE", "SW"))
        with self.assertRaises(ValueError):
            frames_for_plan(None)

    def test_concrete_ports_are_checked_against_declared_abstract_seams(self):
        original = bounded_valence_chain_presentations(3, 1)[0]
        plan = atomic_plumbing_plan(abbreviate_disjoint_rings(original, 1))
        self.assertEqual(len(plan.seams), 2)
        for seam in plan.seams:
            for v4 in range(4):
                result = proposed_seam_puncture_matching(
                    plan, seam.parent, seam.child, v4
                )
                self.assertEqual(result.parent.port, seam.parent)
                self.assertEqual(result.child.port, seam.child)
                self.assertEqual(len({right for _, right in result.pairs()}), 4)
                self.assertEqual(result.mapping.v4, v4)
                self.assertEqual(result.status, "ABSTRACT_ONLY")
        seam = plan.seams[0]
        with self.assertRaises(ValueError):
            proposed_seam_puncture_matching(
                plan, seam.child, seam.parent, 0
            )
        with self.assertRaises(ValueError):
            proposed_seam_puncture_matching(
                plan, Port((), 999), seam.child, 0
            )
        with self.assertRaises(ValueError):
            proposed_seam_puncture_matching(
                plan, seam.parent, seam.child, True
            )
        with self.assertRaises(ValueError):
            proposed_seam_puncture_matching(
                None, seam.parent, seam.child, 0
            )

    def test_invalid_v4_puncture_elements_and_names(self):
        for bad in (True, -1, 4, 1.0, None):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                PuncturePermutation(bad)
        with self.assertRaises(ValueError):
            PuncturePermutation(1).apply("TOP")
        with self.assertRaises(ValueError):
            PuncturePermutation(1).compose("not a permutation")


class TwoBraidDiagramTests(unittest.TestCase):
    def test_crossings_over_under_and_exact_endpoints(self):
        for word in (
            (), (1,), (-1,), (1, -1), (-1, 1),
            (1, 1, -1), (-1, -1, 1, 1),
        ):
            b = Braid2Diagram(word)
            self.assertEqual(len(b.crossing_events()), len(word))
            self.assertEqual(b.exponent, sum(word))
            self.assertEqual(
                b.endpoint_matching(),
                ((("NW", "SE"), ("NE", "SW")) if len(word) % 2
                 else (("NW", "SW"), ("NE", "SE"))),
            )
            paths = b.projection_paths()
            self.assertEqual(paths[0][0], (0.0, 0.0))
            self.assertEqual(paths[1][0], (2.0, 0.0))
            for i, event in enumerate(b.crossing_events()):
                self.assertEqual(event.level, i)
                self.assertEqual(event.sign, word[i])
                self.assertEqual({event.over, event.under}, {"NW", "NE"})
                self.assertEqual(event.over,
                                 event.left_at_top if word[i] == 1
                                 else event.right_at_top)
                self.assertIn((1.0, i+.5), paths[0])
                self.assertIn((1.0, i+.5), paths[1])
            y_end = float(len(word)+1)
            self.assertEqual(paths[0][-1], (
                0.0 if len(word)%2 == 0 else 2.0, y_end))
            self.assertEqual(paths[1][-1], (
                2.0 if len(word)%2 == 0 else 0.0, y_end))
            svg = b.render_svg()
            self.assertTrue(svg.startswith("<svg "))
            self.assertTrue(svg.endswith("</svg>"))
            self.assertEqual(svg.count("<desc>"), len(word))

    def test_every_inverse_cancel_is_replayable_actual_rii(self):
        for length in range(7):
            for word in product((-1, 1), repeat=length):
                b = Braid2Diagram(tuple(word))
                cert = reduce_two_braid(b)
                self.assertTrue(cert.verify())
                self.assertEqual(cert.target.exponent, b.exponent)
                self.assertEqual(cert.target.crossing_count, abs(b.exponent))
                self.assertEqual(len(cert.edits),
                                 (len(word) - abs(b.exponent))//2)
                steps = b
                for action in cert.edits:
                    self.assertEqual(action.kind, "delete")
                    nxt = action.apply(steps)
                    self.assertEqual(action.inverse().apply(nxt), steps)
                    self.assertEqual(steps.exponent, nxt.exponent)
                    steps = nxt
                self.assertEqual(steps, cert.target)

    def test_reidemeister_ii_equivalence_certificate_two_routes(self):
        examples = (
            ((1, -1), ()),
            ((-1, 1), ()),
            ((1, 1, -1, -1, 1), (1,)),
            ((-1, -1, 1, 1, -1), (-1,)),
            ((1, -1, -1, 1), (-1, 1)),
            ((1,)*7, (1, 1, -1, 1, 1, 1, 1, 1, 1)),
        )
        for source, target in examples:
            with self.subTest(source=source, target=target):
                a, b = Braid2Diagram(source), Braid2Diagram(target)
                self.assertEqual(a.exponent, b.exponent)
                cert = prove_two_braid_equivalence(a, b)
                self.assertTrue(cert.verify())
                self.assertEqual(cert.start, a)
                self.assertEqual(cert.target, b)
                self.assertEqual(a.endpoint_matching(), b.endpoint_matching())

    def test_signed_writhe_is_only_a_braid_isotopy_obstruction(self):
        for a, b in (((1,), (-1,)), ((), (1,)),
                     ((1,1), ()), ((-1,-1), (1,1))):
            with self.subTest(a=a, b=b):
                self.assertNotEqual(sum(a), sum(b))
                with self.assertRaisesRegex(ValueError, "exponents differ"):
                    prove_two_braid_equivalence(Braid2Diagram(a), Braid2Diagram(b))
        # Negative result MUST NOT be interpreted as tangle nonisotopy
        # under an arbitrary 3D homeomorphism of its sphere boundary.

    def test_invalid_diagrams_crossings_and_rii_operations_fail_closed(self):
        for bad in ((0,), (True,), (2,), [-1, 1], "++", (1, 1.0)):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                Braid2Diagram(bad)
        with self.assertRaises(ValueError):
            Crossing(0, 1, "NW", "NE", "NE", "NW")
        with self.assertRaises(ValueError):
            Crossing(True, 1, "NW", "NE", "NW", "NE")
        with self.assertRaises(ValueError):
            RIIEdit("delete", 0, True)
        with self.assertRaises(ValueError):
            RIIEdit("delete", -1, 1)
        with self.assertRaises(ValueError):
            RIIEdit("wrong", 0, 1)
        with self.assertRaises(ValueError):
            RIIEdit("delete", 0, 1).apply(Braid2Diagram((1, 1)))
        with self.assertRaises(ValueError):
            RIIEdit("insert", 2, 1).apply(Braid2Diagram(()))
        with self.assertRaises(ValueError):
            reduce_two_braid(None)
        with self.assertRaises(ValueError):
            prove_two_braid_equivalence(None, Braid2Diagram(()))
        self.assertFalse(Braid2IsotopyCertificate(
            Braid2Diagram((1, -1)), Braid2Diagram(()), (),
        ).verify())

    def test_rendered_diagram_has_no_ring_tile_or_geometry_approval(self):
        src = placement_presentations(3, 1)[0]
        p = atomic_plumbing_plan(abbreviate_disjoint_rings(src, 1))
        self.assertEqual(p.topology_status, "UNVERIFIED_GEOMETRY")
        with self.assertRaises(NotImplementedError):
            p.require_verified_isotopy()
        # A two-braid projection is an independently verified 4-ended
        # elementary tangle; it is not any Fig 12.29 ring-atomic tile.
        diagram = Braid2Diagram((1, -1))
        self.assertIn("over=NW", diagram.render_svg())
        self.assertEqual(diagram.exponent, 0)


if __name__ == "__main__":
    unittest.main()
