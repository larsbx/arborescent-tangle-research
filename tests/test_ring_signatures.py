"""Exact source-scoped ring obstruction; NO R± rewrite implementation."""
import unittest

from arborescent.internal_forks import (
    InternalFork, Tip, Port, compile_fork, generate, tangle_screen,
)
from arborescent.partial_moves import Move, apply, explore
from arborescent.ring_signatures import (
    ring_absence_by_unit_weight, ring_signatures, vertex_weights,
)
from arborescent.tangle_labels import (
    tangle_flip_f1, tangle_flip_f2, V4LabeledTree,
)
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


class RingSignatureTests(unittest.TestCase):
    def test_source_necessary_positive_and_negative_motifs(self):
        for sign in (1, -1):
            with self.subTest(sign=sign):
                source = compile_fork(InternalFork(
                    (Tip(2 * sign), Tip(2 * sign)), sign, 0
                ))
                motifs = ring_signatures(source)
                self.assertEqual(len(motifs), 1)
                self.assertEqual(motifs[0].center, 0)
                self.assertEqual(motifs[0].sign, sign)
                self.assertEqual(motifs[0].leaf_vertices, (1, 2))
                self.assertEqual(motifs[0].stem_bond, "f0")
                self.assertFalse(ring_absence_by_unit_weight(source))
                # These particular examples pass the implemented
                # branched tangle screen (not full canonicality).
                self.assertTrue(tangle_screen(source, "+" if sign > 0 else "-").passed)
                self.assertEqual(len(ring_signatures(tangle_flip_f1(source))), 1)
                self.assertEqual(len(ring_signatures(tangle_flip_f2(source, 0))), 1)

    def test_necessary_does_not_mean_sufficient(self):
        # Unit-weight root without both required ±2 terminal leaves
        # cannot instantiate the source's drawn modified-ring pattern.
        missing = compile_fork(InternalFork((Tip(2), Tip(3)), 1))
        wrong_sign = compile_fork(InternalFork((Tip(2), Tip(-2)), 1))
        high_center = compile_fork(InternalFork((Tip(2), Tip(2)), 3))
        for obj in (missing, wrong_sign, high_center):
            self.assertEqual(ring_signatures(obj), ())
        self.assertTrue(ring_absence_by_unit_weight(high_center))
        self.assertFalse(ring_absence_by_unit_weight(missing))

    def test_exact_default_a2_4_grammar_obstruction(self):
        # PR #9 default forks {0,3}; the old A1 grammar uses zero
        # inner forks and tips |w|>=2. Neither admits a ±1 vertex.
        for n in (1, 2):
            family = generate(n, tips=(3,), weights=(0, 3), extra_ports=1)
            self.assertTrue(family)
            for obj in family:
                self.assertTrue(ring_absence_by_unit_weight(obj))
                self.assertEqual(ring_signatures(obj), ())

    def test_individual_vertex_totals_conserved_on_complete_partial_orbit(self):
        # Includes strict F3' with an actual nonzero weighted pivot.
        source = compile_fork(
            InternalFork((InternalFork((Tip(3), Port())), Tip(2)), 3)
        )
        orbit = explore(source)
        self.assertTrue(orbit.complete)
        orig = vertex_weights(source)
        for target in orbit.predecessors:
            self.assertEqual(vertex_weights(target), orig)
            self.assertTrue(ring_absence_by_unit_weight(target))
            self.assertEqual(ring_signatures(target), ())
            self.assertTrue(orbit.trace(target).verify())

    def test_minimal_ring_ready_extension_and_sector_independence(self):
        # Add a +1 root decoration to weights (0,1), tips=(2,),
        # without pretending a ring move has been implemented.
        family = generate(1, tips=(2,), weights=(0, 1), extra_ports=0)
        self.assertEqual(len(family), 4)
        self.assertEqual(sum(bool(ring_signatures(t)) for t in family), 3)
        self.assertEqual(sum(ring_absence_by_unit_weight(t) for t in family), 1)

    def test_pattern_persists_under_index_renaming_and_cyclic_rotation(self):
        src = compile_fork(InternalFork((Tip(2), Tip(2)), 1))
        changed = WeightedPlaneTree(tuple(
            Vertex(
                tuple("e" + str(int(s[1:]) + 47)
                      if s.startswith("e") else s
                      for s in v.slots[1:] + v.slots[:1]),
                v.sectors[1:] + v.sectors[:1],
            ) for v in src.tree.vertices[::-1]
        ))
        self.assertEqual(len(ring_signatures(V4LabeledTree(changed, (("f0", 0),)))), 1)

    def test_input_guards(self):
        with self.assertRaises(ValueError):
            ring_signatures(None)
        with self.assertRaises(ValueError):
            vertex_weights(4)
        with self.assertRaises(ValueError):
            ring_absence_by_unit_weight(object())


if __name__ == "__main__":
    unittest.main()
