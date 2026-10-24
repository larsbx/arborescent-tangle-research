"""A5: geometry is NOT inferred from port/seam instructions."""
import unittest
from dataclasses import replace

from arborescent.atomic_plumbing_contract import (
    AtomicInstruction, AtomicPlumbingPlan, FreeBoundary, Port, Seam,
    atomic_plumbing_plan,
)
from arborescent.abbreviated_flips import (
    abbreviated_flip_f1, abbreviated_flip_f2, flip_abbreviation_square,
)
from arborescent.placement_fibers import (
    placement_presentations, bounded_valence_chain_presentations,
)
from arborescent.ring_abbreviation import Abbreviated, abbreviate_disjoint_rings
from arborescent.ring_grammar import annotate_abbreviated
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree, Fork, Leaf


def permute_storage(obj: Abbreviated) -> Abbreviated:
    """Change vertex indices, interior edge names, and cyclic starting slots."""
    tree = obj.labeled.tree
    new_vertices = []
    for vertex in reversed(tree.vertices):
        slots = tuple("e"+str(100+int(s[1:])) if s.startswith("e") else s
                      for s in vertex.slots)
        new_vertices.append(Vertex(slots[1:]+slots[:1],
                                   vertex.sectors[1:]+vertex.sectors[:1]))
    labeled = V4LabeledTree(WeightedPlaneTree(tuple(new_vertices)),
                            obj.labeled.labels)
    return Abbreviated(labeled, obj.rings[::-1], obj.sign)


class AtomicPlumbingInterfaceTests(unittest.TestCase):
    def test_corolla_ring_count_and_port_bijection(self):
        for sign in (-1, 1):
            for d, rings in ((3, 0), (3, 1), (3, 3), (5, 2)):
                for original in placement_presentations(d, rings, sign):
                    abbreviation = abbreviate_disjoint_rings(original, sign)
                    plan = atomic_plumbing_plan(abbreviation)
                    self.assertEqual(plan.topology_status, "UNVERIFIED_GEOMETRY")
                    self.assertEqual(len(plan.atoms), 1)
                    self.assertEqual(len(plan.seams), 0)
                    self.assertEqual(len(plan.free), d)
                    self.assertEqual(plan.atoms[0].ring_count, rings)
                    self.assertEqual(plan.atoms[0].sector_weights, (0,)*d)
                    self.assertEqual(plan.rooted_instruction_code(),
                                     abbreviation.rooted_code())
                    self.assertEqual(sum(x.marked_root for x in plan.free), 1)
                    with self.assertRaisesRegex(NotImplementedError,
                                                "UNVERIFIED_GEOMETRY"):
                        plan.require_verified_isotopy()

    def test_multi_base_signed_and_v4_data_survive(self):
        for sign in (-1, 1):
            for k in (2, 3, 5):
                original = bounded_valence_chain_presentations(k, sign)[0]
                labeled = V4LabeledTree(
                    original.tree, tuple((port, i % 4)
                                         for i, (port, _) in enumerate(original.labels))
                )
                abbreviation = abbreviate_disjoint_rings(labeled, sign)
                plan = atomic_plumbing_plan(abbreviation)
                self.assertEqual(len(plan.atoms), k)
                self.assertEqual(len(plan.seams), k-1)
                self.assertEqual(len(plan.free), k)
                self.assertEqual(tuple(a.ring_count for a in plan.atoms),
                                 (1,)*k)
                self.assertEqual(plan.rooted_instruction_code(),
                                 abbreviation.rooted_code())
                # Abstract free ports follow the source's rooted
                # cyclic/DFS order, which need not be lexical f0,f1,...
                # Each local V4 location is checked by the exact
                # independent rooted-code round trip above.
                self.assertEqual(
                    sorted(f.v4 for f in plan.free),
                    sorted(i % 4 for i in range(k)),
                )
                self.assertEqual(
                    [f.v4 for f in plan.free if f.marked_root], [0]
                )
                self.assertEqual(sum(a.total_weight for a in plan.atoms), 0)

    def test_full_placement_fiber_is_identical_abstract_plumbing_input(self):
        # Repositioned rings disappear only into ring numbers, not via
        # any assertion of a topological R± equivalence.
        for sign in (-1, 1):
            for d, n in ((3, 2), (4, 2), (5, 3)):
                codes = set()
                for original in placement_presentations(d, n, sign):
                    abbreviation = abbreviate_disjoint_rings(original, sign)
                    codes.add(atomic_plumbing_plan(abbreviation)
                              .rooted_instruction_code())
                self.assertEqual(len(codes), 1)

    def test_storage_changes_do_not_change_rooted_plumbing_input(self):
        for sign in (-1, 1):
            for original in bounded_valence_chain_presentations(3, sign):
                abbreviation = abbreviate_disjoint_rings(original, sign)
                changed = permute_storage(abbreviation)
                plan = atomic_plumbing_plan(abbreviation)
                other = atomic_plumbing_plan(changed)
                self.assertEqual(plan, other)
                self.assertEqual(plan.rooted_instruction_code(),
                                 changed.rooted_code())

    def test_free_v4_and_angular_sector_variants_are_distinguished(self):
        source = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e0", "f1"), (5, 0, 0)),
            Vertex(("e0", "f2", "f3"), (0, 0, 0)),
        )), (("f0", 0), ("f1", 1), ("f2", 2), ("f3", 3)))
        baseline = Abbreviated(source, (2, 0), 1)
        standard = atomic_plumbing_plan(baseline)
        self.assertEqual(standard.rooted_instruction_code(),
                         baseline.rooted_code())
        changed_labels = V4LabeledTree(
            source.tree, (("f0", 0), ("f1", 1), ("f2", 2), ("f3", 0))
        )
        label_plan = atomic_plumbing_plan(Abbreviated(changed_labels, (2, 0), 1))
        self.assertNotEqual(standard, label_plan)
        self.assertNotEqual(standard.rooted_instruction_code(),
                            label_plan.rooted_instruction_code())
        changed_sector = V4LabeledTree(WeightedPlaneTree((
            Vertex(("f0", "e0", "f1"), (0, 5, 0)),
            source.tree.vertices[1],
        )), source.labels)
        sector_plan = atomic_plumbing_plan(Abbreviated(changed_sector, (2, 0), 1))
        self.assertNotEqual(standard.rooted_instruction_code(),
                            sector_plan.rooted_instruction_code())
        changed_ring = atomic_plumbing_plan(Abbreviated(source, (1, 1), 1))
        self.assertNotEqual(standard.rooted_instruction_code(),
                            changed_ring.rooted_instruction_code())

    def test_flip_actions_have_structural_representations_not_topology(self):
        for sign in (-1, 1):
            original = bounded_valence_chain_presentations(3, sign)[0]
            self.assertEqual(flip_abbreviation_square(original, sign),
                             (True, True))
            baseline = abbreviate_disjoint_rings(original, sign)
            variants = (
                baseline, abbreviated_flip_f1(baseline),
                abbreviated_flip_f2(baseline),
            )
            for abbreviated in variants:
                with self.subTest(sign=sign):
                    plan = atomic_plumbing_plan(abbreviated)
                    self.assertEqual(plan.rooted_instruction_code(),
                                     abbreviated.rooted_code())
                    self.assertEqual(tuple(a.ring_count for a in plan.atoms),
                                     (1, 1, 1))
                    self.assertEqual(plan.topology_status, "UNVERIFIED_GEOMETRY")

    def test_guards_refuse_malformed_abstract_sewing(self):
        base = AtomicInstruction((), (0, 0), 1, 1)
        child = AtomicInstruction((1,), (2,), 0, 1)
        root = FreeBoundary(Port((), 0), 0, True)
        with self.assertRaises(ValueError):
            AtomicInstruction((), (0, 0), True, 1)
        with self.assertRaises(ValueError):
            Port((0,), 1)
        with self.assertRaises(ValueError):
            Seam(Port((), 1), Port((2,), 0))
        with self.assertRaises(ValueError):
            FreeBoundary(Port((), 0), True)
        with self.assertRaises(ValueError):
            AtomicPlumbingPlan(
                1, (base, child),
                (Seam(Port((), 1), Port((1,), 0)),),
                (root, root),
            )
        with self.assertRaises(ValueError):
            AtomicPlumbingPlan(1, (base,), (), (root,))
        with self.assertRaises(ValueError):
            AtomicPlumbingPlan(
                1, (base, child),
                (Seam(Port((), 1), Port((1,), 0)),),
                (root,),
                "VERIFIED_TOPOLOGY",
            )

    def test_invalid_source_and_missing_root_fail_closed(self):
        for source, root in ((None, "f0"), ("wrong", "f0")):
            with self.assertRaises(ValueError):
                atomic_plumbing_plan(source, root)
        source = abbreviate_disjoint_rings(placement_presentations(3, 1)[0], 1)
        for root in ("f99", "e0", None, ""):
            with self.subTest(root=root), self.assertRaises(ValueError):
                atomic_plumbing_plan(source, root)

    def test_bounded_syntax_ring_numbers_supported_without_fake_isotopy(self):
        # Independent A3 source grammar uses only 2/3-child forks and
        # a nonzero ring count at each internal fork; not all such
        # presentations are source-canonical tangles.
        spec = Fork((Fork((Leaf(2), Leaf(3))), Leaf(2)))
        abbreviated = annotate_abbreviated(spec, (2, 1))
        plan = atomic_plumbing_plan(abbreviated)
        self.assertEqual(tuple(a.ring_count for a in plan.atoms[:2]), (2, 1))
        self.assertEqual(plan.rooted_instruction_code(), abbreviated.rooted_code())


if __name__ == "__main__":
    unittest.main()
