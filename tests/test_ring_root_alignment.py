"""A4.3: ring normal keys cannot depend on arbitrary cyclic start slots."""
import unittest

from arborescent.placement_fibers import (
    bounded_valence_chain_presentations, placement_presentations,
)
from arborescent.ring_normalizer import (
    align_rooted_cyclic_slots, normalize_ring_positions,
    normalize_rooted_ring_positions, strict_ring_normal_key,
)
from arborescent.ring_relocation import explore_ring_slides
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def rotate_storage(source, shifts):
    """No graph or angular-order change: rotate slot AND sector arrays."""
    vertices = []
    for i, v in enumerate(source.tree.vertices):
        offset = shifts[i] % len(v.slots)
        vertices.append(Vertex(
            v.slots[offset:] + v.slots[:offset],
            v.sectors[offset:] + v.sectors[:offset],
        ))
    return V4LabeledTree(WeightedPlaneTree(tuple(vertices)), source.labels)


def renumber_vertices_and_bonds(source):
    """Presentation isomorphism: reverse vertex storage and rename edges."""
    names = sorted({bond for v in source.tree.vertices for bond in v.slots
                    if bond.startswith("e")})
    mapping = {name: f"e{100+i}" for i, name in enumerate(names)}
    vertices = tuple(
        Vertex(tuple(mapping.get(s, s) for s in v.slots), v.sectors)
        for v in reversed(source.tree.vertices)
    )
    return V4LabeledTree(WeightedPlaneTree(vertices), source.labels)


class RootAlignedRingTests(unittest.TestCase):
    def test_cyclic_start_slots_cannot_change_normal_code(self):
        for sign in (-1, 1):
            for object_ in placement_presentations(4, 2, sign):
                base = strict_ring_normal_key(object_, sign)
                versions = (
                    rotate_storage(object_, [i+1 for i in
                                             range(len(object_.tree.vertices))]),
                    rotate_storage(object_, [2*i+1 for i in
                                             range(len(object_.tree.vertices))]),
                    renumber_vertices_and_bonds(object_),
                )
                for candidate in versions:
                    with self.subTest(sign=sign):
                        self.assertEqual(
                            candidate.labelled_rooted_code(),
                            object_.labelled_rooted_code(),
                        )
                        self.assertEqual(strict_ring_normal_key(candidate, sign), base)
                        cert = normalize_rooted_ring_positions(candidate, sign)
                        self.assertTrue(cert.verify())
                        self.assertEqual(
                            strict_ring_normal_key(cert.normal, sign), base
                        )

    def test_slot_zero_cap_is_only_a_storage_artifact(self):
        src = placement_presentations(3, 1)[0]
        rotated = rotate_storage(src, [1]+[0]*(len(src.tree.vertices)-1))
        with self.assertRaisesRegex(ValueError, "anchored"):
            normalize_ring_positions(rotated, 1)
        cert = normalize_rooted_ring_positions(rotated, 1)
        self.assertTrue(cert.verify())
        self.assertEqual(cert.aligned, align_rooted_cyclic_slots(rotated))
        self.assertEqual(
            strict_ring_normal_key(src, 1),
            strict_ring_normal_key(rotated, 1),
        )
        self.assertEqual(cert.source, rotated)

    def test_bounded_valence_full_orbit_and_cyclic_storage(self):
        for sign in (-1, 1):
            for bases in (2, 3, 4):
                source = bounded_valence_chain_presentations(bases, sign)[0]
                reachable = explore_ring_slides(source, sign)
                self.assertTrue(reachable.complete)
                keys = set()
                for state in reachable.incoming:
                    keys.add(strict_ring_normal_key(state, sign))
                    stored = rotate_storage(
                        state, [i+1 for i in range(len(state.tree.vertices))]
                    )
                    keys.add(strict_ring_normal_key(stored, sign))
                    self.assertTrue(normalize_rooted_ring_positions(stored, sign).verify())
                self.assertEqual(len(keys), 1)

    def test_all_marked_v4_labels_and_sign_remain_in_key(self):
        src = placement_presentations(3, 1)[0]
        original = strict_ring_normal_key(src, 1)
        edited = V4LabeledTree(
            src.tree, tuple((bond, 3 if bond == "f2" else value)
                            for bond, value in src.labels)
        )
        self.assertNotEqual(strict_ring_normal_key(edited, 1), original)
        self.assertNotEqual(
            strict_ring_normal_key(rotate_storage(edited, [1,1,1,1]), 1),
            original,
        )
        with self.assertRaises(ValueError):
            strict_ring_normal_key(src, -1)

    def test_root_validation_and_idempotence(self):
        src = placement_presentations(3, 2)[0]
        for root in ("missing", "e0", "", 12, "f100"):
            with self.subTest(root=root), self.assertRaises(ValueError):
                align_rooted_cyclic_slots(src, root)
        aligned = align_rooted_cyclic_slots(src)
        self.assertEqual(align_rooted_cyclic_slots(aligned), aligned)
        self.assertEqual(src.labelled_rooted_code(), aligned.labelled_rooted_code())
        cert = normalize_rooted_ring_positions(src, 1)
        self.assertEqual(
            normalize_rooted_ring_positions(cert.normal, 1).slide_certificate.moves,
            (),
        )

    def test_nonzero_sector_rotation_preserves_weights_and_code(self):
        base = Vertex(("f0", "e0", "e1"), (3, 0, 0))
        from arborescent.weighted_trees import WeightedPlaneTree
        s = V4LabeledTree(WeightedPlaneTree((
            base, Vertex(("e0",), (2,)), Vertex(("e1",), (3,))
        )), (("f0", 1),))
        shifted = rotate_storage(s, [1, 0, 0])
        self.assertEqual(s.labelled_rooted_code(), shifted.labelled_rooted_code())
        self.assertEqual(align_rooted_cyclic_slots(s),
                         align_rooted_cyclic_slots(shifted))


if __name__ == "__main__":
    unittest.main()
