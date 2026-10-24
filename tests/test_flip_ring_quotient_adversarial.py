"""Independent finite checks of the restricted formal calculus, not isotopy."""
from dataclasses import replace
from itertools import combinations, product
from math import comb
import unittest

from arborescent.flip_ring_quotient import (
    apply_joint, explore_joint, flip_ring_normal_key, joint_successors,
    normalize_flip_ring,
)
from arborescent.placement_fibers import (
    placement_presentations, bounded_valence_chain_presentations,
)
from arborescent.ring_normalizer import normalize_rooted_ring_positions
from arborescent.tangle_labels import V4LabeledTree, tangle_flip_group
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def rotate(source, offsets):
    return V4LabeledTree(WeightedPlaneTree(tuple(
        Vertex(v.slots[n:] + v.slots[:n], v.sectors[n:] + v.sectors[:n])
        for v, n in zip(source.tree.vertices, offsets)
    )), source.labels)


def independent_placements(source, sign):
    """Enumerate binary shuffles directly, without slides, BFS or normalizers.

    In these explicit corolla/path fixtures, ring centers are precisely
    the vertices of total weight sign adjacent to a zero-weight base.
    Keep the anchor and each of the R/R and O/O subsequences fixed.
    """
    vertices = source.tree.vertices
    owners = {}
    for i, v in enumerate(vertices):
        for b in v.slots:
            owners.setdefault(b, []).append(i)
    options = []
    for i, v in enumerate(vertices):
        if any(v.sectors):
            options.append((v,))
            continue
        rings = tuple(b for b in v.slots[1:] if b.startswith('e')
                      and vertices[next(j for j in owners[b] if j != i)].weight == sign)
        ordinary = tuple(b for b in v.slots[1:] if b not in rings)
        choices = []
        for positions in combinations(range(len(v.slots)-1), len(rings)):
            r, o = iter(rings), iter(ordinary)
            slots = (v.slots[0],) + tuple(
                next(r) if j in positions else next(o)
                for j in range(len(v.slots)-1))
            choices.append(Vertex(slots, v.sectors))
        options.append(tuple(choices))
    return {V4LabeledTree(WeightedPlaneTree(vs), source.labels)
            for vs in product(*options)}


class AdversarialJointTests(unittest.TestCase):
    def test_nested_root_and_sign_substitutions_are_rejected(self):
        for sign in (-1, 1):
            for source in placement_presentations(3, 1, sign):
                c = normalize_flip_ring(source, sign)
                self.assertTrue(c.verify())
                other = normalize_rooted_ring_positions(c.flip_gauge, sign, 'f1')
                self.assertFalse(replace(c, slide=other).verify())
                self.assertFalse(replace(c, sign=-sign).verify())
                self.assertFalse(replace(c, slide=replace(
                    c.slide, slide_certificate=replace(
                        c.slide.slide_certificate, sign=-sign))).verify())
                for bad in (True, False, 1.0, 0, None):
                    self.assertFalse(replace(c, sign=bad).verify())
        # Sign binding matters even when there are no caps to detect it.
        c = normalize_flip_ring(placement_presentations(3, 0)[0], 1)
        self.assertFalse(replace(c, sign=-1).verify())

    def test_each_cyclic_start_and_simultaneous_rotations_replay(self):
        for sign in (-1, 1):
            for source in (placement_presentations(3, 1, sign)[0],
                           bounded_valence_chain_presentations(3, sign)[0]):
                n = len(source.tree.vertices)
                offsets = [tuple(1 % len(v.slots) for v in source.tree.vertices)]
                for i, v in enumerate(source.tree.vertices):
                    for shift in range(len(v.slots)):
                        offsets.append(tuple(shift if j == i else 0 for j in range(n)))
                reference = explore_joint(source, sign)
                for shifts in offsets:
                    rotated = rotate(source, shifts)
                    self.assertEqual(flip_ring_normal_key(rotated, sign),
                                     flip_ring_normal_key(source, sign))
                    children = joint_successors(rotated, sign)
                    self.assertEqual({c.labelled_rooted_code() for _, c in children},
                                     {c.labelled_rooted_code() for _, c in
                                      joint_successors(source, sign)})
                    for move, child in children:
                        self.assertEqual(apply_joint(rotated, move), child)
                    orbit = explore_joint(rotated, sign)
                    self.assertTrue(orbit.complete)
                    self.assertEqual(len(orbit.incoming), len(reference.incoming))
                    self.assertEqual({t.labelled_rooted_code() for t in orbit.incoming},
                                     {t.labelled_rooted_code() for t in reference.incoming})
                    self.assertTrue(all(orbit.trace(t).verify() for t in orbit.incoming))

    def test_nondefault_root_replay(self):
        for sign in (-1, 1):
            s = placement_presentations(3, 1, sign)[0]
            s = rotate(s, tuple(1 % len(v.slots) for v in s.tree.vertices))
            orbit = explore_joint(s, sign, root='f1')
            self.assertTrue(orbit.complete)
            self.assertEqual(len(orbit.incoming), 12)
            self.assertTrue(all(orbit.trace(t).verify() for t in orbit.incoming))
            self.assertEqual(len({flip_ring_normal_key(t, sign, 'f1')
                                  for t in orbit.incoming}), 1)

    def test_independent_shuffle_orbits_and_normal_form_claims(self):
        for sign in (-1, 1):
            fixtures = [(placement_presentations(d, r, sign)[0], 4*comb(d+r-1, r))
                        for d, r in ((3, 0), (3, 1), (3, 2), (4, 2))]
            fixtures += [(bounded_valence_chain_presentations(k, sign)[0],
                          16*3**(k-2)) for k in (2, 3, 4)]
            for source, count in fixtures:
                placements = independent_placements(source, sign)
                expected = {q for p in placements for q in tangle_flip_group(p)}
                orbit = explore_joint(source, sign)
                self.assertTrue(orbit.complete)
                self.assertEqual(len(expected), count)
                self.assertEqual(set(orbit.incoming), expected)
                self.assertEqual(len({flip_ring_normal_key(t, sign) for t in expected}), 1)


if __name__ == '__main__':
    unittest.main()
