"""A4.10 Gaussian distance ranks and actual shortest formal replay paths."""
import unittest
from collections import Counter
from itertools import product
from math import comb

from arborescent.ring_geodesics import (
    gaussian_ring_polynomial, path_rank_polynomial,
    formal_slide_distance, formal_slide_geodesic,
)
from arborescent.ring_normalizer import (
    normalize_rooted_ring_positions, ring_position_inversions,
)
from arborescent.ring_relocation import explore_ring_slides
from arborescent.placement_fibers import (
    bounded_valence_chain_presentations, placement_presentations,
    placement_fiber_formula,
)
from arborescent.tangle_labels import V4LabeledTree
from arborescent.weighted_trees import Vertex, WeightedPlaneTree


def rename_and_rotate(source):
    """A presentation-isomorphic copy with same vertex array (new edge IDs)."""
    renamed = []
    for vertex in source.tree.vertices:
        slots = tuple(
            "e" + str(100 + int(b[1:])) if b.startswith("e") else b
            for b in vertex.slots
        )
        renamed.append(Vertex(
            slots[1:] + slots[:1], vertex.sectors[1:] + vertex.sectors[:1]
        ))
    return V4LabeledTree(WeightedPlaneTree(tuple(renamed)), source.labels)


def rotate_cyclic_slots(source):
    return V4LabeledTree(WeightedPlaneTree(tuple(
        Vertex(v.slots[1:]+v.slots[:1], v.sectors[1:]+v.sectors[:1])
        for v in source.tree.vertices
    )), source.labels)


class GaussianRingDistanceTests(unittest.TestCase):
    def test_known_q_binomials_and_boundary_cases(self):
        self.assertEqual(gaussian_ring_polynomial(0, 5), (1,))
        self.assertEqual(gaussian_ring_polynomial(5, 0), (1,))
        self.assertEqual(gaussian_ring_polynomial(1, 2), (1, 1, 1))
        self.assertEqual(gaussian_ring_polynomial(2, 2), (1, 1, 2, 1, 1))
        self.assertEqual(gaussian_ring_polynomial(3, 2),
                         (1, 1, 2, 2, 2, 1, 1))
        for o, r in product(range(5), repeat=2):
            p = gaussian_ring_polynomial(o, r)
            self.assertEqual(sum(p), comb(o+r, r))
            self.assertEqual(p, p[::-1])
            self.assertEqual(len(p), o*r+1 if o and r else 1)

    def test_actual_corolla_rank_histograms_match_q_binomial(self):
        # Completely independent from DP: enumerate stored placements,
        # explicitly count R-before-O inversions with the normalizer.
        for sign in (-1, 1):
            for d, r in ((3, 1), (3, 2), (3, 3), (4, 2), (5, 2)):
                with self.subTest(sign=sign, d=d, rings=r):
                    sources = placement_presentations(d, r, sign)
                    ranks = Counter(ring_position_inversions(s, sign)
                                    for s in sources)
                    actual = tuple(ranks.get(i, 0)
                                   for i in range((d-1)*r+1))
                    self.assertEqual(actual,
                                     gaussian_ring_polynomial(d-1, r))
                    self.assertEqual(len(sources),
                                     placement_fiber_formula(d, r))

    def test_full_bounded_degree_path_rank_polynomial(self):
        for sign in (-1, 1):
            for k in range(2, 7):
                with self.subTest(sign=sign, k=k):
                    seed = bounded_valence_chain_presentations(k, sign)[0]
                    orbit = explore_ring_slides(seed, sign, max_states=4096)
                    self.assertTrue(orbit.complete)
                    self.assertEqual(len(orbit.incoming), 4*3**(k-2))
                    ranks = Counter(
                        ring_position_inversions(s, sign)
                        for s in orbit.incoming
                    )
                    expected = path_rank_polynomial(k)
                    self.assertEqual(
                        tuple(ranks.get(i, 0) for i in range(len(expected))),
                        expected,
                    )
                    self.assertEqual(sum(expected), 4*3**(k-2))
                    self.assertEqual(expected, expected[::-1])

    def test_exact_geodesic_matches_bfs_shortest_paths(self):
        for sign in (-1, 1):
            for d, r in ((3, 1), (3, 2), (4, 2), (4, 3)):
                with self.subTest(sign=sign, d=d, r=r):
                    seed = placement_presentations(d, r, sign)[0]
                    orbit = explore_ring_slides(seed, sign)
                    self.assertTrue(orbit.complete)
                    for state in orbit.incoming:
                        expected = len(orbit.trace(state).moves)
                        self.assertEqual(
                            formal_slide_distance(seed, state, sign),
                            expected,
                        )
                        certificate = formal_slide_geodesic(seed, state, sign)
                        self.assertEqual(certificate.distance, expected)
                        self.assertTrue(certificate.verify())
                        self.assertEqual(certificate.trace.target,
                                         certificate.aligned_target)

    def test_geodesics_between_arbitrary_states_both_directions(self):
        for sign in (-1, 1):
            sources = placement_presentations(4, 2, sign)
            for start in sources:
                for goal in sources:
                    forward = formal_slide_geodesic(start, goal, sign)
                    backward = formal_slide_geodesic(goal, start, sign)
                    self.assertEqual(forward.distance, backward.distance)
                    self.assertTrue(forward.verify())
                    self.assertTrue(backward.verify())

    def test_bounded_path_interior_geodesics_match_bfs(self):
        for sign in (-1, 1):
            for k in (2, 3, 4):
                seed = bounded_valence_chain_presentations(k, sign)[0]
                orbit = explore_ring_slides(seed, sign)
                for state in orbit.incoming:
                    steps = formal_slide_geodesic(seed, state, sign)
                    self.assertEqual(steps.distance,
                                     len(orbit.trace(state).moves))
                    self.assertTrue(steps.verify())

    def test_cyclic_storage_isomorphism_alignment_is_zero_cost(self):
        for sign in (-1, 1):
            source = placement_presentations(3, 2, sign)[0]
            equivalent = rotate_cyclic_slots(source)
            self.assertEqual(formal_slide_distance(source, equivalent, sign), 0)
            cert = formal_slide_geodesic(equivalent, source, sign)
            self.assertEqual(cert.distance, 0)
            self.assertEqual(cert.trace.moves, ())
            self.assertTrue(cert.verify())

    def test_invalid_comparisons_fail_closed(self):
        source = placement_presentations(3, 1)[0]
        other_sign = placement_presentations(3, 1, -1)[0]
        invalids = (
            (source, other_sign, 1),
            (source, source, -1),
            (source, None, 1),
            (source, placement_presentations(3, 2)[0], 1),
        )
        for start, goal, sign in invalids:
            with self.subTest(sign=sign, target=type(goal)):
                with self.assertRaises(ValueError):
                    formal_slide_distance(start, goal, sign)
        for o, rings in ((-1, 1), (1, -1), (True, 1), (2, False)):
            with self.assertRaises(ValueError):
                gaussian_ring_polynomial(o, rings)
        for k in (0, 1, True, 2.4):
            with self.assertRaises(ValueError):
                path_rank_polynomial(k)
        # Formal geodesics are intentionally *same indexed tree*
        # paths, and do not silently identify graph/edge relabellings.
        with self.assertRaises(ValueError):
            formal_slide_distance(source, rename_and_rotate(source), 1)


if __name__ == "__main__":
    unittest.main()
