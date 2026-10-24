"""A5.2e2b4b: a marked K ⊂ ∂B is NOT the whole band boundary.

Bonahon–Siebenmann §12.2 printed p.142 Fig.12.8 describes the band
as a surface 'with boundary CONTAINING the knot'. The model atomic
TANGLE is (M,K), with each Conway sphere in ∂M meeting K at FOUR
points. Thus for n source boundaries, K ∩ M must be 2n properly
embedded ARC components with 4n individual boundary endpoints.

A5.2e2b2–b3b correctly constructs closed annular/Möbius PL band
candidates, but those entire band-boundary circles are NOT the knot:
the source K is a *distinguished subset* of the band boundary,
omitting the newly added transverse/collar rim edges. If a source
comparison equated ∂B and K, it would misidentify a closed link with
a properly embedded tangle. This independent audit enforces the
distinction directly on actual 3D triangle coordinates.

The positive result is only about the *chosen* source-UNVERIFIED
candidate; neither source Fig12.10 framing nor R± is certified.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_return_strands import independent_pl_arc_audit, Arc
from .marked_conway_geometry import MarkedConwayCube
from .spatial_twist_blocks import Point3

Triangle = tuple[Point3,Point3,Point3]
Edge = frozenset[Point3]


def _edge(a: Point3, b: Point3) -> Edge:
    if not isinstance(a,Point3) or not isinstance(b,Point3) or a==b:
        raise ValueError("a PL graph edge has distinct exact rational endpoints")
    return frozenset((a,b))


def _graph_components(graph: dict[Point3,set[Point3]]) -> int:
    remaining=set(graph)
    count=0
    while remaining:
        count+=1
        pending=[remaining.pop()]
        while pending:
            for next_node in graph[pending.pop()]:
                if next_node in remaining:
                    remaining.remove(next_node)
                    pending.append(next_node)
    return count


def _boundary_edges(triangles: tuple[Triangle,...]) -> set[Edge]:
    if not isinstance(triangles,tuple) or not triangles:
        raise ValueError("explicit 3D band triangles required")
    edges: Counter[Edge]=Counter()
    for tri in triangles:
        if (not isinstance(tri,tuple) or len(tri)!=3
                or len(set(tri))!=3
                or any(not isinstance(p,Point3) for p in tri)):
            raise ValueError("invalid PL surface triangle")
        for p,q in zip(tri,tri[1:]+tri[:1]):
            edges[_edge(p,q)]+=1
    if any(k>2 for k in edges.values()):
        raise ValueError("band has nonmanifold edge multiplicity")
    return {edge for edge,count in edges.items() if count==1}


def _edge_graph(edges: set[Edge]) -> dict[Point3,set[Point3]]:
    graph: dict[Point3,set[Point3]]=defaultdict(set)
    for edge in edges:
        a,b=tuple(edge)
        graph[a].add(b)
        graph[b].add(a)
    return graph


@dataclass(frozen=True)
class DistinguishedKnotAudit:
    band_boundary_loops: int
    band_boundary_edges: int
    distinguished_knot_arcs: int
    marked_conway_endpoints: int
    knot_edges: int
    complementary_band_edges: int
    status: Literal["MARKED_KNOT_SOURCE_UNVERIFIED"] = (
        "MARKED_KNOT_SOURCE_UNVERIFIED"
    )


def verify_distinguished_knot(
    triangles: tuple[Triangle,...],
    cubes: tuple[MarkedConwayCube,...],
    knot_arcs: tuple[Arc,...]
) -> DistinguishedKnotAudit:
    """Independent relative-boundary marked K-on-B checker.

    Input contains ONLY physical triangles, source-labelled sphere
    cubes and proposed physical knot arcs. It does NOT consult an
    atom instruction's claimed weights or a stored result status.

    Verifies:
      * ∂B has closed cycles, each vertex degree exactly 2;
      * source K consists of 2n arcs, with all 4n punctures used once;
      * EVERY K segment is an actual exposed edge of the surface ∂B,
        not an interior edge or an invented chord between punctures;
      * its ONLY degree-1 vertices are precisely source punctures;
      * K approaches each Conway sphere from outside with NONZERO
        normal component at that puncture (PL transversality).
    """
    if (not isinstance(cubes,tuple) or not cubes
            or not isinstance(knot_arcs,tuple)):
        raise ValueError("Conway cubes and knot arcs must be exact tuples")
    if len(knot_arcs)!=2*len(cubes):
        raise ValueError("proper source tangle needs exactly 2n K arcs")
    punctures={}
    for cube in cubes:
        if not isinstance(cube,MarkedConwayCube) or cube.port is None:
            raise ValueError("all Conway spheres must be source-port marked")
        for name,p in cube.punctures():
            if p in punctures:
                raise ValueError("overlapping or duplicated marked Conway punctures")
            punctures[p]=(cube,name)
    if len(punctures)!=4*len(cubes):
        raise ValueError("missing source Conway boundary mark")
    # Independent physical arc tester checks disjointness,
    # closed-cube avoidance and complete proper endpoint matching.
    physical_pairings=independent_pl_arc_audit(cubes,knot_arcs)
    if len(physical_pairings)!=2*len(cubes):
        raise ValueError("bad boundary-marked source strand count")

    all_band_boundary=_boundary_edges(triangles)
    whole_graph=_edge_graph(all_band_boundary)
    if not whole_graph or any(len(nbr)!=2 for nbr in whole_graph.values()):
        raise ValueError("completed ribbon must have circular whole-band boundary")
    whole_loops=_graph_components(whole_graph)

    knot_edges=set()
    for arc in knot_arcs:
        for p,q in zip(arc,arc[1:]):
            e=_edge(p,q)
            if e in knot_edges or e not in all_band_boundary:
                raise ValueError("distinguished knot edge is repeated or not on band boundary")
            knot_edges.add(e)
        # Actual geometric transverse contact with each Conway sphere:
        # the first segment must leave its boundary face OUTWARD along
        # x, and the last must arrive from the exterior along x.
        for terminal,near in ((arc[0],arc[1]),(arc[-1],arc[-2])):
            cube,_=punctures[terminal]
            if terminal.x==cube.west:
                if not Fraction(near.x)<Fraction(cube.west):
                    raise ValueError("K tangent or enters the west Conway boundary")
            elif terminal.x==cube.east:
                if not Fraction(near.x)>Fraction(cube.east):
                    raise ValueError("K tangent or enters the east Conway boundary")
            else:
                raise ValueError("K endpoint not on a prescribed puncture face")
    if not knot_edges:
        raise ValueError("distinguished knot has no actual physical edges")
    knot_graph=_edge_graph(knot_edges)
    degree1={p for p,nbr in knot_graph.items() if len(nbr)==1}
    if degree1!=set(punctures):
        raise ValueError("K's true degree-1 endpoints differ from four Conway marks")
    if any(len(neighbors) not in (1,2) for neighbors in knot_graph.values()):
        raise ValueError("K is not a properly embedded PL one-manifold")
    if _graph_components(knot_graph)!=2*len(cubes):
        raise ValueError("K has closed components or wrong proper arc count")

    # The source wording is 'band boundary CONTAINING K',
    # NOT 'K = the entire band boundary'.
    if knot_edges==all_band_boundary:
        raise ValueError("whole closed band boundary cannot be source K proper arcs")
    return DistinguishedKnotAudit(
        whole_loops,len(all_band_boundary),2*len(cubes),
        len(punctures),len(knot_edges),
        len(all_band_boundary)-len(knot_edges)
    )


def audit_source_marking_on_candidate(candidate) -> DistinguishedKnotAudit:
    """Integrate with an existing physically embedded r=0 PL band candidate."""
    # Avoid importing the heavy geometry stack until actually needed.
    from .exterior_conway_collars import ExteriorCollaredBandCandidate
    if not isinstance(candidate,ExteriorCollaredBandCandidate):
        raise ValueError("verified *chosen* exterior band candidate required")
    candidate.audit()
    return verify_distinguished_knot(
        candidate.triangles(),
        candidate.patches.strands.cubes,
        candidate.patches.strands.embedded_arcs()
    )
