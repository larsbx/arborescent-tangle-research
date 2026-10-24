"""A5.2d: exact four-marked PL Conway sphere *patches* and braid corridors.

This constructs a genuine disjoint arrangement of cubical 2-spheres
in Q^3, each with four cyclically marked punctures, together with
spatially embedded two-strand connectors in the complement. The
signed sector list of an unringed A5.1 atom is preserved in order.

This is a CUT-OPEN LINEAR SKELETON of the source Fig. 12.8 / 12.29,
NOT its finished atomic tangle: the left auxiliary sphere is a
bookkeeping anchor, not the source's returning band, and the long
return arcs, twisted band, and 3D plumbing are NOT present.
A marked cubical sphere is not itself an orientation-reversing seam
homeomorphism. Consequently there is NO R+/- or global isotopy proof.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction, Port
from .conway_braid_diagrams import Braid2Diagram
from .spatial_twist_blocks import Point3, SpatialTwoBraid, spatial_two_braid


def _cross(u: tuple[int, int, int], v: tuple[int, int, int]) -> tuple[int, int, int]:
    return (u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0])


def _sub(p: Point3, q: Point3) -> tuple[int, int, int]:
    return int(p.x-q.x), int(p.y-q.y), int(p.z-q.z)


@dataclass(frozen=True)
class MarkedConwayCube:
    """Boundary of a solid axis-aligned cube, with four marked points.

    Coordinate convention: x longitudinal, y across strands and z
    height. Punctures lie inside the west/east square faces and run
    in the cyclic equatorial order NW, NE, SE, SW. This is a CHOSEN
    PL frame, not a source-verified Conway parametrization.
    """
    center_x: int
    port: Port | None = None

    def __post_init__(self) -> None:
        if type(self.center_x) is not int or (self.port is not None
                                             and not isinstance(self.port, Port)):
            raise ValueError("invalid cube position/port")
        self.audit()

    @property
    def west(self) -> int:
        return self.center_x - 2

    @property
    def east(self) -> int:
        return self.center_x + 2

    def punctures(self) -> tuple[tuple[str, Point3], ...]:
        return (("NW", Point3(self.west, -2, 0)),
                ("NE", Point3(self.east, -2, 0)),
                ("SE", Point3(self.east, 2, 0)),
                ("SW", Point3(self.west, 2, 0)))

    def vertices(self) -> tuple[Point3, ...]:
        return tuple(Point3(self.center_x+2*sx, 3*sy, 3*sz)
                     for sx, sy, sz in product((-1, 1), repeat=3))

    def triangles(self) -> tuple[tuple[Point3, Point3, Point3], ...]:
        """A consistently outward oriented 12-triangle cubical sphere."""
        triangles = []
        widths = (2, 3, 3)
        center = (self.center_x, 0, 0)
        # (u,v) basis has cross product +axis on all three axes.
        for axis in range(3):
            u, v = (axis+1)%3, (axis+2)%3
            for sign in (-1, 1):
                corners = []
                for su, sv in ((-1,-1), (1,-1), (1,1), (-1,1)):
                    coords = list(center)
                    coords[axis] += sign*widths[axis]
                    coords[u] += su*widths[u]
                    coords[v] += sv*widths[v]
                    corners.append(Point3(*coords))
                if sign == -1:
                    corners.reverse()
                triangles.extend(((corners[0], corners[1], corners[2]),
                                  (corners[0], corners[2], corners[3])))
        return tuple(triangles)

    def audit(self) -> None:
        """Independent combinatorial surface topology and 3D orientation checks."""
        vertices = set(self.vertices())
        triangles = self.triangles()
        if len(vertices) != 8 or len(triangles) != 12:
            raise ValueError("not a triangulated cube")
        edges: Counter[frozenset[Point3]] = Counter()
        oriented_edges: Counter[tuple[Point3, Point3]] = Counter()
        center = Point3(self.center_x, 0, 0)
        for a,b,c in triangles:
            u,v = _sub(b,a), _sub(c,a)
            normal = _cross(u,v)
            outward = _sub(a,center)
            if sum(normal[k]*outward[k] for k in range(3)) <= 0:
                raise ValueError("surface face orientation is not outward")
            if len({a,b,c}) != 3:
                raise ValueError("degenerate face")
            for start,end in ((a,b),(b,c),(c,a)):
                edges[frozenset((start,end))] += 1
                oriented_edges[(start,end)] += 1
        if (len(edges) != 18 or any(n != 2 for n in edges.values())
                or len(vertices)-len(edges)+len(triangles) != 2
                or any(oriented_edges[(a,b)] != oriented_edges[(b,a)]
                       for a,b in oriented_edges)):
            raise ValueError("invalid closed oriented PL 2-sphere")
        punctures = self.punctures()
        if len({name for name,_ in punctures}) != 4:
            raise ValueError("four named punctures required")
        if len({p for _,p in punctures}) != 4:
            raise ValueError("punctures must be distinct")
        # Each puncture lies strictly inside a distinct west/east face
        # position (not on a mesh edge), with no point at a corner.
        for _,p in punctures:
            if (p.x not in (self.west,self.east)
                    or not (-3 < p.y < 3 and -3 < p.z < 3)):
                raise ValueError("puncture not in the intended boundary face")


@dataclass(frozen=True)
class TwoStrandCorridor:
    """Two exactly embedded arcs between facing marked cube boundaries."""
    left: MarkedConwayCube
    right: MarkedConwayCube
    sector_index: int
    signed_weight: int
    local: SpatialTwoBraid

    def __post_init__(self) -> None:
        self.audit()

    def strands(self) -> tuple[tuple[Point3, ...], tuple[Point3, ...]]:
        def transform(points: tuple[Point3, ...]) -> tuple[Point3, ...]:
            return tuple(Point3(self.left.east + p.y, p.x, p.z) for p in points)
        return transform(self.local.nw), transform(self.local.ne)

    def audit(self) -> None:
        if (not isinstance(self.left, MarkedConwayCube)
                or not isinstance(self.right, MarkedConwayCube)
                or type(self.sector_index) is not int or self.sector_index < 0
                or type(self.signed_weight) is not int
                or not isinstance(self.local, SpatialTwoBraid)):
            raise ValueError("bad Conway corridor data")
        word = ((1,)*self.signed_weight if self.signed_weight >= 0
                else (-1,)*(-self.signed_weight))
        if self.local.word != word:
            raise ValueError("corridor crossings differ from source sector")
        # Constant extra longitudinal length 1 accommodates the
        # zero-crossing case while fixing the boundary-end coordinates.
        if self.right.west - self.left.east != 2*abs(self.signed_weight)+1:
            raise ValueError("gap must fit each exact crossing cell")
        self.local.audit()
        west, east = self.strands()
        if (west[0] != Point3(self.left.east,-2,0)
                or east[0] != Point3(self.left.east,2,0)
                or {west[-1],east[-1]} != {
                    Point3(self.right.west,-2,0),Point3(self.right.west,2,0)
                }):
            raise ValueError("corridor does not meet its four marked face points")
        # Consecutive PL vertices have strictly increasing longitudinal
        # coordinates, so every open segment lies outside both closed
        # cubes; no additional cube can intrude into this exact gap.
        for strand in (west,east):
            if (any(a.x >= b.x for a,b in zip(strand,strand[1:]))
                    or any(abs(p.y)>2 or abs(p.z)>1 for p in strand)):
                raise ValueError("curve escapes its disjoint longitudinal corridor")

    def serialize(self) -> tuple:
        return (self.sector_index, self.signed_weight,
                tuple(tuple(p.serial() for p in arc) for arc in self.strands()))


@dataclass(frozen=True)
class OpenAtomicConwayChain:
    """All unringed sectors, but deliberately no closed outer band."""
    atom: AtomicInstruction
    cubes: tuple[MarkedConwayCube, ...]
    corridors: tuple[TwoStrandCorridor, ...]
    status: Literal["OPEN_CHAIN_NOT_ATOMIC"] = "OPEN_CHAIN_NOT_ATOMIC"

    def __post_init__(self) -> None:
        self.audit()

    def audit(self) -> None:
        if (not isinstance(self.atom, AtomicInstruction)
                or self.atom.ring_count != 0 or self.atom.valence < 2
                or self.status != "OPEN_CHAIN_NOT_ATOMIC"
                or len(self.cubes) != self.atom.valence+1
                or len(self.corridors) != self.atom.valence):
            raise ValueError("not a complete cut-open unringed atomic skeleton")
        expected_centers = [0]
        for w in self.atom.sector_weights:
            expected_centers.append(expected_centers[-1]+5+2*abs(w))
        if tuple(c.center_x for c in self.cubes) != tuple(expected_centers):
            raise ValueError("cube positions cannot overlap or be altered")
        if (self.cubes[0].port is not None
                or tuple(c.port for c in self.cubes[1:]) !=
                tuple(Port(self.atom.address, i) for i in range(self.atom.valence))):
            raise ValueError("missing source port/sector alignment")
        for i, corridor in enumerate(self.corridors):
            if (corridor.left != self.cubes[i]
                    or corridor.right != self.cubes[i+1]
                    or corridor.sector_index != i
                    or corridor.signed_weight != self.atom.sector_weights[i]):
                raise ValueError("corridor disconnected or reordered")
            corridor.audit()
        # The cubes have x half-width 2; consecutive centres are >=5.
        # Thus all closed cubes are disjoint, every corridor interior
        # is in one open gap, and different corridors do not intersect.

    def unused_anchor_endpoints(self) -> tuple[Point3, Point3]:
        return tuple(p for name,p in self.cubes[0].punctures()
                     if name in ("NW","SW"))  # type: ignore[return-value]

    def unused_terminal_endpoints(self) -> tuple[Point3, Point3]:
        return tuple(p for name,p in self.cubes[-1].punctures()
                     if name in ("NE","SE"))  # type: ignore[return-value]

    def canonical_coordinates(self) -> tuple:
        return (self.atom.address, self.atom.sign, self.atom.sector_weights,
                tuple((cube.center_x, tuple((name,p.serial())
                                            for name,p in cube.punctures()))
                      for cube in self.cubes),
                tuple(c.serialize() for c in self.corridors))

    def require_source_atomic_realization(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: missing Fig.12.29 return arcs, complete "
            "twisted band, r-ring embeddings, and verified Conway gluing"
        )


def build_open_unringed_atomic_chain(atom: AtomicInstruction) -> OpenAtomicConwayChain:
    """Embed an exact linear array of holes/twists, with explicit omissions."""
    if not isinstance(atom, AtomicInstruction):
        raise ValueError("A5.1 atom required")
    if atom.ring_count != 0:
        raise NotImplementedError("UNVERIFIED_GEOMETRY: rings cannot be omitted")
    if atom.valence < 2:
        raise NotImplementedError("UNVERIFIED_GEOMETRY: valence-one model differs")
    centers = [0]
    for w in atom.sector_weights:
        centers.append(centers[-1]+5+2*abs(w))
    cubes = tuple(MarkedConwayCube(center, None if i == 0
                                   else Port(atom.address,i-1))
                  for i,center in enumerate(centers))
    corridors = tuple(
        TwoStrandCorridor(
            cubes[i], cubes[i+1], i, w,
            spatial_two_braid(Braid2Diagram(
                (1,)*w if w >= 0 else (-1,)*(-w)
            )),
        ) for i,w in enumerate(atom.sector_weights)
    )
    return OpenAtomicConwayChain(atom,cubes,corridors)
