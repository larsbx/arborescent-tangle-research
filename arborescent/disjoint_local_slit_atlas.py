"""A5.2e2b4c2a: a finite atlas of DISJOINT local ambient 3D slit charts.

Each chart is a real Q^3 tetrahedral product neighborhood of a
chosen two-arrow rectangular patch on a verified candidate band.
All n such local 3-balls are checked for physical disjointness.
Their n simultaneous *local* slit certificates produce 2n cut
3-ball pieces, and exact inverse reglues are independently replayed.

These separate charts are NOT a connected triangulation of the full
ambient S^3 or knot complement. Critically, their cut disks have NO
K-punctures. Thus merely batching local slits CANNOT construct the
source four-punctured Conway spheres. This is a verified intermediate
3D geometric atlas, not the global Bonahon--Siebenmann theorem.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .atomic_plumbing_contract import AtomicInstruction, Port
from .exterior_conway_collars import (
    ExteriorCollaredBandCandidate, build_exterior_collared_band,
)
from .local_plumbing_slit import (
    LocalSlitWitness, TetrahedralProductNeighborhood,
    collar_plumbing_chart, product_neighborhood, slit_product,
)
from .spatial_twist_blocks import Point3


def _bounds(vertices: tuple[Point3, ...]) -> tuple[
    tuple[Fraction, Fraction], tuple[Fraction, Fraction],
    tuple[Fraction, Fraction]
]:
    if not isinstance(vertices,tuple) or not vertices or any(
        not isinstance(p,Point3) for p in vertices
    ):
        raise ValueError("exact nonempty 3D vertex table required")
    return tuple(
        (min(Fraction(getattr(p,axis)) for p in vertices),
         max(Fraction(getattr(p,axis)) for p in vertices))
        for axis in ("x","y","z")
    )  # type: ignore[return-value]


@dataclass(frozen=True)
class LocalSite:
    """One real physical product 3-ball, with two local cut copies."""
    port: Port
    solid: TetrahedralProductNeighborhood
    slit: LocalSlitWitness

    def __post_init__(self) -> None:
        if (not isinstance(self.port,Port)
                or not isinstance(self.solid,TetrahedralProductNeighborhood)
                or not isinstance(self.slit,LocalSlitWitness)
                or self.slit.solid != self.solid
                or not self.slit.verify()):
            raise ValueError("local slit site lacks exact replayable 3D data")

    def spatial_bounds(self):
        return _bounds(self.solid.vertices)

    def canonical(self):
        return (
            self.port.address,self.port.slot,
            tuple(v.serial() for v in self.solid.vertices),
            self.solid.tetrahedra,
            tuple(tuple((v.original_index,v.side) for v in t)
                  for t in self.slit.cut_tetrahedra),
            tuple(tuple((v.original_index,v.side) for v in f)
                  for f in self.slit.lower_cut),
            tuple(tuple((v.original_index,v.side) for v in f)
                  for f in self.slit.upper_cut),
        )


@dataclass(frozen=True)
class LocalSlitAtlasMetrics:
    local_product_components: int
    tetrahedra_before_slit: int
    tetrahedra_after_slit: int
    local_cut_3balls: int
    genuine_conway_boundary_spheres: int
    knot_piercings_of_cut_disks: int
    minimum_longitudinal_gap: Fraction | None
    total_local_volume: Fraction

    def __post_init__(self) -> None:
        values=(self.local_product_components,self.tetrahedra_before_slit,
                self.tetrahedra_after_slit,self.local_cut_3balls,
                self.genuine_conway_boundary_spheres,
                self.knot_piercings_of_cut_disks)
        if any(type(v) is not int or v<0 for v in values):
            raise ValueError("nonnegative exact local atlas metrics required")
        if self.minimum_longitudinal_gap is not None and (
            type(self.minimum_longitudinal_gap) is not Fraction
            or self.minimum_longitudinal_gap<=0
        ):
            raise ValueError("chart separation must be positive and exact")
        if type(self.total_local_volume) is not Fraction or self.total_local_volume<=0:
            raise ValueError("physical local product volume required")


@dataclass(frozen=True)
class DisjointLocalSlitAtlas:
    """n fully checked DISCONNECTED product charts on one chosen band."""
    candidate: ExteriorCollaredBandCandidate
    sites: tuple[LocalSite, ...]
    thickness: Fraction
    status: Literal["DISJOINT_LOCAL_ATLAS_NOT_GLOBAL_MANIFOLD"] = (
        "DISJOINT_LOCAL_ATLAS_NOT_GLOBAL_MANIFOLD"
    )

    def __post_init__(self) -> None:
        self.audit()

    def audit(self) -> LocalSlitAtlasMetrics:
        if (not isinstance(self.candidate,ExteriorCollaredBandCandidate)
                or not isinstance(self.sites,tuple)
                or type(self.thickness) is not Fraction
                or not Fraction(0)<self.thickness<=Fraction(1,4)
                or self.status!="DISJOINT_LOCAL_ATLAS_NOT_GLOBAL_MANIFOLD"):
            raise ValueError("restricted exact local 3D slit atlas required")

        cubes=self.candidate.patches.strands.cubes
        if (len(self.sites)!=len(cubes)
                or len(cubes)<2
                or len(self.candidate.collars)!=len(cubes)):
            raise ValueError("atlas must cover every physical candidate boundary")
        if tuple(site.port for site in self.sites)!=tuple(c.port for c in cubes):
            raise ValueError("site order or source port address differs")

        bounds=[]
        total_volume=Fraction(0)
        for i,site in enumerate(self.sites):
            if not isinstance(site,LocalSite):
                raise ValueError("site must contain its typed cut and 3-ball")
            site.solid.audit()
            if (site.solid.half_thickness!=self.thickness
                    or site.solid.patch!=collar_plumbing_chart(
                        self.candidate.collars[i])):
                raise ValueError("local 3D chart not drawn from correct PL band faces")
            if not site.slit.verify():
                raise ValueError("local 3D cut/reglue cannot replay")
            xyz=site.spatial_bounds()
            cube=cubes[i]
            # Source candidate's closed Conway cube is contained in
            # z<=3; local 3-ball rests strictly above that: z>=3.75.
            if not xyz[2][0]>Fraction(3):
                raise ValueError("product intrudes into Conway ball roof")
            if not xyz[2][1]>=xyz[2][0]:
                raise ValueError("degenerate 3D product height")
            bounds.append(xyz)
            dx=xyz[0][1]-xyz[0][0]
            dy=xyz[1][1]-xyz[1][0]
            dz=xyz[2][1]-xyz[2][0]
            total_volume+=dx*dy*dz

        gap=None
        for left,right in zip(bounds,bounds[1:]):
            # x is the global spatial longitudinal axis, with source
            # ports in precisely the verified A5.2e2b3b cube order.
            if not left[0][1]<right[0][0]:
                raise ValueError("two distinct ambient 3-balls intersect")
            clearance=right[0][0]-left[0][1]
            gap=clearance if gap is None else min(gap,clearance)
        # Check source proper K's ACTUAL coordinate vertices, not an
        # assumed crossing word. All segments interpolate their ends,
        # so a strict bounding z-separation certifies no K piercing.
        k_arcs=self.candidate.patches.strands.embedded_arcs()
        knot_zmax=max(Fraction(p.z) for arc in k_arcs for p in arc)
        knot_zmin=min(Fraction(p.z) for arc in k_arcs for p in arc)
        for xyz in bounds:
            if not (knot_zmax<xyz[2][0] or knot_zmin>xyz[2][1]):
                raise ValueError("K pierces local 3D cut neighborhood")

        n=len(self.sites)
        return LocalSlitAtlasMetrics(
            n,12*n,12*n,2*n,
            # Cut product 3-ball boundaries are unpunctured: neither
            # is the source's four-punctured Conway sphere!
            0,0,gap,total_volume)

    def canonical(self):
        return (self.candidate.patches.strands.atom.address,
                self.candidate.patches.strands.atom.sector_weights,
                str(self.thickness),
                tuple(site.canonical() for site in self.sites))

    def require_global_source_slitting(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: a DISJOINT atlas of local 3-balls "
            "contains no global ambient S3 triangulation and zero "
            "K-piercings, hence does not create 4-punctured source "
            "Conway boundary spheres or verify source isotopy")


def build_disjoint_local_slit_atlas(
    atom: AtomicInstruction,
    thickness: Fraction = Fraction(1,16),
) -> DisjointLocalSlitAtlas:
    if (not isinstance(atom,AtomicInstruction)
            or type(thickness) is not Fraction
            or not Fraction(0)<thickness<=Fraction(1,4)):
        raise ValueError("supported atom and 0<h<=1/4 exact Fraction required")
    candidate=build_exterior_collared_band(atom)
    sites=[]
    for cube,collar in zip(candidate.patches.strands.cubes,candidate.collars):
        patch=collar_plumbing_chart(collar)
        solid=product_neighborhood(patch,thickness)
        sites.append(LocalSite(cube.port,solid,slit_product(solid)))
    return DisjointLocalSlitAtlas(candidate,tuple(sites),thickness)
