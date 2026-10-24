"""A5.1: abstract port/sewing instructions for §12.5 atomic tangles with rings.

Bonahon--Siebenmann §12.5 printed p.158, Figures 12.28--12.29,
defines a *model atomic tangle with rings* at each weighted vertex.
The sectors, cyclic bond positions and ring numbers specify the
intended local model, and internal bonds specify plumbing seams.

THIS MODULE IS NOT THAT TOPOLOGICAL REALIZATION. It constructs a
typed, losslessly auditable, INDEX/FREE-BOND-RENAME-INDEPENDENT
*instruction graph* with an explicit port for every incident bond.
It does not choose geometric strand pairings, crossing conventions,
3D embeddings or V4 parametrization maps across Conway seams.

In particular a sewing instruction is NOT an ambient-isotopy witness.
The public geometry promotion API fails closed until an independent
Conway/plumbing implementation and verified certificates exist.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .ring_abbreviation import Abbreviated
from .tangle_labels import V4LabeledTree


@dataclass(frozen=True, order=True)
class Port:
    """Canonical rooted address and position in the parent-anchored cyclic list."""
    address: tuple[int, ...]
    slot: int

    def __post_init__(self) -> None:
        if (not isinstance(self.address, tuple)
                or any(type(i) is not int or i < 1 for i in self.address)
                or type(self.slot) is not int or self.slot < 0):
            raise ValueError("invalid rooted address or port position")


@dataclass(frozen=True)
class AtomicInstruction:
    """Combinatorial information required at an atomic-ring-tangle vertex.

    This does not assert how to draw its strands in a 3-manifold.
    """
    address: tuple[int, ...]
    sector_weights: tuple[int, ...]
    ring_count: int
    sign: int

    def __post_init__(self) -> None:
        if (not isinstance(self.address, tuple)
                or any(type(i) is not int or i < 1 for i in self.address)
                or not isinstance(self.sector_weights, tuple)
                or len(self.sector_weights) == 0
                or any(type(w) is not int for w in self.sector_weights)
                or type(self.ring_count) is not int or self.ring_count < 0
                or type(self.sign) is not int or self.sign not in (-1, 1)):
            raise ValueError("invalid atomic instruction")

    @property
    def valence(self) -> int:
        return len(self.sector_weights)

    @property
    def total_weight(self) -> int:
        return sum(self.sector_weights)


@dataclass(frozen=True)
class Seam:
    """An UNORIENTED combinatorial plumbing seam; no 3D gluing chosen."""
    parent: Port
    child: Port

    def __post_init__(self) -> None:
        if not (isinstance(self.parent, Port)
                and isinstance(self.child, Port)
                and self.child.address == (*self.parent.address, self.parent.slot)
                and self.child.slot == 0):
            raise ValueError("a seam must connect an ordinary parent port to child 0")


@dataclass(frozen=True)
class FreeBoundary:
    """A free bond with source V4 label, including the distinguished root."""
    port: Port
    v4: int
    marked_root: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.port, Port) or type(self.v4) is not int or self.v4 not in range(4):
            raise ValueError("invalid free-boundary V4 data")
        if type(self.marked_root) is not bool:
            raise ValueError("root indicator must be a bool")


@dataclass(frozen=True)
class AtomicPlumbingPlan:
    """Lossless abstract sewing plan. It has NO topological soundness status."""
    sign: int
    atoms: tuple[AtomicInstruction, ...]
    seams: tuple[Seam, ...]
    free: tuple[FreeBoundary, ...]
    topology_status: Literal["UNVERIFIED_GEOMETRY"] = "UNVERIFIED_GEOMETRY"

    def __post_init__(self) -> None:
        if type(self.sign) is not int or self.sign not in (-1, 1):
            raise ValueError("invalid ring sign")
        if self.topology_status != "UNVERIFIED_GEOMETRY":
            raise ValueError("a sewing plan cannot self-certify 3D topology")
        self.audit()

    def audit(self) -> None:
        """Check every abstract port is used exactly once, and tree wiring."""
        address_map = {a.address: a for a in self.atoms}
        if len(address_map) != len(self.atoms) or () not in address_map:
            raise ValueError("duplicate or missing rooted atomic addresses")
        if any(a.sign != self.sign for a in self.atoms):
            raise ValueError("atom ring signs inconsistent")
        ports = {
            Port(a.address, j) for a in self.atoms
            for j in range(a.valence)
        }
        used: list[Port] = []
        for seam in self.seams:
            if not isinstance(seam, Seam):
                raise ValueError("invalid seam")
            used.extend((seam.parent, seam.child))
        for free in self.free:
            if not isinstance(free, FreeBoundary):
                raise ValueError("invalid free port")
            used.append(free.port)
        if set(used) != ports or len(used) != len(set(used)):
            raise ValueError("plumbing seam/free endpoints do not partition atomic ports")
        if len(self.seams) != len(self.atoms) - 1:
            raise ValueError("an atomic sewing tree requires V-1 seams")
        children = {s.child.address for s in self.seams}
        if children != set(address_map) - {()}:
            raise ValueError("all nonroot vertices must have one parent seam")
        if sum(f.marked_root for f in self.free) != 1:
            raise ValueError("exactly one marked root free port required")
        if not any(f.marked_root and f.port == Port((), 0) for f in self.free):
            raise ValueError("root free port must anchor atomic address () slot 0")

    def rooted_instruction_code(self) -> tuple:
        """Reconstruct full Abbreviated rooted code FROM sewing instructions.

        The procedure deliberately does not read the original tree and
        checks parity of child/free incidence via a second representation.
        """
        address_map = {a.address: a for a in self.atoms}
        free_map = {f.port: f for f in self.free}
        seam_map = {s.parent: s.child for s in self.seams}
        root_label = next(f.v4 for f in self.free if f.marked_root)

        def walk(address: tuple[int, ...]) -> tuple:
            atom = address_map[address]
            descendants = []
            for slot in range(1, atom.valence):
                port = Port(address, slot)
                if port in free_map:
                    descendants.append(("free", free_map[port].v4))
                else:
                    descendants.append(walk(seam_map[port].address))
            return (atom.ring_count, atom.sector_weights, tuple(descendants))

        return (self.sign, root_label, walk(()))

    def require_verified_isotopy(self) -> None:
        """Fail closed: an abstract sewing plan is not a 3D proof."""
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: must independently implement atomic "
            "Conway tangle strands, 3D plumbing maps, and an auditable "
            "boundary-parametrized isotopy certificate"
        )


def atomic_plumbing_plan(source: Abbreviated, root: str = "f0") -> AtomicPlumbingPlan:
    """Convert strict §12.5 data into root-addressed abstract atomic ports.

    Every angular-sector integer and ring number survives unchanged.
    Every internal bond becomes ONE unoriented seam between TWO ports.
    Every free bond becomes a V4-labelled abstract external boundary.
    Names/storage indices/cyclic starting slots do not enter the plan.

    The structural round-trip is checked; no topology is constructed.
    """
    if not isinstance(source, Abbreviated):
        raise ValueError("expected ring-numbered abbreviated presentation")
    if not isinstance(root, str) or not root.startswith("f"):
        raise ValueError("root must be an existing marked free bond")
    tree = source.labeled.tree
    owners: dict[str, list[int]] = {}
    for i, vertex in enumerate(tree.vertices):
        for bond in vertex.slots:
            owners.setdefault(bond, []).append(i)
    if len(owners.get(root, ())) != 1:
        raise ValueError("root free bond must occur exactly once")

    atoms: list[AtomicInstruction] = []
    seams: list[Seam] = []
    free: list[FreeBoundary] = []
    seen: set[int] = set()

    def visit(index: int, incoming: str, address: tuple[int, ...]) -> None:
        if index in seen:
            raise ValueError("non-tree or repeated index in sewing input")
        seen.add(index)
        vertex = tree.vertices[index]
        shift = vertex.slots.index(incoming)
        bonds = vertex.slots[shift:] + vertex.slots[:shift]
        sectors = vertex.sectors[shift:] + vertex.sectors[:shift]
        atoms.append(AtomicInstruction(address, sectors, source.rings[index], source.sign))
        if address == ():
            free.append(FreeBoundary(Port((), 0), source.labeled.label(root), True))
        for j, bond in enumerate(bonds[1:], 1):
            here = Port(address, j)
            if bond.startswith("f"):
                free.append(FreeBoundary(here, source.labeled.label(bond)))
            else:
                ends = owners[bond]
                if len(ends) != 2:
                    raise ValueError("every plumbing seam needs exactly two ends")
                child_idx = next(k for k in ends if k != index)
                child_addr = (*address, j)
                seams.append(Seam(here, Port(child_addr, 0)))
                visit(child_idx, bond, child_addr)

    visit(owners[root][0], root, ())
    if len(seen) != len(tree.vertices):
        raise ValueError("not all weighted vertices appear in the atomic plan")
    plan = AtomicPlumbingPlan(source.sign, tuple(atoms), tuple(seams), tuple(free))
    if plan.rooted_instruction_code() != source.rooted_code(root):
        raise AssertionError("source §12.5 presentation not losslessly encoded")
    return plan
