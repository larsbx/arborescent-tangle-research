"""Weighted planar-tree *presentations* and a finite 2/3-branch grammar.

W/F/S checks below are strict nonexceptional necessary conditions from
Bonahon--Siebenmann sec. 12.4. They are NOT a knot-isotopy classifier.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import factorial
from typing import TypeAlias


@dataclass(frozen=True)
class Vertex:
    """Cyclic bond slots and associated angular-sector integer weights.

    Internal bond 'eN' occurs at two distinct vertices; free bond 'fN'
    occurs at one. sectors[i] is just after slots[i] cyclically.
    """
    slots: tuple[str, ...]
    sectors: tuple[int, ...]

    @property
    def weight(self) -> int:
        return sum(self.sectors)


@dataclass(frozen=True)
class WeightedPlaneTree:
    """Planar rotation system, sector data, and free legs; no topology quotient."""
    vertices: tuple[Vertex, ...]

    def __post_init__(self) -> None:
        if not self.vertices:
            raise ValueError("at least one vertex required")
        occ: dict[str, list[int]] = {}
        for i, v in enumerate(self.vertices):
            if not isinstance(v, Vertex) or not v.slots:
                raise ValueError("vertices require incident bonds")
            if len(v.slots) != len(v.sectors) or len(set(v.slots)) != len(v.slots):
                raise ValueError("distinct slots and matching sectors required")
            if any(type(w) is not int for w in v.sectors):
                raise ValueError("weights must be integers")
            for name in v.slots:
                if (not isinstance(name, str) or len(name) < 2
                        or name[0] not in 'ef' or not name[1:].isdigit()):
                    raise ValueError("bond names must be eN or fN")
                occ.setdefault(name, []).append(i)
        edges = []
        for name, owners in occ.items():
            if name[0] == 'e':
                if len(owners) != 2 or owners[0] == owners[1]:
                    raise ValueError("an edge joins two distinct vertices")
                edges.append(tuple(owners))
            elif len(owners) != 1:
                raise ValueError("a free bond has exactly one endpoint")
        if len(edges) != len(self.vertices) - 1:
            raise ValueError("not a tree: incorrect edge count")
        neighbors = [set() for _ in self.vertices]
        for a, b in edges:
            neighbors[a].add(b)
            neighbors[b].add(a)
        seen, stack = set(), [0]
        while stack:
            v = stack.pop()
            if v not in seen:
                seen.add(v)
                stack.extend(neighbors[v] - seen)
        if len(seen) != len(self.vertices):
            raise ValueError("not a tree: disconnected")

    def edges(self) -> tuple[tuple[int, int], ...]:
        owners: dict[str, list[int]] = {}
        for i, vertex in enumerate(self.vertices):
            for slot in vertex.slots:
                if slot.startswith('e'):
                    owners.setdefault(slot, []).append(i)
        return tuple(tuple(ends) for _, ends in sorted(owners.items()))

    def rooted_plane_code(self, marked_free_bond: str = 'f0') -> tuple:
        """Canonical isomorphism code of *marked presentations*, not isotopy.

        Invariant under vertex/edge renaming and cyclic starting-slot
        choices, but sensitive to child order and sector placement.
        """
        if not marked_free_bond.startswith('f'):
            raise ValueError("root must be a free bond")
        owners: dict[str, list[int]] = {}
        for i, v in enumerate(self.vertices):
            for slot in v.slots:
                owners.setdefault(slot, []).append(i)
        if marked_free_bond not in owners or len(owners[marked_free_bond]) != 1:
            raise ValueError("root must be a present free bond")

        def code(index: int, incoming: str) -> tuple:
            vertex = self.vertices[index]
            shift = vertex.slots.index(incoming)
            slots = vertex.slots[shift:] + vertex.slots[:shift]
            weights = vertex.sectors[shift:] + vertex.sectors[:shift]
            descendants = []
            for slot in slots[1:]:
                if slot.startswith('f'):
                    descendants.append(('free',))
                else:
                    other = next(i for i in owners[slot] if i != index)
                    descendants.append(code(other, slot))
            return (weights, tuple(descendants))

        return code(owners[marked_free_bond][0], marked_free_bond)

    def wfs_failures(self) -> tuple[tuple[str, int], ...]:
        """Strict W, F, S screen; excludes the manuscript's small exceptions.

        W: at most one nonzero angular sector at each vertex.
        F: free bonds only at zero-weight vertices of bond-valence >= 3.
        S: delete valence >= 3 vertices; every remaining stick has
           nonzero alternating vertex weights and endpoints |weight|>=2.
        Passing does NOT establish positive/negative canonicity.
        """
        bad: list[tuple[str, int]] = []
        deg = [len(v.slots) for v in self.vertices]
        for i, v in enumerate(self.vertices):
            if sum(w != 0 for w in v.sectors) > 1:
                bad.append(('W', i))
            if any(s.startswith('f') for s in v.slots):
                if deg[i] < 3 or v.weight != 0:
                    bad.append(('F', i))
        neighbors = [set() for _ in self.vertices]
        for a, b in self.edges():
            neighbors[a].add(b)
            neighbors[b].add(a)
        permitted = {i for i, d in enumerate(deg) if d <= 2}
        unseen = set(permitted)
        while unseen:
            start = min(unseen)
            component, stack = set(), [start]
            while stack:
                i = stack.pop()
                if i not in unseen:
                    continue
                unseen.remove(i)
                component.add(i)
                stack.extend((neighbors[i] & permitted) & unseen)
            ends = sorted(i for i in component if len(neighbors[i] & component) <= 1)
            route, parent, current = [], -1, ends[0]
            while True:
                route.append(current)
                nxt = sorted((neighbors[current] & component) - {parent})
                if not nxt:
                    break
                parent, current = current, nxt[0]
            weights = [self.vertices[i].weight for i in route]
            if any(w == 0 for w in weights):
                bad.append(('S-zero', start))
            if any(a * b >= 0 for a, b in zip(weights, weights[1:])):
                bad.append(('S-sign', start))
            if abs(weights[0]) < 2 or abs(weights[-1]) < 2:
                bad.append(('S-end', start))
        return tuple(sorted(set(bad)))


@dataclass(frozen=True)
class Leaf:
    weight: int


@dataclass(frozen=True)
class Fork:
    children: tuple['Spec', ...]


Spec: TypeAlias = Leaf | Fork


def to_weighted_tree(spec: Spec, labels: tuple[int, ...] = (2, 3)) -> WeightedPlaneTree:
    """Root free port f0; 2/3-child zero-weight forks and weighted leaves."""
    if not isinstance(spec, Fork):
        raise ValueError("root must branch")
    if (not labels or len(set(labels)) != len(labels)
            or any(type(w) is not int or abs(w) < 2 for w in labels)):
        raise ValueError("labels must be distinct integers of |weight|>=2")
    vertices: list[Vertex | None] = []
    count = 0

    def visit(node: Spec, incoming: str) -> None:
        nonlocal count
        index = len(vertices)
        vertices.append(None)
        if isinstance(node, Leaf):
            if node.weight not in labels:
                raise ValueError("leaf weight outside the alphabet")
            vertices[index] = Vertex((incoming,), (node.weight,))
        elif isinstance(node, Fork):
            if len(node.children) not in (2, 3):
                raise ValueError("only 2/3-ary forks")
            outgoing = []
            for child in node.children:
                bond = f'e{count}'
                count += 1
                outgoing.append(bond)
                visit(child, bond)
            slots = (incoming, *outgoing)
            vertices[index] = Vertex(slots, (0,) * len(slots))
        else:
            raise ValueError("unknown syntax node")
    visit(spec, 'f0')
    return WeightedPlaneTree(tuple(v for v in vertices if v is not None))


def _splits(total: int, arity: int):
    if arity == 1:
        yield (total,)
    else:
        for head in range(total + 1):
            for tail in _splits(total - head, arity - 1):
                yield (head, *tail)


@lru_cache(maxsize=None)
def enumerate_specs(n: int, labels: tuple[int, ...] = (2, 3)) -> tuple[Spec, ...]:
    """Exact exhaustive syntax enumeration for small n (exponential work)."""
    if type(n) is not int or n < 0:
        raise ValueError("nonnegative integer size required")
    if not labels or len(set(labels)) != len(labels) or any(type(w) is not int for w in labels):
        raise ValueError("labels must be distinct integers")
    if n == 0:
        return tuple(Leaf(w) for w in labels)
    out = []
    for arity in (2, 3):
        for sizes in _splits(n - 1, arity):
            choices = [enumerate_specs(s, labels) for s in sizes]
            out.extend(Fork(tuple(children)) for children in product(*choices))
    return tuple(out)


def syntax_counts(max_n: int, m: int) -> tuple[int, ...]:
    """Count by internal vertices via F(x)=m+x(F(x)^2+F(x)^3)."""
    if type(max_n) is not int or max_n < 0 or type(m) is not int or m < 1:
        raise ValueError("max_n>=0 and m>=1 integers required")
    a = [m]
    for n in range(1, max_n + 1):
        pair = sum(a[i] * a[n - 1 - i] for i in range(n))
        triple = sum(a[i] * a[j] * a[n - 1 - i - j]
                     for i in range(n) for j in range(n - i))
        a.append(pair + triple)
    return tuple(a)


def count_by_arity(n: int, ternary: int, m: int) -> int:
    """Exact cycle-lemma count with n forks, ternary ternary forks, m leaves."""
    if (any(type(v) is not int for v in (n, ternary, m))
            or not (n >= 0 and 0 <= ternary <= n and m >= 1)):
        raise ValueError("invalid size or alphabet")
    leaves = n + ternary + 1
    return (factorial(2 * n + ternary) //
            (factorial(n - ternary) * factorial(ternary) * factorial(leaves))
            * m ** leaves)


def syntax_count_closed(n: int, m: int) -> int:
    return sum(count_by_arity(n, j, m) for j in range(n + 1))
