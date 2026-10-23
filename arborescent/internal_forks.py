"""A2.4 binary internal-fork grammar; presentation rewrites, not isotopy."""
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from .weighted_trees import Vertex, WeightedPlaneTree
from .tangle_labels import label_free_ports
from .partial_moves import explore

@dataclass(frozen=True)
class Tip:
    weight: int

@dataclass(frozen=True)
class Port:
    pass

@dataclass(frozen=True)
class InternalFork:
    children: tuple
    weight: int = 0
    sector: int = 0


def compile_fork(spec):
    """Compile ordered binary forks, weighted tips and free-port terminals."""
    if not isinstance(spec, InternalFork):
        raise ValueError('root must fork')
    vertices = []
    edge, free = 0, 1
    def visit(node, incoming):
        nonlocal edge, free
        if isinstance(node, Tip):
            if type(node.weight) is not int or abs(node.weight) < 2:
                raise ValueError('tips require integer absolute weight >=2')
            vertices.append(Vertex((incoming,), (node.weight,)))
            return
        if not isinstance(node, InternalFork) or len(node.children) != 2:
            raise ValueError('binary fork required')
        if (type(node.weight) is not int or type(node.sector) is not int
                or node.sector not in range(3) or (node.weight == 0 and node.sector != 0)):
            raise ValueError('invalid decoration')
        index = len(vertices)
        vertices.append(None)
        slots = [incoming]
        for child in node.children:
            if isinstance(child, Port):
                slots.append(f'f{free}')
                free += 1
            else:
                bond = f'e{edge}'
                edge += 1
                slots.append(bond)
                visit(child, bond)
        sectors = [0]*3
        sectors[node.sector] = node.weight
        vertices[index] = Vertex(tuple(slots), tuple(sectors))
    visit(spec, 'f0')
    return label_free_ports(WeightedPlaneTree(tuple(vertices)))


def generate(forks, tips=(2, 3), weights=(0, 3), extra_ports=1):
    """F -> (T|port|F)(T|port|F), with one sector decoration per F.

    Exact fork and extra-port counts; identity V4 labels. Zero decorations
    have one encoding. Size counters bound generation, not local productions.
    """
    if (type(forks) is not int or forks < 1 or type(extra_ports) is not int
            or not 0 <= extra_ports <= forks+1):
        raise ValueError('invalid bounds')
    for alphabet, tip in ((tips, True), (weights, False)):
        if (not isinstance(alphabet, tuple) or not alphabet or len(set(alphabet)) != len(alphabet)
                or any(type(w) is not int or (tip and abs(w)<2) for w in alphabet)):
            raise ValueError('invalid alphabet')
    decorations = [(w,s) for w in weights for s in (range(3) if w else (0,))]
    @lru_cache(None)
    def sub(n,p):
        if n == 0:
            return (Port(),) if p == 1 else tuple(Tip(w) for w in tips) if p == 0 else ()
        if not 0 <= p <= n+1:
            return ()
        return tuple(InternalFork((a,b),w,s)
                     for left in range(n) for lp in range(p+1)
                     for a,b in product(sub(left,lp), sub(n-1-left,p-lp))
                     for w,s in decorations)
    return tuple(compile_fork(s) for s in sub(forks,extra_ports))


@dataclass(frozen=True)
class Screen:
    failures: tuple
    unresolved: tuple = ()
    @property
    def passed(self):
        return not self.failures and not self.unresolved


def tangle_screen(source, sign='+'):
    """§12.7.2 W/S/P or N on branched trees; unbranched exceptions unknown.

    P/N excludes singleton -2/+2 sticks without free bonds (p.153
    diagrams). Isolated ±1 whole-tree exceptions remain outside this screen.
    F is absent. Passing is structural, not an isotopy certificate.
    """
    if sign not in ('+','-'):
        raise ValueError('invalid sign')
    tree = source.tree
    bad = [('W',i) for i,v in enumerate(tree.vertices) if sum(w != 0 for w in v.sectors)>1]
    if all(len(v.slots)<=2 for v in tree.vertices):
        return Screen(tuple(bad), ('unbranched whole-tree exceptions',))
    adj = [set() for _ in tree.vertices]
    for a,b in tree.edges():
        adj[a].add(b)
        adj[b].add(a)
    unseen = {i for i,v in enumerate(tree.vertices) if len(v.slots)<=2}
    while unseen:
        component, stack = set(), [min(unseen)]
        while stack:
            i = stack.pop()
            if i not in unseen:
                continue
            unseen.remove(i)
            component.add(i)
            stack.extend(adj[i] & unseen)
        ends = sorted(i for i in component if len(adj[i]&component)<=1)
        route, prev, cur = [], -1, ends[0]
        while True:
            route.append(cur)
            nxt = sorted((adj[cur]&component)-{prev})
            if not nxt:
                break
            prev,cur = cur,nxt[0]
        free = {i for i in component if any(s.startswith('f') for s in tree.vertices[i].slots)}
        for i in route:
            w = tree.vertices[i].weight
            if w == 0 and not (i in ends and i in free):
                bad.append(('S-zero',i))
            if i in ends and abs(w)==1 and i not in free:
                bad.append(('S-end',i))
        nz = [tree.vertices[i].weight for i in route if tree.vertices[i].weight]
        if any(a*b >= 0 for a,b in zip(nz,nz[1:])):
            bad.append(('S-sign',route[0]))
        if len(component)==1 and not free and tree.vertices[route[0]].weight == (-2 if sign=='+' else 2):
            bad.append(('P' if sign=='+' else 'N',route[0]))
    return Screen(tuple(sorted(bad)))


def elementary_states(source):
    """All indexed F1/F2 states, including both pivot parities."""
    from .tangle_labels import tangle_flip_f1 as f1, tangle_flip_f2 as f2
    neighbor = next((b if a==0 else a for a,b in source.tree.edges() if 0 in (a,b)),0)
    even,odd = f2(source,0),f2(source,neighbor)
    base = (source,even,odd,f2(even,neighbor))
    return set(base) | {f1(x) for x in base}


def compare_partial_orbits(seeds, max_states=32768):
    """Count seed intersections of complete indexed orbits; replay new merges.

    Paths may leave grammar slices; no isomorphism or isotopy quotient.
    Incomplete searches raise, never produce a misleading count.
    """
    pending = set(seeds)
    todo = set(pending)
    elementary_count = partial_count = 0
    traces = []
    while todo:
        s = min(todo,key=repr)
        todo -= elementary_states(s)
        elementary_count += 1
    while pending:
        s = min(pending,key=repr)
        orbit = explore(s,max_states)
        if not orbit.complete:
            raise ValueError('incomplete orbit')
        reached = pending & orbit.predecessors.keys()
        traces.extend(orbit.trace(t) for t in sorted(reached-elementary_states(s),key=repr))
        pending -= reached
        partial_count += 1
    return elementary_count,partial_count,tuple(traces)
