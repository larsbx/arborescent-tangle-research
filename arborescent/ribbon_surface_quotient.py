"""A5.2e2b1: exact finite 2-complex gluing of already audited PL ribbon disks.

A5.2e2a constructs n embedded ribbon *patches*, not the one global
source twisted band in Bonahon--Siebenmann Fig.12.10. This module
makes a missing attachment choice an explicit finite surface object.
It computes Euler characteristic, boundary circuits, orientability,
connected components and vertex-link manifold checks independently
from the proposed cap identifications.

CRITICAL: cap identifications are ABSTRACT quotient operations, and
DO NOT construct a geometric identification in Q^3, a ribbon framing,
the source's band, or any Conway 3D sewing homeomorphism. Identical
puncture sets may support distinct surface topologies. All results
are SOURCE_UNVERIFIED; R+/- cannot be certified here.
"""
from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from typing import Literal

from .pl_ribbon_patches import TriangulatedRibbon


@dataclass(frozen=True)
class RibbonCap:
    patch: int
    side: Literal["start", "end"]

    def __post_init__(self) -> None:
        if type(self.patch) is not int or self.patch < 0 or self.side not in ("start","end"):
            raise ValueError("cap needs patch index and start/end side")


@dataclass(frozen=True)
class CapAttachment:
    first: RibbonCap
    second: RibbonCap
    reverse_width: bool = False
    status: Literal["ABSTRACT_CAP_PAIRING"] = "ABSTRACT_CAP_PAIRING"

    def __post_init__(self) -> None:
        if (not isinstance(self.first,RibbonCap)
                or not isinstance(self.second,RibbonCap)
                or self.first == self.second
                or type(self.reverse_width) is not bool
                or self.status != "ABSTRACT_CAP_PAIRING"):
            raise ValueError("invalid abstract ribbon cap identification")


class _DSU:
    def __init__(self, nodes):
        self.parent={a:a for a in nodes}

    def find(self, a):
        root=a
        while root!=self.parent[root]:
            root=self.parent[root]
        while a!=root:
            nxt=self.parent[a]
            self.parent[a]=root
            a=nxt
        return root

    def join(self, a, b):
        x,y=self.find(a),self.find(b)
        if x!=y:
            self.parent[y]=x


def _cap(ribbon: TriangulatedRibbon, side: str) -> tuple[int,int,int]:
    if side not in ("start","end"):
        raise ValueError("unknown cap")
    j=0 if side=="start" else len(ribbon.first)-1
    return (3*j,3*j+1,3*j+2)


def _components(adjacency: dict, nodes: set) -> tuple[frozenset, ...]:
    remaining=set(nodes)
    components=[]
    while remaining:
        start=min(remaining)
        queue=[start]
        remaining.remove(start)
        component={start}
        while queue:
            for dest in adjacency.get(queue.pop(),()):
                if dest in remaining:
                    remaining.remove(dest)
                    queue.append(dest)
                    component.add(dest)
        components.append(frozenset(component))
    return tuple(sorted(components,key=lambda x:min(x)))


@dataclass(frozen=True)
class SurfaceTopology:
    """Computed result for an abstract 2-manifold with boundary."""
    components: int
    euler: int
    boundary_components: int
    orientable: bool
    vertex_classes: int
    edge_classes: int
    face_count: int


@dataclass(frozen=True)
class AbstractRibbonSurface:
    patches: tuple[TriangulatedRibbon, ...]
    attachments: tuple[CapAttachment, ...]
    status: Literal["ABSTRACT_SURFACE_SOURCE_UNVERIFIED"] = "ABSTRACT_SURFACE_SOURCE_UNVERIFIED"

    def __post_init__(self) -> None:
        if (not isinstance(self.patches,tuple) or not self.patches
                or any(not isinstance(p,TriangulatedRibbon) for p in self.patches)
                or not isinstance(self.attachments,tuple)
                or any(not isinstance(a,CapAttachment) for a in self.attachments)
                or self.status != "ABSTRACT_SURFACE_SOURCE_UNVERIFIED"):
            raise ValueError("only typed PL patches, abstract attachments supported")
        self.analyze()

    def analyze(self) -> SurfaceTopology:
        """Independently compute quotient topology, checking vertex links.

        Each mesh face remains a distinct 2-simplex. Cap edges become
        paired *abstract* edges with an explicit endpoint map. Both
        edge incidence and corner-link connectivity are checked, to
        reject pinched vertices or triple/nonmanifold edge gluings.
        """
        for ribbon in self.patches:
            ribbon.audit()
        ids={(p,i) for p,patch in enumerate(self.patches)
             for i in range(len(patch.vertices))}
        vertex=_DSU(ids)
        edge_occ: dict[tuple[int,int,int],list[tuple[int,tuple[int,int]]]]=defaultdict(list)
        faces=[]
        for p,patch in enumerate(self.patches):
            for tri in patch.faces:
                fid=len(faces)
                faces.append((p,tri))
                for i,j in zip(tri,tri[1:]+tri[:1]):
                    edge_occ[(p,min(i,j),max(i,j))].append((fid,(i,j)))
        if any(len(v) not in (1,2) for v in edge_occ.values()):
            raise ValueError("invalid triangle-edge incidence")
        edge=_DSU(edge_occ)
        cap_seen=set()
        # For an abstract paired edge, explicitly record the bijection
        # of its ORIGINAL vertex ends (needed for Möbius twists).
        seam_map={}
        for attach in self.attachments:
            c1,c2=attach.first,attach.second
            if (c1.patch>=len(self.patches) or c2.patch>=len(self.patches)
                    or c1 in cap_seen or c2 in cap_seen):
                raise ValueError("nonexistent or reused ribbon cap")
            cap_seen.update((c1,c2))
            v1=_cap(self.patches[c1.patch],c1.side)
            v2=_cap(self.patches[c2.patch],c2.side)
            if attach.reverse_width:
                v2=v2[::-1]
            for i,j in zip(v1,v2):
                vertex.join((c1.patch,i),(c2.patch,j))
            for k in range(2):
                s,t=v1[k:k+2],v2[k:k+2]
                e1=(c1.patch,min(s),max(s))
                e2=(c2.patch,min(t),max(t))
                if (len(edge_occ[e1])!=1 or len(edge_occ[e2])!=1
                        or e1==e2 or e1 in seam_map or e2 in seam_map):
                    raise ValueError("only unused cap boundary edges can be paired")
                seam_map[e1]=(e2,{s[0]:t[0],s[1]:t[1]})
                seam_map[e2]=(e1,{t[0]:s[0],t[1]:s[1]})
                edge.join(e1,e2)

        classes: dict[tuple[int,int,int],list[tuple[int,tuple[int,int],tuple[int,int,int]]]]=defaultdict(list)
        for e,occurrences in edge_occ.items():
            for fid,direction in occurrences:
                classes[edge.find(e)].append((fid,direction,e))
        if any(len(group) not in (1,2) for group in classes.values()):
            raise ValueError("edge has nonmanifold face incidence")

        # Face orientation constraints and corner-link adjacency.
        face_adj=defaultdict(list)
        corner_adj=defaultdict(list)
        for group in classes.values():
            if len(group)==1:
                continue
            (f,d,e),(g,h,k)=group
            pf=faces[f][0]
            pg=faces[g][0]
            if e==k:
                # Two original triangles on one physical patch mesh edge.
                matching={e[1]:e[1],e[2]:e[2]}
            elif k==seam_map.get(e,(None,None))[0]:
                matching=seam_map[e][1]
            else:
                raise ValueError("unexplained edge identification")
            if any((vertex.find((pf,v))!=vertex.find((pg,matching[v])))
                   for v in d):
                raise ValueError("identified edge endpoints disagree")
            align=1 if (matching[d[0]],matching[d[1]])==h else -1
            if align==-1 and (matching[d[0]],matching[d[1]])!=h[::-1]:
                raise ValueError("edge directions incompatible")
            face_adj[f].append((g,-align))
            face_adj[g].append((f,-align))
            for local in d:
                other=matching[local]
                a=(f,(pf,local))
                b=(g,(pg,other))
                corner_adj[a].append(b)
                corner_adj[b].append(a)

        # At each abstract vertex the local link must be ONE path
        # (boundary point) or ONE cycle (interior point).
        vertex_nodes: dict[tuple[int,int],set]=defaultdict(set)
        for fid,(p,tri) in enumerate(faces):
            if len({vertex.find((p,i)) for i in tri})!=3:
                raise ValueError("collapsed triangle under vertex identification")
            for i in tri:
                vertex_nodes[vertex.find((p,i))].add((fid,(p,i)))
        for v,nodes in vertex_nodes.items():
            graph={n:corner_adj.get(n,[]) for n in nodes}
            if any(len(neighbors)>2 or any(other not in nodes for other in neighbors)
                   for neighbors in graph.values()):
                raise ValueError("nonmanifold corner incidence")
            if len(_components(graph,nodes))!=1:
                raise ValueError("pinched or disconnected vertex link")
            degrees=sorted(len(graph[n]) for n in nodes)
            if not (all(d==2 for d in degrees)
                    or degrees==[0]
                    or degrees.count(1)==2 and degrees.count(0)==0
                    and all(d in (1,2) for d in degrees)):
                raise ValueError("vertex link is neither circle nor interval")

        # Boundary is a degree-two multigraph, including loop edges.
        boundary_graph=defaultdict(list)
        boundary_vertices=set()
        for group in classes.values():
            if len(group)!=1:
                continue
            _,(i,j),edge_key=group[0]
            p=edge_key[0]
            a,b=vertex.find((p,i)),vertex.find((p,j))
            boundary_graph[a].append(b)
            boundary_graph[b].append(a)
            boundary_vertices.update((a,b))
        if any(len(nbr)!=2 for nbr in boundary_graph.values()):
            raise ValueError("branching or dangling boundary graph")
        boundary_count=len(_components(boundary_graph,boundary_vertices))

        # Whole surface connected components and signed face orientation.
        face_nodes=set(range(len(faces)))
        face_graph={i:[other for other,_ in face_adj[i]] for i in face_nodes}
        component_count=len(_components(face_graph,face_nodes))
        orientation={}
        orientable=True
        for base in sorted(face_nodes):
            if base in orientation:
                continue
            orientation[base]=1
            queue=[base]
            while queue:
                current=queue.pop()
                for neighbor,relative in face_adj[current]:
                    sign=orientation[current]*relative
                    if neighbor in orientation:
                        if orientation[neighbor]!=sign:
                            orientable=False
                    else:
                        orientation[neighbor]=sign
                        queue.append(neighbor)
        return SurfaceTopology(
            component_count,len({vertex.find(i) for i in ids})
            -len(classes)+len(faces),boundary_count,orientable,
            len({vertex.find(i) for i in ids}),len(classes),len(faces))

    def require_source_band_certificate(self) -> None:
        raise NotImplementedError(
            "UNVERIFIED_GEOMETRY: abstract cap gluings are not an "
            "embedded source Fig.12.10 twisted band, framing or PL isotopy"
        )
