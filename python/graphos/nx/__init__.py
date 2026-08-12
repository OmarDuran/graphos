"""graphos.nx — NetworkX integration.

Level 2: read-only graph VIEWS over graphos complexes, duck-typing the
NetworkX protocol directly on the CSR arrays. Nodes are (dim, index)
tuples for the Hasse and incidence views (nfempy's convention), plain
cell indices for the adjacency view. Views are structural: computed
attributes only ("dim" on nodes, "sign" on Hasse arcs); mutation is not
supported — build a new complex instead (out-of-place epochs, see
THEORY.md P3).

Level 3: the "graphos" NetworkX backend. Views carry
__networkx_backend__ = "graphos"; once registered (automatically via the
packaging entry point, or explicitly with graphos.nx.register() when
running from a source tree), plain calls like

    nx.connected_components(adjacency_graph(c, 2, 1))

dispatch to graphos's native C++ implementations, and algorithms without
a native implementation fall back to NetworkX through to_networkx()
conversion (enable nx.config.fallback_to_nx).
"""

from collections.abc import Mapping

from .. import _core

__all__ = [
    "hasse_graph",
    "incidence_graph",
    "adjacency_graph",
    "HasseGraph",
    "IncidenceGraph",
    "AdjacencyGraph",
    "register",
]


class _Row(Mapping):
    """neighbor -> edge-attribute mapping for one CSR row"""

    def __init__(self, targets, attrs):
        self._targets = targets
        self._attrs = attrs

    def __getitem__(self, key):
        for t, a in zip(self._targets, self._attrs):
            if t == key:
                return a
        raise KeyError(key)

    def __iter__(self):
        return iter(self._targets)

    def __len__(self):
        return len(self._targets)


class _NodeView(Mapping):
    def __init__(self, nodes, attr_of):
        self._nodes = nodes
        self._attr_of = attr_of

    def __getitem__(self, key):
        if key not in self._nodes:
            raise KeyError(key)
        return self._attr_of(key)

    def __iter__(self):
        return iter(self._nodes)

    def __len__(self):
        return len(self._nodes)

    def __call__(self, data=False):
        if not data:
            return iter(self._nodes)
        return ((n, self._attr_of(n)) for n in self._nodes)


class _GraphBase:
    """the read-only slice of the NetworkX protocol the views share"""

    __networkx_backend__ = "graphos"

    graph = {}
    name = ""

    def is_multigraph(self):
        return False

    def __len__(self):
        return self.number_of_nodes()

    def __iter__(self):
        return iter(self.nodes)

    def __contains__(self, n):
        return n in self.nodes

    def __getitem__(self, n):
        return self.adj[n]

    def order(self):
        return self.number_of_nodes()

    def size(self):
        return self.number_of_edges()

    # networkx internals reach for the private mappings directly
    @property
    def _adj(self):
        return self.adj

    @property
    def _node(self):
        return self.nodes


class HasseGraph(_GraphBase):
    """Directed view of the Hasse diagram of a frozen complex: arcs run
    from each k-cell to its (k-1)-faces, carrying the orientation sign."""

    def __init__(self, frozen):
        self._f = frozen
        self._b = {k: frozen.boundary(k) for k in range(1, frozen.dim + 1)}
        self._c = {k: frozen.coboundary(k) for k in range(0, frozen.dim)}
        node_list = [
            (k, int(i)) for k in range(frozen.dim + 1) for i in range(frozen.count(k))
        ]
        self._node_set = set(node_list)
        self.nodes = _NodeView(node_list, lambda n: {"dim": n[0]})

    def is_directed(self):
        return True

    def number_of_nodes(self):
        return sum(self._f.count(k) for k in range(self._f.dim + 1))

    def number_of_edges(self):
        return sum(len(self._b[k][1]) for k in self._b)

    def successors(self, n):
        k, i = n
        if k == 0:
            return iter(())
        off, idx, _ = self._b[k]
        return ((k - 1, int(f)) for f in idx[off[i]:off[i + 1]])

    def predecessors(self, n):
        k, i = n
        if k == self._f.dim:
            return iter(())
        off, idx, _ = self._c[k]
        return ((k + 1, int(e)) for e in idx[off[i]:off[i + 1]])

    def _succ_row(self, n):
        k, i = n
        if k == 0:
            return _Row([], [])
        off, idx, sgn = self._b[k]
        targets = [(k - 1, int(f)) for f in idx[off[i]:off[i + 1]]]
        attrs = [{"sign": int(s)} for s in sgn[off[i]:off[i + 1]]]
        return _Row(targets, attrs)

    @property
    def adj(self):
        return _NodeView(list(self.nodes), self._succ_row)

    @property
    def pred(self):
        def row(n):
            k, i = n
            if k == self._f.dim:
                return _Row([], [])
            off, idx, sgn = self._c[k]
            targets = [(k + 1, int(e)) for e in idx[off[i]:off[i + 1]]]
            attrs = [{"sign": int(s)} for s in sgn[off[i]:off[i + 1]]]
            return _Row(targets, attrs)

        return _NodeView(list(self.nodes), row)

    @property
    def succ(self):
        return _NodeView(list(self.nodes), self._succ_row)

    _succ = succ
    _pred = pred

    def edges(self, data=False):
        for k, (off, idx, sgn) in self._b.items():
            for e in range(self._f.count(k)):
                for m in range(int(off[e]), int(off[e + 1])):
                    u, v = (k, e), (k - 1, int(idx[m]))
                    yield (u, v, {"sign": int(sgn[m])}) if data else (u, v)

    def to_networkx(self):
        import networkx as nx

        g = nx.DiGraph()
        g.add_nodes_from(self.nodes(data=True))
        g.add_edges_from(self.edges(data=True))
        return g


class IncidenceGraph(_GraphBase):
    """Directed bipartite view of the transitive incidence I(k, j):
    arcs from each k-cell to the j-cells of its closure/star — the
    nfempy build_graph(dim, codim) shape."""

    def __init__(self, complex_, k, j):
        self._c = complex_
        self._k, self._j = k, j
        self._off, self._idx = _core.incidence(complex_, k, j)
        node_list = [(k, int(i)) for i in range(complex_.count(k))] + [
            (j, int(i)) for i in range(complex_.count(j))
        ]
        self._node_set = set(node_list)
        self.nodes = _NodeView(node_list, lambda n: {"dim": n[0]})

    def is_directed(self):
        return True

    def number_of_nodes(self):
        return self._c.count(self._k) + self._c.count(self._j)

    def number_of_edges(self):
        return len(self._idx)

    def successors(self, n):
        k, i = n
        if k != self._k:
            return iter(())
        return ((self._j, int(f)) for f in self._idx[self._off[i]:self._off[i + 1]])

    def _row(self, n):
        k, i = n
        if k != self._k:
            return _Row([], [])
        targets = [(self._j, int(f)) for f in self._idx[self._off[i]:self._off[i + 1]]]
        return _Row(targets, [{} for _ in targets])

    @property
    def succ(self):
        return _NodeView(list(self.nodes), self._row)

    adj = succ

    def edges(self, data=False):
        for i in range(self._c.count(self._k)):
            for f in self._idx[self._off[i]:self._off[i + 1]]:
                u, v = (self._k, i), (self._j, int(f))
                yield (u, v, {}) if data else (u, v)

    def to_networkx(self):
        import networkx as nx

        g = nx.DiGraph()
        g.add_nodes_from(self.nodes(data=True))
        g.add_edges_from(self.edges())
        return g


class AdjacencyGraph(_GraphBase):
    """Undirected view of the generalized adjacency of k-cells through
    shared via-cells (via = k-1 is the mesh dual graph). Nodes are plain
    cell indices. Carries its construction data so backend dispatch can
    run graphos-native algorithms."""

    def __init__(self, complex_, k, via, exclude=None):
        self._c = complex_
        self._k, self._via, self._exclude = k, via, exclude
        self._off, self._idx = _core.adjacency(complex_, k, via)
        node_list = list(range(complex_.count(k)))
        self._node_set = set(node_list)
        self.nodes = _NodeView(node_list, lambda n: {"dim": k})

    def is_directed(self):
        return False

    def number_of_nodes(self):
        return self._c.count(self._k)

    def number_of_edges(self):
        return len(self._idx) // 2

    def neighbors(self, n):
        return (int(v) for v in self._idx[self._off[n]:self._off[n + 1]])

    def _row(self, n):
        targets = [int(v) for v in self._idx[self._off[n]:self._off[n + 1]]]
        return _Row(targets, [{} for _ in targets])

    @property
    def adj(self):
        return _NodeView(list(self.nodes), self._row)

    def degree(self, n=None):
        if n is not None:
            return int(self._off[n + 1] - self._off[n])
        return ((u, int(self._off[u + 1] - self._off[u])) for u in self.nodes)

    def edges(self, data=False):
        for u in range(self._c.count(self._k)):
            for v in self._idx[self._off[u]:self._off[u + 1]]:
                if u < int(v):
                    yield (u, int(v), {}) if data else (u, int(v))

    def to_networkx(self):
        import networkx as nx

        g = nx.Graph()
        g.add_nodes_from(self.nodes)
        g.add_edges_from(self.edges())
        return g


def hasse_graph(frozen):
    """Directed Hasse-diagram view of a FrozenComplex (zero-copy CSR)."""
    return HasseGraph(frozen)


def incidence_graph(complex_, k, j):
    """Bipartite transitive-incidence view I(k, j) of a Complex."""
    return IncidenceGraph(complex_, k, j)


def adjacency_graph(complex_, k, via, exclude=None):
    """Generalized-adjacency view of the k-cells of a Complex."""
    return AdjacencyGraph(complex_, k, via, exclude)


def register():
    """Register the 'graphos' backend with NetworkX at runtime (needed
    when running from a source tree; installed packages register through
    the packaging entry point)."""
    from importlib.metadata import EntryPoint

    from networkx.utils import backends as _nxb

    ep = EntryPoint(
        name="graphos", value="graphos.nx.backend:GraphosBackend", group="networkx.backends"
    )
    _nxb.backends["graphos"] = ep
    if hasattr(_nxb, "backend_info"):
        _nxb.backend_info.setdefault("graphos", {})
