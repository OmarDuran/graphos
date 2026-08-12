"""Level 2/3 NetworkX integration tests: view protocol over the CSR
arrays, the nfempy-style predecessors/successors pattern, backend
dispatch of nx algorithms to the native C++ engine, and conversion
fallback for algorithms without a native implementation."""

import sys

import graphos as g
import graphos.nx as gnx

try:
    import networkx as nx
except ImportError:  # pragma: no cover
    print("SKIP: networkx not installed")
    sys.exit(0)


def make_fan():
    c = g.Complex(dim=2)
    c.attach_vertices(5)
    for i, (a, b) in enumerate([(0, 4), (1, 4), (2, 4), (3, 4), (0, 1), (1, 2), (2, 3), (3, 0)]):
        c.attach_cell(1, [a, b], [-1, +1])
    c.attach_cell(2, [4, 1, 0], [+1, +1, -1])
    c.attach_cell(2, [5, 2, 1], [+1, +1, -1])
    c.attach_cell(2, [6, 3, 2], [+1, +1, -1])
    c.attach_cell(2, [7, 0, 3], [+1, +1, -1])
    return c


def make_disk():
    return g.from_polygons(4, [[0, 1, 2], [1, 0, 3]])


def test_hasse_view_protocol():
    f = g.freeze(make_fan())
    h = gnx.hasse_graph(f)
    assert h.is_directed() and not h.is_multigraph()
    assert len(h) == 5 + 8 + 4
    assert h.number_of_edges() == 8 * 2 + 4 * 3
    assert (2, 0) in h and (3, 7) not in h
    # nfempy migration pattern: successors/predecessors with (dim, id) nodes
    assert sorted(h.successors((2, 0))) == [(1, 0), (1, 1), (1, 4)]
    assert sorted(h.predecessors((1, 0))) == [(2, 0), (2, 3)]  # spoke s0
    # signs ride as edge attributes
    assert h.adj[(2, 0)][(1, 0)] == {"sign": -1}
    assert h.nodes[(1, 4)] == {"dim": 1}


def test_incidence_view_is_build_graph():
    c = make_fan()
    iv = gnx.incidence_graph(c, 2, 0)  # nfempy build_graph(dim, codim=2)
    assert sorted(iv.successors((2, 0))) == [(0, 0), (0, 1), (0, 4)]
    assert iv.number_of_edges() == 4 * 3


def test_adjacency_view():
    a = gnx.adjacency_graph(make_fan(), 2, 1)  # the dual graph
    assert not a.is_directed()
    assert len(a) == 4
    assert sorted(a.neighbors(0)) == [1, 3]
    assert a.degree(0) == 2
    assert a.number_of_edges() == 4  # the 4-cycle of triangles


def test_to_networkx_materialization():
    f = g.freeze(make_fan())
    h = gnx.hasse_graph(f).to_networkx()
    assert isinstance(h, nx.DiGraph)
    assert h.number_of_nodes() == 17 and h.number_of_edges() == 28
    assert h[(2, 0)][(1, 0)]["sign"] == -1
    a = gnx.adjacency_graph(make_fan(), 2, 1).to_networkx()
    assert nx.is_connected(a)


def test_backend_dispatch_runs_native():
    gnx.register()
    disk = make_disk()
    interface = g.Marker(disk)
    interface.mark_where(1, lambda e: e == 0)
    cut = g.cut_along(disk, interface)

    # plain networkx call on a graphos view: dispatches to the C++ engine
    a = gnx.adjacency_graph(cut.complex, 2, 1)
    comps = list(nx.connected_components(a))
    assert len(comps) == 2
    assert sorted(map(sorted, comps)) == [[0], [1]]
    assert nx.number_connected_components(a) == 2
    assert not nx.is_connected(a)

    # excluded connectors: the sides-of-a-cut computation through nx syntax
    b = gnx.adjacency_graph(disk, 2, 1, exclude=interface)
    assert nx.number_connected_components(b) == 2


def test_backend_fallback_converts():
    gnx.register()
    a = gnx.adjacency_graph(make_fan(), 2, 1)
    old = nx.config.fallback_to_nx
    nx.config.fallback_to_nx = True
    try:
        # no native implementation: falls back through convert_to_nx
        assert nx.diameter(a) == 2
        tree = nx.bfs_tree(a, 0)
        assert tree.number_of_nodes() == 4
    finally:
        nx.config.fallback_to_nx = old


if __name__ == "__main__":
    import traceback

    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"[PASS] {name}")
            except Exception:
                failures += 1
                print(f"[FAIL] {name}")
                traceback.print_exc()
    sys.exit(1 if failures else 0)
