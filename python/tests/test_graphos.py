"""graphos Python binding tests: the same mathematical certificates the
C++ suites assert, re-asserted through the binding layer — plus the
binding-specific contracts (zero-copy frozen views, predicate markers,
exception mapping). Runnable directly (python test_graphos.py) or under
pytest."""

import numpy as np

import graphos as g


def make_triangle():
    c = g.Complex(dim=2)
    c.attach_vertices(3)
    c.attach_cell(1, [0, 1], [-1, +1])
    c.attach_cell(1, [1, 2], [-1, +1])
    c.attach_cell(1, [2, 0], [-1, +1])
    c.attach_cell(2, [0, 1, 2], [+1, +1, +1])
    return c


def make_circle():
    c = g.Complex(dim=1)
    c.attach_vertices(3)
    c.attach_cell(1, [0, 1], [-1, +1])
    c.attach_cell(1, [1, 2], [-1, +1])
    c.attach_cell(1, [2, 0], [-1, +1])
    return c


def make_disk():
    # two triangles sharing edge (0,1), opposite windings: consistent
    return g.from_polygons(4, [[0, 1, 2], [1, 0, 3]])


def test_build_and_invariants():
    c = make_triangle()
    c.validate()
    assert c.counts() == [3, 3, 1]
    assert g.d_squared_is_zero(c)
    assert g.euler_characteristic(c) == 1
    tet = g.from_simplices(3, 4, [[0, 1, 2, 3]])
    assert tet.counts() == [4, 6, 4, 1]
    assert g.d_squared_is_zero(tet)


def test_homology_of_constructed_spaces():
    torus = g.product(make_circle(), make_circle()).complex
    assert g.betti_numbers_z2(torus) == [1, 2, 1]
    assert g.euler_characteristic(torus) == 0
    sphere = g.join(g.Complex(dim=0), make_circle())  # needs two points
    two_points = g.Complex(dim=0)
    two_points.attach_vertices(2)
    s2 = g.join(two_points, make_circle()).complex
    assert g.betti_numbers_z2(s2) == [1, 0, 1]
    # relative homology: H(D², S¹) = {0, 0, 1}
    disk = make_disk()
    assert g.betti_numbers_z2(disk, g.classify_facets(disk).boundary) == [0, 0, 1]


def test_marker_predicates_from_python():
    c = make_disk()
    interface = g.Marker(c)
    interface.mark_where(1, lambda e: e == 0)  # a Python callable as predicate
    assert interface.marked_count(1) == 1
    cut = g.cut_along(c, interface)
    cut.complex.validate()
    assert g.d_squared_is_zero(cut.complex)
    assert g.euler_characteristic(cut.complex) == 3  # two disks + a segment
    assert g.connected_components(cut.complex).count == 3
    # ancestor map: the copies descend from the interface edge
    anc = cut.ancestor.index(1)
    assert int(anc[5]) == 0 and int(anc[6]) == 0


def test_frozen_views_are_zero_copy():
    f = g.freeze(make_triangle())
    offsets, indices, signs = f.boundary(2)
    assert offsets.tolist() == [0, 3]
    assert indices.tolist() == [0, 1, 2]
    assert signs.tolist() == [1, 1, 1]
    # the arrays are views: their base keeps the frozen complex alive
    assert offsets.base is not None
    o2, _, _ = f.boundary(2)
    assert np.shares_memory(offsets, o2)
    # queries
    assert f.star(0, 0) == [[0], [0, 2], [0]]
    assert f.link(0, 0) == [[1, 2], [1], []]


def test_fracture_pipeline_end_to_end():
    # quad + triangle polytopal mesh, cut, classify, extract
    c = g.from_polygons(5, [[0, 1, 2, 3], [1, 4, 2]])
    assert g.check_manifold(c).manifold_like
    interface = g.Marker(c)
    interface.mark_where(1, lambda e: e == 1)  # the shared edge (1,2)
    cut = g.cut_along(c, interface)
    cls = g.classify_facets(cut.complex)
    assert cls.maximal.marked(1, 1)  # the detached interface original
    fracture = g.subcomplex(cut.complex, g.Marker.from_cells(cut.complex, [[], [1]]))
    assert fracture.complex.counts() == [2, 1, 0]
    sides = g.connected_components(cut.complex, 2, 1)
    assert sides.count == 2


def test_coarsening_and_orientation():
    disk = make_disk()
    o = g.orient(disk)
    assert o.orientable
    agg = g.agglomerate(o.complex, [0, 0])
    assert agg.complex.counts() == [4, 4, 1]
    assert g.d_squared_is_zero(agg.complex)
    # amalgamation safety queries
    assert g.amalgamates_to_cell(disk, 2, 0, 1)
    cb = g.common_boundary(disk, 2, 0, 1)
    assert cb.n_facets == 1 and cb.acyclic


def test_subdivision_signed_carrier():
    sd = g.barycentric_subdivision(make_triangle())
    assert sd.complex.counts() == [7, 12, 6]
    assert g.d_squared_is_zero(sd.complex)
    signs = sd.carrier_sign(2)
    assert sorted(np.abs(signs).tolist()) == [1] * 6  # top flags all ±1
    assert sd.carrier_sign(0).tolist() == [1, 1, 1, 0, 0, 0, 0]


def test_collapse_preserves_homotopy():
    disk = make_disk()
    res = g.collapse(disk)
    assert res.complex.counts() == [1, 0, 0]  # collapsible to a point
    assert g.betti_numbers_z2(res.complex) == [1, 0, 0]


def test_exceptions_map_to_python():
    c = make_triangle()
    try:
        c.attach_cell(1, [0, 99], [-1, +1])
        assert False, "expected IndexError"
    except IndexError:
        pass
    try:
        g.agglomerate(c, [0, 0])  # wrong label count
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_chain_map_composition():
    disk = make_disk()
    m = g.Marker(disk)
    m.mark(0, 3)
    sd = g.star_deletion(disk, m)
    assert sd.complex.counts() == [3, 3, 1]
    idx = sd.map.index(0)
    assert int(idx[3]) == g.INVALID_INDEX
    assert [int(v) for v in idx[:3]] == [0, 1, 2]


if __name__ == "__main__":
    import sys
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
