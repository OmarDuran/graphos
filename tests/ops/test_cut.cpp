#include <vector>

#include "fixtures.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos_test.hpp"

// Witnesses the through-cut: the interface reaches ∂K at both endpoints, so
// the 1-cell and both vertices split two ways. The bulk separates into two
// disjoint disks and cl(S) survives as a detached 1-dimensional subcomplex.
GRAPHOS_TEST(through_cut_separates_bulk_and_detaches_interface) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  // the interface, marked by predicate: the rank-independent selection form
  graphos::Marker interface(c);
  interface.mark_where(1, [](graphos::Index i) { return i == 0; });
  const auto cut = graphos::cut_along(c, interface);
  cut.complex.validate();

  CHECK(cut.complex.count(0) == 8);  // 4 originals + 2 copies each endpoint
  CHECK(cut.complex.count(1) == 7);  // 5 originals + 2 copies of e0
  CHECK(cut.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(cut.complex));
  // two disjoint disks and one detached segment
  CHECK(graphos::euler_characteristic(cut.complex) == 3);

  // ancestry: a copy descends from cl(S), an original from itself
  CHECK(cut.ancestor.index[1][5] == 0);
  CHECK(cut.ancestor.index[1][6] == 0);
  CHECK(cut.ancestor.index[0][4] == 0);
  CHECK(cut.ancestor.index[0][5] == 0);
  CHECK(cut.ancestor.index[0][6] == 1);
  CHECK(cut.ancestor.index[0][7] == 1);
  CHECK(cut.ancestor.index[2][0] == 0);

  // the two 2-cells reference distinct copies of the interface 1-cell, and
  // nothing references the original
  const graphos::BoundaryOperator& f = cut.complex.boundary(2);
  bool original_referenced = false;
  std::vector<graphos::Index> interface_refs;
  for (graphos::Index e = 0; e < cut.complex.count(2); ++e) {
    for (graphos::Index m = f.offsets[e]; m < f.offsets[e + 1]; ++m) {
      if (f.indices[m] == 0) original_referenced = true;
      if (cut.ancestor.index[1][f.indices[m]] == 0) interface_refs.push_back(f.indices[m]);
    }
  }
  CHECK(!original_referenced);
  CHECK(interface_refs.size() == 2);
  CHECK(interface_refs[0] != interface_refs[1]);
}

// Witnesses the crack front: the interface terminates at the interior vertex,
// which therefore has one side and is not copied, while its outer endpoint on
// ∂K splits. The bulk stays connected around the tip.
GRAPHOS_TEST(tip_vertex_is_not_duplicated) {
  const graphos::Complex c = graphos_test::make_fan();
  c.validate();
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);

  const auto cut = graphos::cut_along(c, graphos::Marker::from_cells(c, {{}, {0}}));
  cut.complex.validate();
  CHECK(cut.complex.count(0) == 7);   // corner 0 splits; center 4 does NOT
  CHECK(cut.complex.count(1) == 10);  // s0 gets two copies
  CHECK(cut.complex.count(2) == 4);
  CHECK(graphos::d_squared_is_zero(cut.complex));
  // a disk with a whisker attached at the tip: contractible
  CHECK(graphos::euler_characteristic(cut.complex) == 1);

  // both copies of the spoke still reach the original interior vertex
  const graphos::BoundaryOperator& e = cut.complex.boundary(1);
  for (graphos::Index copy = 8; copy <= 9; ++copy) {
    CHECK(cut.ancestor.index[1][copy] == 0);
    bool touches_center = false;
    for (graphos::Index m = e.offsets[copy]; m < e.offsets[copy + 1]; ++m) {
      if (e.indices[m] == 4) touches_center = true;
    }
    CHECK(touches_center);
  }
  // the bulk references exactly two distinct copies of the spoke
  const graphos::BoundaryOperator& f = cut.complex.boundary(2);
  std::vector<graphos::Index> spoke_refs;
  for (graphos::Index t = 0; t < cut.complex.count(2); ++t) {
    for (graphos::Index m = f.offsets[t]; m < f.offsets[t + 1]; ++m) {
      if (cut.ancestor.index[1][f.indices[m]] == 0) spoke_refs.push_back(f.indices[m]);
    }
  }
  CHECK(spoke_refs.size() == 2);
  CHECK(spoke_refs[0] != spoke_refs[1]);
  CHECK(spoke_refs[0] != 0 && spoke_refs[1] != 0);
}

GRAPHOS_TEST(rejects_marks_outside_interface_dimension) {
  const graphos::Complex c = graphos_test::make_triangle();
  graphos::Marker m(c);
  m.mark(0, 1);  // a vertex is not an (n-1)-cell
  CHECK_THROWS(graphos::cut_along(c, m));
  const graphos::Complex other = graphos_test::make_fan();
  CHECK_THROWS(graphos::cut_along(other, graphos::Marker(c)));  // wrong complex
}

// A copy carries the ∂ of the cell it descends from, so the ancestor map is a
// chain map and cochain transport through it commutes with d.
GRAPHOS_TEST(the_ancestor_map_is_a_chain_map) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::Marker iface(disk);
  iface.mark(1, 1);
  const graphos::CutResult cut = graphos::cut_along(disk, iface);
  CHECK(graphos::commutes_with_boundary(cut.complex, disk, cut.ancestor));

  // an empty cut is the identity, trivially a chain map
  const graphos::Marker none(disk);
  const graphos::CutResult id = graphos::cut_along(disk, none);
  CHECK(graphos::commutes_with_boundary(id.complex, disk, id.ancestor));
}

GRAPHOS_TEST_MAIN()
