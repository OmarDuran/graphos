#include <vector>

#include "fixtures.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos_test.hpp"

// Cut a two-triangle disk along the shared edge. The interface spans the
// whole domain (both endpoints reach the boundary), so the edge AND both
// endpoints split two ways: the bulk separates into two disjoint triangles
// and the original interface cells survive as a detached segment.
GRAPHOS_TEST(through_cut_separates_bulk_and_detaches_interface) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  // predicate-marked interface: the distribution-stable selection form
  graphos::Marker interface(c);
  interface.mark_where(1, [](graphos::Index i) { return i == 0; });
  const auto cut = graphos::cut_along(c, interface);
  cut.complex.validate();

  CHECK(cut.complex.count(0) == 8);  // 4 originals + 2 copies each endpoint
  CHECK(cut.complex.count(1) == 7);  // 5 originals + 2 copies of e0
  CHECK(cut.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(cut.complex));
  // two disjoint triangle disks + one detached segment
  CHECK(graphos::euler_characteristic(cut.complex) == 3);

  // ancestry: copies descend from the interface closure, originals from
  // themselves
  CHECK(cut.ancestor.index[1][5] == 0);
  CHECK(cut.ancestor.index[1][6] == 0);
  CHECK(cut.ancestor.index[0][4] == 0);
  CHECK(cut.ancestor.index[0][5] == 0);
  CHECK(cut.ancestor.index[0][6] == 1);
  CHECK(cut.ancestor.index[0][7] == 1);
  CHECK(cut.ancestor.index[2][0] == 0);

  // the two 2-cells reference DIFFERENT copies of the interface edge, and
  // nothing references the original anymore
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

// Cut one spoke of a 4-triangle fan around a center vertex. The interface
// ends at the center, so the center is a tip: the bulk stays connected
// around it and it must NOT be duplicated, while the outer endpoint (on the
// domain boundary) splits.
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
  // disk with an attached whisker (interface hangs on the tip): contractible
  CHECK(graphos::euler_characteristic(cut.complex) == 1);

  // both copies of the spoke still reach the ORIGINAL center vertex
  const graphos::BoundaryOperator& e = cut.complex.boundary(1);
  for (graphos::Index copy = 8; copy <= 9; ++copy) {
    CHECK(cut.ancestor.index[1][copy] == 0);
    bool touches_center = false;
    for (graphos::Index m = e.offsets[copy]; m < e.offsets[copy + 1]; ++m) {
      if (e.indices[m] == 4) touches_center = true;
    }
    CHECK(touches_center);
  }
  // exactly two distinct spoke copies are referenced by the bulk
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

GRAPHOS_TEST_MAIN()
