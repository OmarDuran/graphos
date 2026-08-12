#include "graphos/queries/components.hpp"

#include "graphos/ops/cut.hpp"
#include "graphos/ops/disjoint_union.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(disk_top_cells_are_one_component) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto cc = graphos::connected_components(c, 2, 1);
  CHECK(cc.count == 1);
  CHECK(cc.label[0] == 0);
  CHECK(cc.label[1] == 0);
}

// The "sides of a cut" computation: excluding the interface facet as a
// connector separates the two triangles.
GRAPHOS_TEST(excluding_the_interface_separates_sides) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cc = graphos::connected_components(c, 2, 1, interface);
  CHECK(cc.count == 2);
  // deterministic: components numbered by lowest-indexed cell
  CHECK(cc.label[0] == 0);
  CHECK(cc.label[1] == 1);
}

GRAPHOS_TEST(cut_result_has_two_bulk_components) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cut = graphos::cut_along(c, interface);
  const auto cc = graphos::connected_components(cut.complex, 2, 1);
  CHECK(cc.count == 2);
}

GRAPHOS_TEST(vertices_via_edges) {
  const graphos::Complex fan = graphos_test::make_fan();
  CHECK(graphos::connected_components(fan, 0, 1).count == 1);
}

GRAPHOS_TEST(whole_complex_components_count_beta_zero) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  CHECK(graphos::connected_components(c).count == 1);

  // disjoint union: two pieces
  const auto du = graphos::disjoint_union(c, graphos_test::make_segment());
  const auto cc = graphos::connected_components(du.complex);
  CHECK(cc.count == 2);

  // after a through-cut: two triangles + the detached interface segment
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cut = graphos::cut_along(c, interface);
  const auto cut_cc = graphos::connected_components(cut.complex);
  CHECK(cut_cc.count == 3);
  // the fracture original and its endpoints share a component
  CHECK(cut_cc.label[1][0] == cut_cc.label[0][0]);
  CHECK(cut_cc.label[1][0] == cut_cc.label[0][1]);
}

GRAPHOS_TEST(rejects_bad_dimensions) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK_THROWS(graphos::connected_components(c, 2, 2));
  CHECK_THROWS(graphos::connected_components(c, 3, 1));
}

GRAPHOS_TEST_MAIN()
