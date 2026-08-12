#include "graphos/ops/quotient.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

GRAPHOS_TEST(identification_chains_resolve_to_root) {
  graphos::Complex c(0);
  c.attach_vertices(3);
  std::vector<std::vector<graphos::Identification>> ids(1);
  ids[0] = {{2, 1, +1}, {1, 0, +1}};  // 2 ~ 1 ~ 0
  const auto q = graphos::quotient(c, ids);
  CHECK(q.complex.count(0) == 1);
  CHECK(q.map.index[0][0] == 0);
  CHECK(q.map.index[0][1] == 0);
  CHECK(q.map.index[0][2] == 0);
}

GRAPHOS_TEST(cyclic_identifications_throw) {
  graphos::Complex c(0);
  c.attach_vertices(2);
  std::vector<std::vector<graphos::Identification>> ids(1);
  ids[0] = {{0, 1, +1}, {1, 0, +1}};
  CHECK_THROWS(graphos::quotient(c, ids));
}

GRAPHOS_TEST(quotient_rewrites_boundary_through_vertex_gluing) {
  // glue segment [0,1] and segment [2,3] end to end: 1 ~ 2
  graphos::Complex c(1);
  c.attach_vertices(4);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {2, 3}, {-1, +1});
  std::vector<std::vector<graphos::Identification>> ids(2);
  ids[0] = {{2, 1, +1}};
  const auto q = graphos::quotient(c, ids);
  q.complex.validate();
  CHECK(q.complex.count(0) == 3);
  CHECK(q.complex.count(1) == 2);
  const graphos::BoundaryOperator& e = q.complex.boundary(1);
  CHECK(e.indices[e.offsets[1]] == 1);  // second edge now starts at glued vertex
}

GRAPHOS_TEST(find_parallel_cells_reports_orientation) {
  graphos::Complex c(1);
  c.attach_vertices(2);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {0, 1}, {+1, -1});  // same faces, reversed
  c.attach_cell(1, {0, 1}, {-1, +1});  // same faces, same orientation
  const auto dups = graphos::find_parallel_cells(c, 1);
  CHECK(dups.size() == 2);
  CHECK(dups[0].from == 1 && dups[0].to == 0 && dups[0].rel_sign == -1);
  CHECK(dups[1].from == 2 && dups[1].to == 0 && dups[1].rel_sign == +1);
}

GRAPHOS_TEST_MAIN()
