#include "graphos/core/complex.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

GRAPHOS_TEST(triangle_invariants) {
  const graphos::Complex c = graphos_test::make_triangle();
  c.validate();
  CHECK(c.count(0) == 3);
  CHECK(c.count(1) == 3);
  CHECK(c.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);  // disk
}

GRAPHOS_TEST(segment_invariants) {
  const graphos::Complex c = graphos_test::make_segment();
  c.validate();
  CHECK(graphos::euler_characteristic(c) == 1);
}

GRAPHOS_TEST(attach_cell_rejects_bad_input) {
  graphos::Complex c(2);
  c.attach_vertices(2);
  CHECK_THROWS(c.attach_cell(1, {0, 7}, {-1, +1}));  // face out of range
  CHECK_THROWS(c.attach_cell(1, {0, 1}, {-1, +2}));  // sign not +/-1
  CHECK_THROWS(c.attach_cell(3, {0, 1}, {-1, +1}));  // dimension out of range
  CHECK(c.count(1) == 0);
}

GRAPHOS_TEST(d_squared_detects_broken_orientation) {
  graphos::Complex c(2);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  // e0 + e1 - e2 does not close: ∂∂ = 2v2 - 2v0 != 0
  c.attach_cell(2, {0, 1, 2}, {+1, +1, -1});
  c.validate();  // structurally fine...
  CHECK(!graphos::d_squared_is_zero(c));  // ...but not a chain complex
}

GRAPHOS_TEST(bulk_constructor_checks_stratification) {
  std::vector<graphos::BoundaryOperator> strata(2);
  strata[1].append_row(std::vector<graphos::Index>{0, 1}, std::vector<graphos::Sign>{-1, +1});
  // counts claim two edges but the operator has one row
  CHECK_THROWS(graphos::Complex({2, 2}, std::move(strata)));
}

GRAPHOS_TEST_MAIN()
