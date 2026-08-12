#include "graphos/ops/disjoint_union.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

GRAPHOS_TEST(mixed_dimensional_coproduct) {
  const graphos::Complex a = graphos_test::make_triangle();
  const graphos::Complex b = graphos_test::make_segment();
  const auto du = graphos::disjoint_union(a, b);
  du.complex.validate();
  CHECK(du.complex.dim() == 2);
  CHECK(du.complex.count(0) == 5);
  CHECK(du.complex.count(1) == 4);
  CHECK(du.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(du.complex));
  CHECK(graphos::euler_characteristic(du.complex) == 2);  // disk + segment
}

GRAPHOS_TEST(chain_maps_shift_b_after_a) {
  const graphos::Complex a = graphos_test::make_triangle();
  const graphos::Complex b = graphos_test::make_segment();
  const auto du = graphos::disjoint_union(a, b);
  for (graphos::Index i = 0; i < a.count(0); ++i) CHECK(du.a_map.index[0][i] == i);
  CHECK(du.b_map.index[0][0] == 3);
  CHECK(du.b_map.index[0][1] == 4);
  CHECK(du.b_map.index[1][0] == 3);
  // B's edge references B's shifted vertices
  const graphos::BoundaryOperator& e = du.complex.boundary(1);
  CHECK(e.indices[e.offsets[3]] == 3);
  CHECK(e.indices[e.offsets[3] + 1] == 4);
}

GRAPHOS_TEST_MAIN()
