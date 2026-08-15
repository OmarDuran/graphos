#include "fixtures.hpp"
#include "graphos/ops/join.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

graphos::Complex make_point() {
  graphos::Complex c(0);
  c.attach_vertices(1);
  return c;
}

graphos::Complex make_two_points() {
  graphos::Complex c(0);
  c.attach_vertices(2);  // S⁰
  return c;
}

graphos::Complex make_circle() {
  graphos::Complex c(1);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  return c;
}

}  // namespace

GRAPHOS_TEST(point_joined_with_segment_is_a_triangle) {
  const auto j = graphos::join(make_point(), graphos_test::make_segment());
  j.complex.validate();
  CHECK(j.complex.dim() == 2);
  CHECK(j.complex.count(0) == 3);
  CHECK(j.complex.count(1) == 3);
  CHECK(j.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(j.complex));  // the join signs
  CHECK(graphos::euler_characteristic(j.complex) == 1);
}

GRAPHOS_TEST(cone_over_a_circle_is_a_disk) {
  const auto j = graphos::join(make_point(), make_circle());
  j.complex.validate();
  CHECK(j.complex.count(0) == 4);
  CHECK(j.complex.count(1) == 6);
  CHECK(j.complex.count(2) == 3);
  CHECK(graphos::d_squared_is_zero(j.complex));
  CHECK(graphos::betti_numbers_z2(j.complex) == (std::vector<Index>{1, 0, 0}));
}

GRAPHOS_TEST(suspension_of_a_circle_is_a_sphere) {
  const auto j = graphos::join(make_two_points(), make_circle());
  j.complex.validate();
  CHECK(j.complex.count(0) == 5);
  CHECK(j.complex.count(1) == 9);
  CHECK(j.complex.count(2) == 6);
  CHECK(graphos::d_squared_is_zero(j.complex));
  CHECK(graphos::euler_characteristic(j.complex) == 2);  // χ(S²)
  CHECK(graphos::betti_numbers_z2(j.complex) == (std::vector<Index>{1, 0, 1}));
}

GRAPHOS_TEST(factors_embed_as_subcomplexes) {
  const graphos::Complex a = graphos_test::make_segment();
  const graphos::Complex b = make_circle();
  const auto j = graphos::join(a, b);
  // A occupies the leading indices, B follows
  CHECK(j.a_map.index[1][0] == 0);
  CHECK(j.b_map.index[0][0] == a.count(0));
  CHECK(j.b_map.index[1][0] == a.count(1));
  // B's embedded edge references B's embedded vertices
  const graphos::BoundaryOperator& e = j.complex.boundary(1);
  const Index be = j.b_map.index[1][0];
  CHECK(e.indices[e.offsets[be]] == j.b_map.index[0][0]);
}

GRAPHOS_TEST_MAIN()
