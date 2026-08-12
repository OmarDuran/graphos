#include "graphos/ops/product.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

// circle: three vertices, three edges, no faces (χ = 0)
graphos::Complex make_circle() {
  graphos::Complex c(1);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  return c;
}

}  // namespace

GRAPHOS_TEST(interval_times_interval_is_a_quad) {
  const graphos::Complex seg = graphos_test::make_segment();
  const auto prod = graphos::product(seg, seg);
  prod.complex.validate();
  CHECK(prod.complex.dim() == 2);
  CHECK(prod.complex.count(0) == 4);
  CHECK(prod.complex.count(1) == 4);
  CHECK(prod.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(prod.complex));  // the Leibniz signs
  CHECK(graphos::euler_characteristic(prod.complex) == 1);
}

GRAPHOS_TEST(triangle_extruded_is_a_prism) {
  const auto prod = graphos::product(graphos_test::make_triangle(), graphos_test::make_segment());
  prod.complex.validate();
  CHECK(prod.complex.dim() == 3);
  CHECK(prod.complex.count(0) == 6);
  CHECK(prod.complex.count(1) == 9);   // 3 edges x 2 layers + 3 verticals
  CHECK(prod.complex.count(2) == 5);   // 2 lids + 3 side quads
  CHECK(prod.complex.count(3) == 1);
  CHECK(graphos::d_squared_is_zero(prod.complex));
  CHECK(graphos::euler_characteristic(prod.complex) == 1);

  // block layout: ascending p, i-major
  CHECK(prod.blocks[2].size() == 2);
  CHECK(prod.blocks[2][0].p == 1);  // side quads: edge x edge
  CHECK(prod.blocks[2][1].p == 2);  // lids: face x vertex
  CHECK(prod.blocks[2][1].offset == 3);
}

GRAPHOS_TEST(circle_times_circle_is_a_torus) {
  const graphos::Complex circle = make_circle();
  const auto prod = graphos::product(circle, circle);
  prod.complex.validate();
  CHECK(prod.complex.count(0) == 9);
  CHECK(prod.complex.count(1) == 18);
  CHECK(prod.complex.count(2) == 9);
  CHECK(graphos::d_squared_is_zero(prod.complex));
  CHECK(graphos::euler_characteristic(prod.complex) == 0);  // χ(T²) = 0
}

GRAPHOS_TEST(euler_characteristic_multiplies) {
  const auto prod = graphos::product(graphos_test::make_two_triangle_disk(), make_circle());
  prod.complex.validate();
  CHECK(graphos::d_squared_is_zero(prod.complex));
  CHECK(graphos::euler_characteristic(prod.complex) == 0);  // 1 x 0
}

GRAPHOS_TEST_MAIN()
