#include "fixtures.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

// the square as two 2-cells: b=[0,1](0), r=[1,2](1), t=[3,2](2), l=[0,3](3),
// diagonal d=[0,2](4)
graphos::Complex make_square() {
  graphos::Complex c(2);
  c.attach_vertices(4);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {3, 2}, {-1, +1});
  c.attach_cell(1, {0, 3}, {-1, +1});
  c.attach_cell(1, {0, 2}, {-1, +1});
  c.attach_cell(2, {0, 1, 4}, {+1, +1, -1});
  c.attach_cell(2, {4, 2, 3}, {+1, -1, -1});
  c.validate();
  return c;
}

}  // namespace

// Witnesses that nothing lifts without fully mapped faces: identifying the
// ends of a chain of 1-cells leaves the interior vertex unmapped, so no 1-cell
// lifts and the quotient is a circle.
GRAPHOS_TEST(chain_to_circle_lifts_nothing_above_vertices) {
  graphos::Complex c(1);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});

  const auto ids = graphos::lift_identifications(c, {{2, 0, +1}});
  CHECK(ids[0].size() == 1);
  CHECK(ids[1].empty());

  const auto q = graphos::quotient(c, ids);
  CHECK(q.complex.count(0) == 2);
  CHECK(q.complex.count(1) == 2);
  CHECK(graphos::euler_characteristic(q.complex) == 0);  // a circle
}

// Witnesses the periodic identification: pairing the vertices of the left
// 1-cell with those of the right lifts the 1-cells themselves, and the
// quotient is a cylinder.
GRAPHOS_TEST(periodic_square_becomes_a_cylinder) {
  const graphos::Complex c = make_square();
  const auto ids = graphos::lift_identifications(c, {{0, 1, +1}, {3, 2, +1}});
  CHECK(ids[0].size() == 2);
  CHECK(ids[1].size() == 1);
  CHECK(ids[1][0].from == 3);       // l
  CHECK(ids[1][0].to == 1);         // r
  CHECK(ids[1][0].rel_sign == +1);  // same orientation
  CHECK(ids[2].empty());

  const auto q = graphos::quotient(c, ids);
  q.complex.validate();
  CHECK(q.complex.count(0) == 2);
  CHECK(q.complex.count(1) == 4);
  CHECK(q.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(q.complex));
  CHECK(graphos::euler_characteristic(q.complex) == 0);  // cylinder
}

// Witnesses the reversed lift: the twisted pairing carries the left 1-cell
// onto the right with rel_sign = −1, and the quotient is a Möbius band —
// non-orientable, which orient() then reports.
GRAPHOS_TEST(twisted_pairing_becomes_a_moebius_band) {
  const graphos::Complex c = make_square();
  const auto ids = graphos::lift_identifications(c, {{0, 2, +1}, {3, 1, +1}});
  CHECK(ids[1].size() == 1);
  CHECK(ids[1][0].from == 3);
  CHECK(ids[1][0].to == 1);
  CHECK(ids[1][0].rel_sign == -1);  // the twist

  const auto q = graphos::quotient(c, ids);
  q.complex.validate();
  CHECK(q.complex.count(0) == 2);
  CHECK(q.complex.count(1) == 4);
  CHECK(q.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(q.complex));
  CHECK(graphos::euler_characteristic(q.complex) == 0);  // Möbius band
}

GRAPHOS_TEST(rejects_bad_vertex_pairs) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK_THROWS(graphos::lift_identifications(c, {{0, 99, +1}}));
  CHECK_THROWS(graphos::lift_identifications(c, {{0, 1, +2}}));
}

GRAPHOS_TEST_MAIN()
