#include "fixtures.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos_test.hpp"

using graphos::Index;

// Witnesses that subcomplex extracts cl(S): restricting the disk to one 2-cell
// yields that cell with its closure.
GRAPHOS_TEST(one_face_restriction_is_a_triangle) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker m(c);
  m.mark(2, 0);  // face A
  const auto sub = graphos::subcomplex(c, m);
  sub.complex.validate();

  CHECK(sub.complex.count(0) == 3);
  CHECK(sub.complex.count(1) == 3);
  CHECK(sub.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(sub.complex));
  CHECK(graphos::euler_characteristic(sub.complex) == 1);

  // B's private cells lie outside cl(S)
  CHECK(sub.map.index[0][3] == graphos::invalid_index);
  CHECK(sub.map.index[1][3] == graphos::invalid_index);
  CHECK(sub.map.index[1][4] == graphos::invalid_index);
  CHECK(sub.map.index[2][1] == graphos::invalid_index);
  // the embedding inverts the map on survivors
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < sub.embedding.index[k].size(); ++i) {
      const Index parent = sub.embedding.index[k][i];
      CHECK(sub.map.index[k][static_cast<std::size_t>(parent)] == static_cast<Index>(i));
    }
  }
}

// Witnesses skeleton extraction: marking every 1-cell yields the 1-skeleton,
// a circle, with χ = 0.
GRAPHOS_TEST(one_skeleton_of_a_triangle_is_a_circle) {
  const graphos::Complex c = graphos_test::make_triangle();
  graphos::Marker m(c);
  m.mark_where(1, [](Index) { return true; });
  const auto sub = graphos::subcomplex(c, m);
  sub.complex.validate();
  CHECK(sub.complex.count(0) == 3);
  CHECK(sub.complex.count(1) == 3);
  CHECK(sub.complex.count(2) == 0);
  CHECK(graphos::euler_characteristic(sub.complex) == 0);  // circle
}

// Witnesses that cut and subcomplex compose: cutting the disk along the shared
// 1-cell and extracting the detached interface gives it as its own complex.
GRAPHOS_TEST(extracts_fracture_domain_after_a_cut) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cut = graphos::cut_along(c, interface);

  // the original interface survives the cut at its old index
  graphos::Marker fracture(cut.complex);
  fracture.mark(1, 0);
  const auto sub = graphos::subcomplex(cut.complex, fracture);
  sub.complex.validate();

  // a segment: the interface 1-cell with its two original endpoints
  CHECK(sub.complex.count(0) == 2);
  CHECK(sub.complex.count(1) == 1);
  CHECK(sub.complex.count(2) == 0);
  CHECK(graphos::euler_characteristic(sub.complex) == 1);
  // embedded at the parent's original cells
  CHECK(sub.embedding.index[1][0] == 0);
  CHECK(sub.embedding.index[0][0] == 0);
  CHECK(sub.embedding.index[0][1] == 1);
}

GRAPHOS_TEST(empty_marker_yields_empty_complex) {
  const graphos::Complex c = graphos_test::make_triangle();
  const auto sub = graphos::subcomplex(c, graphos::Marker(c));
  CHECK(sub.complex.count(0) == 0);
  CHECK(sub.complex.count(1) == 0);
  CHECK(sub.complex.count(2) == 0);
}

GRAPHOS_TEST(rejects_marker_for_wrong_complex) {
  const graphos::Complex tri = graphos_test::make_triangle();
  const graphos::Complex fan = graphos_test::make_fan();
  CHECK_THROWS(graphos::subcomplex(fan, graphos::Marker(tri)));
}

GRAPHOS_TEST_MAIN()
