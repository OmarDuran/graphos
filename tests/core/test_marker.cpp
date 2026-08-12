#include "graphos/core/marker.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

using graphos::Index;
using graphos::Marker;

GRAPHOS_TEST(mark_and_query) {
  const graphos::Complex c = graphos_test::make_triangle();
  Marker m(c);
  CHECK(m.dim() == 2);
  CHECK(m.marked_count(1) == 0);
  m.mark(1, 0).mark(1, 2);
  CHECK(m.marked(1, 0));
  CHECK(!m.marked(1, 1));
  CHECK(m.marked(1, 2));
  CHECK(m.marked_count(1) == 2);
  CHECK(m.marked_count(0) == 0);
}

GRAPHOS_TEST(mark_where_predicate) {
  const graphos::Complex c = graphos_test::make_fan();
  Marker m(c);
  // the distribution-stable form: a locally evaluated predicate
  m.mark_where(1, [](Index i) { return i < 4; });  // the four spokes
  CHECK(m.marked_count(1) == 4);
  CHECK(m.marked(1, 3));
  CHECK(!m.marked(1, 4));
}

GRAPHOS_TEST(from_cells_factory) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const Marker m = Marker::from_cells(c, {{3}, {0, 4}});
  CHECK(m.marked(0, 3));
  CHECK(m.marked(1, 0));
  CHECK(m.marked(1, 4));
  CHECK(m.marked_count(2) == 0);
}

GRAPHOS_TEST(rejects_bad_input) {
  const graphos::Complex c = graphos_test::make_triangle();
  Marker m(c);
  CHECK_THROWS(m.mark(3, 0));
  CHECK_THROWS(m.mark(0, 99));
  CHECK_THROWS(Marker::from_cells(c, {{99}}));
}

GRAPHOS_TEST(validate_for_detects_wrong_complex) {
  const graphos::Complex tri = graphos_test::make_triangle();
  const graphos::Complex fan = graphos_test::make_fan();
  const Marker m(tri);
  m.validate_for(tri);  // fine
  CHECK_THROWS(m.validate_for(fan));
}

GRAPHOS_TEST_MAIN()
