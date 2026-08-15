#include "fixtures.hpp"
#include "graphos/queries/neighborhood.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(star_of_the_fan_center) {
  const graphos::Complex c = graphos_test::make_fan();
  graphos::Marker s(c);
  s.mark(0, 4);
  const graphos::Marker st = graphos::star_of(c, s);
  CHECK(st.marked_count(0) == 1);
  CHECK(st.marked_count(1) == 4);  // the spokes
  CHECK(st.marked_count(2) == 4);  // all faces
  CHECK(!st.marked(1, 4));         // boundary edges excluded
}

GRAPHOS_TEST(closure_of_a_face) {
  const graphos::Complex c = graphos_test::make_fan();
  graphos::Marker s(c);
  s.mark(2, 0);
  const graphos::Marker cl = graphos::closure_of(c, s);
  CHECK(cl.marked_count(2) == 1);
  CHECK(cl.marked_count(1) == 3);
  CHECK(cl.marked_count(0) == 3);
}

GRAPHOS_TEST(link_and_frontier_agree_for_a_vertex) {
  const graphos::Complex c = graphos_test::make_fan();
  graphos::Marker s(c);
  s.mark(0, 4);
  const graphos::Marker lk = graphos::link_of(c, s);
  const graphos::Marker fr = graphos::frontier_of(c, s);
  // both are the boundary square: 4 corners + 4 boundary edges
  for (int k = 0; k <= 2; ++k) {
    for (Index i = 0; i < c.count(k); ++i) {
      CHECK(lk.marked(k, i) == fr.marked(k, i));
    }
  }
  CHECK(lk.marked_count(0) == 4);
  CHECK(lk.marked_count(1) == 4);
  CHECK(lk.marked_count(2) == 0);
}

GRAPHOS_TEST(link_and_frontier_differ_for_a_top_cell) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker s(c);
  s.mark(2, 0);  // face A
  // the link of a top cell is empty...
  const graphos::Marker lk = graphos::link_of(c, s);
  for (int k = 0; k <= 2; ++k) CHECK(lk.marked_count(k) == 0);
  // ...but its excision frontier is its boundary
  const graphos::Marker fr = graphos::frontier_of(c, s);
  CHECK(fr.marked_count(0) == 3);
  CHECK(fr.marked_count(1) == 3);
  CHECK(fr.marked_count(2) == 0);
}

GRAPHOS_TEST_MAIN()
