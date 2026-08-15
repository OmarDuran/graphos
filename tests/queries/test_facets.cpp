#include "fixtures.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos/queries/facets.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(disk_facets_split_into_boundary_and_interior) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto cls = graphos::classify_facets(c);
  CHECK(cls.interior.marked(1, 0));  // the shared edge
  CHECK(cls.interior.marked_count(1) == 1);
  CHECK(cls.boundary.marked_count(1) == 4);
  CHECK(cls.maximal.marked_count(1) == 0);
  CHECK(cls.nonmanifold.marked_count(1) == 0);
}

// Witnesses the classification after a through-cut: the detached interface has
// no coface and is maximal, while every remaining bulk facet has exactly
// one.
GRAPHOS_TEST(cut_turns_interface_maximal_and_copies_into_boundary) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cut = graphos::cut_along(c, interface);

  const auto cls = graphos::classify_facets(cut.complex);
  CHECK(cls.maximal.marked(1, 0));  // the fracture domain cell
  CHECK(cls.maximal.marked_count(1) == 1);
  CHECK(cls.boundary.marked_count(1) == 6);  // both triangles fully exposed
  CHECK(cls.interior.marked_count(1) == 0);
}

// Witnesses the non-manifold case: three 2-cells over one 1-cell give it three
// cofaces, a book junction.
GRAPHOS_TEST(book_spine_is_nonmanifold) {
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 1}, {-1, +1});  // the spine
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {1, 3}, {-1, +1});
  c.attach_cell(1, {3, 0}, {-1, +1});
  c.attach_cell(1, {1, 4}, {-1, +1});
  c.attach_cell(1, {4, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  c.attach_cell(2, {0, 3, 4}, {+1, +1, +1});
  c.attach_cell(2, {0, 5, 6}, {+1, +1, +1});
  c.validate();
  CHECK(graphos::d_squared_is_zero(c));

  const auto cls = graphos::classify_facets(c);
  CHECK(cls.nonmanifold.marked(1, 0));
  CHECK(cls.nonmanifold.marked_count(1) == 1);
  CHECK(cls.boundary.marked_count(1) == 6);
  CHECK(cls.interior.marked_count(1) == 0);
}

// Witnesses the composition the classification exists for: extracting ∂K as
// its own complex.
GRAPHOS_TEST(boundary_complex_of_a_disk_is_a_circle) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto sub = graphos::subcomplex(c, graphos::classify_facets(c).boundary);
  sub.complex.validate();
  CHECK(sub.complex.count(0) == 4);
  CHECK(sub.complex.count(1) == 4);
  CHECK(sub.complex.count(2) == 0);
  CHECK(graphos::d_squared_is_zero(sub.complex));
  CHECK(graphos::euler_characteristic(sub.complex) == 0);  // a circle
}

GRAPHOS_TEST(closedness_and_boundary_of_boundary) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  CHECK(!graphos::is_closed(disk));
  // the extracted boundary is itself closed: ∂∂K = ∅
  const auto bd = graphos::subcomplex(disk, graphos::classify_facets(disk).boundary);
  CHECK(graphos::is_closed(bd.complex));

  graphos::Complex circle(1);
  circle.attach_vertices(3);
  circle.attach_cell(1, {0, 1}, {-1, +1});
  circle.attach_cell(1, {1, 2}, {-1, +1});
  circle.attach_cell(1, {2, 0}, {-1, +1});
  CHECK(graphos::is_closed(circle));
}

GRAPHOS_TEST(rejects_zero_dimensional_complex) {
  graphos::Complex c(0);
  c.attach_vertices(3);
  CHECK_THROWS(graphos::classify_facets(c));
}

GRAPHOS_TEST_MAIN()
