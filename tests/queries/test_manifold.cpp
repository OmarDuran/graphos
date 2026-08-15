#include "fixtures.hpp"
#include "graphos/core/build.hpp"
#include "graphos/queries/manifold.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(disk_and_fan_are_manifold_like) {
  const auto disk = graphos::check_manifold(graphos_test::make_two_triangle_disk());
  CHECK(disk.manifold_like);
  const auto fan = graphos::check_manifold(graphos_test::make_fan());
  CHECK(fan.manifold_like);
}

GRAPHOS_TEST(book_junction_fails_the_facet_check) {
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 1}, {-1, +1});  // spine
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {1, 3}, {-1, +1});
  c.attach_cell(1, {3, 0}, {-1, +1});
  c.attach_cell(1, {1, 4}, {-1, +1});
  c.attach_cell(1, {4, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  c.attach_cell(2, {0, 3, 4}, {+1, +1, +1});
  c.attach_cell(2, {0, 5, 6}, {+1, +1, +1});

  const auto rep = graphos::check_manifold(c);
  CHECK(!rep.manifold_like);
  CHECK(!rep.facet_condition);
  CHECK(rep.pure);
  CHECK(rep.offending.marked(1, 0));  // the spine
  CHECK(rep.offending.marked_count(1) == 1);
}

GRAPHOS_TEST(pinched_vertex_fails_link_connectivity) {
  // two 2-cells meeting in a single vertex: a pinch point
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {0, 3}, {-1, +1});
  c.attach_cell(1, {3, 4}, {-1, +1});
  c.attach_cell(1, {4, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  c.attach_cell(2, {3, 4, 5}, {+1, +1, +1});

  const auto rep = graphos::check_manifold(c);
  CHECK(!rep.manifold_like);
  CHECK(!rep.links_connected);
  CHECK(rep.pure);
  CHECK(rep.facet_condition);
  CHECK(rep.offending.marked(0, 0));  // the pinch point
  CHECK(rep.offending.marked_count(0) == 1);
}

GRAPHOS_TEST(mixed_dimensional_complex_reports_impure) {
  // a 2-cell with a maximal 1-cell attached: not pure, and reported as such
  graphos::Complex c(2);
  c.attach_vertices(4);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {0, 3}, {-1, +1});  // whisker
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});

  const auto rep = graphos::check_manifold(c);
  CHECK(!rep.manifold_like);
  CHECK(!rep.pure);
  CHECK(rep.offending.marked(1, 3));
}

GRAPHOS_TEST(vertex_links_classify_sphere_and_ball) {
  // in the fan, lk(interior vertex) is a circle and lk(rim vertex) an
  // interval
  const auto fan = graphos::classify_vertex_links(graphos_test::make_fan());
  CHECK(fan.sphere.marked(0, 4));
  CHECK(fan.sphere.marked_count(0) == 1);
  CHECK(fan.ball.marked_count(0) == 4);
  CHECK(fan.other.marked_count(0) == 0);

  // in the path, lk(interior vertex) = S⁰ and the endpoints are balls
  graphos::Complex path(1);
  path.attach_vertices(3);
  path.attach_cell(1, {0, 1}, {-1, +1});
  path.attach_cell(1, {1, 2}, {-1, +1});
  const auto pl = graphos::classify_vertex_links(path);
  CHECK(pl.sphere.marked(0, 1));
  CHECK(pl.ball.marked(0, 0));
  CHECK(pl.ball.marked(0, 2));
}

GRAPHOS_TEST(pinch_vertex_link_is_a_defect) {
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {0, 3}, {-1, +1});
  c.attach_cell(1, {3, 4}, {-1, +1});
  c.attach_cell(1, {4, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  c.attach_cell(2, {3, 4, 5}, {+1, +1, +1});
  const auto lc = graphos::classify_vertex_links(c);
  CHECK(lc.other.marked(0, 0));  // the pinch: link is two disjoint intervals
  CHECK(lc.other.marked_count(0) == 1);
  CHECK(lc.ball.marked_count(0) == 4);
}

GRAPHOS_TEST(tetrahedron_vertex_links_are_disks) {
  const graphos::Complex tet = graphos::from_simplices(3, 4, {{0, 1, 2, 3}});
  const auto lc = graphos::classify_vertex_links(tet);
  CHECK(lc.ball.marked_count(0) == 4);  // each link is one triangle: a disk
  CHECK(lc.sphere.marked_count(0) == 0);
}

GRAPHOS_TEST(points_are_a_zero_manifold) {
  graphos::Complex c(0);
  c.attach_vertices(7);
  CHECK(graphos::check_manifold(c).manifold_like);
}

GRAPHOS_TEST_MAIN()
