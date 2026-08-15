#include "fixtures.hpp"
#include "graphos/ops/collapse.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(free_faces_of_the_fan_are_its_boundary_edges) {
  const graphos::Marker ff = graphos::free_faces(graphos_test::make_fan());
  CHECK(ff.marked_count(0) == 0);  // every vertex has several cofaces
  CHECK(ff.marked_count(1) == 4);  // exactly the boundary edges
  CHECK(ff.marked(1, 4));
  CHECK(!ff.marked(1, 0));  // spokes have two cofaces
}

GRAPHOS_TEST(disk_collapses_to_a_point) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto res = graphos::collapse(c);
  CHECK(res.complex.count(0) == 1);
  CHECK(res.complex.count(1) == 0);
  CHECK(res.complex.count(2) == 0);
  // (4 + 5 + 2 cells − 1 survivor) / 2 elementary collapses
  CHECK(res.removed_pairs == 5);
  // one cell survives, and it is a vertex: the complex is collapsible
  Index survivors = 0;
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < res.map.index[k].size(); ++i) {
      if (res.map.index[k][i] != graphos::invalid_index) ++survivors;
    }
  }
  CHECK(survivors == 1);
}

GRAPHOS_TEST(collapse_preserves_the_homotopy_type) {
  const graphos::Complex fan = graphos_test::make_fan();
  const auto before = graphos::betti_numbers_z2(fan);
  const auto res = graphos::collapse(fan);
  const auto after = graphos::betti_numbers_z2(res.complex);
  CHECK(before == after);  // simple homotopy equivalence
}

GRAPHOS_TEST(moebius_band_collapses_to_its_core_circle) {
  graphos::Complex square(2);
  square.attach_vertices(4);
  square.attach_cell(1, {0, 1}, {-1, +1});
  square.attach_cell(1, {1, 2}, {-1, +1});
  square.attach_cell(1, {3, 2}, {-1, +1});
  square.attach_cell(1, {0, 3}, {-1, +1});
  square.attach_cell(1, {0, 2}, {-1, +1});
  square.attach_cell(2, {0, 1, 4}, {+1, +1, -1});
  square.attach_cell(2, {4, 2, 3}, {+1, -1, -1});
  const auto ids = graphos::lift_identifications(square, {{0, 2, +1}, {3, 1, +1}});
  const auto moebius = graphos::quotient(square, ids);

  const auto res = graphos::collapse(moebius.complex);
  CHECK(res.complex.count(2) == 0);  // the band is gone
  CHECK(res.complex.count(1) > 0);   // the core circle is not
  CHECK(graphos::betti_numbers_z2(res.complex) == (std::vector<Index>{1, 1, 0}));
  // the endpoint carries no free face
  const graphos::Marker ff = graphos::free_faces(res.complex);
  for (int k = 0; k <= 2; ++k) CHECK(ff.marked_count(k) == 0);
}

GRAPHOS_TEST_MAIN()
