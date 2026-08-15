#include "fixtures.hpp"
#include "graphos/ops/agglomerate.hpp"
#include "graphos/ops/orient.hpp"
#include "graphos_test.hpp"

using graphos::Index;

// Witnesses cancellation: merging both 2-cells of the disk gives one polytopal
// cell whose shared 1-cell cancels from Σ ∂(members), leaving four faces.
GRAPHOS_TEST(disk_agglomerates_into_one_polytopal_cell) {
  const auto oriented = graphos::orient(graphos_test::make_two_triangle_disk());
  const auto agg = graphos::agglomerate(oriented.complex, {0, 0});
  agg.complex.validate();

  CHECK(agg.complex.count(0) == 4);
  CHECK(agg.complex.count(1) == 4);  // the shared edge is gone
  CHECK(agg.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(agg.complex));
  CHECK(graphos::euler_characteristic(agg.complex) == 1);

  // the interior facet maps to 0, the top cells to the aggregate
  CHECK(agg.map.index[1][0] == graphos::invalid_index);
  CHECK(agg.map.index[2][0] == 0);
  CHECK(agg.map.index[2][1] == 0);
}

GRAPHOS_TEST(fan_agglomerates_into_two_cells) {
  const graphos::Complex c = graphos_test::make_fan();  // already consistent
  const auto agg = graphos::agglomerate(c, {0, 0, 1, 1});
  agg.complex.validate();

  // the interfaces s0 and s2 survive; s1 and s3 are interior
  CHECK(agg.complex.count(0) == 5);  // center survives on the interface
  CHECK(agg.complex.count(1) == 6);  // 4 boundary + 2 interface spokes
  CHECK(agg.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(agg.complex));
  CHECK(graphos::euler_characteristic(agg.complex) == 1);
  CHECK(agg.map.index[1][1] == graphos::invalid_index);  // s1 sent to zero
  CHECK(agg.map.index[1][3] == graphos::invalid_index);  // s3 sent to zero
  CHECK(agg.map.index[1][0] != graphos::invalid_index);  // s0 survives
}

// Witnesses mixed-dimensional passthrough: a maximal lower-dimensional
// stratum survives agglomeration of the bulk untouched.
GRAPHOS_TEST(preserves_maximal_lower_dimensional_cells) {
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {3, 4}, {-1, +1});  // detached segment
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});

  const auto agg = graphos::agglomerate(c, {0});
  agg.complex.validate();
  CHECK(agg.complex.count(0) == 5);
  CHECK(agg.complex.count(1) == 4);
  CHECK(agg.complex.count(2) == 1);
  CHECK(agg.map.index[1][3] != graphos::invalid_index);
}

GRAPHOS_TEST(detects_inconsistent_orientation) {
  // the raw disk references the shared 1-cell with +1 from both cofaces, so
  // the coefficient would be 2 and the merge is rejected
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  CHECK_THROWS(graphos::agglomerate(c, {0, 0}));
}

GRAPHOS_TEST(rejects_bad_labels) {
  const graphos::Complex c = graphos_test::make_fan();
  CHECK_THROWS(graphos::agglomerate(c, {0, 0, 1}));      // wrong size
  CHECK_THROWS(graphos::agglomerate(c, {0, 0, 2, 2}));   // id 1 unused
  CHECK_THROWS(graphos::agglomerate(c, {0, 0, -1, 1}));  // negative
}

GRAPHOS_TEST_MAIN()
