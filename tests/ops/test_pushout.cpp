#include "fixtures.hpp"
#include "graphos/ops/pushout.hpp"
#include "graphos_test.hpp"

// Witnesses that deduplication carries relative orientation: B is attached
// through reversed vertex identifications, so the shared 1-cell is discovered
// with rel_sign = −1 and every reference to it flips.
GRAPHOS_TEST(pushout_with_orientation_flip) {
  const graphos::Complex a = graphos_test::make_triangle();
  const graphos::Complex b = graphos_test::make_triangle();
  const std::vector<graphos::Identification> vertex_ids = {{0, 1, +1}, {1, 0, +1}};
  const auto po = graphos::pushout(a, b, vertex_ids, /*deduplicate=*/true);
  po.complex.validate();

  CHECK(po.complex.count(0) == 4);
  CHECK(po.complex.count(1) == 5);
  CHECK(po.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(po.complex));
  CHECK(graphos::euler_characteristic(po.complex) == 1);  // still a disk

  // every cell of A survives with its orientation
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < po.a_map.index[k].size(); ++i) {
      CHECK(po.a_map.index[k][i] != graphos::invalid_index);
      CHECK(po.a_map.sign[k][i] == 1);
    }
  }
  // B's e0 is identified with A's through a flip; B's other cells survive
  CHECK(po.b_map.index[1][0] == po.a_map.index[1][0]);
  CHECK(po.b_map.sign[1][0] == -1);
  CHECK(po.b_map.index[0][0] == po.a_map.index[0][1]);
  CHECK(po.b_map.index[0][1] == po.a_map.index[0][0]);
  CHECK(po.b_map.index[0][2] != graphos::invalid_index);
  CHECK(po.b_map.index[2][0] != graphos::invalid_index);
  CHECK(po.b_map.index[2][0] != po.a_map.index[2][0]);

  // B's 2-cell references the shared 1-cell with the flipped incidence
  // number
  const graphos::Index shared_edge = po.a_map.index[1][0];
  const graphos::Index b_face = po.b_map.index[2][0];
  const graphos::BoundaryOperator& f = po.complex.boundary(2);
  bool found = false;
  for (graphos::Index m = f.offsets[b_face]; m < f.offsets[b_face + 1]; ++m) {
    if (f.indices[m] == shared_edge) {
      found = true;
      CHECK(f.signs[m] == -1);
    }
  }
  CHECK(found);
}

// Witnesses the un-deduplicated pushout: both copies of the shared 1-cell
// survive as parallel cells over one vertex pair.
GRAPHOS_TEST(pushout_without_dedupe_keeps_parallel_cells) {
  const graphos::Complex a = graphos_test::make_triangle();
  const graphos::Complex b = graphos_test::make_triangle();
  const std::vector<graphos::Identification> vertex_ids = {{0, 1, +1}, {1, 0, +1}};
  const auto po = graphos::pushout(a, b, vertex_ids, /*deduplicate=*/false);
  po.complex.validate();
  CHECK(po.complex.count(0) == 4);
  CHECK(po.complex.count(1) == 6);
  CHECK(po.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(po.complex));
}

GRAPHOS_TEST_MAIN()
