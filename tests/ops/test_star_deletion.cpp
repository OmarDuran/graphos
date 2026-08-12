#include "graphos/ops/star_deletion.hpp"

#include "graphos/ops/pushout.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

// Deleting the star of one apex of a two-triangle disk sweeps out its two
// private edges and its 2-cell, leaving a single triangle.
GRAPHOS_TEST(apex_star_deletion_leaves_one_triangle) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker apex(c);
  apex.mark(0, 3);  // B's apex
  const auto sd = graphos::star_deletion(c, apex);
  sd.complex.validate();

  CHECK(sd.complex.count(0) == 3);
  CHECK(sd.complex.count(1) == 3);
  CHECK(sd.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(sd.complex));
  CHECK(graphos::euler_characteristic(sd.complex) == 1);

  // casualties: the apex, its two edges, B's face; everything else survives
  CHECK(sd.map.index[0][3] == graphos::invalid_index);
  CHECK(sd.map.index[1][3] == graphos::invalid_index);
  CHECK(sd.map.index[1][4] == graphos::invalid_index);
  CHECK(sd.map.index[2][1] == graphos::invalid_index);
  for (graphos::Index i = 0; i < 3; ++i) {
    CHECK(sd.map.index[0][static_cast<std::size_t>(i)] == i);
    CHECK(sd.map.index[1][static_cast<std::size_t>(i)] == i);
  }
  CHECK(sd.map.index[2][0] == 0);
}

// Same deletion built through a pushout: chain maps must compose so that all
// of A survives and B's share of the interface persists with its flip.
GRAPHOS_TEST(chain_maps_compose_across_pushout_and_deletion) {
  const graphos::Complex a = graphos_test::make_triangle();
  const graphos::Complex b = graphos_test::make_triangle();
  const std::vector<graphos::Identification> vertex_ids = {{0, 1, +1}, {1, 0, +1}};
  const auto po = graphos::pushout(a, b, vertex_ids, /*deduplicate=*/true);

  graphos::Marker apex(po.complex);
  apex.mark(0, po.b_map.index[0][2]);  // B's apex in the pushout
  const auto sd = graphos::star_deletion(po.complex, apex);
  sd.complex.validate();
  CHECK(graphos::d_squared_is_zero(sd.complex));

  const auto a_final = graphos::compose(po.a_map, sd.map);
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < a_final.index[k].size(); ++i) {
      CHECK(a_final.index[k][i] != graphos::invalid_index);
    }
  }
  const auto b_final = graphos::compose(po.b_map, sd.map);
  CHECK(b_final.index[0][2] == graphos::invalid_index);
  CHECK(b_final.index[1][1] == graphos::invalid_index);
  CHECK(b_final.index[1][2] == graphos::invalid_index);
  CHECK(b_final.index[2][0] == graphos::invalid_index);
  // ...but B's share of the interface persists, with the flip preserved
  CHECK(b_final.index[0][0] != graphos::invalid_index);
  CHECK(b_final.index[1][0] != graphos::invalid_index);
  CHECK(b_final.sign[1][0] == -1);
}

GRAPHOS_TEST(deleting_a_maximal_lower_cell_keeps_its_faces) {
  // mixed-dimensional: a bare segment; deleting the edge keeps its vertices
  const graphos::Complex c = graphos_test::make_segment();
  const auto sd = graphos::star_deletion(c, graphos::Marker::from_cells(c, {{}, {0}}));
  CHECK(sd.complex.count(0) == 2);  // orphan pruning is a caller decision
  CHECK(sd.complex.count(1) == 0);
}

GRAPHOS_TEST_MAIN()
