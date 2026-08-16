#include <iterator>
#include <map>

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

// THE MAP IS A LABELLING, NOT A CHAIN MAP, and the direction is why.
//
// A coarse cell is a CHAIN of fine cells, so the aggregation ι runs
// coarse → fine and is multi-valued; ChainMap is single-valued and carries its
// transpose. That transpose does not commute with ∂: for one member σ of
// aggregate A, ∂(πσ) is the boundary of ALL of A, while π(∂σ) is σ's share.
GRAPHOS_TEST(the_labelling_is_not_a_chain_map_but_the_aggregation_is) {
  const auto oriented = graphos::orient(graphos_test::make_two_triangle_disk());
  const graphos::Complex& fine = oriented.complex;
  const auto ag = graphos::agglomerate(fine, {0, 0});
  CHECK(!graphos::commutes_with_boundary(fine, ag.complex, ag.map));

  // a 2-cell fibre, so no single-valued map inverts it
  int members = 0;
  for (graphos::Index i = 0; i < fine.count(2); ++i) {
    if (ag.map.index[2][static_cast<std::size_t>(i)] != graphos::invalid_index) ++members;
  }
  CHECK(members == 2);

  // What does hold, and what a coarse operator rests on: ∂(Σ members) is the
  // coarse boundary, the interior facet cancelling. The chain-map condition
  // for ι, stated on chains because ι is not a map of generators.
  std::map<graphos::Index, int> lhs;  // Σ ∂(members), in fine 1-cells
  const graphos::BoundaryOperator& bf = fine.boundary(2);
  for (graphos::Index i = 0; i < fine.count(2); ++i) {
    const int sg = ag.map.sign[2][static_cast<std::size_t>(i)];
    for (graphos::Index m = bf.offsets[i]; m < bf.offsets[i + 1]; ++m) {
      lhs[bf.indices[m]] += sg * bf.signs[m];
    }
  }
  for (auto it = lhs.begin(); it != lhs.end();) {
    it = (it->second == 0) ? lhs.erase(it) : std::next(it);  // the interior cancels
  }

  std::vector<graphos::Index> fine_of(static_cast<std::size_t>(ag.complex.count(1)),
                                      graphos::invalid_index);
  for (graphos::Index e = 0; e < fine.count(1); ++e) {
    const graphos::Index t = ag.map.index[1][static_cast<std::size_t>(e)];
    if (t != graphos::invalid_index) fine_of[static_cast<std::size_t>(t)] = e;
  }
  std::map<graphos::Index, int> rhs;  // ∂(A), pulled back
  const graphos::BoundaryOperator& bc = ag.complex.boundary(2);
  for (graphos::Index m = bc.offsets[0]; m < bc.offsets[1]; ++m) {
    const graphos::Index e = fine_of[static_cast<std::size_t>(bc.indices[m])];
    CHECK(e != graphos::invalid_index);
    rhs[e] += bc.signs[m] * ag.map.sign[1][static_cast<std::size_t>(e)];
  }
  CHECK(lhs == rhs);
}

GRAPHOS_TEST_MAIN()
