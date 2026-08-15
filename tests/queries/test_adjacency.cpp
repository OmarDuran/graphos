#include <vector>

#include "fixtures.hpp"
#include "graphos/queries/adjacency.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

std::vector<Index> row(const graphos::Adjacency& a, Index i) {
  return {a.indices.begin() + a.offsets[i], a.indices.begin() + a.offsets[i + 1]};
}

}  // namespace

GRAPHOS_TEST(facet_adjacency_is_the_dual_graph) {
  const graphos::Complex fan = graphos_test::make_fan();
  const graphos::Adjacency adj = graphos::adjacency(fan, 2, 1);
  // T0 shares spokes with T1 (s1) and T3 (s0); boundary edges connect nothing
  CHECK(row(adj, 0) == (std::vector<Index>{1, 3}));
  CHECK(row(adj, 1) == (std::vector<Index>{0, 2}));
  CHECK(row(adj, 2) == (std::vector<Index>{1, 3}));
  CHECK(row(adj, 3) == (std::vector<Index>{0, 2}));
}

GRAPHOS_TEST(node_adjacency_is_coarser_than_facet_adjacency) {
  const graphos::Complex fan = graphos_test::make_fan();
  const graphos::Adjacency adj = graphos::adjacency(fan, 2, 0);
  // every triangle shares the center vertex with every other
  for (Index t = 0; t < 4; ++t) CHECK(row(adj, t).size() == 3);
}

GRAPHOS_TEST(vertex_adjacency_through_edges) {
  // path 0-1-2: the middle vertex neighbors both ends
  graphos::Complex c(1);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  const graphos::Adjacency adj = graphos::adjacency(c, 0, 1);
  CHECK(row(adj, 0) == (std::vector<Index>{1}));
  CHECK(row(adj, 1) == (std::vector<Index>{0, 2}));
  CHECK(row(adj, 2) == (std::vector<Index>{1}));
}

GRAPHOS_TEST(rejects_bad_dimensions) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK_THROWS(graphos::adjacency(c, 2, 2));
  CHECK_THROWS(graphos::adjacency(c, 3, 0));
}

GRAPHOS_TEST_MAIN()
