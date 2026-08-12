#include "graphos/core/coboundary.hpp"

#include <vector>

#include "fixtures.hpp"
#include "graphos_test.hpp"

namespace {

std::vector<graphos::Index> row(const graphos::Adjacency& a, graphos::Index i) {
  return {a.indices.begin() + a.offsets[i], a.indices.begin() + a.offsets[i + 1]};
}

}  // namespace

GRAPHOS_TEST(vertex_to_edge_cofaces_of_triangle) {
  const graphos::Complex c = graphos_test::make_triangle();
  const graphos::Adjacency a = graphos::coboundary(c, 0);
  CHECK(row(a, 0) == (std::vector<graphos::Index>{0, 2}));
  CHECK(row(a, 1) == (std::vector<graphos::Index>{0, 1}));
  CHECK(row(a, 2) == (std::vector<graphos::Index>{1, 2}));
}

GRAPHOS_TEST(edge_to_face_cofaces_of_two_triangle_disk) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const graphos::Adjacency a = graphos::coboundary(c, 1);
  CHECK(row(a, 0) == (std::vector<graphos::Index>{0, 1}));  // shared edge: both faces
  CHECK(row(a, 1) == (std::vector<graphos::Index>{0}));
  CHECK(row(a, 4) == (std::vector<graphos::Index>{1}));
}

GRAPHOS_TEST(rejects_top_dimension) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK_THROWS(graphos::coboundary(c, 2));
  CHECK_THROWS(graphos::coboundary(c, -1));
}

GRAPHOS_TEST_MAIN()
