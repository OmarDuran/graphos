#include "graphos/core/coboundary.hpp"

#include <vector>

#include "fixtures.hpp"
#include "graphos_test.hpp"

namespace {

std::vector<graphos::Index> row(const graphos::CoboundaryOperator& a, graphos::Index i) {
  return {a.indices.begin() + a.offsets[i], a.indices.begin() + a.offsets[i + 1]};
}

std::vector<int> row_signs(const graphos::CoboundaryOperator& a, graphos::Index i) {
  return {a.signs.begin() + a.offsets[i], a.signs.begin() + a.offsets[i + 1]};
}

}  // namespace

GRAPHOS_TEST(vertex_to_edge_cofaces_of_triangle) {
  const graphos::Complex c = graphos_test::make_triangle();
  const graphos::CoboundaryOperator a = graphos::coboundary(c, 0);
  CHECK(row(a, 0) == (std::vector<graphos::Index>{0, 2}));
  CHECK(row(a, 1) == (std::vector<graphos::Index>{0, 1}));
  CHECK(row(a, 2) == (std::vector<graphos::Index>{1, 2}));
}

GRAPHOS_TEST(coboundary_carries_the_transposed_signs) {
  // edges [a,b] enter with -1 at a and +1 at b, so δ_0 rows mirror that:
  // v0 is the tail of e0 (-1) and the head of e2 (+1)
  const graphos::Complex c = graphos_test::make_triangle();
  const graphos::CoboundaryOperator a = graphos::coboundary(c, 0);
  CHECK(row_signs(a, 0) == (std::vector<int>{-1, +1}));
  CHECK(row_signs(a, 1) == (std::vector<int>{+1, -1}));
  CHECK(row_signs(a, 2) == (std::vector<int>{+1, -1}));
}

GRAPHOS_TEST(edge_to_face_cofaces_of_two_triangle_disk) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const graphos::CoboundaryOperator a = graphos::coboundary(c, 1);
  CHECK(row(a, 0) == (std::vector<graphos::Index>{0, 1}));  // shared edge: both faces
  CHECK(row_signs(a, 0) == (std::vector<int>{+1, +1}));
  CHECK(row(a, 1) == (std::vector<graphos::Index>{0}));
  CHECK(row(a, 4) == (std::vector<graphos::Index>{1}));
}

GRAPHOS_TEST(rejects_top_dimension) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK_THROWS(graphos::coboundary(c, 2));
  CHECK_THROWS(graphos::coboundary(c, -1));
}

GRAPHOS_TEST_MAIN()
