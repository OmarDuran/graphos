#include <vector>

#include "fixtures.hpp"
#include "graphos/core/incidence.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

std::vector<Index> row(const graphos::Adjacency& a, Index i) {
  return {a.indices.begin() + a.offsets[i], a.indices.begin() + a.offsets[i + 1]};
}

}  // namespace

GRAPHOS_TEST(downward_closure_levels) {
  const graphos::Complex c = graphos_test::make_triangle();
  const auto f2v = graphos::incidence(c, 2, 0);  // face -> vertices
  CHECK(row(f2v, 0) == (std::vector<Index>{0, 1, 2}));
  const auto f2e = graphos::incidence(c, 2, 1);  // face -> edges
  CHECK(row(f2e, 0) == (std::vector<Index>{0, 1, 2}));
  const auto e2v = graphos::incidence(c, 1, 0);  // one level: matches ∂ pattern
  CHECK(row(e2v, 0) == (std::vector<Index>{0, 1}));
  CHECK(row(e2v, 2) == (std::vector<Index>{0, 2}));
}

GRAPHOS_TEST(upward_star_levels) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto v2f = graphos::incidence(c, 0, 2);      // vertex -> faces
  CHECK(row(v2f, 0) == (std::vector<Index>{0, 1}));  // shared vertex: both
  CHECK(row(v2f, 2) == (std::vector<Index>{0}));     // A's apex
  CHECK(row(v2f, 3) == (std::vector<Index>{1}));     // B's apex
}

GRAPHOS_TEST(same_dimension_is_identity) {
  const graphos::Complex c = graphos_test::make_triangle();
  const auto id = graphos::incidence(c, 1, 1);
  for (Index e = 0; e < 3; ++e) CHECK(row(id, e) == (std::vector<Index>{e}));
}

GRAPHOS_TEST(fan_cell_to_vertices) {
  const graphos::Complex c = graphos_test::make_fan();
  const auto f2v = graphos::incidence(c, 2, 0);
  CHECK(row(f2v, 0) == (std::vector<Index>{0, 1, 4}));  // T0
  CHECK(row(f2v, 3) == (std::vector<Index>{0, 3, 4}));  // T3
}

GRAPHOS_TEST(rejects_out_of_range_dimensions) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK_THROWS(graphos::incidence(c, 3, 0));
  CHECK_THROWS(graphos::incidence(c, 0, -1));
}

GRAPHOS_TEST_MAIN()
