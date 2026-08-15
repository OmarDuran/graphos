// Witnesses that agglomeration realizes the multilevel hierarchy on three
// levels of tensor-refined complexes (products of path complexes), in
// dimension 2 and 3:
//
//  - agglomerating the fine level by parent labels reproduces the middle
//    level's top stratum — counts, χ, ∂∘∂ = 0, boundary arity — over the
//    refined middle skeleton;
//  - every surviving interface facet joins the aggregates that are exactly
//    the parents of its fine cofaces;
//  - the composite of the level-by-level chain maps equals the one-step map,
//    entry by entry.

#include <algorithm>
#include <vector>

#include "graphos/graphos.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

graphos::Complex path(int n) {
  std::vector<std::vector<Index>> segs;
  for (Index i = 0; i < n; ++i) segs.push_back({i, i + 1});
  return graphos::from_edges(n + 1, segs);
}

graphos::Complex quad_mesh(int n) {
  const graphos::Complex p = path(n);
  return graphos::product(p, p).complex;
}

graphos::Complex hex_mesh(int n) {
  const graphos::Complex p = path(n);
  return graphos::product(quad_mesh(n), p).complex;
}

// the product layout is i-major, so a 2-cell (i,j) has index in+j and a
// 3-cell (i,j,k) has (in+j)n+k; the parent halves each coordinate
std::vector<Index> quad_parents(int n) {
  std::vector<Index> labels(static_cast<std::size_t>(n) * n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      labels[static_cast<std::size_t>(i * n + j)] = static_cast<Index>((i / 2) * (n / 2) + (j / 2));
  return labels;
}

std::vector<Index> hex_parents(int n) {
  std::vector<Index> labels(static_cast<std::size_t>(n) * n * n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        labels[static_cast<std::size_t>((i * n + j) * n + k)] =
            static_cast<Index>(((i / 2) * (n / 2) + (j / 2)) * (n / 2) + (k / 2));
  return labels;
}

void check_level(const graphos::Complex& c) {
  c.validate();
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);
}

// each surviving facet's aggregated cofaces are exactly the parents of its
// fine cofaces
void check_seams(const graphos::Complex& fine, const graphos::AgglomerationResult& agg,
                 const std::vector<Index>& parents) {
  const int n = fine.dim();
  const graphos::CoboundaryOperator cob_f = graphos::coboundary(fine, n - 1);
  const graphos::CoboundaryOperator cob_a = graphos::coboundary(agg.complex, n - 1);
  for (Index f = 0; f < fine.count(n - 1); ++f) {
    const Index fa = agg.map.index[static_cast<std::size_t>(n) - 1][static_cast<std::size_t>(f)];
    if (fa == graphos::invalid_index) continue;
    std::vector<Index> from_fine;
    for (Index m = cob_f.offsets[static_cast<std::size_t>(f)];
         m < cob_f.offsets[static_cast<std::size_t>(f) + 1]; ++m) {
      from_fine.push_back(
          parents[static_cast<std::size_t>(cob_f.indices[static_cast<std::size_t>(m)])]);
    }
    std::sort(from_fine.begin(), from_fine.end());
    from_fine.erase(std::unique(from_fine.begin(), from_fine.end()), from_fine.end());
    std::vector<Index> from_agg;
    for (Index m = cob_a.offsets[static_cast<std::size_t>(fa)];
         m < cob_a.offsets[static_cast<std::size_t>(fa) + 1]; ++m) {
      from_agg.push_back(cob_a.indices[static_cast<std::size_t>(m)]);
    }
    std::sort(from_agg.begin(), from_agg.end());
    from_agg.erase(std::unique(from_agg.begin(), from_agg.end()), from_agg.end());
    CHECK(from_fine == from_agg);
  }
}

void check_maps_equal(const graphos::ChainMap& a, const graphos::ChainMap& b) {
  CHECK(a.index.size() == b.index.size());
  for (std::size_t k = 0; k < a.index.size(); ++k) CHECK(a.index[k] == b.index[k]);
}

}  // namespace

GRAPHOS_TEST(quad_hierarchy_three_levels) {
  const graphos::Complex coarse = quad_mesh(1);
  const graphos::Complex middle = quad_mesh(2);
  const graphos::Complex fine = quad_mesh(4);
  check_level(coarse);
  check_level(middle);
  check_level(fine);
  CHECK(coarse.count(2) == 1);
  CHECK(middle.count(2) == 4);
  CHECK(fine.count(2) == 16);

  // fine → middle: four aggregates over the refined middle skeleton
  const auto parents = quad_parents(4);
  const auto a1 = graphos::agglomerate(fine, parents);
  check_level(a1.complex);
  CHECK(a1.complex.count(2) == middle.count(2));
  CHECK(a1.complex.count(1) == 24);  // 12 middle edges, each two fine edges
  CHECK(a1.complex.count(0) == 21);  // 25 fine vertices minus 4 cell centers
  // each aggregated 2-cell is bounded by eight fine 1-cells, two per side
  const graphos::BoundaryOperator& b1 = a1.complex.boundary(2);
  for (Index e = 0; e < a1.complex.count(2); ++e) {
    CHECK(b1.offsets[e + 1] - b1.offsets[e] == 8);
  }
  check_seams(fine, a1, parents);

  // middle → coarse: the four aggregates merge into one
  const std::vector<Index> to_one(4, 0);
  const auto a2 = graphos::agglomerate(a1.complex, to_one);
  check_level(a2.complex);
  CHECK(a2.complex.count(2) == 1);
  const graphos::BoundaryOperator& b2 = a2.complex.boundary(2);
  CHECK(b2.offsets[1] - b2.offsets[0] == 16);  // the fine domain boundary

  // level-by-level equals one step: the chain maps agree
  const std::vector<Index> all_zero(16, 0);
  const auto direct = graphos::agglomerate(fine, all_zero);
  CHECK(direct.complex.count(0) == a2.complex.count(0));
  CHECK(direct.complex.count(1) == a2.complex.count(1));
  CHECK(direct.complex.count(2) == a2.complex.count(2));
  check_maps_equal(graphos::compose(a1.map, a2.map), direct.map);
}

GRAPHOS_TEST(hex_hierarchy_three_levels) {
  const graphos::Complex coarse = hex_mesh(1);
  const graphos::Complex middle = hex_mesh(2);
  const graphos::Complex fine = hex_mesh(4);
  check_level(coarse);
  check_level(middle);
  check_level(fine);
  CHECK(coarse.count(3) == 1);
  CHECK(coarse.count(2) == 6);
  CHECK(middle.count(3) == 8);
  CHECK(fine.count(3) == 64);
  CHECK(fine.count(0) == 125);
  CHECK(fine.count(1) == 300);
  CHECK(fine.count(2) == 240);

  // fine → middle: eight aggregates over the refined middle skeleton
  const auto parents = hex_parents(4);
  const auto a1 = graphos::agglomerate(fine, parents);
  check_level(a1.complex);
  CHECK(a1.complex.count(3) == middle.count(3));
  CHECK(a1.complex.count(2) == 144);  // 36 middle faces, each four fine faces
  CHECK(a1.complex.count(1) == 252);  // refined middle skeleton edges
  CHECK(a1.complex.count(0) == 117);  // 125 minus 8 cell centers
  // each aggregated 3-cell is bounded by 24 fine 2-cells, four per side
  const graphos::BoundaryOperator& b1 = a1.complex.boundary(3);
  for (Index e = 0; e < a1.complex.count(3); ++e) {
    CHECK(b1.offsets[e + 1] - b1.offsets[e] == 24);
  }
  check_seams(fine, a1, parents);

  // → coarse, and level-by-level equals one step
  const std::vector<Index> to_one(8, 0);
  const auto a2 = graphos::agglomerate(a1.complex, to_one);
  check_level(a2.complex);
  CHECK(a2.complex.count(3) == 1);
  const graphos::BoundaryOperator& b2 = a2.complex.boundary(3);
  CHECK(b2.offsets[1] - b2.offsets[0] == 96);  // 6 sides x 16 fine faces

  const std::vector<Index> all_zero(64, 0);
  const auto direct = graphos::agglomerate(fine, all_zero);
  for (int k = 0; k <= 3; ++k) CHECK(direct.complex.count(k) == a2.complex.count(k));
  check_maps_equal(graphos::compose(a1.map, a2.map), direct.map);
}

GRAPHOS_TEST_MAIN()
