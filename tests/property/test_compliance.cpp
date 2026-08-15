// Compliance: each case restates a claim the README makes about the surface
// and checks the code satisfies it, recomputing the expected answer from the
// stated definition rather than from the implementation under test.
//
// Cross-cutting by design, so it lives beside the other property suites rather
// than mirroring a single header.

#include <set>
#include <utility>
#include <vector>

#include "fixtures.hpp"
#include "graphos/graphos.hpp"
#include "graphos_test.hpp"

namespace {

using graphos::Complex;
using graphos::Index;
using graphos::Marker;

using CellSet = std::set<std::pair<int, Index>>;

CellSet to_set(const Complex& c, const Marker& m) {
  CellSet out;
  for (int k = 0; k <= c.dim(); ++k) {
    for (Index i = 0; i < c.count(k); ++i) {
      if (m.marked(k, i)) out.insert({k, i});
    }
  }
  return out;
}

Marker of_set(const Complex& c, const CellSet& s) {
  Marker m(c);
  for (const auto& [k, i] : s) m.mark(k, i);
  return m;
}

Marker one(const Complex& c, int k, Index cell) {
  Marker m(c);
  m.mark(k, cell);
  return m;
}

CellSet difference(const CellSet& a, const CellSet& b) {
  CellSet out;
  for (const auto& e : a) {
    if (b.find(e) == b.end()) out.insert(e);
  }
  return out;
}

// the j-slice of a marked set
std::vector<Index> slice(const CellSet& s, int j) {
  std::vector<Index> out;
  for (const auto& [k, i] : s) {
    if (k == j) out.push_back(i);
  }
  return out;
}

std::vector<Complex> corpus() {
  std::vector<Complex> out;
  out.push_back(graphos_test::make_triangle());
  out.push_back(graphos_test::make_two_triangle_disk());
  out.push_back(graphos_test::make_fan());
  out.push_back(graphos::from_simplices(3, 4, {{0, 1, 2, 3}}));
  return out;
}

}  // namespace

// README: "lk(S) = cl(st S) ∖ st(cl S)"
GRAPHOS_TEST(link_is_closed_star_minus_star_closure) {
  for (const Complex& c : corpus()) {
    for (int k = 0; k <= c.dim(); ++k) {
      for (Index cell = 0; cell < c.count(k); ++cell) {
        const Marker s = one(c, k, cell);
        const CellSet cl_st = to_set(c, closure_of(c, of_set(c, to_set(c, star_of(c, s)))));
        const CellSet st_cl = to_set(c, star_of(c, of_set(c, to_set(c, closure_of(c, s)))));
        CHECK(to_set(c, link_of(c, s)) == difference(cl_st, st_cl));
      }
    }
  }
}

// README: "the frontier cl(st S) ∖ st(S)"
GRAPHOS_TEST(frontier_is_closed_star_minus_open_star) {
  for (const Complex& c : corpus()) {
    for (int k = 0; k <= c.dim(); ++k) {
      for (Index cell = 0; cell < c.count(k); ++cell) {
        const Marker s = one(c, k, cell);
        const CellSet st = to_set(c, star_of(c, s));
        const CellSet cl_st = to_set(c, closure_of(c, of_set(c, st)));
        CHECK(to_set(c, frontier_of(c, s)) == difference(cl_st, st));
      }
    }
  }
}

// README: "per k-cell σ, the j-cells of cl(σ) (j < k), of st(σ) (j > k), or σ
// itself"
GRAPHOS_TEST(incidence_is_closure_below_and_star_above) {
  for (const Complex& c : corpus()) {
    for (int k = 0; k <= c.dim(); ++k) {
      for (int j = 0; j <= c.dim(); ++j) {
        const graphos::Adjacency inc = incidence(c, k, j);
        for (Index sigma = 0; sigma < c.count(k); ++sigma) {
          std::vector<Index> got(
              inc.indices.begin() + inc.offsets[static_cast<std::size_t>(sigma)],
              inc.indices.begin() + inc.offsets[static_cast<std::size_t>(sigma) + 1]);
          std::sort(got.begin(), got.end());

          const Marker s = one(c, k, sigma);
          std::vector<Index> want;
          if (j < k) {
            want = slice(to_set(c, closure_of(c, s)), j);
          } else if (j > k) {
            want = slice(to_set(c, star_of(c, s)), j);
          } else {
            want = {sigma};
          }
          std::sort(want.begin(), want.end());
          CHECK(got == want);
        }
      }
    }
  }
}

// README: "free_faces(c) — τ with one proper coface, occurring once in ∂σ"
GRAPHOS_TEST(free_faces_is_exactly_the_whitehead_condition) {
  for (const Complex& c : corpus()) {
    CellSet want;
    for (int k = 0; k < c.dim(); ++k) {
      const graphos::CoboundaryOperator cob = coboundary(c, k);
      for (Index tau = 0; tau < c.count(k); ++tau) {
        const Index lo = cob.offsets[static_cast<std::size_t>(tau)];
        const Index hi = cob.offsets[static_cast<std::size_t>(tau) + 1];
        if (hi - lo != 1) continue;  // needs exactly one proper coface
        const Index sigma = cob.indices[static_cast<std::size_t>(lo)];
        const graphos::BoundaryOperator& bnd = c.boundary(k + 1);
        int occurrences = 0;
        for (Index m = bnd.offsets[static_cast<std::size_t>(sigma)];
             m < bnd.offsets[static_cast<std::size_t>(sigma) + 1]; ++m) {
          if (bnd.indices[static_cast<std::size_t>(m)] == tau) ++occurrences;
        }
        if (occurrences == 1) want.insert({k, tau});  // and occurring once
      }
    }
    CHECK(to_set(c, free_faces(c)) == want);
  }
}

// README: "partitions the facets by top-coface count: maximal (0), boundary
// (1), interior (2), nonmanifold (3+)"
GRAPHOS_TEST(facet_classes_partition_by_coface_count) {
  for (const Complex& c : corpus()) {
    const int n = c.dim();
    const graphos::FacetClassification fc = classify_facets(c);
    const graphos::CoboundaryOperator cob = coboundary(c, n - 1);
    Index counted = 0;
    for (Index f = 0; f < c.count(n - 1); ++f) {
      const Index degree =
          cob.offsets[static_cast<std::size_t>(f) + 1] - cob.offsets[static_cast<std::size_t>(f)];
      const int hits = static_cast<int>(fc.maximal.marked(n - 1, f)) +
                       static_cast<int>(fc.boundary.marked(n - 1, f)) +
                       static_cast<int>(fc.interior.marked(n - 1, f)) +
                       static_cast<int>(fc.nonmanifold.marked(n - 1, f));
      CHECK(hits == 1);  // exactly one class
      CHECK(fc.maximal.marked(n - 1, f) == (degree == 0));
      CHECK(fc.boundary.marked(n - 1, f) == (degree == 1));
      CHECK(fc.interior.marked(n - 1, f) == (degree == 2));
      CHECK(fc.nonmanifold.marked(n - 1, f) == (degree >= 3));
      ++counted;
    }
    CHECK(fc.maximal.marked_count(n - 1) + fc.boundary.marked_count(n - 1) +
              fc.interior.marked_count(n - 1) + fc.nonmanifold.marked_count(n - 1) ==
          counted);
  }
}

// README: "Closedness ∂K = ∅" — the predicate is exactly "no facet with one
// coface"
GRAPHOS_TEST(is_closed_agrees_with_the_boundary_class) {
  for (const Complex& c : corpus()) {
    const graphos::FacetClassification fc = classify_facets(c);
    CHECK(is_closed(c) == (fc.boundary.marked_count(c.dim() - 1) == 0));
  }
  // a closed surface: ∂Δ³ has no boundary facet
  const Complex tet = graphos::from_simplices(3, 4, {{0, 1, 2, 3}});
  const auto shell = subcomplex(tet, classify_facets(tet).boundary);
  CHECK(is_closed(shell.complex));
}

// README: "connected_components(c) is β₀ of the whole complex"
GRAPHOS_TEST(connected_components_of_the_complex_is_betti_zero) {
  const Complex tri = graphos_test::make_triangle();
  const Complex two = disjoint_union(tri, graphos_test::make_fan()).complex;
  const Complex three = disjoint_union(two, graphos_test::make_segment()).complex;
  for (const Complex& c : {tri, two, three}) {
    CHECK(connected_components(c).count == betti_numbers_z2(c)[0]);
  }
}

// README: "k-cells whose boundary chains agree as signed sets up to a uniform
// flip"
GRAPHOS_TEST(parallel_cells_require_a_uniform_flip) {
  Complex c(2);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  c.attach_cell(2, {0, 1, 2}, {-1, -1, -1});  // the same chain, uniformly flipped
  const auto dups = find_parallel_cells(c, 2);
  CHECK(dups.size() == 1);
  if (dups.size() == 1) CHECK(dups[0].rel_sign == -1);

  // distinct boundary chains are not parallel
  const Complex disk = graphos_test::make_two_triangle_disk();
  CHECK(find_parallel_cells(disk, 2).empty());
}

// README: "orient(c, k) — coherently orients the maximal cells of stratum k";
// the classes are returned
GRAPHOS_TEST(orient_takes_a_stratum_and_returns_its_classes) {
  const Complex c = graphos_test::make_two_triangle_disk();
  for (int k = 1; k <= c.dim(); ++k) {
    const auto r = orient(c, k);
    CHECK(r.stratum == k);
    CHECK(r.class_of.size() == static_cast<std::size_t>(c.count(k)) || r.classes == 0);
  }
  // the no-stratum form is the top stratum
  const auto top = orient(c);
  CHECK(top.stratum == c.dim());
  CHECK(top.orientable == orient(c, c.dim()).orientable);
}

GRAPHOS_TEST_MAIN()
