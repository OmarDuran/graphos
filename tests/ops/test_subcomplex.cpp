#include "fixtures.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos_test.hpp"

using graphos::Index;

// Witnesses that subcomplex extracts cl(S): restricting the disk to one 2-cell
// yields that cell with its closure.
GRAPHOS_TEST(one_face_restriction_is_a_triangle) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker m(c);
  m.mark(2, 0);  // face A
  const auto sub = graphos::subcomplex(c, m);
  sub.complex.validate();

  CHECK(sub.complex.count(0) == 3);
  CHECK(sub.complex.count(1) == 3);
  CHECK(sub.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(sub.complex));
  CHECK(graphos::euler_characteristic(sub.complex) == 1);

  // B's private cells lie outside cl(S)
  CHECK(sub.map.index[0][3] == graphos::invalid_index);
  CHECK(sub.map.index[1][3] == graphos::invalid_index);
  CHECK(sub.map.index[1][4] == graphos::invalid_index);
  CHECK(sub.map.index[2][1] == graphos::invalid_index);
  // the embedding inverts the map on survivors
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < sub.embedding.index[k].size(); ++i) {
      const Index parent = sub.embedding.index[k][i];
      CHECK(sub.map.index[k][static_cast<std::size_t>(parent)] == static_cast<Index>(i));
    }
  }
}

// Witnesses skeleton extraction: marking every 1-cell yields the 1-skeleton,
// a circle, with χ = 0.
GRAPHOS_TEST(one_skeleton_of_a_triangle_is_a_circle) {
  const graphos::Complex c = graphos_test::make_triangle();
  graphos::Marker m(c);
  m.mark_where(1, [](Index) { return true; });
  const auto sub = graphos::subcomplex(c, m);
  sub.complex.validate();
  CHECK(sub.complex.count(0) == 3);
  CHECK(sub.complex.count(1) == 3);
  CHECK(sub.complex.count(2) == 0);
  CHECK(graphos::euler_characteristic(sub.complex) == 0);  // circle
}

// Witnesses that cut and subcomplex compose: cutting the disk along the shared
// 1-cell and extracting the detached interface gives it as its own complex.
GRAPHOS_TEST(extracts_fracture_domain_after_a_cut) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cut = graphos::cut_along(c, interface);

  // the original interface survives the cut at its old index
  graphos::Marker fracture(cut.complex);
  fracture.mark(1, 0);
  const auto sub = graphos::subcomplex(cut.complex, fracture);
  sub.complex.validate();

  // a segment: the interface 1-cell with its two original endpoints
  CHECK(sub.complex.count(0) == 2);
  CHECK(sub.complex.count(1) == 1);
  CHECK(sub.complex.count(2) == 0);
  CHECK(graphos::euler_characteristic(sub.complex) == 1);
  // embedded at the parent's original cells
  CHECK(sub.embedding.index[1][0] == 0);
  CHECK(sub.embedding.index[0][0] == 0);
  CHECK(sub.embedding.index[0][1] == 1);
}

GRAPHOS_TEST(empty_marker_yields_empty_complex) {
  const graphos::Complex c = graphos_test::make_triangle();
  const auto sub = graphos::subcomplex(c, graphos::Marker(c));
  CHECK(sub.complex.count(0) == 0);
  CHECK(sub.complex.count(1) == 0);
  CHECK(sub.complex.count(2) == 0);
}

GRAPHOS_TEST(rejects_marker_for_wrong_complex) {
  const graphos::Complex tri = graphos_test::make_triangle();
  const graphos::Complex fan = graphos_test::make_fan();
  CHECK_THROWS(graphos::subcomplex(fan, graphos::Marker(tri)));
}

// -- subcomplex_from: cl(S) by descent -------------------------------------
//
// The specification is agreement with subcomplex() cell for cell, so it is
// asserted differentially over every S. A closure with its own numbering would
// be a second dialect.

namespace {

// S as a Marker
graphos::Marker marker_of(const graphos::Complex& c,
                          const std::vector<std::vector<graphos::Index>>& seeds) {
  graphos::Marker mk(c);
  for (std::size_t k = 0; k < seeds.size() && static_cast<int>(k) <= c.dim(); ++k) {
    for (const graphos::Index s : seeds[k]) mk.mark(static_cast<int>(k), s);
  }
  return mk;
}

// the two extractions of cl(S), compared on every observable
bool agrees(const graphos::Complex& c, const std::vector<std::vector<graphos::Index>>& seeds,
            graphos::SubcomplexWorkspace& ws) {
  const graphos::SubcomplexResult dense = graphos::subcomplex(c, marker_of(c, seeds));
  const graphos::SparseSubcomplexResult sparse = graphos::subcomplex_from(c, seeds, ws);

  if (dense.complex.dim() != sparse.complex.dim()) return false;
  for (int k = 0; k <= c.dim(); ++k) {
    if (dense.complex.count(k) != sparse.complex.count(k)) return false;
    // sub -> parent
    if (dense.embedding.index[static_cast<std::size_t>(k)] !=
        sparse.embedding.index[static_cast<std::size_t>(k)]) {
      return false;
    }
    // strictly increasing, which is what to_local() binary-searches
    const auto& e = sparse.embedding.index[static_cast<std::size_t>(k)];
    for (std::size_t i = 1; i < e.size(); ++i) {
      if (!(e[i - 1] < e[i])) return false;
    }
    // parent -> sub, against the dense ChainMap
    for (graphos::Index p = 0; p < c.count(k); ++p) {
      if (graphos::to_local(sparse, k, p) != dense.map.index[static_cast<std::size_t>(k)][
              static_cast<std::size_t>(p)]) {
        return false;
      }
    }
  }
  for (int k = 1; k <= c.dim(); ++k) {
    const graphos::BoundaryOperator& a = dense.complex.boundary(k);
    const graphos::BoundaryOperator& b = sparse.complex.boundary(k);
    if (a.offsets != b.offsets || a.indices != b.indices || a.signs != b.signs) return false;
  }
  return true;
}

// every S ⊆ stratum k
bool agrees_on_all_subsets(const graphos::Complex& c, int k) {
  graphos::SubcomplexWorkspace ws(c);
  const graphos::Index n = c.count(k);
  if (n > 12) return true;  // keep the sweep finite
  for (int mask = 0; mask < (1 << n); ++mask) {
    std::vector<std::vector<graphos::Index>> seeds(static_cast<std::size_t>(c.dim()) + 1);
    for (graphos::Index i = 0; i < n; ++i) {
      if (mask & (1 << i)) seeds[static_cast<std::size_t>(k)].push_back(i);
    }
    if (!agrees(c, seeds, ws)) return false;
  }
  return true;
}

}  // namespace

GRAPHOS_TEST(sparse_agrees_with_dense_on_every_subset) {
  const graphos::Complex tri = graphos_test::make_triangle();
  for (int k = 0; k <= tri.dim(); ++k) CHECK(agrees_on_all_subsets(tri, k));
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  for (int k = 0; k <= disk.dim(); ++k) CHECK(agrees_on_all_subsets(disk, k));
  const graphos::Complex fan = graphos_test::make_fan();
  for (int k = 0; k <= fan.dim(); ++k) CHECK(agrees_on_all_subsets(fan, k));
  const graphos::Complex seg = graphos_test::make_segment();
  for (int k = 0; k <= seg.dim(); ++k) CHECK(agrees_on_all_subsets(seg, k));
}

GRAPHOS_TEST(sparse_agrees_on_mixed_dimension_selections) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::SubcomplexWorkspace ws(disk);
  // S spanning several strata, as the Marker contract allows
  CHECK(agrees(disk, {{0}, {2}, {0}}, ws));
  CHECK(agrees(disk, {{}, {0, 3}, {1}}, ws));
  CHECK(agrees(disk, {{1, 2}, {}, {}}, ws));
  // the 1-skeleton
  std::vector<std::vector<graphos::Index>> skel(3);
  for (graphos::Index i = 0; i < disk.count(0); ++i) skel[0].push_back(i);
  for (graphos::Index i = 0; i < disk.count(1); ++i) skel[1].push_back(i);
  CHECK(agrees(disk, skel, ws));
}

GRAPHOS_TEST(sparse_tolerates_unsorted_and_repeated_seeds) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::SubcomplexWorkspace ws(disk);
  // the seed list denotes a set: order and multiplicity are inert
  const graphos::SparseSubcomplexResult a = graphos::subcomplex_from(disk, {{}, {}, {1, 0}}, ws);
  const graphos::SparseSubcomplexResult b =
      graphos::subcomplex_from(disk, {{}, {}, {0, 1, 1, 0, 1}}, ws);
  CHECK(a.embedding.index == b.embedding.index);
  for (int k = 1; k <= disk.dim(); ++k) {
    CHECK(a.complex.boundary(k).indices == b.complex.boundary(k).indices);
  }
  CHECK(agrees(disk, {{}, {}, {1, 0}}, ws));
}

// A slot left dirty corrupts the NEXT extraction, not this one, so the check
// is a repeated sequence rather than a single call.
GRAPHOS_TEST(sparse_workspace_is_reusable_across_extractions) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::SubcomplexWorkspace ws(disk);
  const std::vector<std::vector<graphos::Index>> a = {{}, {}, {0}};
  const std::vector<std::vector<graphos::Index>> b = {{}, {}, {1}};
  const std::vector<std::vector<graphos::Index>> both = {{}, {}, {0, 1}};

  const auto ref_a = graphos::subcomplex_from(disk, a, ws).embedding.index;
  // interleaved; each S must reproduce its own cl(S)
  for (int r = 0; r < 4; ++r) {
    CHECK(agrees(disk, a, ws));
    CHECK(agrees(disk, both, ws));
    CHECK(agrees(disk, b, ws));
    CHECK(graphos::subcomplex_from(disk, a, ws).embedding.index == ref_a);
  }
  // S = ∅ must leave the slots clean too
  CHECK(agrees(disk, {}, ws));
  CHECK(graphos::subcomplex_from(disk, a, ws).embedding.index == ref_a);
}

GRAPHOS_TEST(sparse_empty_selection_is_the_empty_complex) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::SubcomplexWorkspace ws(disk);
  const graphos::SparseSubcomplexResult r = graphos::subcomplex_from(disk, {}, ws);
  CHECK(r.complex.dim() == disk.dim());  // dimension is the parent's
  for (int k = 0; k <= disk.dim(); ++k) CHECK(r.complex.count(k) == 0);
  CHECK(graphos::to_local(r, 2, 0) == graphos::invalid_index);
  // and with explicit empty strata
  CHECK(agrees(disk, {{}, {}, {}}, ws));
}

GRAPHOS_TEST(sparse_rejects_a_foreign_workspace_and_bad_seeds) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  const graphos::Complex tri = graphos_test::make_triangle();
  graphos::SubcomplexWorkspace ws_tri(tri);
  CHECK_THROWS(graphos::subcomplex_from(disk, {{}, {}, {0}}, ws_tri));

  graphos::SubcomplexWorkspace ws(disk);
  CHECK_THROWS(graphos::subcomplex_from(disk, {{}, {}, {disk.count(2)}}, ws));
  CHECK_THROWS(graphos::subcomplex_from(disk, {{}, {}, {-1}}, ws));
  // the workspace survives a rejected call
  CHECK(agrees(disk, {{}, {}, {0}}, ws));
}

GRAPHOS_TEST(sparse_inherits_the_chain_complex_law) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::SubcomplexWorkspace ws(disk);
  const graphos::SparseSubcomplexResult r = graphos::subcomplex_from(disk, {{}, {}, {0}}, ws);
  r.complex.validate();
  CHECK(graphos::d_squared_is_zero(r.complex));
  // cl(S) is closed: every face reference lands inside it
  for (int k = 1; k <= r.complex.dim(); ++k) {
    const graphos::BoundaryOperator& b = r.complex.boundary(k);
    for (graphos::Index e = 0; e < r.complex.count(k); ++e) {
      for (graphos::Index m = b.offsets[e]; m < b.offsets[e + 1]; ++m) {
        CHECK(b.indices[m] >= 0 && b.indices[m] < r.complex.count(k - 1));
      }
    }
  }
}

// The inclusion cl(S) ↪ C is a chain map: ∂ is inherited, so it commutes.
// Both extractions must satisfy it, being the same inclusion.
GRAPHOS_TEST(the_embedding_is_a_chain_map) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::Marker mk(disk);
  mk.mark(2, 0);
  const graphos::SubcomplexResult dense = graphos::subcomplex(disk, mk);
  CHECK(graphos::commutes_with_boundary(dense.complex, disk, dense.embedding));

  graphos::SubcomplexWorkspace ws(disk);
  const graphos::SparseSubcomplexResult sparse = graphos::subcomplex_from(disk, {{}, {}, {0}}, ws);
  CHECK(graphos::commutes_with_boundary(sparse.complex, disk, sparse.embedding));

  // and for the 1-skeleton, where the inclusion is not top-dimensional
  graphos::Marker skel(disk);
  for (graphos::Index e = 0; e < disk.count(1); ++e) skel.mark(1, e);
  const graphos::SubcomplexResult k1 = graphos::subcomplex(disk, skel);
  CHECK(graphos::commutes_with_boundary(k1.complex, disk, k1.embedding));
}

GRAPHOS_TEST_MAIN()
