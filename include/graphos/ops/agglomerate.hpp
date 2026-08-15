#pragma once

#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"

namespace graphos {

struct AgglomerationResult {
  // One coarse top cell per aggregate, its boundary Σ ∂(members) with
  // interior facets cancelling. The surviving lower strata are the closure of
  // the interfaces, ∂K, and the maximal lower-dimensional cells.
  Complex complex;
  // fine → coarse: a top cell maps to its aggregate, a cell interior to an
  // aggregate to 0. Restriction and prolongation key off this.
  ChainMap map;
};

// Merges the top cells into the aggregates named by `labels`, one label in
// [0, n_aggregates) per top cell. The direction opposite refinement, and the
// constructor of polytopal coarse spaces.
//
// The coarse boundary is Σ ∂(members). Interior facets cancel exactly when
// the complex is coherently oriented, so an incoherent one is detected — a
// facet with coefficient ±2 — and reported; run orient() first. Aggregates
// should be facet-connected; a disconnected one yields a valid chain but not
// a cell-like one.
inline AgglomerationResult agglomerate(const Complex& c, const std::vector<Index>& labels) {
  const int n = c.dim();
  if (n < 1) throw std::invalid_argument("agglomerate: complex must have dimension >= 1");
  if (static_cast<Index>(labels.size()) != c.count(n)) {
    throw std::invalid_argument("agglomerate: one label per top cell required");
  }
  Index n_agg = 0;
  for (const Index a : labels) {
    if (a < 0) throw std::invalid_argument("agglomerate: labels must be non-negative");
    n_agg = std::max(n_agg, a + 1);
  }
  std::vector<char> used(static_cast<std::size_t>(n_agg), 0);
  for (const Index a : labels) used[static_cast<std::size_t>(a)] = 1;
  for (const char u : used) {
    if (!u) throw std::invalid_argument("agglomerate: aggregate ids must be contiguous");
  }

  // coarse top boundaries: the signed facet sum per aggregate
  const BoundaryOperator& bnd = c.boundary(n);
  std::vector<std::map<Index, int>> acc(static_cast<std::size_t>(n_agg));
  for (Index e = 0; e < c.count(n); ++e) {
    std::map<Index, int>& row = acc[static_cast<std::size_t>(labels[static_cast<std::size_t>(e)])];
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      row[bnd.indices[m]] += bnd.signs[m];
    }
  }
  for (const auto& row : acc) {
    for (const auto& [f, coeff] : row) {
      (void)f;
      if (coeff < -1 || coeff > 1) {
        throw std::invalid_argument(
            "agglomerate: facet coefficient exceeds 1 - orientation is inconsistent, "
            "run orient() first");
      }
    }
  }

  // survivors below the top: cl({facets with nonzero coefficient} ∪ {maximal
  // lower-dimensional cells})
  std::vector<std::vector<char>> keep(static_cast<std::size_t>(n));
  for (int k = 0; k < n; ++k) {
    keep[static_cast<std::size_t>(k)].assign(static_cast<std::size_t>(c.count(k)), 0);
  }
  for (const auto& row : acc) {
    for (const auto& [f, coeff] : row) {
      if (coeff != 0) keep[static_cast<std::size_t>(n) - 1][static_cast<std::size_t>(f)] = 1;
    }
  }
  for (int k = 0; k < n; ++k) {
    const CoboundaryOperator cob = coboundary(c, k);
    for (Index i = 0; i < c.count(k); ++i) {
      if (cob.offsets[static_cast<std::size_t>(i) + 1] ==
          cob.offsets[static_cast<std::size_t>(i)]) {
        keep[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = 1;  // maximal lower cell
      }
    }
  }
  for (int k = n - 1; k >= 1; --k) {
    const BoundaryOperator& b = c.boundary(k);
    for (Index e = 0; e < c.count(k); ++e) {
      if (!keep[static_cast<std::size_t>(k)][static_cast<std::size_t>(e)]) continue;
      for (Index m = b.offsets[e]; m < b.offsets[e + 1]; ++m) {
        keep[static_cast<std::size_t>(k) - 1][static_cast<std::size_t>(b.indices[m])] = 1;
      }
    }
  }

  // compact, and record the chain map
  ChainMap map = ChainMap::sized(c.counts());
  std::vector<Index> coarse_counts(static_cast<std::size_t>(n) + 1, 0);
  for (int k = 0; k < n; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    Index nn = 0;
    for (Index i = 0; i < c.count(k); ++i) {
      map.index[sk][static_cast<std::size_t>(i)] =
          keep[sk][static_cast<std::size_t>(i)] ? nn++ : invalid_index;
    }
    coarse_counts[sk] = nn;
  }
  coarse_counts[static_cast<std::size_t>(n)] = n_agg;
  for (Index e = 0; e < c.count(n); ++e) {
    map.index[static_cast<std::size_t>(n)][static_cast<std::size_t>(e)] =
        labels[static_cast<std::size_t>(e)];
  }

  // assemble the coarse strata
  std::vector<BoundaryOperator> strata(static_cast<std::size_t>(n) + 1);
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (int k = 1; k < n; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const BoundaryOperator& b = c.boundary(k);
    for (Index e = 0; e < c.count(k); ++e) {
      if (!keep[sk][static_cast<std::size_t>(e)]) continue;
      row_idx.clear();
      row_sg.clear();
      for (Index m = b.offsets[e]; m < b.offsets[e + 1]; ++m) {
        row_idx.push_back(map.index[sk - 1][static_cast<std::size_t>(b.indices[m])]);
        row_sg.push_back(b.signs[m]);
      }
      strata[sk].append_row(row_idx, row_sg);
    }
  }
  for (Index a = 0; a < n_agg; ++a) {
    row_idx.clear();
    row_sg.clear();
    for (const auto& [f, coeff] : acc[static_cast<std::size_t>(a)]) {
      if (coeff == 0) continue;
      row_idx.push_back(map.index[static_cast<std::size_t>(n) - 1][static_cast<std::size_t>(f)]);
      row_sg.push_back(static_cast<Sign>(coeff));
    }
    strata[static_cast<std::size_t>(n)].append_row(row_idx, row_sg);
  }

  return AgglomerationResult{Complex(std::move(coarse_counts), std::move(strata)), std::move(map)};
}

}  // namespace graphos
