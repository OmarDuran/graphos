#pragma once

#include <algorithm>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

struct DisjointUnionResult {
  Complex complex;
  ChainMap a_map;  // cells of A into the union
  ChainMap b_map;  // cells of B into the union
};

// The coproduct A ⊔ B: strata are concatenated with B's indices shifted.
// Complexes of different dimension are allowed (mixed-dimensional geometry);
// the result has the larger dimension.
//
// Collective: every rank calls with its local partitions of A and B; the
// result is the distributed coproduct. P=1 today.
inline DisjointUnionResult disjoint_union(const Complex& a, const Complex& b) {
  const int dim = std::max(a.dim(), b.dim());
  Complex out(dim);
  out.attach_vertices(a.count(0) + b.count(0));

  for (int k = 1; k <= dim; ++k) {
    const Index shift = a.count(k - 1);
    if (k <= a.dim()) {
      const BoundaryOperator& bnd = a.boundary(k);
      for (Index e = 0; e < a.count(k); ++e) {
        std::span<const Index> idx(bnd.indices.data() + bnd.offsets[e],
                                   static_cast<std::size_t>(bnd.offsets[e + 1] - bnd.offsets[e]));
        std::span<const Sign> sg(bnd.signs.data() + bnd.offsets[e],
                                 static_cast<std::size_t>(bnd.offsets[e + 1] - bnd.offsets[e]));
        out.attach_cell(k, idx, sg);
      }
    }
    if (k <= b.dim()) {
      const BoundaryOperator& bnd = b.boundary(k);
      std::vector<Index> shifted;
      for (Index e = 0; e < b.count(k); ++e) {
        shifted.assign(bnd.indices.begin() + bnd.offsets[e],
                       bnd.indices.begin() + bnd.offsets[e + 1]);
        for (Index& i : shifted) i += shift;
        std::span<const Sign> sg(bnd.signs.data() + bnd.offsets[e],
                                 static_cast<std::size_t>(bnd.offsets[e + 1] - bnd.offsets[e]));
        out.attach_cell(k, std::span<const Index>(shifted), sg);
      }
    }
  }

  DisjointUnionResult res{std::move(out), {}, {}};
  std::vector<Index> a_counts(static_cast<std::size_t>(dim) + 1, 0);
  std::vector<Index> b_counts(static_cast<std::size_t>(dim) + 1, 0);
  for (int k = 0; k <= dim; ++k) {
    a_counts[static_cast<std::size_t>(k)] = a.count(k);
    b_counts[static_cast<std::size_t>(k)] = b.count(k);
  }
  res.a_map = ChainMap::sized(a_counts);
  res.b_map = ChainMap::sized(b_counts);
  for (int k = 0; k <= dim; ++k) {
    for (Index i = 0; i < a.count(k); ++i) res.a_map.index[static_cast<std::size_t>(k)][i] = i;
    for (Index i = 0; i < b.count(k); ++i) {
      res.b_map.index[static_cast<std::size_t>(k)][i] = a.count(k) + i;
    }
  }
  return res;
}

}  // namespace graphos
