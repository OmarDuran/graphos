#pragma once

#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

// Unsigned upward adjacency: for each k-cell, the (k+1)-cells having it on
// their boundary. The sparsity transpose of ∂_{k+1}; signs are not carried
// because coface traversal is a structural query.
struct Adjacency {
  std::vector<Index> offsets;
  std::vector<Index> indices;
};

inline Adjacency coboundary(const Complex& c, int k) {
  if (k < 0 || k >= c.dim()) {
    throw std::invalid_argument("coboundary: dimension out of range");
  }
  const BoundaryOperator& bnd = c.boundary(k + 1);
  const Index n_lo = c.count(k);
  const Index n_hi = c.count(k + 1);

  Adjacency adj;
  adj.offsets.assign(static_cast<std::size_t>(n_lo) + 1, 0);
  for (const Index f : bnd.indices) ++adj.offsets[static_cast<std::size_t>(f) + 1];
  for (Index i = 0; i < n_lo; ++i) {
    adj.offsets[static_cast<std::size_t>(i) + 1] += adj.offsets[static_cast<std::size_t>(i)];
  }
  adj.indices.resize(bnd.indices.size());
  std::vector<Index> cursor(adj.offsets.begin(), adj.offsets.end() - 1);
  for (Index e = 0; e < n_hi; ++e) {
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      adj.indices[static_cast<std::size_t>(cursor[static_cast<std::size_t>(bnd.indices[m])]++)] = e;
    }
  }
  return adj;
}

}  // namespace graphos
