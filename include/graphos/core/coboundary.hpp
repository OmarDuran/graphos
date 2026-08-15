#pragma once

#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

// δ_k = ∂_{k+1}^T, the coboundary operator, in CSR. Row τ lists the cofaces σ
// of τ with [σ : τ]. Applied to a k-cochain it is the discrete exterior
// derivative.
struct CoboundaryOperator {
  std::vector<Index> offsets;
  std::vector<Index> indices;
  std::vector<Sign> signs;
};

inline CoboundaryOperator coboundary(const Complex& c, int k) {
  if (k < 0 || k >= c.dim()) {
    throw std::invalid_argument("coboundary: dimension out of range");
  }
  const BoundaryOperator& bnd = c.boundary(k + 1);
  const Index n_lo = c.count(k);
  const Index n_hi = c.count(k + 1);

  CoboundaryOperator cob;
  cob.offsets.assign(static_cast<std::size_t>(n_lo) + 1, 0);
  for (const Index f : bnd.indices) ++cob.offsets[static_cast<std::size_t>(f) + 1];
  for (Index i = 0; i < n_lo; ++i) {
    cob.offsets[static_cast<std::size_t>(i) + 1] += cob.offsets[static_cast<std::size_t>(i)];
  }
  cob.indices.resize(bnd.indices.size());
  cob.signs.resize(bnd.signs.size());
  std::vector<Index> cursor(cob.offsets.begin(), cob.offsets.end() - 1);
  for (Index e = 0; e < n_hi; ++e) {
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      const std::size_t w =
          static_cast<std::size_t>(cursor[static_cast<std::size_t>(bnd.indices[m])]++);
      cob.indices[w] = e;
      cob.signs[w] = bnd.signs[m];
    }
  }
  return cob;
}

}  // namespace graphos
