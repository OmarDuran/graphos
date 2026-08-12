#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"

namespace graphos {

// Unsigned sparse adjacency between two strata, CSR layout.
struct Adjacency {
  std::vector<Index> offsets{0};
  std::vector<Index> indices;
};

// The transitive incidence I(k, j): for each k-cell, the j-cells of its
// closure (j < k), its star (j > k), or itself (j == k). This is the
// operator DoF gathers, NetworkX incidence views, and nfempy's
// build_graph(dim, codim) sit on.
//
// Deliberately UNSIGNED: orientation coefficients compose only across
// single levels (multi-level compositions telescope and cancel — that is
// d∘d = 0); the signed operators are the one-level ∂_k and δ_k. Rows are
// sorted and duplicate-free.
//
// Local: computed per rank over its own cells (plus ghosts under
// distribution); no communication is implied.
inline Adjacency incidence(const Complex& c, int k, int j) {
  if (k < 0 || k > c.dim() || j < 0 || j > c.dim()) {
    throw std::invalid_argument("incidence: dimension out of range");
  }

  Adjacency out;
  const Index n = c.count(k);
  out.offsets.reserve(static_cast<std::size_t>(n) + 1);

  if (k == j) {
    out.indices.resize(static_cast<std::size_t>(n));
    for (Index e = 0; e < n; ++e) {
      out.indices[static_cast<std::size_t>(e)] = e;
      out.offsets.push_back(e + 1);
    }
    return out;
  }

  // one-level operators for the traversal direction, gathered once
  std::vector<CoboundaryOperator> up;
  if (j > k) {
    up.reserve(static_cast<std::size_t>(j - k));
    for (int level = k; level < j; ++level) up.push_back(coboundary(c, level));
  }

  std::vector<Index> cur, nxt;
  for (Index e = 0; e < n; ++e) {
    cur.assign(1, e);
    if (j < k) {
      for (int level = k; level > j; --level) {
        const BoundaryOperator& bnd = c.boundary(level);
        nxt.clear();
        for (const Index x : cur) {
          nxt.insert(nxt.end(), bnd.indices.begin() + bnd.offsets[x],
                     bnd.indices.begin() + bnd.offsets[x + 1]);
        }
        std::sort(nxt.begin(), nxt.end());
        nxt.erase(std::unique(nxt.begin(), nxt.end()), nxt.end());
        cur.swap(nxt);
      }
    } else {
      for (int level = k; level < j; ++level) {
        const CoboundaryOperator& cob = up[static_cast<std::size_t>(level - k)];
        nxt.clear();
        for (const Index x : cur) {
          nxt.insert(nxt.end(), cob.indices.begin() + cob.offsets[x],
                     cob.indices.begin() + cob.offsets[x + 1]);
        }
        std::sort(nxt.begin(), nxt.end());
        nxt.erase(std::unique(nxt.begin(), nxt.end()), nxt.end());
        cur.swap(nxt);
      }
    }
    out.indices.insert(out.indices.end(), cur.begin(), cur.end());
    out.offsets.push_back(static_cast<Index>(out.indices.size()));
  }
  return out;
}

}  // namespace graphos
