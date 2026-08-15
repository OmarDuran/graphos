#pragma once

#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/ops/star_deletion.hpp"

namespace graphos {

// Marks the free faces: τ with exactly one proper coface σ, appearing in ∂σ
// exactly once. A τ glued into σ twice — the vertex of a loop edge — is not
// free. A free pair (τ, σ) is the site of an elementary Whitehead collapse:
// star_deletion at τ removes exactly {τ, σ} and is a simple homotopy
// equivalence.
inline Marker free_faces(const Complex& c) {
  Marker out(c);
  const int dim = c.dim();
  for (int k = 0; k < dim; ++k) {
    // direct cofaces with multiplicity: δ keeps duplicates
    const CoboundaryOperator cob = coboundary(c, k);
    // any higher coface disqualifies
    std::vector<Index> higher(static_cast<std::size_t>(c.count(k)), 0);
    for (int j = k + 2; j <= dim; ++j) {
      const Adjacency inc = incidence(c, k, j);
      for (Index i = 0; i < c.count(k); ++i) {
        higher[static_cast<std::size_t>(i)] +=
            inc.offsets[static_cast<std::size_t>(i) + 1] - inc.offsets[static_cast<std::size_t>(i)];
      }
    }
    for (Index i = 0; i < c.count(k); ++i) {
      const Index direct =
          cob.offsets[static_cast<std::size_t>(i) + 1] - cob.offsets[static_cast<std::size_t>(i)];
      if (direct == 1 && higher[static_cast<std::size_t>(i)] == 0) out.mark(k, i);
    }
  }
  return out;
}

struct CollapseResult {
  Complex complex;
  // the composite of the elementary-collapse chain maps
  ChainMap map;
  Index removed_pairs{0};
};

// Removes free pairs until none remains, taking the lowest-dimension,
// lowest-index free face each step. Every step is a simple homotopy
// equivalence, so the homotopy type is preserved: a collapsible complex
// reduces to a point, a Möbius band to its core circle. The endpoint has no
// free face; it need not be a minimal model.
//
// Serial by design — each collapse changes the free-face set — at
// O(pairs × N). A diagnostic, not an assembly kernel.
inline CollapseResult collapse(const Complex& c) {
  Complex cur = c;
  ChainMap total = ChainMap::sized(c.counts());
  for (int k = 0; k <= c.dim(); ++k) {
    for (Index i = 0; i < c.count(k); ++i) {
      total.index[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = i;
    }
  }

  Index removed = 0;
  for (;;) {
    const Marker ff = free_faces(cur);
    int k_pick = -1;
    Index i_pick = invalid_index;
    for (int k = 0; k <= cur.dim() && k_pick < 0; ++k) {
      for (Index i = 0; i < cur.count(k); ++i) {
        if (ff.marked(k, i)) {
          k_pick = k;
          i_pick = i;
          break;
        }
      }
    }
    if (k_pick < 0) break;
    Marker one(cur);
    one.mark(k_pick, i_pick);
    StarDeletionResult sd = star_deletion(cur, one);
    cur = std::move(sd.complex);
    total = compose(total, sd.map);
    ++removed;
  }

  return CollapseResult{std::move(cur), std::move(total), removed};
}

}  // namespace graphos
