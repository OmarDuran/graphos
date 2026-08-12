#pragma once

#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/ops/star_deletion.hpp"

namespace graphos {

// Marks every free face: a cell τ with exactly one proper coface σ, in
// whose boundary τ appears exactly ONCE (a τ glued into σ twice — e.g. the
// vertex of a loop edge — is not free: removing it would not be a
// collapse). A free pair (τ, σ) is the site of an elementary (Whitehead)
// collapse — star_deletion at τ removes exactly {τ, σ} and is a simple
// homotopy equivalence.
//
// Collective. P=1 today.
inline Marker free_faces(const Complex& c) {
  Marker out(c);
  const int dim = c.dim();
  for (int k = 0; k < dim; ++k) {
    // direct cofaces WITH multiplicity (the coboundary keeps duplicates)
    const CoboundaryOperator cob = coboundary(c, k);
    // cofaces further up (any single one disqualifies)
    std::vector<Index> higher(static_cast<std::size_t>(c.count(k)), 0);
    for (int j = k + 2; j <= dim; ++j) {
      const Adjacency inc = incidence(c, k, j);
      for (Index i = 0; i < c.count(k); ++i) {
        higher[static_cast<std::size_t>(i)] +=
            inc.offsets[static_cast<std::size_t>(i) + 1] - inc.offsets[static_cast<std::size_t>(i)];
      }
    }
    for (Index i = 0; i < c.count(k); ++i) {
      const Index direct = cob.offsets[static_cast<std::size_t>(i) + 1] -
                           cob.offsets[static_cast<std::size_t>(i)];
      if (direct == 1 && higher[static_cast<std::size_t>(i)] == 0) out.mark(k, i);
    }
  }
  return out;
}

struct CollapseResult {
  Complex complex;
  // composition of all elementary-collapse chain maps: removed cells go to
  // zero, survivors keep track of their compacted indices
  ChainMap map;
  Index removed_pairs{0};
};

// Collapses the complex as far as elementary collapses reach: repeatedly
// remove a free pair (deterministically the lowest-dimension,
// lowest-index free face) until none remains. Every step is a simple
// homotopy equivalence, so the result has the same homotopy type — a
// collapsible complex (any disk) reduces to a point, a Möbius band to its
// core circle. The endpoint is a complex with no free faces, not
// necessarily a minimal model.
//
// Serial by design (each collapse changes the free-face set); cost is
// O(pairs × complex size) — a diagnostic tool, not an assembly kernel.
//
// Collective. P=1 today.
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
