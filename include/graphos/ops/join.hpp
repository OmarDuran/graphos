#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

struct JoinResult {
  Complex complex;

  // Both factors embed as subcomplexes of the join; these are the
  // embeddings (all signs +1).
  ChainMap a_map;
  ChainMap b_map;

  // Join cells of dimension k: A's k-cells first, then B's k-cells, then
  // pair blocks (α, β) with dim α + dim β + 1 = k, ascending in p = dim α,
  // i-major: index = offset + i * b_count + j. Empty blocks omitted.
  struct Block {
    int p;
    int q;
    Index offset;
    Index a_count;
    Index b_count;
  };
  std::vector<std::vector<Block>> pair_blocks;
};

// The join A ∗ B: every cell of A, every cell of B, and one cell α∗β of
// dimension dim α + dim β + 1 for every pair, with the boundary
//   ∂(α∗β) = ∂α∗β + (−1)^{dim α + 1} α∗∂β,
// where ∂(vertex)∗β degenerates to β itself (and symmetrically). This is
// the constructor of the PL toolkit: join(point, X) is the CONE over X,
// join(S⁰, X) the SUSPENSION, and the regular-neighborhood factorization
// cl(st σ) ≅ σ ∗ lk(σ) becomes computable on both sides.
//
// Collective. P=1 today.
inline JoinResult join(const Complex& a, const Complex& b) {
  const int dim = a.dim() + b.dim() + 1;

  JoinResult res{Complex(dim), ChainMap::sized(a.counts()), ChainMap::sized(b.counts()), {}};
  res.pair_blocks.resize(static_cast<std::size_t>(dim) + 1);
  std::vector<Index> counts(static_cast<std::size_t>(dim) + 1, 0);
  for (int k = 0; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    counts[sk] = a.count(k) + b.count(k);
    for (int p = std::max(0, k - 1 - b.dim()); p <= std::min(k - 1, a.dim()); ++p) {
      const int q = k - 1 - p;
      if (a.count(p) == 0 || b.count(q) == 0) continue;
      const long long block = static_cast<long long>(a.count(p)) * b.count(q);
      if (block + counts[sk] > static_cast<long long>(std::numeric_limits<Index>::max())) {
        throw std::overflow_error("join: cell count exceeds Index capacity");
      }
      res.pair_blocks[sk].push_back(JoinResult::Block{p, q, counts[sk], a.count(p), b.count(q)});
      counts[sk] += static_cast<Index>(block);
    }
  }

  const auto pair_offset = [&res](int k, int p) -> Index {
    for (const JoinResult::Block& blk : res.pair_blocks[static_cast<std::size_t>(k)]) {
      if (blk.p == p) return blk.offset;
    }
    throw std::logic_error("join: missing pair block");
  };

  Complex out(dim);
  out.attach_vertices(counts[0]);
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (int k = 1; k <= dim; ++k) {
    // A's cells, verbatim (A occupies the leading indices of every stratum)
    if (k <= a.dim()) {
      const BoundaryOperator& bnd = a.boundary(k);
      for (Index e = 0; e < a.count(k); ++e) {
        row_idx.assign(bnd.indices.begin() + bnd.offsets[e],
                       bnd.indices.begin() + bnd.offsets[e + 1]);
        row_sg.assign(bnd.signs.begin() + bnd.offsets[e], bnd.signs.begin() + bnd.offsets[e + 1]);
        out.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
      }
    }
    // B's cells, shifted past A's
    if (k <= b.dim()) {
      const BoundaryOperator& bnd = b.boundary(k);
      for (Index e = 0; e < b.count(k); ++e) {
        row_idx.clear();
        row_sg.clear();
        for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
          row_idx.push_back(a.count(k - 1) + bnd.indices[m]);
          row_sg.push_back(bnd.signs[m]);
        }
        out.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
      }
    }
    // pair cells
    for (const JoinResult::Block& blk : res.pair_blocks[static_cast<std::size_t>(k)]) {
      const Sign flip = (blk.p % 2 == 0) ? Sign{-1} : Sign{1};  // (−1)^{p+1}
      for (Index i = 0; i < blk.a_count; ++i) {
        for (Index j = 0; j < blk.b_count; ++j) {
          row_idx.clear();
          row_sg.clear();
          if (blk.p >= 1) {
            const BoundaryOperator& ba = a.boundary(blk.p);
            const Index off = pair_offset(k - 1, blk.p - 1);
            for (Index m = ba.offsets[i]; m < ba.offsets[i + 1]; ++m) {
              row_idx.push_back(off + ba.indices[m] * blk.b_count + j);
              row_sg.push_back(ba.signs[m]);
            }
          } else {
            // ∂(vertex)∗β degenerates to β itself
            row_idx.push_back(a.count(k - 1) + j);
            row_sg.push_back(Sign{1});
          }
          if (blk.q >= 1) {
            const BoundaryOperator& bb = b.boundary(blk.q);
            const Index off = pair_offset(k - 1, blk.p);
            const Index b_faces = b.count(blk.q - 1);
            for (Index m = bb.offsets[j]; m < bb.offsets[j + 1]; ++m) {
              row_idx.push_back(off + i * b_faces + bb.indices[m]);
              row_sg.push_back(static_cast<Sign>(flip * bb.signs[m]));
            }
          } else {
            // α∗∂(vertex) degenerates to α itself
            row_idx.push_back(i);
            row_sg.push_back(flip);
          }
          out.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
        }
      }
    }
  }

  for (int k = 0; k <= a.dim(); ++k) {
    for (Index i = 0; i < a.count(k); ++i) {
      res.a_map.index[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = i;
    }
  }
  for (int k = 0; k <= b.dim(); ++k) {
    for (Index i = 0; i < b.count(k); ++i) {
      res.b_map.index[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = a.count(k) + i;
    }
  }
  res.complex = std::move(out);
  return res;
}

}  // namespace graphos
