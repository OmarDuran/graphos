#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

struct ProductResult {
  Complex complex;

  // Product cells are pairs (α, β) with dim(α×β) = dim α + dim β. For each
  // product dimension k, cells are laid out in blocks of fixed factor
  // dimensions (p, q = k−p), ascending in p, i-major within a block:
  //   index(α_i × β_j) = offset + i * b_count + j.
  // Empty blocks (either factor stratum empty) are omitted.
  struct Block {
    int p;
    int q;
    Index offset;
    Index a_count;
    Index b_count;
  };
  std::vector<std::vector<Block>> blocks;  // per product dimension
};

// The Cartesian product A × B, with boundary by the Leibniz rule:
//   ∂(α×β) = ∂α×β + (−1)^{dim α} α×∂β.
// product(mesh, segment) is mesh EXTRUSION (triangles become prisms, one
// layer per segment edge); products of intervals build tensor cells. The
// factor structure is preserved in `blocks`, and the Euler characteristic
// multiplies: χ(A×B) = χ(A)·χ(B).
//
// Collective. P=1 today.
inline ProductResult product(const Complex& a, const Complex& b) {
  const int dim = a.dim() + b.dim();

  ProductResult res{Complex(dim), {}};
  res.blocks.resize(static_cast<std::size_t>(dim) + 1);
  std::vector<Index> counts(static_cast<std::size_t>(dim) + 1, 0);
  for (int k = 0; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    for (int p = std::max(0, k - b.dim()); p <= std::min(k, a.dim()); ++p) {
      const int q = k - p;
      if (a.count(p) == 0 || b.count(q) == 0) continue;
      const long long block = static_cast<long long>(a.count(p)) * b.count(q);
      if (block + counts[sk] > static_cast<long long>(std::numeric_limits<Index>::max())) {
        throw std::overflow_error("product: cell count exceeds Index capacity");
      }
      res.blocks[sk].push_back(ProductResult::Block{p, q, counts[sk], a.count(p), b.count(q)});
      counts[sk] += static_cast<Index>(block);
    }
  }

  const auto block_offset = [&res](int k, int p) -> Index {
    for (const ProductResult::Block& blk : res.blocks[static_cast<std::size_t>(k)]) {
      if (blk.p == p) return blk.offset;
    }
    throw std::logic_error("product: missing face block");
  };

  Complex out(dim);
  out.attach_vertices(counts[0]);
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (int k = 1; k <= dim; ++k) {
    for (const ProductResult::Block& blk : res.blocks[static_cast<std::size_t>(k)]) {
      const BoundaryOperator* ba = blk.p >= 1 ? &a.boundary(blk.p) : nullptr;
      const BoundaryOperator* bb = blk.q >= 1 ? &b.boundary(blk.q) : nullptr;
      const Sign flip = (blk.p % 2 == 0) ? Sign{1} : Sign{-1};
      for (Index i = 0; i < blk.a_count; ++i) {
        for (Index j = 0; j < blk.b_count; ++j) {
          row_idx.clear();
          row_sg.clear();
          if (ba != nullptr) {
            const Index off = block_offset(k - 1, blk.p - 1);
            for (Index m = ba->offsets[i]; m < ba->offsets[i + 1]; ++m) {
              row_idx.push_back(off + ba->indices[m] * blk.b_count + j);
              row_sg.push_back(ba->signs[m]);
            }
          }
          if (bb != nullptr) {
            const Index off = block_offset(k - 1, blk.p);
            const Index b_faces = b.count(blk.q - 1);
            for (Index m = bb->offsets[j]; m < bb->offsets[j + 1]; ++m) {
              row_idx.push_back(off + i * b_faces + bb->indices[m]);
              row_sg.push_back(static_cast<Sign>(flip * bb->signs[m]));
            }
          }
          out.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
        }
      }
    }
  }

  res.complex = std::move(out);
  return res;
}

}  // namespace graphos
