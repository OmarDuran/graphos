#pragma once

#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "graphos/core/types.hpp"

namespace graphos {

// ∂_k ∈ {−1, 0, +1}^{N_{k−1} × N_k}, stored row-major by k-cell: row σ holds
// the (k−1)-faces τ of σ with incidence numbers [σ : τ] ∈ {±1}. Absence
// encodes 0, so the CSR rows are the columns of ∂_k.
//
// Variable-arity CSR is the general layout; a fixed-arity strided form is
// planned for complexes of a single cell type.
struct BoundaryOperator {
  std::vector<Index> offsets{0};
  std::vector<Index> indices;
  std::vector<Sign> signs;

  Index rows() const { return static_cast<Index>(offsets.size()) - 1; }

  void append_row(std::span<const Index> row_indices, std::span<const Sign> row_signs) {
    if (row_indices.size() != row_signs.size()) {
      throw std::invalid_argument("BoundaryOperator::append_row: indices/signs size mismatch");
    }
    // offsets are Index-typed, so nnz ≤ max(Index) ≈ 2.1e9
    if (indices.size() + row_indices.size() >
        static_cast<std::size_t>(std::numeric_limits<Index>::max())) {
      throw std::overflow_error("BoundaryOperator: nnz exceeds Index capacity");
    }
    indices.insert(indices.end(), row_indices.begin(), row_indices.end());
    signs.insert(signs.end(), row_signs.begin(), row_signs.end());
    offsets.push_back(static_cast<Index>(indices.size()));
  }
};

// A finite stratified signed incidence structure C = (N₀, …, N_n ; ∂₁, …, ∂_n).
// A k-cell is an index; no realization, coordinate or metric is stored.
//
// C is a chain complex exactly when ∂_{k−1} ∘ ∂_k = 0; the representation does
// not impose it (see d_squared_is_zero). Deliberately admitted: mixed-
// dimensional (maximal cells in any stratum), polytopal, non-regular,
// non-manifold, non-orientable. Validity is queried, not presupposed.
class Complex {
 public:
  explicit Complex(int dim)
      : counts_(static_cast<std::size_t>(dim) + 1, 0),
        boundary_(static_cast<std::size_t>(dim) + 1) {
    if (dim < 0) throw std::invalid_argument("Complex: dimension must be >= 0");
  }

  // Stratum-wise constructor, for operations that assemble ∂_k in kernels.
  Complex(std::vector<Index> cell_counts, std::vector<BoundaryOperator> boundary)
      : counts_(std::move(cell_counts)), boundary_(std::move(boundary)) {
    if (counts_.empty() || boundary_.size() != counts_.size()) {
      throw std::invalid_argument("Complex: counts/boundary stratification mismatch");
    }
    for (int k = 1; k <= dim(); ++k) {
      if (boundary_[static_cast<std::size_t>(k)].rows() != counts_[static_cast<std::size_t>(k)]) {
        throw std::invalid_argument("Complex: boundary operator row count mismatch");
      }
    }
  }

  int dim() const { return static_cast<int>(counts_.size()) - 1; }

  // N_k, global and rank-invariant (P = 1: the local count).
  Index count(int k) const {
    return (k >= 0 && k <= dim()) ? counts_[static_cast<std::size_t>(k)] : Index{0};
  }

  std::vector<Index> counts() const { return counts_; }

  Index attach_vertices(Index n) {
    if (n < 0) throw std::invalid_argument("Complex::attach_vertices: negative count");
    counts_[0] += n;
    return counts_[0];
  }

  // Attaches a k-cell along the (k−1)-chain ∂σ — the CW paradigm. Returns σ.
  Index attach_cell(int k, std::span<const Index> bnd, std::span<const Sign> sg) {
    if (k < 1 || k > dim()) {
      throw std::invalid_argument("Complex::attach_cell: dimension out of range");
    }
    for (const Index b : bnd) {
      if (b < 0 || b >= count(k - 1)) {
        throw std::out_of_range("Complex::attach_cell: boundary index out of range");
      }
    }
    for (const Sign s : sg) {
      if (s != 1 && s != -1) {
        throw std::invalid_argument("Complex::attach_cell: signs must be +/-1");
      }
    }
    boundary_[static_cast<std::size_t>(k)].append_row(bnd, sg);
    return counts_[static_cast<std::size_t>(k)]++;
  }

  Index attach_cell(int k, std::initializer_list<Index> bnd, std::initializer_list<int> sg) {
    std::vector<Index> bv(bnd);
    std::vector<Sign> sv;
    sv.reserve(sg.size());
    for (const int s : sg) sv.push_back(static_cast<Sign>(s));
    return attach_cell(k, std::span<const Index>(bv), std::span<const Sign>(sv));
  }

  const BoundaryOperator& boundary(int k) const {
    if (k < 1 || k > dim()) {
      throw std::invalid_argument("Complex::boundary: dimension out of range");
    }
    return boundary_[static_cast<std::size_t>(k)];
  }

  // Structural invariants: row counts equal N_k, face indices in [0, N_{k−1}),
  // incidence numbers strictly ±1. Not ∂∘∂ = 0, which is d_squared_is_zero.
  void validate() const {
    for (int k = 1; k <= dim(); ++k) {
      const BoundaryOperator& bnd = boundary_[static_cast<std::size_t>(k)];
      if (bnd.rows() != count(k)) {
        throw std::runtime_error("Complex::validate: row count mismatch in dimension " +
                                 std::to_string(k));
      }
      for (std::size_t m = 0; m < bnd.indices.size(); ++m) {
        if (bnd.indices[m] < 0 || bnd.indices[m] >= count(k - 1)) {
          throw std::runtime_error("Complex::validate: boundary index out of range in dimension " +
                                   std::to_string(k));
        }
        if (bnd.signs[m] != 1 && bnd.signs[m] != -1) {
          throw std::runtime_error("Complex::validate: sign not +/-1 in dimension " +
                                   std::to_string(k));
        }
      }
    }
  }

 private:
  std::vector<Index> counts_;
  std::vector<BoundaryOperator> boundary_;  // boundary_[k] valid for k >= 1
};

// ∂_{k−1} ∘ ∂_k = 0, the defining identity of a chain complex: every
// (k−2)-face of σ is reached twice with opposite incidence. Rank-invariant.
inline bool d_squared_is_zero(const Complex& c) {
  for (int k = 2; k <= c.dim(); ++k) {
    const BoundaryOperator& hi = c.boundary(k);
    const BoundaryOperator& lo = c.boundary(k - 1);
    for (Index e = 0; e < c.count(k); ++e) {
      std::map<Index, int> acc;
      for (Index m = hi.offsets[e]; m < hi.offsets[e + 1]; ++m) {
        const Index f = hi.indices[m];
        const int s1 = hi.signs[m];
        for (Index p = lo.offsets[f]; p < lo.offsets[f + 1]; ++p) {
          acc[lo.indices[p]] += s1 * lo.signs[p];
        }
      }
      for (const auto& [g, v] : acc) {
        (void)g;
        if (v != 0) return false;
      }
    }
  }
  return true;
}

// ∂π = π∂ on every generator: π is a chain map C_*(a) → C_*(b), degree 0.
// It then descends to homology, and its dual π* commutes with d — which is
// what a level transition, a subcomplex embedding or a quotient must satisfy
// for the coarse operator to discretize the same complex. invalid_index sends
// a generator to 0; the zero map qualifies.
inline bool commutes_with_boundary(const Complex& a, const Complex& b, const ChainMap& pi) {
  if (static_cast<int>(pi.index.size()) <= a.dim()) return false;
  for (int k = 0; k <= a.dim(); ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    if (pi.index[sk].size() != static_cast<std::size_t>(a.count(k))) return false;
    if (pi.sign[sk].size() != pi.index[sk].size()) return false;
    for (const Index t : pi.index[sk]) {
      if (t == invalid_index) continue;
      if (k > b.dim() || t < 0 || t >= b.count(k)) return false;  // off the target complex
    }
  }
  for (int k = 1; k <= a.dim(); ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const BoundaryOperator& ba = a.boundary(k);
    for (Index e = 0; e < a.count(k); ++e) {
      std::map<Index, int> acc;  // ∂(πe) − π(∂e), as a chain in b
      const Index im = pi.index[sk][static_cast<std::size_t>(e)];
      if (im != invalid_index) {
        const int sg = pi.sign[sk][static_cast<std::size_t>(e)];
        const BoundaryOperator& bb = b.boundary(k);
        for (Index m = bb.offsets[im]; m < bb.offsets[im + 1]; ++m) {
          acc[bb.indices[m]] += sg * bb.signs[m];
        }
      }
      for (Index m = ba.offsets[e]; m < ba.offsets[e + 1]; ++m) {
        const Index f = ba.indices[m];
        const Index jm = pi.index[sk - 1][static_cast<std::size_t>(f)];
        if (jm == invalid_index) continue;
        acc[jm] -= ba.signs[m] * pi.sign[sk - 1][static_cast<std::size_t>(f)];
      }
      for (const auto& [g, v] : acc) {
        (void)g;
        if (v != 0) return false;
      }
    }
  }
  return true;
}

// χ(C) = Σ_k (−1)^k N_k. Rank-invariant.
inline long long euler_characteristic(const Complex& c) {
  long long chi = 0;
  for (int k = 0; k <= c.dim(); ++k) {
    chi += (k % 2 == 0 ? 1LL : -1LL) * static_cast<long long>(c.count(k));
  }
  return chi;
}

}  // namespace graphos
