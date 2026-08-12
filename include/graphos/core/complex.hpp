#pragma once

#include <cstdlib>
#include <initializer_list>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "graphos/core/types.hpp"

namespace graphos {

// The boundary operator ∂_k as a signed sparse matrix from k-cells to
// (k-1)-cells, in CSR layout. Row e lists the (k-1)-cells on the boundary of
// cell e with their orientation coefficients (+/-1 only).
//
// Variable-arity CSR is the general layout; a fixed-arity strided
// specialization (no offsets) is planned for single-cell-type complexes.
struct BoundaryOperator {
  std::vector<Index> offsets{0};
  std::vector<Index> indices;
  std::vector<Sign> signs;

  Index rows() const { return static_cast<Index>(offsets.size()) - 1; }

  void append_row(std::span<const Index> row_indices, std::span<const Sign> row_signs) {
    if (row_indices.size() != row_signs.size()) {
      throw std::invalid_argument("BoundaryOperator::append_row: indices/signs size mismatch");
    }
    indices.insert(indices.end(), row_indices.begin(), row_indices.end());
    signs.insert(signs.end(), row_signs.begin(), row_signs.end());
    offsets.push_back(static_cast<Index>(indices.size()));
  }
};

// A finite, metric-free cell complex, stratified by dimension: cell counts
// per dimension plus the boundary operators ∂_k. Cells are nothing but
// indices; all geometric meaning (coordinates, metrics) lives outside
// graphos.
//
// The complex may be mixed-dimensional: any k-skeleton may contain maximal
// cells (cells that are not faces of anything above), which is how
// lower-dimensional subdomains (fractures, interfaces, wells) coexist with
// the bulk.
class Complex {
 public:
  explicit Complex(int dim)
      : counts_(static_cast<std::size_t>(dim) + 1, 0),
        boundary_(static_cast<std::size_t>(dim) + 1) {
    if (dim < 0) throw std::invalid_argument("Complex: dimension must be >= 0");
  }

  // Bulk constructor for operations that assemble whole strata in kernels.
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

  // The GLOBAL number of k-cells: a collective value, identical on every
  // rank (P=1: trivially the local count).
  Index count(int k) const {
    return (k >= 0 && k <= dim()) ? counts_[static_cast<std::size_t>(k)] : Index{0};
  }

  std::vector<Index> counts() const { return counts_; }

  Index attach_vertices(Index n) {
    if (n < 0) throw std::invalid_argument("Complex::attach_vertices: negative count");
    counts_[0] += n;
    return counts_[0];
  }

  // Attaches a k-cell along its boundary chain. Returns the new cell's index.
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

  // Structural invariants: row counts match cell counts, boundary indices in
  // range, orientation coefficients strictly +/-1.
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

// The fundamental identity of a chain complex: ∂_{k-1} ∘ ∂_k = 0. This is
// the check that catches every orientation-convention mistake.
// Collective: the verdict is global and identical on every rank (P=1 today).
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

// Collective: the global alternating sum, identical on every rank.
inline long long euler_characteristic(const Complex& c) {
  long long chi = 0;
  for (int k = 0; k <= c.dim(); ++k) {
    chi += (k % 2 == 0 ? 1LL : -1LL) * static_cast<long long>(c.count(k));
  }
  return chi;
}

}  // namespace graphos
