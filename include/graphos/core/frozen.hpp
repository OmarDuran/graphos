#pragma once

#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/exec/array.hpp"

namespace graphos {

// Signed CSR incidence in device-capable storage (exec::Array, the CHAI
// seam), with the POD view kernels capture by value.
struct FrozenCsr {
  exec::Array<Index> offsets;
  exec::Array<Index> indices;
  exec::Array<Sign> signs;
};

struct CsrView {
  exec::Array<Index>::view_type offsets;
  exec::Array<Index>::view_type indices;
  exec::Array<Sign>::view_type signs;
};

// The immutable half of the builder/frozen split: Complex is the mutable
// builder, freezing yields a FrozenComplex carrying ∂_k and the derived δ_k
// in device-capable arrays. Queries, kernels and the NetworkX views sit on
// it; immutability is what makes zero-copy views sound.
//
// Operations build new complexes and freeze them again, so freezing is the
// epoch boundary of the bulk-edit model.
class FrozenComplex {
 public:
  // Under distribution, freezing is where partitioning, ownership and halo
  // construction happen; halo_depth bounds the neighbourhood st/cl/lk can
  // answer. At P = 1 the partition is the whole complex and it has no effect.
  explicit FrozenComplex(const Complex& c, int halo_depth = 1)
      : dim_(c.dim()),
        halo_depth_(halo_depth),
        counts_(c.counts()),
        boundary_(static_cast<std::size_t>(c.dim()) + 1),
        coboundary_(static_cast<std::size_t>(c.dim()) + 1) {
    if (halo_depth < 1) throw std::invalid_argument("FrozenComplex: halo_depth must be >= 1");
    c.validate();
    for (int k = 1; k <= dim_; ++k) {
      const BoundaryOperator& bnd = c.boundary(k);
      boundary_[static_cast<std::size_t>(k)] = FrozenCsr{
          exec::Array<Index>(bnd.offsets),
          exec::Array<Index>(bnd.indices),
          exec::Array<Sign>(bnd.signs),
      };
    }
    for (int k = 0; k < dim_; ++k) {
      const CoboundaryOperator cob = coboundary(c, k);
      coboundary_[static_cast<std::size_t>(k)] = FrozenCsr{
          exec::Array<Index>(cob.offsets),
          exec::Array<Index>(cob.indices),
          exec::Array<Sign>(cob.signs),
      };
    }
  }

  // move-only: exec::Array owns device-capable memory. Stated explicitly so
  // type traits stay honest for the bindings.
  FrozenComplex(const FrozenComplex&) = delete;
  FrozenComplex& operator=(const FrozenComplex&) = delete;
  FrozenComplex(FrozenComplex&&) = default;
  FrozenComplex& operator=(FrozenComplex&&) = default;

  int dim() const { return dim_; }

  int halo_depth() const { return halo_depth_; }

  // N_k, global and rank-invariant (P = 1: the local count).
  Index count(int k) const {
    return (k >= 0 && k <= dim_) ? counts_[static_cast<std::size_t>(k)] : Index{0};
  }

  std::vector<Index> counts() const { return counts_; }

  // host-side rows of ∂_k and δ_k
  struct Row {
    const Index* indices;
    const Sign* signs;
    Index size;
  };

  Row boundary_row(int k, Index cell) const {
    check_cell(k, cell, 1, dim_);
    return row(boundary_[static_cast<std::size_t>(k)], cell);
  }

  Row coboundary_row(int k, Index cell) const {
    check_cell(k, cell, 0, dim_ - 1);
    return row(coboundary_[static_cast<std::size_t>(k)], cell);
  }

  // kernel-facing views: capture CsrView by value in RAJA lambdas
  CsrView boundary_view(int k) {
    if (k < 1 || k > dim_) throw std::invalid_argument("boundary_view: dimension out of range");
    FrozenCsr& m = boundary_[static_cast<std::size_t>(k)];
    return CsrView{m.offsets.view(), m.indices.view(), m.signs.view()};
  }

  CsrView coboundary_view(int k) {
    if (k < 0 || k >= dim_) throw std::invalid_argument("coboundary_view: dimension out of range");
    FrozenCsr& m = coboundary_[static_cast<std::size_t>(k)];
    return CsrView{m.offsets.view(), m.indices.view(), m.signs.view()};
  }

  // --- closure, star, link -------------------------------------------------
  // Results are per-stratum sorted cell lists, sized dim()+1.
  //
  // Local: evaluated on the partition plus its ghost ring, correct when the
  // neighbourhood lies within freeze()'s halo_depth. No communication.

  // cl(σ): σ with all its faces, transitively.
  std::vector<std::vector<Index>> closure(int k, Index cell) const {
    check_cell(k, cell, 0, dim_);
    std::vector<std::set<Index>> acc(static_cast<std::size_t>(dim_) + 1);
    descend(k, cell, acc);
    return to_lists(acc);
  }

  // st(σ): σ with all its cofaces, transitively.
  std::vector<std::vector<Index>> star(int k, Index cell) const {
    check_cell(k, cell, 0, dim_);
    std::vector<std::set<Index>> acc(static_cast<std::size_t>(dim_) + 1);
    ascend(k, cell, acc);
    return to_lists(acc);
  }

  // lk(σ) = cl(st(σ)) \ st(cl(σ)). The link of an interior vertex of a
  // 2-manifold complex is a cycle.
  std::vector<std::vector<Index>> link(int k, Index cell) const {
    check_cell(k, cell, 0, dim_);
    std::vector<std::set<Index>> st(static_cast<std::size_t>(dim_) + 1);
    ascend(k, cell, st);
    std::vector<std::set<Index>> cl_st(static_cast<std::size_t>(dim_) + 1);
    for (int d = 0; d <= dim_; ++d) {
      for (const Index i : st[static_cast<std::size_t>(d)]) descend(d, i, cl_st);
    }
    std::vector<std::set<Index>> cl(static_cast<std::size_t>(dim_) + 1);
    descend(k, cell, cl);
    std::vector<std::set<Index>> st_cl(static_cast<std::size_t>(dim_) + 1);
    for (int d = 0; d <= dim_; ++d) {
      for (const Index i : cl[static_cast<std::size_t>(d)]) ascend(d, i, st_cl);
    }
    for (int d = 0; d <= dim_; ++d) {
      for (const Index i : st_cl[static_cast<std::size_t>(d)]) {
        cl_st[static_cast<std::size_t>(d)].erase(i);
      }
    }
    return to_lists(cl_st);
  }

 private:
  static Row row(const FrozenCsr& m, Index cell) {
    const Index lo = m.offsets[static_cast<std::size_t>(cell)];
    const Index hi = m.offsets[static_cast<std::size_t>(cell) + 1];
    return Row{m.indices.data() + lo, m.signs.data() + lo, hi - lo};
  }

  void check_cell(int k, Index cell, int k_min, int k_max) const {
    if (k < k_min || k > k_max) {
      throw std::invalid_argument("FrozenComplex: dimension out of range");
    }
    if (cell < 0 || cell >= count(k)) {
      throw std::out_of_range("FrozenComplex: cell index out of range");
    }
  }

  void descend(int k, Index cell, std::vector<std::set<Index>>& acc) const {
    if (!acc[static_cast<std::size_t>(k)].insert(cell).second) return;
    if (k == 0) return;
    const Row r = row(boundary_[static_cast<std::size_t>(k)], cell);
    for (Index m = 0; m < r.size; ++m) descend(k - 1, r.indices[m], acc);
  }

  void ascend(int k, Index cell, std::vector<std::set<Index>>& acc) const {
    if (!acc[static_cast<std::size_t>(k)].insert(cell).second) return;
    if (k == dim_) return;
    const Row r = row(coboundary_[static_cast<std::size_t>(k)], cell);
    for (Index m = 0; m < r.size; ++m) ascend(k + 1, r.indices[m], acc);
  }

  std::vector<std::vector<Index>> to_lists(const std::vector<std::set<Index>>& acc) const {
    std::vector<std::vector<Index>> out(acc.size());
    for (std::size_t d = 0; d < acc.size(); ++d) out[d].assign(acc[d].begin(), acc[d].end());
    return out;
  }

  int dim_;
  int halo_depth_;
  std::vector<Index> counts_;
  std::vector<FrozenCsr> boundary_;    // [k] valid for k >= 1
  std::vector<FrozenCsr> coboundary_;  // [k] valid for k < dim
};

inline FrozenComplex freeze(const Complex& c, int halo_depth = 1) {
  return FrozenComplex(c, halo_depth);
}

// The face poset with its order reversed. The dual k-cell is the primal
// (n−k)-cell — the same index — and ∂^dual_k = δ_{n−k}, so the view is
// zero-copy. This is the combinatorial half of the DEC dual mesh: the metric
// layer realizes dual cells (by barycentric subdivision), puts measures on
// them, and the discrete Hodge star maps primal k-cochains to dual
// (n−k)-cochains.
//
// Signs are the transpose coefficients. DEC conventions inserting
// (−1)^{k(n−k)} factors belong to the metric layer, as does the truncation of
// the dual of a boundary cell.
//
// A view; no communication implied.
class DualView {
 public:
  explicit DualView(const FrozenComplex& f) : f_(&f) {}

  int dim() const { return f_->dim(); }

  Index count(int k) const { return f_->count(f_->dim() - k); }

  // ∂^dual_k = δ_{n−k}
  FrozenComplex::Row boundary_row(int k, Index cell) const {
    return f_->coboundary_row(f_->dim() - k, cell);
  }

  // δ^dual_k = ∂_{n−k}
  FrozenComplex::Row coboundary_row(int k, Index cell) const {
    return f_->boundary_row(f_->dim() - k, cell);
  }

 private:
  const FrozenComplex* f_;
};

inline DualView dual(const FrozenComplex& f) { return DualView(f); }

}  // namespace graphos
