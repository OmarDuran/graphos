#pragma once

#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/exec/array.hpp"

namespace graphos {

// Signed CSR incidence in persistent device-capable storage (exec::Array,
// the CHAI seam), plus the POD view kernels capture by value.
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

// The frozen half of the builder/frozen split. A Complex is the mutable
// host-side builder (attach_cell, and the ops that assemble new complexes);
// freezing it yields an immutable FrozenComplex whose boundary operators ∂_k
// AND derived signed coboundary operators δ_k live in device-capable arrays.
// This is the object queries, device kernels, and the NetworkX views sit on:
// its immutability is what makes zero-copy views and CHAI's copy-on-capture
// semantics sound.
//
// Topology-changing operations produce new complexes (via the builder path),
// which are then frozen again — freezing is the epoch boundary of the
// frozen-complexes/bulk-edits model.
class FrozenComplex {
 public:
  explicit FrozenComplex(const Complex& c)
      : dim_(c.dim()),
        counts_(c.counts()),
        boundary_(static_cast<std::size_t>(c.dim()) + 1),
        coboundary_(static_cast<std::size_t>(c.dim()) + 1) {
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

  int dim() const { return dim_; }

  Index count(int k) const {
    return (k >= 0 && k <= dim_) ? counts_[static_cast<std::size_t>(k)] : Index{0};
  }

  std::vector<Index> counts() const { return counts_; }

  // host-side row access into ∂_k / δ_k
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

  // kernel-facing views (capture the CsrView by value in RAJA lambdas)
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

  // --- combinatorial topology queries -------------------------------------
  // Results are per-dimension sorted cell lists, sized dim()+1.

  // cl(σ): σ and all its faces, recursively.
  std::vector<std::vector<Index>> closure(int k, Index cell) const {
    check_cell(k, cell, 0, dim_);
    std::vector<std::set<Index>> acc(static_cast<std::size_t>(dim_) + 1);
    descend(k, cell, acc);
    return to_lists(acc);
  }

  // st(σ): σ and all cells having σ as a face.
  std::vector<std::vector<Index>> star(int k, Index cell) const {
    check_cell(k, cell, 0, dim_);
    std::vector<std::set<Index>> acc(static_cast<std::size_t>(dim_) + 1);
    ascend(k, cell, acc);
    return to_lists(acc);
  }

  // lk(σ) = cl(st(σ)) \ st(cl(σ)): the boundary of a neighborhood of σ,
  // e.g. the link of an interior vertex of a 2-complex is a cycle.
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
  std::vector<Index> counts_;
  std::vector<FrozenCsr> boundary_;    // [k] valid for k >= 1
  std::vector<FrozenCsr> coboundary_;  // [k] valid for k < dim
};

inline FrozenComplex freeze(const Complex& c) { return FrozenComplex(c); }

}  // namespace graphos
