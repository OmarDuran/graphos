#pragma once

#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

// A collective selection of cells — the argument form for every
// topology-changing operation that acts on a chosen set of cells.
//
// Distribution contract: a Marker is a LOCAL object with COLLECTIVE meaning.
// Each rank marks cells of its own partition (mark() takes local indices,
// predicates evaluate locally); the union across ranks is the logical
// selection. Under P=1 — the current implementation — the local partition is
// the whole complex. Because selection is expressed by marking rather than
// by naming global index lists, the same program text is meaningful at any
// rank count.
class Marker {
 public:
  explicit Marker(const Complex& c) : counts_(c.counts()) {
    flags_.resize(counts_.size());
    for (std::size_t k = 0; k < counts_.size(); ++k) {
      flags_[k].assign(static_cast<std::size_t>(counts_[k]), 0);
    }
  }

  int dim() const { return static_cast<int>(counts_.size()) - 1; }

  Marker& mark(int k, Index cell) {
    check(k, cell);
    flags_[static_cast<std::size_t>(k)][static_cast<std::size_t>(cell)] = 1;
    return *this;
  }

  // Marks every local k-cell satisfying pred(index). This is the
  // distribution-stable selection form: the predicate runs on each rank
  // over its own cells.
  template <typename Pred>
  Marker& mark_where(int k, Pred&& pred) {
    if (k < 0 || k > dim()) throw std::invalid_argument("Marker: dimension out of range");
    for (Index i = 0; i < counts_[static_cast<std::size_t>(k)]; ++i) {
      if (pred(i)) flags_[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = 1;
    }
    return *this;
  }

  bool marked(int k, Index cell) const {
    check(k, cell);
    return flags_[static_cast<std::size_t>(k)][static_cast<std::size_t>(cell)] != 0;
  }

  // number of locally marked k-cells
  Index marked_count(int k) const {
    if (k < 0 || k > dim()) throw std::invalid_argument("Marker: dimension out of range");
    Index n = 0;
    for (const char f : flags_[static_cast<std::size_t>(k)]) n += (f != 0);
    return n;
  }

  const std::vector<char>& flags(int k) const {
    if (k < 0 || k > dim()) throw std::invalid_argument("Marker: dimension out of range");
    return flags_[static_cast<std::size_t>(k)];
  }

  // Ops call this to guarantee the marker was built for the complex at hand.
  void validate_for(const Complex& c) const {
    if (counts_ != c.counts()) {
      throw std::invalid_argument("Marker: built for a different complex (cell counts differ)");
    }
  }

  static Marker from_cells(const Complex& c, const std::vector<std::vector<Index>>& cells) {
    Marker m(c);
    for (std::size_t k = 0; k < cells.size(); ++k) {
      if (static_cast<int>(k) > m.dim() && !cells[k].empty()) {
        throw std::invalid_argument("Marker::from_cells: dimension out of range");
      }
      for (const Index i : cells[k]) m.mark(static_cast<int>(k), i);
    }
    return m;
  }

 private:
  void check(int k, Index cell) const {
    if (k < 0 || k > dim()) throw std::invalid_argument("Marker: dimension out of range");
    if (cell < 0 || cell >= counts_[static_cast<std::size_t>(k)]) {
      throw std::out_of_range("Marker: cell index out of range");
    }
  }

  std::vector<Index> counts_;
  std::vector<std::vector<char>> flags_;
};

}  // namespace graphos
