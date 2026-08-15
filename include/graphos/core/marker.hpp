#pragma once

#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

// A selection of cells: the argument form of every operation that acts on a
// chosen subcomplex.
//
// A Marker is local data with collective meaning. Each rank marks cells of its
// own partition and the union across ranks is the selection; at P = 1 the
// partition is the whole complex. Selection by predicate rather than by global
// index list is what makes the same program text rank-independent.
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

  // Marks every local k-cell satisfying pred. Rank-independent: the predicate
  // runs on each rank over its own cells.
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

  // locally marked k-cells
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

  // Guards that the marker was built against the complex being operated on.
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
