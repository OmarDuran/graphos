#pragma once

#include <numeric>
#include <vector>

#include "graphos/core/types.hpp"

namespace graphos {

// Path-compressing union-find over the index range [0, n): the workhorse of
// component labeling. Serial by design — the device/kernel form of
// component labeling is iterative label propagation, which replaces this on
// GPU rather than porting it.
class UnionFind {
 public:
  explicit UnionFind(Index n) : parent_(static_cast<std::size_t>(n)) {
    std::iota(parent_.begin(), parent_.end(), Index{0});
  }

  Index find(Index a) {
    while (parent_[static_cast<std::size_t>(a)] != a) {
      parent_[static_cast<std::size_t>(a)] =
          parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(a)])];
      a = parent_[static_cast<std::size_t>(a)];
    }
    return a;
  }

  void unite(Index a, Index b) {
    a = find(a);
    b = find(b);
    if (a != b) parent_[static_cast<std::size_t>(b)] = a;
  }

  Index size() const { return static_cast<Index>(parent_.size()); }

 private:
  std::vector<Index> parent_;
};

}  // namespace graphos
