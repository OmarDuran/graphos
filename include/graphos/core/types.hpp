#pragma once

#include <cstdint>
#include <vector>

namespace graphos {

// Index of a cell within one stratum of one complex. Global 64-bit IDs exist
// only in the distributed indexing layer, never in kernels.
using Index = std::int32_t;

// Incidence number [σ : τ] ∈ {−1, +1}. Zero is never stored; absence encodes
// it.
using Sign = std::int8_t;

inline constexpr Index invalid_index = -1;

// The relation σ ~ τ with relative orientation rel_sign, generating the
// equivalence a quotient collapses. Which cells are identified is decided
// outside graphos; this is how that decision enters the complex.
struct Identification {
  Index from{invalid_index};
  Index to{invalid_index};
  Sign rel_sign{1};
};

// The chain map f_* : C_k(C) → C_k(C′) induced on generators by an operation:
// index[k][σ] is the image cell, sign[k][σ] its orientation coefficient. A
// generator sent to 0 carries invalid_index. Cochains are transported by
// gathering through f_*.
struct ChainMap {
  std::vector<std::vector<Index>> index;
  std::vector<std::vector<Sign>> sign;

  static ChainMap sized(const std::vector<Index>& cell_counts) {
    ChainMap m;
    m.index.resize(cell_counts.size());
    m.sign.resize(cell_counts.size());
    for (std::size_t k = 0; k < cell_counts.size(); ++k) {
      m.index[k].resize(static_cast<std::size_t>(cell_counts[k]));
      m.sign[k].assign(static_cast<std::size_t>(cell_counts[k]), Sign{1});
    }
    return m;
  }
};

// (second ∘ first)_*, multiplying orientation coefficients. A generator sent
// to 0 by either factor is sent to 0 by the composite.
inline ChainMap compose(const ChainMap& first, const ChainMap& second) {
  ChainMap out;
  const std::size_t nk = first.index.size();
  out.index.resize(nk);
  out.sign.resize(nk);
  for (std::size_t k = 0; k < nk; ++k) {
    const std::size_t n = first.index[k].size();
    out.index[k].assign(n, invalid_index);
    out.sign[k].assign(n, Sign{1});
    for (std::size_t i = 0; i < n; ++i) {
      const Index mid = first.index[k][i];
      if (mid == invalid_index || k >= second.index.size() ||
          static_cast<std::size_t>(mid) >= second.index[k].size()) {
        continue;
      }
      const Index fin = second.index[k][static_cast<std::size_t>(mid)];
      out.index[k][i] = fin;
      if (fin != invalid_index) {
        out.sign[k][i] =
            static_cast<Sign>(first.sign[k][i] * second.sign[k][static_cast<std::size_t>(mid)]);
      }
    }
  }
  return out;
}

}  // namespace graphos
