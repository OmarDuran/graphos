#pragma once

#include <cstdint>
#include <vector>

namespace graphos {

// Local cell index within one dimension of one complex. Global IDs (64-bit)
// exist only in the distributed indexing layer, never in kernels.
using Index = std::int32_t;

// Orientation coefficient. Stored entries are strictly -1 or +1; 0 never
// appears in CSR storage (absence of an arc encodes 0).
using Sign = std::int8_t;

inline constexpr Index invalid_index = -1;

// Declares that cell `from` is the same cell as `to`, with `from` oriented
// `rel_sign` relative to `to`. Which cells are "the same" is a geometric or
// user-level decision made outside graphos; this struct is how it enters the
// quotient.
struct Identification {
  Index from{invalid_index};
  Index to{invalid_index};
  Sign rel_sign{1};
};

// The cellular chain map induced by a topology-changing operation: for each
// dimension k, index[k][old] is the cell's index in the target complex, and
// sign[k][old] the orientation coefficient picked up along the way. A cell
// sent to zero (deleted, or merged away into its representative's star) has
// index invalid_index. Cochains (fields, DoFs) attached to cells must be
// transported through this map.
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

// Composition of chain maps: follows `first` then `second`. A cell sent to
// zero at either step is sent to zero by the composition.
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
