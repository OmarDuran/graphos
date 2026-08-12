#pragma once

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"

namespace graphos {

// Generalized adjacency between k-cells: two k-cells are adjacent when
// they share a via-cell. via = k-1 is facet-adjacency (the mesh dual
// graph); via = 0 is node-adjacency; any stratum in between works. Rows
// are sorted, duplicate-free, and exclude the cell itself.
//
// Local: computed per rank over its own cells; no communication implied.
inline Adjacency adjacency(const Complex& c, int k, int via) {
  if (k < 0 || k > c.dim() || via < 0 || via > c.dim()) {
    throw std::invalid_argument("adjacency: dimension out of range");
  }
  if (via == k) throw std::invalid_argument("adjacency: via dimension must differ from k");

  const Adjacency down = incidence(c, k, via);
  const Adjacency back = incidence(c, via, k);

  Adjacency out;
  out.offsets.reserve(static_cast<std::size_t>(c.count(k)) + 1);
  std::vector<Index> row;
  for (Index e = 0; e < c.count(k); ++e) {
    row.clear();
    for (Index m = down.offsets[static_cast<std::size_t>(e)];
         m < down.offsets[static_cast<std::size_t>(e) + 1]; ++m) {
      const Index v = down.indices[static_cast<std::size_t>(m)];
      row.insert(row.end(), back.indices.begin() + back.offsets[static_cast<std::size_t>(v)],
                 back.indices.begin() + back.offsets[static_cast<std::size_t>(v) + 1]);
    }
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
    row.erase(std::remove(row.begin(), row.end(), e), row.end());
    out.indices.insert(out.indices.end(), row.begin(), row.end());
    out.offsets.push_back(static_cast<Index>(out.indices.size()));
  }
  return out;
}

}  // namespace graphos
