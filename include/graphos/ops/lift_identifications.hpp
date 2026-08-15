#pragma once

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

// Extends a vertex-level identification upward through the strata: a k-cell
// x is identified with a k-cell y when every face of x is (already) mapped
// and the image of x's boundary chain equals y's boundary chain up to a
// uniform orientation flip (which is reported). Processing runs bottom-up,
// so edge identifications derive from the vertex pairs, face
// identifications from the edges, and so on.
//
// This is the missing piece for PERIODIC boundary conditions: pair the
// slave boundary vertices with their master images (a geometric decision),
// lift, and quotient — a square becomes a cylinder, a twisted pairing a
// Möbius band. It is the cross-level generalization of find_parallel_cells,
// the same matching pushout's deduplication performs implicitly.
//
// Preconditions: vertex pairs must point at final representatives (no
// chains a→b→c); cells whose faces are only partially mapped simply do not
// lift, which is the correct behavior at the rim of a periodic patch.
//
// Logically collective: every rank supplies the vertex pairs it owns; the
// result feeds quotient(). P=1 today.
inline std::vector<std::vector<Identification>> lift_identifications(
    const Complex& c, const std::vector<Identification>& vertex_identifications) {
  const int dim = c.dim();
  std::vector<std::vector<Identification>> out(static_cast<std::size_t>(dim) + 1);

  std::vector<std::map<Index, std::pair<Index, Sign>>> lifted(static_cast<std::size_t>(dim) + 1);
  for (const Identification& id : vertex_identifications) {
    if (id.from < 0 || id.from >= c.count(0) || id.to < 0 || id.to >= c.count(0)) {
      throw std::out_of_range("lift_identifications: vertex index out of range");
    }
    if (id.rel_sign != 1 && id.rel_sign != -1) {
      throw std::invalid_argument("lift_identifications: rel_sign must be +/-1");
    }
    if (id.from == id.to) continue;
    out[0].push_back(id);
    lifted[0][id.from] = {id.to, id.rel_sign};
  }

  std::vector<std::pair<Index, Sign>> row;
  for (int k = 1; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const BoundaryOperator& bnd = c.boundary(k);

    // targets, keyed by sorted boundary index list (degenerate rows skipped)
    std::map<std::vector<Index>, std::pair<Index, std::vector<Sign>>> targets;
    for (Index y = 0; y < c.count(k); ++y) {
      row.clear();
      for (Index m = bnd.offsets[y]; m < bnd.offsets[y + 1]; ++m) {
        row.emplace_back(bnd.indices[m], bnd.signs[m]);
      }
      if (row.empty()) continue;
      std::sort(row.begin(), row.end(),
                [](const auto& u, const auto& v) { return u.first < v.first; });
      bool degenerate = false;
      for (std::size_t m = 1; m < row.size(); ++m) {
        if (row[m].first == row[m - 1].first) degenerate = true;
      }
      if (degenerate) continue;
      std::vector<Index> key;
      std::vector<Sign> sg;
      for (const auto& [i, s] : row) {
        key.push_back(i);
        sg.push_back(s);
      }
      targets.emplace(std::move(key), std::make_pair(y, std::move(sg)));
    }

    // lift each fully mapped cell through the image of its boundary chain
    for (Index x = 0; x < c.count(k); ++x) {
      row.clear();
      bool fully_mapped = bnd.offsets[x] < bnd.offsets[x + 1];
      for (Index m = bnd.offsets[x]; m < bnd.offsets[x + 1] && fully_mapped; ++m) {
        const auto it = lifted[sk - 1].find(bnd.indices[m]);
        if (it == lifted[sk - 1].end()) {
          fully_mapped = false;
          break;
        }
        row.emplace_back(it->second.first, static_cast<Sign>(bnd.signs[m] * it->second.second));
      }
      if (!fully_mapped) continue;
      std::sort(row.begin(), row.end(),
                [](const auto& u, const auto& v) { return u.first < v.first; });
      bool degenerate = false;
      for (std::size_t m = 1; m < row.size(); ++m) {
        if (row[m].first == row[m - 1].first) degenerate = true;
      }
      if (degenerate) continue;
      std::vector<Index> key;
      std::vector<Sign> sg;
      for (const auto& [i, s] : row) {
        key.push_back(i);
        sg.push_back(s);
      }
      const auto it = targets.find(key);
      if (it == targets.end()) continue;
      const auto& [y, tg] = it->second;
      if (y == x) continue;
      const int r = static_cast<int>(sg[0]) * static_cast<int>(tg[0]);
      bool consistent = true;
      for (std::size_t m = 0; m < sg.size(); ++m) {
        if (static_cast<int>(sg[m]) != r * static_cast<int>(tg[m])) {
          consistent = false;
          break;
        }
      }
      if (!consistent) continue;
      out[sk].push_back({x, y, static_cast<Sign>(r)});
      lifted[sk][x] = {y, static_cast<Sign>(r)};
    }
  }
  return out;
}

}  // namespace graphos
