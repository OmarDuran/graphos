#pragma once

#include <vector>

#include "graphos/core/types.hpp"

namespace graphos_bench {

// Freudenthal–Kuhn triangulation of an n×n×n grid: six conforming 3-simplices
// per cube over (n+1)³ vertices. Deterministic and file-free at any size.
inline graphos::Index kuhn_vertex_count(int n) {
  return static_cast<graphos::Index>((n + 1) * (n + 1) * (n + 1));
}

inline std::vector<std::vector<graphos::Index>> kuhn_tets(int n) {
  using graphos::Index;
  const auto vid = [n](int i, int j, int k) {
    return static_cast<Index>(i * (n + 1) * (n + 1) + j * (n + 1) + k);
  };
  static constexpr int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                                      {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  std::vector<std::vector<Index>> tets;
  tets.reserve(static_cast<std::size_t>(n) * n * n * 6);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        for (const auto& p : perms) {
          int c[3] = {i, j, k};
          std::vector<Index> t{vid(c[0], c[1], c[2])};
          for (int s = 0; s < 3; ++s) {
            c[p[s]] += 1;
            t.push_back(vid(c[0], c[1], c[2]));
          }
          tets.push_back(std::move(t));
        }
      }
    }
  }
  return tets;
}

// a structured 2-complex: n×n quads split into 2n² 2-simplices
inline graphos::Index tri_grid_vertex_count(int n) {
  return static_cast<graphos::Index>((n + 1) * (n + 1));
}

inline std::vector<std::vector<graphos::Index>> tri_grid(int n) {
  using graphos::Index;
  const auto vid = [n](int i, int j) { return static_cast<Index>(i * (n + 1) + j); };
  std::vector<std::vector<Index>> tris;
  tris.reserve(static_cast<std::size_t>(n) * n * 2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      tris.push_back({vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)});
      tris.push_back({vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)});
    }
  }
  return tris;
}

}  // namespace graphos_bench
