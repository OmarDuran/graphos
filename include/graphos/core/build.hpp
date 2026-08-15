#pragma once

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

// Constructors from mesh connectivity: intermediate strata are derived (each
// edge, each face once) under a deterministic orientation convention, so
// ∂∘∂ = 0 holds by construction.
//
// Polytopal, so the general forms are primitive: a 2-cell is an arbitrary
// vertex cycle, a 3-cell an arbitrary set of polygonal faces (the
// VTK-polyhedron description). Simplices are a wrapper, not a special case.
//
// Conventions, independent of insertion order:
//  - a 1-cell runs from the smaller to the larger vertex index;
//  - a derived 2-cell carries its canonical cycle (rotated to start at its
//    smallest vertex, directed toward that vertex's smaller neighbour); a
//    coface sets [σ : τ] by comparing its traversal against that cycle;
//  - a 2-cell is keyed by its vertex set; two cells giving that set
//    cyclically incompatible orderings is an error.
// Orientability of the result is not asserted here; orient() normalizes it.
//
// Collective (each rank builds its local partition; P = 1 today), serial
// host code.

namespace detail {

// keys for the derivation pools: vertex pairs pack into one word, tuples hash
// by FNV-1a
struct PairHash {
  std::size_t operator()(const std::pair<Index, Index>& p) const {
    const std::uint64_t w = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(p.first)) << 32) |
                            static_cast<std::uint32_t>(p.second);
    return std::hash<std::uint64_t>{}(w);
  }
};

struct VertsHash {
  std::size_t operator()(const std::vector<Index>& v) const {
    std::size_t h = 1469598103934665603ull;
    for (const Index x : v) {
      h ^= static_cast<std::size_t>(static_cast<std::uint32_t>(x));
      h *= 1099511628211ull;
    }
    return h;
  }
};


inline void check_cycle(const std::vector<Index>& cyc, Index n_vertices, const char* what) {
  if (cyc.size() < 3) {
    throw std::invalid_argument(std::string(what) + ": cycle needs at least 3 vertices");
  }
  for (const Index v : cyc) {
    if (v < 0 || v >= n_vertices) {
      throw std::out_of_range(std::string(what) + ": vertex index out of range");
    }
  }
  std::vector<Index> s = cyc;
  std::sort(s.begin(), s.end());
  for (std::size_t i = 1; i < s.size(); ++i) {
    if (s[i] == s[i - 1]) {
      throw std::invalid_argument(std::string(what) + ": repeated vertex in cycle");
    }
  }
}

// canonical form of a vertex cycle, and the input's orientation relative to
// it (±1)
inline std::pair<std::vector<Index>, Sign> canonical_cycle(const std::vector<Index>& cyc) {
  const std::size_t n = cyc.size();
  std::size_t m = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (cyc[i] < cyc[m]) m = i;
  }
  const Index next = cyc[(m + 1) % n];
  const Index prev = cyc[(m + n - 1) % n];
  std::vector<Index> canon(n);
  if (next < prev) {
    for (std::size_t i = 0; i < n; ++i) canon[i] = cyc[(m + i) % n];
    return {std::move(canon), Sign{1}};
  }
  for (std::size_t i = 0; i < n; ++i) canon[i] = cyc[(m + n - i) % n];
  return {std::move(canon), Sign{-1}};
}

// the directed vertex cycle of a 2-cell, reconstructed from its signed edge
// chain (a 1-cell with +1 runs low → high). Cells built by cycle_row store
// their canonical cycle, so canonicalizing the reconstruction recovers it.
inline std::vector<Index> stored_face_cycle(const Complex& c, Index face) {
  const BoundaryOperator& b2 = c.boundary(2);
  const BoundaryOperator& b1 = c.boundary(1);
  std::vector<std::pair<Index, Index>> arcs;  // directed edges
  for (Index m = b2.offsets[face]; m < b2.offsets[face + 1]; ++m) {
    const Index e = b2.indices[m];
    const Index a = b1.indices[b1.offsets[e]];
    const Index b = b1.indices[b1.offsets[e] + 1];
    arcs.emplace_back(b2.signs[m] > 0 ? a : b, b2.signs[m] > 0 ? b : a);
  }
  std::vector<Index> cyc{arcs.front().first};
  Index cur = arcs.front().second;
  while (cur != cyc.front()) {
    cyc.push_back(cur);
    for (const auto& [from, to] : arcs) {
      if (from == cur) {
        cur = to;
        break;
      }
    }
  }
  return canonical_cycle(cyc).first;
}

class EdgePool {
 public:
  explicit EdgePool(Complex& c) : c_(&c) {}

  void reserve(std::size_t n) { edges_.reserve(n); }

  Index get(Index u, Index v) {
    const auto key = std::minmax(u, v);
    const auto it = edges_.find(key);
    if (it != edges_.end()) return it->second;
    const Index idx =
        c_->attach_cell(1, {key.first, key.second}, {-1, +1});  // low -> high
    edges_.emplace(key, idx);
    return idx;
  }

 private:
  Complex* c_;
  std::unordered_map<std::pair<Index, Index>, Index, PairHash> edges_;
};

// the ∂₂ row of a cycle: consecutive 1-cells signed by traversal direction
inline void cycle_row(EdgePool& edges, const std::vector<Index>& cyc, std::vector<Index>& row_idx,
                      std::vector<Sign>& row_sg) {
  row_idx.clear();
  row_sg.clear();
  for (std::size_t i = 0; i < cyc.size(); ++i) {
    const Index u = cyc[i];
    const Index v = cyc[(i + 1) % cyc.size()];
    row_idx.push_back(edges.get(u, v));
    row_sg.push_back(u < v ? Sign{1} : Sign{-1});
  }
}

}  // namespace detail

// 1-complex: one 1-cell per pair, not deduplicated.
inline Complex from_edges(Index n_vertices, const std::vector<std::vector<Index>>& segments) {
  Complex c(1);
  c.attach_vertices(n_vertices);
  for (const std::vector<Index>& s : segments) {
    if (s.size() != 2 || s[0] == s[1]) {
      throw std::invalid_argument("from_edges: each segment is two distinct vertices");
    }
    c.attach_cell(1, {s[0], s[1]}, {-1, +1});
  }
  return c;
}

// 2-complex: 2-cells as vertex cycles, of any arity and any mix.
inline Complex from_polygons(Index n_vertices, const std::vector<std::vector<Index>>& polygons) {
  Complex c(2);
  c.attach_vertices(n_vertices);
  detail::EdgePool edges(c);
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (const std::vector<Index>& poly : polygons) {
    detail::check_cycle(poly, n_vertices, "from_polygons");
    detail::cycle_row(edges, poly, row_idx, row_sg);
    c.attach_cell(2, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
  }
  return c;
}

// 3-complex: each 3-cell a set of 2-cells given as vertex cycles. A face
// shared by two cofaces is derived once, each referencing it with the
// incidence number of its own winding.
//
// Flat CSR: cell_offsets indexes face_offsets, which indexes face_vertices.
// Large meshes enter here — no nested vectors on the input side.
inline Complex from_polyhedra(Index n_vertices, std::span<const Index> cell_offsets,
                              std::span<const Index> face_offsets,
                              std::span<const Index> face_vertices) {
  Complex c(3);
  c.attach_vertices(n_vertices);
  detail::EdgePool edges(c);
  edges.reserve(face_vertices.size() / 2);  // arcs ~ 2x edges on conformal meshes
  // sorted vertex set → face index. The canonical cycle is recoverable from
  // the complex and so is not stored.
  std::unordered_map<std::vector<Index>, Index, detail::VertsHash> faces;
  faces.reserve(face_offsets.empty() ? 0 : face_offsets.size() - 1);  // upper bound
  std::vector<Index> row_idx, cell_idx, cyc, key;
  std::vector<Sign> row_sg, cell_sg;
  for (std::size_t e = 0; e + 1 < cell_offsets.size(); ++e) {
    if (cell_offsets[e] == cell_offsets[e + 1]) {
      throw std::invalid_argument("from_polyhedra: cell without faces");
    }
    cell_idx.clear();
    cell_sg.clear();
    for (Index fc = cell_offsets[e]; fc < cell_offsets[e + 1]; ++fc) {
      cyc.assign(face_vertices.begin() + face_offsets[fc],
                 face_vertices.begin() + face_offsets[fc + 1]);
      detail::check_cycle(cyc, n_vertices, "from_polyhedra");
      auto [canon, rel] = detail::canonical_cycle(cyc);
      key = canon;
      std::sort(key.begin(), key.end());
      const auto it = faces.find(key);
      Index fi;
      if (it == faces.end()) {
        detail::cycle_row(edges, canon, row_idx, row_sg);
        fi = c.attach_cell(2, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
        faces.emplace(std::move(key), fi);
        key = {};
      } else {
        fi = it->second;
        if (detail::stored_face_cycle(c, fi) != canon) {
          throw std::invalid_argument(
              "from_polyhedra: two cells describe the same vertex set with "
              "cyclically incompatible face orderings");
        }
      }
      cell_idx.push_back(fi);
      cell_sg.push_back(rel);
    }
    c.attach_cell(3, std::span<const Index>(cell_idx), std::span<const Sign>(cell_sg));
  }
  return c;
}

inline Complex from_polyhedra(Index n_vertices,
                              const std::vector<std::vector<std::vector<Index>>>& cells) {
  std::vector<Index> cell_offsets{0}, face_offsets{0}, face_vertices;
  for (const auto& cell : cells) {
    for (const auto& cyc : cell) {
      face_vertices.insert(face_vertices.end(), cyc.begin(), cyc.end());
      face_offsets.push_back(static_cast<Index>(face_vertices.size()));
    }
    cell_offsets.push_back(static_cast<Index>(face_offsets.size()) - 1);
  }
  return from_polyhedra(n_vertices, cell_offsets, face_offsets, face_vertices);
}

namespace detail {

// sorted (k+1)-tuple → k-cell, built on demand with the simplex boundary
// ∂[v₀…v_k] = Σ_i (−1)^i [v₀…v̂_i…v_k]
inline Index get_simplex(Complex& c,
                         std::vector<std::unordered_map<std::vector<Index>, Index, VertsHash>>& pool,
                         int k, const std::vector<Index>& sorted_verts) {
  if (k == 0) return sorted_verts[0];
  const auto it = pool[static_cast<std::size_t>(k)].find(sorted_verts);
  if (it != pool[static_cast<std::size_t>(k)].end()) return it->second;
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  std::vector<Index> face;
  for (std::size_t i = 0; i < sorted_verts.size(); ++i) {
    face = sorted_verts;
    face.erase(face.begin() + static_cast<std::ptrdiff_t>(i));
    row_idx.push_back(get_simplex(c, pool, k - 1, face));
    row_sg.push_back(i % 2 == 0 ? Sign{1} : Sign{-1});
  }
  const Index idx =
      c.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
  pool[static_cast<std::size_t>(k)].emplace(sorted_verts, idx);
  return idx;
}

// parity of the permutation taking `sorted` to `input` (±1)
inline Sign permutation_parity(const std::vector<Index>& input) {
  std::vector<Index> v = input;
  int swaps = 0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    std::size_t m = i;
    for (std::size_t j = i + 1; j < v.size(); ++j) {
      if (v[j] < v[m]) m = j;
    }
    if (m != i) {
      std::swap(v[i], v[m]);
      ++swaps;
    }
  }
  return swaps % 2 == 0 ? Sign{1} : Sign{-1};
}

}  // namespace detail

// Simplicial connectivity in any dimension d ≥ 1: a d-cell is a (d+1)-tuple.
// Every (k+1)-subset becomes one k-cell, oriented by increasing vertex order
// under the simplex boundary; a d-cell carries the parity of its input
// ordering, so consistently ordered simplices come out coherently oriented.
// ∂∘∂ = 0 by construction in every stratum.
inline Complex from_simplices(int dim, Index n_vertices,
                              const std::vector<std::vector<Index>>& cells) {
  if (dim < 1) throw std::invalid_argument("from_simplices: dim must be >= 1");
  Complex c(dim);
  c.attach_vertices(n_vertices);
  std::vector<std::unordered_map<std::vector<Index>, Index, detail::VertsHash>> pool(
      static_cast<std::size_t>(dim));

  std::vector<Index> sorted;
  std::vector<Index> face;
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (const std::vector<Index>& cell : cells) {
    if (static_cast<int>(cell.size()) != dim + 1) {
      throw std::invalid_argument("from_simplices: each cell needs dim+1 vertices");
    }
    for (const Index v : cell) {
      if (v < 0 || v >= n_vertices) {
        throw std::out_of_range("from_simplices: vertex index out of range");
      }
    }
    sorted = cell;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 1; i < sorted.size(); ++i) {
      if (sorted[i] == sorted[i - 1]) {
        throw std::invalid_argument("from_simplices: repeated vertex in cell");
      }
    }
    const Sign parity = detail::permutation_parity(cell);
    row_idx.clear();
    row_sg.clear();
    for (std::size_t i = 0; i < sorted.size(); ++i) {
      face = sorted;
      face.erase(face.begin() + static_cast<std::ptrdiff_t>(i));
      row_idx.push_back(detail::get_simplex(c, pool, dim - 1, face));
      row_sg.push_back(static_cast<Sign>(parity * (i % 2 == 0 ? 1 : -1)));
    }
    c.attach_cell(dim, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
  }
  return c;
}

}  // namespace graphos
