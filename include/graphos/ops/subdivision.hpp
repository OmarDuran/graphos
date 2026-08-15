#pragma once

#include <map>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"

namespace graphos {

struct SubdivisionResult {
  // The order complex of the face poset: simplicial, of the same dimension.
  // A k-cell is a strict chain c₀ < c₁ < … < c_k under the simplex boundary
  // ∂[c₀…c_k] = Σ_m (−1)^m [c₀…ĉ_m…c_k].
  Complex complex;

  // The vertex carried by original cell (d, i) has index
  // vertex_offset[d] + i.
  std::vector<Index> vertex_offset;

  // The carrier of each k-cell: the original cell whose interior it
  // subdivides, the maximal element of its chain. This is the refinement
  // relation; prolongation and dual-cell assembly key off it.
  std::vector<std::vector<int>> carrier_dim;
  std::vector<std::vector<Index>> carrier_index;

  // The signed carrier ∏ᵢ [cᵢ : cᵢ₊₁] along the flag: the subdivision
  // chain-map coefficient, up to a fixed dimension-dependent factor. It is
  // ±1 exactly on full flags — chain dimensions 0, 1, …, k — so that
  // sd(σᵏ) = Σ carrier_sign · flag over the full flags through σ, and 0 on
  // dimension-jumping chains, on the vertices of positive-dimensional cells,
  // and on degenerate incidences.
  //
  // Multiplying each top cell by its carrier_sign carries a coherent
  // orientation of the parent to one of the subdivision, so refinement
  // preserves orientation through this map rather than by convention.
  std::vector<std::vector<Sign>> carrier_sign;
};

// sd(C): one vertex per cell of C — barycentre is a name only, since graphos
// places no coordinates — and one k-simplex per strict chain of k+1 cells in
// the face poset. Defined on any complex, polytopal or mixed-dimensional, and
// always simplicial, which is what makes the combinatorial dual realizable on
// the metric side. The carrier relation is also a prolongation pattern.
//
// Deterministic: chains are enumerated in lexicographic order of flattened
// cell ids, so indices are reproducible.
inline SubdivisionResult barycentric_subdivision(const Complex& c) {
  const int dim = c.dim();

  std::vector<Index> vertex_offset(static_cast<std::size_t>(dim) + 1, 0);
  for (int d = 1; d <= dim; ++d) {
    vertex_offset[static_cast<std::size_t>(d)] =
        vertex_offset[static_cast<std::size_t>(d) - 1] + c.count(d - 1);
  }
  const Index n_sd_vertices = vertex_offset[static_cast<std::size_t>(dim)] + c.count(dim);

  // upward incidence for chain extension: up[dx][dy] lists the dy-cells with
  // a given dx-cell in their closure
  std::vector<std::vector<Adjacency>> up(static_cast<std::size_t>(dim) + 1);
  for (int dx = 0; dx < dim; ++dx) {
    up[static_cast<std::size_t>(dx)].resize(static_cast<std::size_t>(dim) + 1);
    for (int dy = dx + 1; dy <= dim; ++dy) {
      up[static_cast<std::size_t>(dx)][static_cast<std::size_t>(dy)] = incidence(c, dx, dy);
    }
  }

  // chains per stratum as flattened-id vectors, ascending since flat ids grow
  // with dimension, with their (dim, index) elements
  struct Chain {
    std::vector<Index> flat;
    std::vector<std::pair<int, Index>> elems;
    bool full{false};  // dimensions are exactly 0..k
    Sign csign{0};     // ∏ incidence signs along the flag (0 unless full)
  };
  std::vector<std::vector<Chain>> chains(static_cast<std::size_t>(dim) + 1);
  std::vector<std::map<std::vector<Index>, Index>> chain_index(static_cast<std::size_t>(dim) + 1);

  for (int d = 0; d <= dim; ++d) {
    for (Index i = 0; i < c.count(d); ++i) {
      const Index flat = vertex_offset[static_cast<std::size_t>(d)] + i;
      chain_index[0][{flat}] = static_cast<Index>(chains[0].size());
      chains[0].push_back(Chain{{flat}, {{d, i}}, d == 0, d == 0 ? Sign{1} : Sign{0}});
    }
  }
  for (int k = 1; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    for (const Chain& base : chains[sk - 1]) {
      const auto [d_top, i_top] = base.elems.back();
      for (int dy = d_top + 1; dy <= dim; ++dy) {
        const Adjacency& a = up[static_cast<std::size_t>(d_top)][static_cast<std::size_t>(dy)];
        for (Index m = a.offsets[i_top]; m < a.offsets[i_top + 1]; ++m) {
          const Index y = a.indices[static_cast<std::size_t>(m)];
          Chain ext = base;
          ext.flat.push_back(vertex_offset[static_cast<std::size_t>(dy)] + y);
          ext.elems.emplace_back(dy, y);
          ext.full = false;
          ext.csign = 0;
          if (base.full && dy == d_top + 1) {
            // [c_top : y], the summed incidence of the extended-over cell in
            // ∂y; ±1 when the complex is regular
            const BoundaryOperator& by = c.boundary(dy);
            int step = 0;
            for (Index q = by.offsets[y]; q < by.offsets[y + 1]; ++q) {
              if (by.indices[q] == i_top) step += by.signs[q];
            }
            if (step == 1 || step == -1) {
              ext.full = true;
              ext.csign = static_cast<Sign>(base.csign * step);
            }
          }
          chain_index[sk][ext.flat] = static_cast<Index>(chains[sk].size());
          chains[sk].push_back(std::move(ext));
        }
      }
    }
  }

  Complex out(dim);
  out.attach_vertices(n_sd_vertices);
  std::vector<std::vector<int>> carrier_dim(static_cast<std::size_t>(dim) + 1);
  std::vector<std::vector<Index>> carrier_index(static_cast<std::size_t>(dim) + 1);
  std::vector<std::vector<Sign>> carrier_sign(static_cast<std::size_t>(dim) + 1);
  for (const Chain& s : chains[0]) {
    carrier_dim[0].push_back(s.elems.back().first);
    carrier_index[0].push_back(s.elems.back().second);
    carrier_sign[0].push_back(s.csign);
  }

  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  std::vector<Index> face;
  for (int k = 1; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    for (const Chain& s : chains[sk]) {
      row_idx.clear();
      row_sg.clear();
      for (std::size_t m = 0; m < s.flat.size(); ++m) {
        face = s.flat;
        face.erase(face.begin() + static_cast<std::ptrdiff_t>(m));
        row_idx.push_back(chain_index[sk - 1].at(face));
        row_sg.push_back(m % 2 == 0 ? Sign{1} : Sign{-1});
      }
      out.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
      carrier_dim[sk].push_back(s.elems.back().first);
      carrier_index[sk].push_back(s.elems.back().second);
      carrier_sign[sk].push_back(s.csign);
    }
  }

  return SubdivisionResult{std::move(out), std::move(vertex_offset), std::move(carrier_dim),
                           std::move(carrier_index), std::move(carrier_sign)};
}

}  // namespace graphos
