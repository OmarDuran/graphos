#pragma once

#include <map>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"

namespace graphos {

struct SubdivisionResult {
  // Simplicial, same dimension as the input: the order complex of the face
  // poset. k-cells are strict chains c_0 < c_1 < ... < c_k, with the
  // standard simplex boundary (remove the m-th element, sign (−1)^m).
  Complex complex;

  // The sd-vertex sitting at the barycenter of original cell (d, i) has
  // index vertex_offset[d] + i.
  std::vector<Index> vertex_offset;

  // The carrier of each sd k-cell: the original cell whose interior it
  // subdivides (the maximal element of its chain). This is the refinement
  // relation — cochain prolongation and dual-cell assembly key off it.
  std::vector<std::vector<int>> carrier_dim;
  std::vector<std::vector<Index>> carrier_index;

  // The SIGNED carrier: the incidence product ∏ᵢ [cᵢ : cᵢ₊₁] along the
  // cell's flag — the subdivision chain-map coefficient (up to a fixed
  // dimension-dependent convention factor). Nonzero (±1) exactly on FULL
  // flags (chain dimensions 0, 1, …, k), i.e. on the sd cells that carry
  // a cell of their own dimension: sd(σᵏ) = Σ carrier_sign · flag over
  // full flags through σ. Zero on dimension-jumping chains and on
  // barycenter vertices of positive-dimensional cells (they lie outside
  // the chain-map image), and on degenerate (non-regular) incidences.
  //
  // ORIENTATION TRANSFER: multiplying each top sd cell by its
  // carrier_sign turns a consistent orientation of the parent into a
  // consistent orientation of the subdivision — refinement preserves
  // orientation through this map rather than by convention.
  std::vector<std::vector<Sign>> carrier_sign;
};

// Barycentric subdivision sd(C): one new vertex per cell of C (its
// "barycenter" — a name only; graphos places no coordinates), one k-simplex
// per strict chain of k+1 cells in the face partial order. Works on ANY
// complex — polytopal, mixed-dimensional — and always yields a simplicial
// one; this is what makes the combinatorial dual geometrically realizable
// on the metric side, and the carrier relation doubles as a multigrid
// prolongation pattern.
//
// Deterministic: chains are enumerated in lexicographic order of their
// flattened cell ids, so sd indices are reproducible.
//
// Collective. P=1 today.
inline SubdivisionResult barycentric_subdivision(const Complex& c) {
  const int dim = c.dim();

  std::vector<Index> vertex_offset(static_cast<std::size_t>(dim) + 1, 0);
  for (int d = 1; d <= dim; ++d) {
    vertex_offset[static_cast<std::size_t>(d)] =
        vertex_offset[static_cast<std::size_t>(d) - 1] + c.count(d - 1);
  }
  const Index n_sd_vertices = vertex_offset[static_cast<std::size_t>(dim)] + c.count(dim);

  // upward incidence for chain extension: up[dx][dy] lists the dy-cells
  // having a given dx-cell in their closure
  std::vector<std::vector<Adjacency>> up(static_cast<std::size_t>(dim) + 1);
  for (int dx = 0; dx < dim; ++dx) {
    up[static_cast<std::size_t>(dx)].resize(static_cast<std::size_t>(dim) + 1);
    for (int dy = dx + 1; dy <= dim; ++dy) {
      up[static_cast<std::size_t>(dx)][static_cast<std::size_t>(dy)] = incidence(c, dx, dy);
    }
  }

  // chains per sd dimension, as flattened-id vectors (ascending, since flat
  // ids grow with cell dimension), plus their (dim, index) elements
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
            // incidence coefficient [c_top : y]: the (summed) sign of the
            // extended-over cell in ∂y — ±1 in regular complexes
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
