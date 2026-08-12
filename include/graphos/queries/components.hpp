#pragma once

#include <stdexcept>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/core/union_find.hpp"

namespace graphos {

// Component labels of the k-cells: label[i] in [0, count), deterministic
// (components numbered by their lowest-indexed cell), so results are
// independent of traversal and — later — of rank count.
struct ComponentLabels {
  std::vector<Index> label;
  Index count{0};
};

namespace detail {

inline ComponentLabels relabel(UnionFind& uf) {
  const Index n = uf.size();
  ComponentLabels out;
  out.label.assign(static_cast<std::size_t>(n), invalid_index);
  std::vector<Index> root_id(static_cast<std::size_t>(n), invalid_index);
  for (Index i = 0; i < n; ++i) {
    const Index r = uf.find(i);
    if (root_id[static_cast<std::size_t>(r)] == invalid_index) {
      root_id[static_cast<std::size_t>(r)] = out.count++;
    }
    out.label[static_cast<std::size_t>(i)] = root_id[static_cast<std::size_t>(r)];
  }
  return out;
}

inline ComponentLabels components_via(const Complex& c, int k, int via,
                                      const Marker* exclude_via) {
  if (k < 0 || k > c.dim() || via < 0 || via > c.dim()) {
    throw std::invalid_argument("connected_components: dimension out of range");
  }
  if (via == k) {
    throw std::invalid_argument("connected_components: via dimension must differ from k");
  }
  if (exclude_via != nullptr) exclude_via->validate_for(c);

  const Adjacency inc = incidence(c, via, k);  // via-cell -> incident k-cells
  UnionFind uf(c.count(k));
  for (Index v = 0; v < c.count(via); ++v) {
    if (exclude_via != nullptr && exclude_via->marked(via, v)) continue;
    const Index lo = inc.offsets[static_cast<std::size_t>(v)];
    const Index hi = inc.offsets[static_cast<std::size_t>(v) + 1];
    for (Index m = lo + 1; m < hi; ++m) {
      uf.unite(inc.indices[static_cast<std::size_t>(lo)], inc.indices[static_cast<std::size_t>(m)]);
    }
  }
  return relabel(uf);
}

}  // namespace detail

// Connected components of the k-cells, two cells adjacent iff incident to a
// common via-cell (via < k: shared faces, the usual case — e.g. top cells
// through facets identifies subdomains; via > k: shared cofaces — e.g.
// vertices through edges).
//
// Collective: labels are globally consistent (under distribution this is a
// distributed label propagation; P=1 today).
inline ComponentLabels connected_components(const Complex& c, int k, int via) {
  return detail::components_via(c, k, via, nullptr);
}

// As above, with marked via-cells excluded as connectors: the "sides of a
// cut" computation — exclude the interface facets and the components of the
// top cells are the subdomains the interface separates. Only marks of
// dimension `via` are consulted.
inline ComponentLabels connected_components(const Complex& c, int k, int via,
                                            const Marker& exclude_via) {
  return detail::components_via(c, k, via, &exclude_via);
}

// Components of the whole complex across all dimensions, cells connected by
// incidence: β0 of the (possibly mixed-dimensional) complex. label[k][i] is
// the component of the i-th k-cell; detached lower-dimensional domains
// (e.g. extracted fracture networks) count as their own components.
//
// Collective.
struct ComplexComponents {
  std::vector<std::vector<Index>> label;
  Index count{0};
};

inline ComplexComponents connected_components(const Complex& c) {
  const int dim = c.dim();
  std::vector<Index> base(static_cast<std::size_t>(dim) + 2, 0);
  for (int k = 0; k <= dim; ++k) {
    base[static_cast<std::size_t>(k) + 1] = base[static_cast<std::size_t>(k)] + c.count(k);
  }
  UnionFind uf(base[static_cast<std::size_t>(dim) + 1]);
  for (int k = 1; k <= dim; ++k) {
    const BoundaryOperator& bnd = c.boundary(k);
    for (Index e = 0; e < c.count(k); ++e) {
      for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
        uf.unite(base[static_cast<std::size_t>(k)] + e,
                 base[static_cast<std::size_t>(k) - 1] + bnd.indices[m]);
      }
    }
  }
  const ComponentLabels flat = detail::relabel(uf);

  ComplexComponents out;
  out.count = flat.count;
  out.label.resize(static_cast<std::size_t>(dim) + 1);
  for (int k = 0; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    out.label[sk].assign(flat.label.begin() + base[sk], flat.label.begin() + base[sk + 1]);
  }
  return out;
}

}  // namespace graphos
