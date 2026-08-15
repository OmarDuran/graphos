#pragma once

#include <stdexcept>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos/queries/components.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos/queries/neighborhood.hpp"

namespace graphos {

// The PL amalgamation lemma: σ ∪ τ glued along a common (k−1)-ball in their
// boundaries is again a k-cell. The hypothesis splits into two combinatorial
// tests on (σ, τ):
//
//   1. Is ∂σ ∩ ∂τ ball-like? Decided as connected and Z₂-acyclic. Exact when
//      the common boundary is 1-dimensional (a connected acyclic graph is a
//      tree, an interval under the right co-degrees); in higher dimension
//      acyclicity is necessary but not sufficient, since H_* does not see π₁.
//   2. Do the closures meet properly, cl σ ∩ cl τ = cl(∂σ ∩ ∂τ)? Excess
//      intersection away from the common boundary pinches σ ∪ τ into a
//      non-manifold wedge.

// ∂σ ∩ ∂τ with the topology of its closure; `acyclic` is connected and
// Z₂-acyclic. Local.
struct CommonBoundary {
  Marker facets;
  Index n_facets{0};
  Index components{0};
  std::vector<Index> betti;
  bool acyclic{false};
};

inline CommonBoundary common_boundary(const Complex& c, int k, Index a, Index b) {
  if (k < 1 || k > c.dim()) {
    throw std::invalid_argument("common_boundary: dimension out of range");
  }
  if (a < 0 || a >= c.count(k) || b < 0 || b >= c.count(k) || a == b) {
    throw std::invalid_argument("common_boundary: invalid cell pair");
  }
  CommonBoundary out{Marker(c), 0, 0, {}, false};
  const BoundaryOperator& bnd = c.boundary(k);
  for (Index m = bnd.offsets[a]; m < bnd.offsets[a + 1]; ++m) {
    for (Index p = bnd.offsets[b]; p < bnd.offsets[b + 1]; ++p) {
      if (bnd.indices[m] == bnd.indices[p] && !out.facets.marked(k - 1, bnd.indices[m])) {
        out.facets.mark(k - 1, bnd.indices[m]);
        ++out.n_facets;
      }
    }
  }
  if (out.n_facets == 0) return out;
  const Complex shared = subcomplex(c, out.facets).complex;
  out.components = connected_components(shared).count;
  out.betti = betti_numbers_z2(shared);
  out.acyclic = out.components == 1;
  for (std::size_t q = 1; q < out.betti.size(); ++q) {
    out.acyclic = out.acyclic && out.betti[q] == 0;
  }
  return out;
}

// cl(σ) ∩ cl(τ) ∖ cl(∂σ ∩ ∂τ). Nonempty means the closures meet improperly:
// ∂(σ ∪ τ) would touch itself at exactly these cells. Local.
inline Marker excess_intersection(const Complex& c, int k, Index a, Index b) {
  const CommonBoundary shared = common_boundary(c, k, a, b);
  Marker ma(c), mb(c);
  ma.mark(k, a);
  mb.mark(k, b);
  const Marker cl_a = closure_of(c, ma);
  const Marker cl_b = closure_of(c, mb);
  const Marker cl_i = closure_of(c, shared.facets);
  Marker out(c);
  for (int d = 0; d < k; ++d) {
    for (Index i = 0; i < c.count(d); ++i) {
      if (cl_a.marked(d, i) && cl_b.marked(d, i) && !cl_i.marked(d, i)) out.mark(d, i);
    }
  }
  return out;
}

// σ ∪ τ is a k-cell iff ∂σ ∩ ∂τ is ball-like and the closures meet
// properly.
inline bool amalgamates_to_cell(const Complex& c, int k, Index a, Index b) {
  if (!common_boundary(c, k, a, b).acyclic) return false;
  const Marker excess = excess_intersection(c, k, a, b);
  for (int d = 0; d < k; ++d) {
    if (excess.marked_count(d) != 0) return false;
  }
  return true;
}

}  // namespace graphos
