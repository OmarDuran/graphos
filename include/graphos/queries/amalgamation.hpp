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

// Cellular amalgamation criteria — the PL lemma behind sound coarsening:
// the union of two k-cells glued along a common (k−1)-ball in their
// boundaries is again a k-cell. Whether a pair (σ, τ) satisfies the
// hypothesis decomposes into two combinatorial questions:
//
//   1. Is the COMMON BOUNDARY ∂σ ∩ ∂τ ball-like? Decided here as
//      connected + Z₂-acyclic: exact for one-dimensional common boundaries
//      (an acyclic connected graph is a tree, hence an interval when the
//      co-degrees are right), the standard practical criterion above —
//      homology cannot see π₁, so acyclicity is necessary, not sufficient,
//      in higher dimension.
//   2. Do the closures intersect PROPERLY — cl σ ∩ cl τ = cl(∂σ ∩ ∂τ)?
//      Any EXCESS INTERSECTION (an isolated vertex or edge contact away
//      from the common boundary) would pinch the amalgamated cell into a
//      non-manifold wedge.

// The common boundary ∂σ ∩ ∂τ of two k-cells, with the topology of its
// closure. `acyclic` = connected and Z₂-acyclic (see caveat above).
//
// Local (a query about two given cells). P=1 today.
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

// The excess intersection of two k-cells: cl(σ) ∩ cl(τ) ∖ cl(∂σ ∩ ∂τ).
// Nonempty means the closures meet IMPROPERLY — the amalgamated cell's
// boundary would touch itself at exactly these cells (pinch points, in
// the non-manifold sense).
//
// Local. P=1 today.
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

// The amalgamation verdict: σ ∪ τ is again a k-cell when the common
// boundary is ball-like and the intersection is proper.
inline bool amalgamates_to_cell(const Complex& c, int k, Index a, Index b) {
  if (!common_boundary(c, k, a, b).acyclic) return false;
  const Marker excess = excess_intersection(c, k, a, b);
  for (int d = 0; d < k; ++d) {
    if (excess.marked_count(d) != 0) return false;
  }
  return true;
}

}  // namespace graphos
