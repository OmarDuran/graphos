#pragma once

#include <stdexcept>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"

namespace graphos {

// Classification of the (n-1)-cells (facets) of an n-complex by their
// number of top-dimensional cofaces. Each facet lands in exactly one
// marker; all four are marked over dimension n-1 and compose directly with
// subcomplex() and cut_along().
struct FacetClassification {
  Marker maximal;      // 0 cofaces: maximal (n-1)-cells — no proper coface
                       // (e.g. a detached interface domain after a cut).
                       // NOT "free" in the Whitehead sense: a free face has
                       // exactly ONE coface (see free_faces in collapse.hpp)
  Marker boundary;     // 1 coface: the topological boundary ∂Ω — the facets
                       // boundary conditions live on
  Marker interior;     // 2 cofaces: manifold interior facets
  Marker nonmanifold;  // 3+ cofaces: junction facets (DFN intersections,
                       // T-junctions) — also a mesh sanity signal
};

// Collective: the classification is computed from global coface counts
// (facets on partition boundaries count cofaces across ranks). P=1 today.
//
// subcomplex(c, classify_facets(c).boundary) extracts ∂Ω as its own
// complex with the embedding chain map.
//
// See also is_closed() below: ∂K = ∅ as a single predicate.
inline FacetClassification classify_facets(const Complex& c) {
  const int n = c.dim();
  if (n < 1) {
    throw std::invalid_argument("classify_facets: complex must have dimension >= 1");
  }
  FacetClassification out{Marker(c), Marker(c), Marker(c), Marker(c)};
  const CoboundaryOperator cob = coboundary(c, n - 1);
  for (Index f = 0; f < c.count(n - 1); ++f) {
    const Index degree =
        cob.offsets[static_cast<std::size_t>(f) + 1] - cob.offsets[static_cast<std::size_t>(f)];
    Marker& target = degree == 0   ? out.maximal
                     : degree == 1 ? out.boundary
                     : degree == 2 ? out.interior
                                   : out.nonmanifold;
    target.mark(n - 1, f);
  }
  return out;
}

// Is the complex CLOSED — is its boundary empty (∂K = ∅)? Evaluated at the
// effective top dimension (largest nonempty stratum): closed iff no facet
// there has exactly one coface. Note ∂(∂K) = ∅ always: the boundary
// subcomplex extracted via classify_facets + subcomplex is itself closed,
// checkable with this predicate.
//
// Collective. P=1 today.
inline bool is_closed(const Complex& c) {
  int d = c.dim();
  while (d > 0 && c.count(d) == 0) --d;
  if (d == 0) return true;
  const CoboundaryOperator cob = coboundary(c, d - 1);
  for (Index f = 0; f < c.count(d - 1); ++f) {
    if (cob.offsets[static_cast<std::size_t>(f) + 1] - cob.offsets[static_cast<std::size_t>(f)] ==
        1) {
      return false;
    }
  }
  return true;
}

}  // namespace graphos
