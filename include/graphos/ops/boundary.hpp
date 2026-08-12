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
  Marker free;         // 0 cofaces: maximal lower-dimensional cells (e.g. a
                       // detached interface/fracture domain after a cut)
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
inline FacetClassification classify_facets(const Complex& c) {
  const int n = c.dim();
  if (n < 1) {
    throw std::invalid_argument("classify_facets: complex must have dimension >= 1");
  }
  FacetClassification out{Marker(c), Marker(c), Marker(c), Marker(c)};
  const CoboundaryOperator cob = coboundary(c, n - 1);
  for (Index f = 0; f < c.count(n - 1); ++f) {
    const Index degree = cob.offsets[static_cast<std::size_t>(f) + 1] -
                         cob.offsets[static_cast<std::size_t>(f)];
    Marker& target = degree == 0   ? out.free
                     : degree == 1 ? out.boundary
                     : degree == 2 ? out.interior
                                   : out.nonmanifold;
    target.mark(n - 1, f);
  }
  return out;
}

}  // namespace graphos
