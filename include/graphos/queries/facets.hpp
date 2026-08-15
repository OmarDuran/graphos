#pragma once

#include <stdexcept>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"

namespace graphos {

// Partition of the facets of an n-complex by top-coface count. Each facet
// lands in exactly one marker, all marked over dimension n−1, composable
// with subcomplex() and cut_along().
struct FacetClassification {
  Marker maximal;      // 0 cofaces: maximal (n-1)-cells — no proper coface
                       // Not free in the Whitehead sense, which requires
                       // exactly one coface (see free_faces in collapse.hpp).
  Marker boundary;     // 1 coface: the topological boundary ∂Ω — the facets
                       // of ∂K
  Marker interior;     // 2 cofaces: manifold interior facets
  Marker nonmanifold;  // 3+ cofaces: junction facets (DFN intersections,
                       // T-junctions)
};

// Computed from global coface counts, so facets on partition boundaries
// count cofaces across ranks.
//
// subcomplex(c, classify_facets(c).boundary) extracts ∂K with its embedding
// chain map; is_closed() below is the predicate ∂K = ∅.
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

// ∂K = ∅, evaluated in the largest nonempty stratum: closed iff no facet
// there has exactly one coface. Since ∂∂K = ∅, the extracted boundary
// subcomplex satisfies this predicate.
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
