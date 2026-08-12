#pragma once

#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/ops/disjoint_union.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos/ops/star_deletion.hpp"

namespace graphos {

struct ReplaceResult {
  Complex complex;
  // host -> result: cells of the excised region go to zero, everything
  // else survives into the repaired complex
  ChainMap map;
  // patch -> result: where each patch cell landed
  ChainMap patch_map;
};

// Local surgery, the primitive of INCREMENTAL IMPRINTING: excise the open
// star of the marked region and amalgamate a replacement patch along the
// frontier — star_deletion followed by the pushout, fused so the chain
// maps compose into single host→result and patch→result maps. This is the
// excision + amalgamation pattern executed as one operation: cutting a
// fault through a cell replaces it (and its split boundary faces) by the
// imprinted polyhedral pieces without rebuilding the rest of the complex.
//
// The glue is a vertex correspondence (patch vertex -> host vertex, host
// indices as in the ORIGINAL complex), lifted through the strata by
// boundary-chain matching — shared edges and faces of patch and frontier
// are found combinatorially and need not be enumerated. Glue targets must
// survive the excision (lie outside the closed star of the region);
// otherwise the surgery is ill-posed and throws.
//
// Logically collective: region and glue are supplied per-rank. P=1 today.
inline ReplaceResult replace(const Complex& c, const Marker& region, const Complex& patch,
                             const std::vector<Identification>& vertex_glue) {
  StarDeletionResult excised = star_deletion(c, region);
  DisjointUnionResult du = disjoint_union(excised.complex, patch);

  std::vector<Identification> pairs;
  pairs.reserve(vertex_glue.size());
  for (const Identification& g : vertex_glue) {
    if (g.from < 0 || g.from >= patch.count(0) || g.to < 0 || g.to >= c.count(0)) {
      throw std::out_of_range("replace: glue vertex out of range");
    }
    const Index host = excised.map.index[0][static_cast<std::size_t>(g.to)];
    if (host == invalid_index) {
      throw std::invalid_argument("replace: glue target was removed by the excision");
    }
    pairs.push_back({du.b_map.index[0][static_cast<std::size_t>(g.from)],
                     du.a_map.index[0][static_cast<std::size_t>(host)], g.rel_sign});
  }

  const auto lifted = lift_identifications(du.complex, pairs);
  QuotientResult q = quotient(du.complex, lifted);

  ReplaceResult res{std::move(q.complex), {}, {}};
  res.map = compose(excised.map, compose(du.a_map, q.map));
  res.patch_map = compose(du.b_map, q.map);
  return res;
}

}  // namespace graphos
