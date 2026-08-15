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
  // host → result: excised cells go to 0, the rest survive
  ChainMap map;
  // patch → result
  ChainMap patch_map;
};

// Surgery: excise st(S) and amalgamate a patch along the frontier —
// star_deletion followed by pushout, fused so the chain maps compose into
// single host→result and patch→result maps. Cutting a fault through a cell
// replaces it and its split faces without rebuilding the complex.
//
// The glue is a vertex correspondence (patch vertex → host vertex, indices as
// in the original complex), lifted through the strata by boundary-chain
// matching, so shared edges and faces need not be enumerated. Glue targets
// must survive the excision — lie outside cl(st S) — or the surgery is
// ill-posed and throws.
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
