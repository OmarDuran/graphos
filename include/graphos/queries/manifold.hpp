#pragma once

#include <map>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/core/incidence.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/core/union_find.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos/queries/components.hpp"
#include "graphos/queries/facets.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos/queries/neighborhood.hpp"

namespace graphos {

// Necessary conditions for a combinatorial n-manifold with boundary.
// Recognition is undecidable in general; these are the decidable local ones:
//
//   pure            every cell is a face of a top cell (a mixed-dimensional
//                   complex fails by construction; test extracted pure
//                   subcomplexes instead)
//   facet_condition every facet has one or two top cofaces (no book
//                   junctions)
//   links_connected for every σ with dim σ ≤ n−2, st(σ) is connected through
//                   the facets containing σ (no pinch points)
//
// `offending` marks each cell that failed.
struct ManifoldReport {
  bool manifold_like{true};  // all of the below
  bool pure{true};
  bool facet_condition{true};
  bool links_connected{true};
  Marker offending;
};

// Collective. P=1 today.
inline ManifoldReport check_manifold(const Complex& c) {
  const int n = c.dim();
  ManifoldReport rep{true, true, true, true, Marker(c)};
  if (n < 1) return rep;  // a set of points is a 0-manifold

  // facets: coface counts
  const CoboundaryOperator cob = coboundary(c, n - 1);
  for (Index f = 0; f < c.count(n - 1); ++f) {
    const Index deg =
        cob.offsets[static_cast<std::size_t>(f) + 1] - cob.offsets[static_cast<std::size_t>(f)];
    if (deg == 0) {
      rep.pure = false;
      rep.offending.mark(n - 1, f);
    } else if (deg > 2) {
      rep.facet_condition = false;
      rep.offending.mark(n - 1, f);
    }
  }

  // lower cells: purity and link connectivity
  for (int k = 0; k <= n - 2; ++k) {
    const Adjacency tops = incidence(c, k, n);
    const Adjacency facets = incidence(c, k, n - 1);
    for (Index x = 0; x < c.count(k); ++x) {
      const Index t_lo = tops.offsets[static_cast<std::size_t>(x)];
      const Index t_hi = tops.offsets[static_cast<std::size_t>(x) + 1];
      if (t_lo == t_hi) {
        rep.pure = false;
        rep.offending.mark(k, x);
        continue;
      }
      std::map<Index, Index> slot;
      for (Index m = t_lo; m < t_hi; ++m) {
        slot[tops.indices[static_cast<std::size_t>(m)]] = m - t_lo;
      }
      UnionFind uf(t_hi - t_lo);
      for (Index m = facets.offsets[static_cast<std::size_t>(x)];
           m < facets.offsets[static_cast<std::size_t>(x) + 1]; ++m) {
        const Index f = facets.indices[static_cast<std::size_t>(m)];
        const Index lo = cob.offsets[static_cast<std::size_t>(f)];
        const Index hi = cob.offsets[static_cast<std::size_t>(f) + 1];
        for (Index p = lo + 1; p < hi; ++p) {
          uf.unite(slot.at(cob.indices[static_cast<std::size_t>(lo)]),
                   slot.at(cob.indices[static_cast<std::size_t>(p)]));
        }
      }
      Index n_comp = 0;
      for (Index t = 0; t < t_hi - t_lo; ++t) {
        if (uf.find(t) == t) ++n_comp;
      }
      if (n_comp > 1) {
        rep.links_connected = false;
        rep.offending.mark(k, x);
      }
    }
  }

  rep.manifold_like = rep.pure && rep.facet_condition && rep.links_connected;
  return rep;
}

// The LINK SPHERE PROPERTY, per vertex: an interior vertex of a
// combinatorial d-manifold has lk(v) ≅ S^{d-1}, a boundary vertex has
// lk(v) ≅ B^{d-1}; anything else (pinch fans, junction stars, isolated
// vertices) is a defect. Classification is by the link's effective
// dimension ℓ using closedness, connectedness, and χ:
//   ℓ = 0: sphere ⇔ two points (S⁰), ball ⇔ one point;
//   ℓ = 1: sphere ⇔ closed ∧ connected ∧ χ = 0 (a circle),
//          ball   ⇔ open ∧ connected ∧ χ = 1 (an interval);
//   ℓ = 2: sphere ⇔ closed ∧ connected ∧ χ = 2 (S², by surface
//          classification), ball ⇔ open ∧ connected ∧ χ = 1 (a disk).
// EXACT for complexes of dimension ≤ 3 (links are curves/surfaces, which
// these invariants classify completely). For ℓ ≥ 3 the test degrades to
// the Z₂-homology profile — necessary, not sufficient (sphere recognition
// is undecidable in high dimension) — and is honest about it.
//
// Collective. P=1 today.
struct LinkClassification {
  Marker sphere;  // interior manifold vertices
  Marker ball;    // boundary manifold vertices
  Marker other;   // defects
};

inline LinkClassification classify_vertex_links(const Complex& c) {
  LinkClassification out{Marker(c), Marker(c), Marker(c)};
  for (Index v = 0; v < c.count(0); ++v) {
    Marker mv(c);
    mv.mark(0, v);
    const Marker lk = link_of(c, mv);
    bool empty = true;
    for (int k = 0; k <= c.dim() && empty; ++k) empty = lk.marked_count(k) == 0;
    if (empty) {
      out.other.mark(0, v);
      continue;
    }
    const Complex L = subcomplex(c, lk).complex;
    int ell = L.dim();
    while (ell > 0 && L.count(ell) == 0) --ell;
    const Index comps = connected_components(L).count;
    const bool closed = is_closed(L);
    const long long chi = euler_characteristic(L);

    bool sphere = false, ball = false;
    if (ell == 0) {
      sphere = L.count(0) == 2;
      ball = L.count(0) == 1;
    } else if (ell == 1) {
      sphere = closed && comps == 1 && chi == 0;
      ball = !closed && comps == 1 && chi == 1;
    } else if (ell == 2) {
      sphere = closed && comps == 1 && chi == 2;
      ball = !closed && comps == 1 && chi == 1;
    } else {
      const auto b = betti_numbers_z2(L);
      bool profile_sphere = comps == 1 && closed;
      bool profile_ball = comps == 1 && !closed;
      for (int k = 1; k <= ell; ++k) {
        const Index expect_s = (k == ell) ? 1 : 0;
        profile_sphere = profile_sphere && b[static_cast<std::size_t>(k)] == expect_s;
        profile_ball = profile_ball && b[static_cast<std::size_t>(k)] == 0;
      }
      sphere = profile_sphere;  // homological only: see caveat above
      ball = profile_ball;
    }
    (sphere ? out.sphere : ball ? out.ball : out.other).mark(0, v);
  }
  return out;
}

}  // namespace graphos
