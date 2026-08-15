#pragma once

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"

namespace graphos {

// Set-level neighbourhoods, the bulk counterparts of the per-cell st/cl/lk on
// FrozenComplex, returned as Markers so they compose with subcomplex(),
// star_deletion() and the relative Betti numbers. Together they make the
// excision decomposition computable:
//   st(S)       = star_of(c, S)                    the open star
//   cl(st(S))   = closure_of(c, star_of(c, S))     the regular neighbourhood
//   lk(S)       = link_of(c, S)                    cl(st S) ∖ st(cl S)
//   frontier(S) = frontier_of(c, S)                cl(st S) ∖ st(S)
// The frontier bounds the void left by star_deletion(c, S) and is what a patch
// is glued to in the pushout. It equals lk(S) for a set of vertices but not in
// general: the link of a top cell is empty, its frontier its boundary sphere.

// st(S): the marked cells with all their cofaces.
inline Marker star_of(const Complex& c, const Marker& cells) {
  cells.validate_for(c);
  Marker out(c);
  const int dim = c.dim();
  for (int k = 0; k <= dim; ++k) {
    const std::vector<char>& f = cells.flags(k);
    for (Index i = 0; i < c.count(k); ++i) {
      if (f[static_cast<std::size_t>(i)]) out.mark(k, i);
    }
  }
  for (int k = 0; k < dim; ++k) {
    const CoboundaryOperator cob = coboundary(c, k);
    for (Index i = 0; i < c.count(k); ++i) {
      if (!out.marked(k, i)) continue;
      for (Index m = cob.offsets[static_cast<std::size_t>(i)];
           m < cob.offsets[static_cast<std::size_t>(i) + 1]; ++m) {
        out.mark(k + 1, cob.indices[static_cast<std::size_t>(m)]);
      }
    }
  }
  return out;
}

// cl(S): the marked cells with all their faces.
inline Marker closure_of(const Complex& c, const Marker& cells) {
  cells.validate_for(c);
  Marker out(c);
  const int dim = c.dim();
  for (int k = 0; k <= dim; ++k) {
    const std::vector<char>& f = cells.flags(k);
    for (Index i = 0; i < c.count(k); ++i) {
      if (f[static_cast<std::size_t>(i)]) out.mark(k, i);
    }
  }
  for (int k = dim; k >= 1; --k) {
    const BoundaryOperator& bnd = c.boundary(k);
    for (Index e = 0; e < c.count(k); ++e) {
      if (!out.marked(k, e)) continue;
      for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
        out.mark(k - 1, bnd.indices[m]);
      }
    }
  }
  return out;
}

// lk(S) = cl(st(S)) ∖ st(cl(S)).
inline Marker link_of(const Complex& c, const Marker& cells) {
  const Marker cl_st = closure_of(c, star_of(c, cells));
  const Marker st_cl = star_of(c, closure_of(c, cells));
  Marker out(c);
  for (int k = 0; k <= c.dim(); ++k) {
    for (Index i = 0; i < c.count(k); ++i) {
      if (cl_st.marked(k, i) && !st_cl.marked(k, i)) out.mark(k, i);
    }
  }
  return out;
}

// cl(st(S)) ∖ st(S): the cells surviving star_deletion(c, S) that meet the
// removed region.
inline Marker frontier_of(const Complex& c, const Marker& cells) {
  const Marker st = star_of(c, cells);
  const Marker cl_st = closure_of(c, st);
  Marker out(c);
  for (int k = 0; k <= c.dim(); ++k) {
    for (Index i = 0; i < c.count(k); ++i) {
      if (cl_st.marked(k, i) && !st.marked(k, i)) out.mark(k, i);
    }
  }
  return out;
}

}  // namespace graphos
