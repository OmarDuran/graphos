#pragma once

#include <queue>
#include <utility>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"

namespace graphos {

struct OrientationResult {
  // Cell boundary rows negated where a flip was needed, in the stratum asked for.
  Complex complex;
  // Identity on indices; sign[dim][i] records each cell's flip, which
  // is exactly how cochains attached to flipped cells transport.
  ChainMap map;
  // False when no consistent orientation exists (e.g. a Möbius band). The
  // returned complex then carries the spanning-tree best effort.
  bool orientable{true};
  // THE ORIENTATION CLASSES of the stratum: class_of[i] is the index of the
  // component cell i was reached in, and classes counts them.
  //
  // This is the half of orientation that IS topological, reported rather than
  // discarded. Coherence fixes every cell of a class RELATIVE to its
  // neighbours and leaves the class as a whole free to be inside-out; no
  // amount of boundary algebra can tell which way is out, because "out" is not
  // a property of a chain complex. So the free signs are exactly one per
  // class, and a caller that does have a notion of out -- a geometric one --
  // needs to decide that many, not one per cell. Handing back the classes is
  // what lets it: it settles each class with a single determinant and flips
  // the class with flip_cells, and no per-cell patch ever has to exist.
  std::vector<Index> class_of;
  Index classes{0};
  // the stratum this result is about
  int stratum{0};
};

// NEGATE THE ORIENTATION of a set of cells in stratum k: their boundary rows
// change sign, which is precisely what reversing a cell's orientation means to
// every operator built on the complex. The inverse of the flips `orient`
// applies, exposed so a caller can apply its own.
inline Complex flip_cells(const Complex& c, int k, const std::vector<Index>& cells) {
  const int n = c.dim();
  if (k < 1 || k > n) return c;
  std::vector<BoundaryOperator> strata(static_cast<std::size_t>(n) + 1);
  for (int j = 1; j <= n; ++j) strata[static_cast<std::size_t>(j)] = c.boundary(j);

  std::vector<bool> turned(static_cast<std::size_t>(c.count(k)), false);
  BoundaryOperator& b = strata[static_cast<std::size_t>(k)];
  for (const Index e : cells) {
    if (turned[static_cast<std::size_t>(e)]) continue;  // an involution, not a toggle
    turned[static_cast<std::size_t>(e)] = true;
    for (Index m = b.offsets[e]; m < b.offsets[e + 1]; ++m) {
      b.signs[m] = static_cast<Sign>(-b.signs[m]);
    }
  }

  // AND THE COLUMN ABOVE, or the complex stops being one. Reversing sigma
  // negates its row of d_k; the term [tau : sigma][sigma : rho] of d_k d_{k+1}
  // then survives only if [tau : sigma] is negated with it. Both factors flip,
  // their product does not, and d.d = 0 is preserved -- which is the whole
  // content of "sigma has an orientation" as opposed to "sigma's row has a
  // sign". For the top stratum there is no column above and this loop is empty.
  if (k < n) {
    BoundaryOperator& up = strata[static_cast<std::size_t>(k) + 1];
    for (std::size_t m = 0; m < up.indices.size(); ++m) {
      if (turned[static_cast<std::size_t>(up.indices[m])]) up.signs[m] = static_cast<Sign>(-up.signs[m]);
    }
  }
  return Complex(c.counts(), std::move(strata));
}

// Constructs a consistent global orientation of the cells of stratum `k`:
// flips are propagated through interior facets (exactly two cofaces within the
// stratum) so that every such facet is induced with opposite signs by its two
// sides — the condition for the sum of the cells to be a fundamental chain,
// and the convention mixed methods, flux continuity, and DEC silently assume.
// Meshes imported with arbitrary per-cell orientations are normalized by one
// call; non-orientability is detected and reported rather than silently
// mis-signed.
//
// ANY STRATUM, not only the top one. A fracture network is a 2-stratum inside
// a three-dimensional complex and its flux continuity is the same statement
// about the same boundary operator one dimension down; a k-form on it needs
// its cells coherently wound just as much. Orienting only the top cells left
// that case to be handled by whoever remembered to, which is to say by nobody.
// The default is still the top stratum, so existing callers are unchanged.
//
// Facets with one coface (boundary), none (detached interfaces), or three
// or more (nonmanifold junctions) impose no constraint and do not
// propagate.
//
// Deterministic: each adjacency component is seeded at its lowest-indexed
// cell with flip +1, and the components are reported in `class_of` so a caller
// with a geometric notion of "outward" can settle the one sign per component
// that coherence necessarily leaves free.
//
// Collective. P=1 today.
inline OrientationResult orient(const Complex& c, int k) {
  const int n = c.dim();
  OrientationResult res{c, ChainMap::sized(c.counts()), true, {}, 0, k};
  for (int j = 0; j <= n; ++j) {
    for (Index i = 0; i < c.count(j); ++i) {
      res.map.index[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = i;
    }
  }
  if (k < 1 || k > n) return res;

  const BoundaryOperator& bnd = c.boundary(k);

  // THE MAXIMAL CELLS OF THE STRATUM, and only those. Coherence is a statement
  // about a stratum that is a manifold in its own right: the top cells of a
  // volume mesh, or a fracture network sitting inside one. A k-cell that is a
  // FACET of something is a different object -- its orientation is the shared
  // reference the cells on either side agree on, not a free choice, and turning
  // it would only relabel which of them calls it positive. So the faces of a
  // three-dimensional mesh are left exactly as they are, which is also why
  // nothing here can disturb d.d = 0.
  std::vector<bool> maximal(static_cast<std::size_t>(c.count(k)), true);
  if (k < n) {
    const BoundaryOperator& up = c.boundary(k + 1);
    for (const Index i : up.indices) maximal[static_cast<std::size_t>(i)] = false;
  }

  const auto signs_of = [&bnd](Index e, Index f) {
    std::vector<Sign> out;
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      if (bnd.indices[m] == f) out.push_back(bnd.signs[m]);
    }
    return out;
  };

  // COFACES WITHIN THE STRATUM. On the top stratum every coface of a facet is
  // a k-cell and this is what coboundary already returns; on a lower one it is
  // not, so the count that decides "interior" has to be taken over the k-cells
  // alone. A fracture edge shared by two fracture faces is interior to the
  // fracture whatever the surrounding volume does.
  std::vector<std::vector<Index>> upper(static_cast<std::size_t>(c.count(k - 1)));
  for (Index e = 0; e < c.count(k); ++e) {
    if (!maximal[static_cast<std::size_t>(e)]) continue;
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      upper[static_cast<std::size_t>(bnd.indices[m])].push_back(e);
    }
  }

  // adjacency through interior facets: (neighbor, own sign, neighbor sign)
  struct Arc {
    Index nb;
    Sign s_own;
    Sign s_nb;
  };
  std::vector<std::vector<Arc>> adj(static_cast<std::size_t>(c.count(k)));
  for (Index f = 0; f < c.count(k - 1); ++f) {
    const std::vector<Index>& up = upper[static_cast<std::size_t>(f)];
    if (up.size() != 2) continue;
    if (!maximal[static_cast<std::size_t>(up[0])] || !maximal[static_cast<std::size_t>(up[1])]) continue;
    const Index e1 = up[0];
    const Index e2 = up[1];
    if (e1 == e2) {
      // the facet is glued to itself within one cell: consistent only if
      // the two occurrences carry opposite signs
      const std::vector<Sign> ss = signs_of(e1, f);
      if (ss.size() == 2 && ss[0] == ss[1]) res.orientable = false;
      continue;
    }
    const Sign s1 = signs_of(e1, f).front();
    const Sign s2 = signs_of(e2, f).front();
    adj[static_cast<std::size_t>(e1)].push_back({e2, s1, s2});
    adj[static_cast<std::size_t>(e2)].push_back({e1, s2, s1});
  }

  std::vector<Sign> flip(static_cast<std::size_t>(c.count(k)), 0);
  res.class_of.assign(static_cast<std::size_t>(c.count(k)), -1);
  for (Index seed = 0; seed < c.count(k); ++seed) {
    if (flip[static_cast<std::size_t>(seed)] != 0) continue;
    if (!maximal[static_cast<std::size_t>(seed)]) continue;
    const Index cls = res.classes++;
    flip[static_cast<std::size_t>(seed)] = 1;
    res.class_of[static_cast<std::size_t>(seed)] = cls;
    std::queue<Index> queue;
    queue.push(seed);
    while (!queue.empty()) {
      const Index cur = queue.front();
      queue.pop();
      for (const Arc& arc : adj[static_cast<std::size_t>(cur)]) {
        // opposite induced signs: flip_cur * s_own = -(flip_nb * s_nb)
        const Sign required = static_cast<Sign>(-flip[static_cast<std::size_t>(cur)] *
                                                arc.s_own * arc.s_nb);
        Sign& nb = flip[static_cast<std::size_t>(arc.nb)];
        if (nb == 0) {
          nb = required;
          res.class_of[static_cast<std::size_t>(arc.nb)] = cls;
          queue.push(arc.nb);
        } else if (nb != required) {
          res.orientable = false;
        }
      }
    }
  }

  // rebuild with flipped rows and record the flips in the chain map
  std::vector<BoundaryOperator> strata(static_cast<std::size_t>(n) + 1);
  for (int j = 1; j <= n; ++j) strata[static_cast<std::size_t>(j)] = c.boundary(j);
  BoundaryOperator& top = strata[static_cast<std::size_t>(k)];
  for (Index e = 0; e < c.count(k); ++e) {
    res.map.sign[static_cast<std::size_t>(k)][static_cast<std::size_t>(e)] =
        flip[static_cast<std::size_t>(e)];
    if (flip[static_cast<std::size_t>(e)] == -1) {
      for (Index m = top.offsets[e]; m < top.offsets[e + 1]; ++m) {
        top.signs[m] = static_cast<Sign>(-top.signs[m]);
      }
    }
  }
  res.complex = Complex(c.counts(), std::move(strata));
  return res;
}

// the top stratum, which is what a caller that does not say means
inline OrientationResult orient(const Complex& c) { return orient(c, c.dim()); }

}  // namespace graphos
