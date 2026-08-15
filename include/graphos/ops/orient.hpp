#pragma once

#include <queue>
#include <utility>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"

namespace graphos {

struct OrientationResult {
  // ∂_k rows negated where a flip was needed, in the stratum asked for.
  Complex complex;
  // Identity on indices; sign[dim][i] is each cell's flip, which is how
  // cochains on flipped cells transport.
  ChainMap map;
  // False when the stratum is non-orientable. The complex then carries the
  // spanning-tree best effort.
  bool orientable{true};
  // The orientation classes: class_of[i] is the adjacency component cell i
  // was reached in, and `classes` counts them.
  //
  // Coherence fixes each cell relative to its neighbours and leaves the class
  // as a whole free to be inside-out. No boundary algebra can decide which
  // way is out — "out" is not a property of a chain complex — so the free
  // signs number exactly one per class. A caller with a geometric notion of
  // outward settles each class with one determinant and flip_cells.
  std::vector<Index> class_of;
  Index classes{0};
  // the stratum this result concerns
  int stratum{0};
};

// Reverses the orientation of cells in stratum k: their ∂_k rows change sign,
// which is what reversal means to every operator over the complex. The
// operation `orient` applies internally, exposed so a caller can apply its
// own.
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

  // and the column above, or ∂∘∂ = 0 fails. Reversing σ negates its row of
  // ∂_k; the term [τ : σ][σ : ρ] of ∂_k ∂_{k+1} then cancels only if [τ : σ]
  // is negated with it. Both factors flip and their product does not — the
  // content of "σ carries an orientation" rather than "σ's row carries a
  // sign". The top stratum has no column above and this loop is empty.
  if (k < n) {
    BoundaryOperator& up = strata[static_cast<std::size_t>(k) + 1];
    for (std::size_t m = 0; m < up.indices.size(); ++m) {
      if (turned[static_cast<std::size_t>(up.indices[m])])
        up.signs[m] = static_cast<Sign>(-up.signs[m]);
    }
  }
  return Complex(c.counts(), std::move(strata));
}

// Coherently orients stratum k: flips propagate through interior facets —
// those with exactly two cofaces within the stratum — until each is induced
// with opposite signs by its two sides. That is the condition for Σ σ to be a
// fundamental chain, and the convention mixed methods, flux continuity and
// DEC assume. Non-orientability is detected, not silently mis-signed.
//
// Any stratum, not only the top. A fracture network is a 2-stratum of a
// 3-complex, and its flux continuity is the same statement about ∂ one
// dimension down. The default remains the top stratum.
//
// Facets with one coface (∂K), none (detached interfaces), or three or more
// (non-manifold junctions) impose no constraint and do not propagate.
//
// Deterministic: each adjacency component is seeded at its lowest-indexed
// cell with flip +1, and reported in `class_of` so a caller can settle the
// one sign per component that coherence leaves free.
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

  // The maximal cells of the stratum, and only those. Coherence concerns a
  // stratum that is a manifold in its own right. A k-cell that is a facet of
  // something is a different object: its orientation is the shared reference
  // its two cofaces agree on, not a free choice, and reversing it would only
  // relabel which calls it positive. So the facets of a 3-complex are left
  // alone, which is also why nothing here can disturb ∂∘∂ = 0.
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

  // Cofaces within the stratum. On the top stratum every coface of a facet is
  // a k-cell, which is what δ returns; lower down it is not, so the count
  // deciding "interior" must be taken over the k-cells alone. A fracture edge
  // shared by two fracture faces is interior to the fracture whatever the
  // surrounding volume does.
  std::vector<std::vector<Index>> upper(static_cast<std::size_t>(c.count(k - 1)));
  for (Index e = 0; e < c.count(k); ++e) {
    if (!maximal[static_cast<std::size_t>(e)]) continue;
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      upper[static_cast<std::size_t>(bnd.indices[m])].push_back(e);
    }
  }

  // adjacency through interior facets: (neighbour, own sign, neighbour sign)
  struct Arc {
    Index nb;
    Sign s_own;
    Sign s_nb;
  };
  std::vector<std::vector<Arc>> adj(static_cast<std::size_t>(c.count(k)));
  for (Index f = 0; f < c.count(k - 1); ++f) {
    const std::vector<Index>& up = upper[static_cast<std::size_t>(f)];
    if (up.size() != 2) continue;
    if (!maximal[static_cast<std::size_t>(up[0])] || !maximal[static_cast<std::size_t>(up[1])])
      continue;
    const Index e1 = up[0];
    const Index e2 = up[1];
    if (e1 == e2) {
      // the facet is glued to itself within one cell: coherent only if the
      // two occurrences carry opposite incidence numbers
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
        // opposite induced signs: flip_cur · s_own = −(flip_nb · s_nb)
        const Sign required =
            static_cast<Sign>(-flip[static_cast<std::size_t>(cur)] * arc.s_own * arc.s_nb);
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

  // rebuild with flipped rows, recording them in the chain map
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

// the top stratum, the default
inline OrientationResult orient(const Complex& c) { return orient(c, c.dim()); }

}  // namespace graphos
