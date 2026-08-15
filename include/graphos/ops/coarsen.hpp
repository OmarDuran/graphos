#pragma once

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/core/union_find.hpp"

namespace graphos {

struct CoarsenResult {
  Complex complex;
  // fine → coarse: index[k][σ] is the coarse cell σ joined (invalid_index if
  // sent to 0), sign[k][σ] its orientation within that cell. Cochain
  // restriction gathers through this map.
  ChainMap map;
};

// Multilevel coarsening by dimensional descent of agglomerate. The top
// stratum is agglomerated by the given labels; descending one stratum at a
// time, the surviving interface cells are agglomerated in turn — cells with
// equal sets of coarse cofaces, connected through faces incident to exactly
// two of them, merge under the same cancellation rule, interior faces
// cancelling and the rim surviving. The descent stops at the vertices, which
// are only relabelled; merging vertices is quotient's operation.
//
// A patch merges only if every condition holds; otherwise its members remain
// individual coarse cells rather than raising:
//   - the patch admits a coherent orientation (sign propagation),
//   - no boundary coefficient exceeds 1 in magnitude,
//   - the merged boundary is nonempty (no closed cell is created),
//   - the patch sits consistently inside each of its coarse cofaces.
//
// Being metric-free, graphos cannot locate a geometric feature: an L-shaped
// coarse edge is combinatorially legitimate. Feature cells are therefore
// declared through `protected_cells`. Protection is a barrier — patches one
// stratum up neither connect nor cancel across a protected cell, which stays
// on the rim — while the protected cell may still merge in its own descent.
// Protecting the coarse lattice frame (corner vertices, then lattice edges,
// and so on by dimension) recovers a tensor hierarchy exactly. Unprotected,
// the merge is maximal.
inline CoarsenResult coarsen(const Complex& c, const std::vector<Index>& top_labels,
                             const Marker* protected_cells = nullptr) {
  const int n = c.dim();
  if (n < 1) throw std::invalid_argument("coarsen: complex must have dimension >= 1");
  if (static_cast<Index>(top_labels.size()) != c.count(n)) {
    throw std::invalid_argument("coarsen: one label per top cell required");
  }
  if (protected_cells != nullptr) protected_cells->validate_for(c);
  const auto guarded = [&](int k, Index i) {
    return protected_cells != nullptr && protected_cells->marked(k, i);
  };

  Index n_agg = 0;
  for (const Index a : top_labels) {
    if (a < 0) throw std::invalid_argument("coarsen: labels must be non-negative");
    n_agg = std::max(n_agg, a + 1);
  }
  std::vector<char> used(static_cast<std::size_t>(n_agg), 0);
  for (const Index a : top_labels) used[static_cast<std::size_t>(a)] = 1;
  for (const char u : used) {
    if (!u) throw std::invalid_argument("coarsen: aggregate ids must be contiguous");
  }

  // per stratum: the fine → coarse assignment, member orientations, and
  // coarse rows in the fine basis of the stratum below
  std::vector<std::vector<Index>> group(static_cast<std::size_t>(n) + 1);
  std::vector<std::vector<Sign>> rel(static_cast<std::size_t>(n) + 1);
  std::vector<std::vector<std::map<Index, int>>> rows(static_cast<std::size_t>(n) + 1);
  std::vector<Index> coarse_counts(static_cast<std::size_t>(n) + 1, 0);

  // ---- top stratum: agglomeration by the given labels ------------------
  {
    const BoundaryOperator& bnd = c.boundary(n);
    group[static_cast<std::size_t>(n)] = top_labels;
    rel[static_cast<std::size_t>(n)].assign(static_cast<std::size_t>(c.count(n)), Sign{1});
    rows[static_cast<std::size_t>(n)].resize(static_cast<std::size_t>(n_agg));
    for (Index e = 0; e < c.count(n); ++e) {
      std::map<Index, int>& row =
          rows[static_cast<std::size_t>(n)]
              [static_cast<std::size_t>(top_labels[static_cast<std::size_t>(e)])];
      for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
        row[bnd.indices[m]] += bnd.signs[m];
      }
    }
    for (auto& row : rows[static_cast<std::size_t>(n)]) {
      for (auto it = row.begin(); it != row.end();) {
        if (it->second < -1 || it->second > 1) {
          throw std::invalid_argument(
              "coarsen: facet coefficient exceeds 1 - orientation is inconsistent, "
              "run orient() first");
        }
        it = (it->second == 0) ? row.erase(it) : std::next(it);
      }
    }
    coarse_counts[static_cast<std::size_t>(n)] = n_agg;
  }

  // ---- dimensional descent ---------------------------------------------
  for (int k = n - 1; k >= 0; --k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const Index nk = c.count(k);

    // signature: a cell's coarse cofaces with their coefficients
    std::vector<std::vector<std::pair<Index, int>>> sig(static_cast<std::size_t>(nk));
    for (Index C = 0; C < coarse_counts[sk + 1]; ++C) {
      for (const auto& [f, a] : rows[sk + 1][static_cast<std::size_t>(C)]) {
        sig[static_cast<std::size_t>(f)].emplace_back(C, a);
      }
    }
    const CoboundaryOperator cob = coboundary(c, k);

    // candidate patches: equal coface sets, unprotected, dim ≥ 1
    struct Patch {
      std::vector<Index> members;
      std::vector<Sign> coeff;
      std::map<Index, int> chain;
    };
    std::vector<Patch> patches;
    std::vector<Index> patch_of(static_cast<std::size_t>(nk), invalid_index);

    if (k >= 1) {
      std::map<std::vector<Index>, std::vector<Index>> groups;
      for (Index f = 0; f < nk; ++f) {
        if (sig[static_cast<std::size_t>(f)].empty()) continue;
        std::vector<Index> key;
        for (const auto& [C, a] : sig[static_cast<std::size_t>(f)]) {
          (void)a;
          key.push_back(C);
        }
        groups[std::move(key)].push_back(f);
      }
      const BoundaryOperator& bnd = c.boundary(k);
      for (const auto& [key, members] : groups) {
        (void)key;
        if (members.size() < 2) continue;
        // unprotected faces incident to exactly two members connect and
        // cancel; the rest is rim
        std::map<Index, std::vector<int>> via;
        for (std::size_t p = 0; p < members.size(); ++p) {
          for (Index m = bnd.offsets[members[p]]; m < bnd.offsets[members[p] + 1]; ++m) {
            if (!guarded(k - 1, bnd.indices[m])) via[bnd.indices[m]].push_back(static_cast<int>(p));
          }
        }
        UnionFind uf(static_cast<Index>(members.size()));
        for (const auto& [f, lst] : via) {
          (void)f;
          if (lst.size() == 2) uf.unite(lst[0], lst[1]);
        }
        std::map<Index, std::vector<int>> comps;
        for (std::size_t p = 0; p < members.size(); ++p) {
          comps[uf.find(static_cast<Index>(p))].push_back(static_cast<int>(p));
        }
        for (const auto& [root, positions] : comps) {
          (void)root;
          if (positions.size() < 2) continue;

          // propagate orientation across two-member faces
          std::vector<int> local(static_cast<std::size_t>(members.size()), -1);
          for (std::size_t i = 0; i < positions.size(); ++i) {
            local[static_cast<std::size_t>(positions[i])] = static_cast<int>(i);
          }
          std::vector<Sign> coeff(positions.size(), 0);
          coeff[0] = 1;
          bool ok = true;
          std::vector<int> stack{positions[0]};
          std::vector<std::vector<std::pair<int, Sign>>> adj(positions.size());
          for (const auto& [f, lst] : via) {
            if (lst.size() != 2) continue;
            const int la = local[static_cast<std::size_t>(lst[0])];
            const int lb = local[static_cast<std::size_t>(lst[1])];
            if (la < 0 || lb < 0) continue;
            const auto sign_in = [&](Index mem) {
              for (Index m = bnd.offsets[mem]; m < bnd.offsets[mem + 1]; ++m) {
                if (bnd.indices[m] == f) return bnd.signs[m];
              }
              return Sign{0};
            };
            const Sign sa = sign_in(members[static_cast<std::size_t>(lst[0])]);
            const Sign sb = sign_in(members[static_cast<std::size_t>(lst[1])]);
            const Sign r = static_cast<Sign>(-sa * sb);
            adj[static_cast<std::size_t>(la)].push_back({lb, r});
            adj[static_cast<std::size_t>(lb)].push_back({la, r});
          }
          std::vector<int> todo{0};
          while (!todo.empty() && ok) {
            const int cur = todo.back();
            todo.pop_back();
            for (const auto& [nb, r] : adj[static_cast<std::size_t>(cur)]) {
              const Sign want = static_cast<Sign>(coeff[static_cast<std::size_t>(cur)] * r);
              if (coeff[static_cast<std::size_t>(nb)] == 0) {
                coeff[static_cast<std::size_t>(nb)] = want;
                todo.push_back(nb);
              } else if (coeff[static_cast<std::size_t>(nb)] != want) {
                ok = false;
                break;
              }
            }
          }
          for (const Sign s : coeff) ok = ok && (s != 0);

          // cancellation: a nonempty merged boundary with coefficients in
          // {−1, 0, 1}
          std::map<Index, int> chain;
          if (ok) {
            for (std::size_t i = 0; i < positions.size() && ok; ++i) {
              const Index mem = members[static_cast<std::size_t>(positions[i])];
              for (Index m = bnd.offsets[mem]; m < bnd.offsets[mem + 1]; ++m) {
                chain[bnd.indices[m]] += coeff[i] * bnd.signs[m];
              }
            }
            for (auto it = chain.begin(); it != chain.end() && ok;) {
              if (it->second < -1 || it->second > 1) ok = false;
              it = (it->second == 0) ? chain.erase(it) : std::next(it);
            }
            ok = ok && !chain.empty();
          }

          // the patch must sit consistently inside every coarse coface
          if (ok) {
            const Index anchor = members[static_cast<std::size_t>(positions[0])];
            for (std::size_t i = 1; i < positions.size() && ok; ++i) {
              const Index mem = members[static_cast<std::size_t>(positions[i])];
              for (std::size_t s = 0; s < sig[static_cast<std::size_t>(anchor)].size() && ok; ++s) {
                const int a0 = sig[static_cast<std::size_t>(anchor)][s].second;
                const int am = sig[static_cast<std::size_t>(mem)][s].second;
                ok = (am * coeff[i] == a0 * coeff[0]);
              }
            }
          }

          if (ok) {
            Patch p;
            for (std::size_t i = 0; i < positions.size(); ++i) {
              p.members.push_back(members[static_cast<std::size_t>(positions[i])]);
              p.coeff.push_back(coeff[i]);
            }
            p.chain = std::move(chain);
            for (const Index mem : p.members) {
              patch_of[static_cast<std::size_t>(mem)] = static_cast<Index>(patches.size());
            }
            patches.push_back(std::move(p));
          }
        }
      }
    }

    // allocate in fine-index order, a patch at its first member; survivors
    // are the cells with a signature or with no cofaces
    group[sk].assign(static_cast<std::size_t>(nk), invalid_index);
    rel[sk].assign(static_cast<std::size_t>(nk), Sign{1});
    std::vector<char> patch_done(patches.size(), 0);
    for (Index f = 0; f < nk; ++f) {
      const std::size_t sf = static_cast<std::size_t>(f);
      const bool maximal = cob.offsets[sf + 1] == cob.offsets[sf];
      if (sig[sf].empty() && !maximal) continue;  // interior: sent to zero
      const Index pid = patch_of[sf];
      if (pid != invalid_index) {
        if (!patch_done[static_cast<std::size_t>(pid)]) {
          patch_done[static_cast<std::size_t>(pid)] = 1;
          const Patch& p = patches[static_cast<std::size_t>(pid)];
          const Index id = coarse_counts[sk]++;
          for (std::size_t i = 0; i < p.members.size(); ++i) {
            group[sk][static_cast<std::size_t>(p.members[i])] = id;
            rel[sk][static_cast<std::size_t>(p.members[i])] = p.coeff[i];
          }
          if (k >= 1) rows[sk].push_back(p.chain);
        }
        continue;
      }
      const Index id = coarse_counts[sk]++;
      group[sk][sf] = id;
      if (k >= 1) {
        const BoundaryOperator& bnd = c.boundary(k);
        std::map<Index, int> row;
        for (Index m = bnd.offsets[f]; m < bnd.offsets[f + 1]; ++m) {
          row[bnd.indices[m]] += bnd.signs[m];
        }
        for (auto it = row.begin(); it != row.end();) {
          it = (it->second == 0) ? row.erase(it) : std::next(it);
        }
        rows[sk].push_back(std::move(row));
      }
    }
  }

  // ---- assemble in the coarse basis ------------------------------------
  std::vector<BoundaryOperator> strata(static_cast<std::size_t>(n) + 1);
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (int k = 1; k <= n; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    for (Index K = 0; K < coarse_counts[sk]; ++K) {
      std::map<Index, int> crow;
      for (const auto& [f, a] : rows[sk][static_cast<std::size_t>(K)]) {
        const Index Kf = group[sk - 1][static_cast<std::size_t>(f)];
        if (Kf == invalid_index) throw std::logic_error("coarsen: zero-mapped face referenced");
        const int b = a * rel[sk - 1][static_cast<std::size_t>(f)];
        auto [it, fresh] = crow.try_emplace(Kf, b);
        if (!fresh && it->second != b) throw std::logic_error("coarsen: inconsistent coarse row");
      }
      row_idx.clear();
      row_sg.clear();
      for (const auto& [Kf, b] : crow) {
        row_idx.push_back(Kf);
        row_sg.push_back(static_cast<Sign>(b));
      }
      strata[sk].append_row(row_idx, row_sg);
    }
  }

  CoarsenResult res{Complex(coarse_counts, std::move(strata)), ChainMap::sized(c.counts())};
  for (int k = 0; k <= n; ++k) {
    res.map.index[static_cast<std::size_t>(k)] = group[static_cast<std::size_t>(k)];
    res.map.sign[static_cast<std::size_t>(k)] = rel[static_cast<std::size_t>(k)];
  }
  return res;
}

inline CoarsenResult coarsen(const Complex& c, const std::vector<Index>& top_labels,
                             const Marker& protected_cells) {
  return coarsen(c, top_labels, &protected_cells);
}

}  // namespace graphos
