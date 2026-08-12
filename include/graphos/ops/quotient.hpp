#pragma once

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"

namespace graphos {

struct QuotientResult {
  Complex complex;
  ChainMap map;
};

// The quotient complex under the given per-dimension identifications: each
// `from` cell is identified with its `to` representative, boundary chains are
// rewritten through the quotient map (orientation flips propagate into the
// star of every reversed cell), and surviving cells are compacted.
//
// Logically collective: every rank participates, supplying identifications
// for cells of its own partition; chains spanning ranks resolve through a
// distributed union-find. P=1 today: all cells are local.
//
// Precondition (caller's responsibility, typically via find_parallel_cells
// or an external geometric decision): identified cells have equal boundary
// chains up to rel_sign. The identified cell's own boundary row is
// discarded, not checked.
//
// Identification chains (a~b, b~c) resolve with sign composition; cycles are
// an error.
inline QuotientResult quotient(const Complex& c,
                               const std::vector<std::vector<Identification>>& identifications) {
  const int dim = c.dim();

  std::vector<std::vector<Index>> rep(static_cast<std::size_t>(dim) + 1);
  std::vector<std::vector<Sign>> rep_sign(static_cast<std::size_t>(dim) + 1);
  for (int k = 0; k <= dim; ++k) {
    const Index n = c.count(k);
    rep[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(n));
    rep_sign[static_cast<std::size_t>(k)].assign(static_cast<std::size_t>(n), Sign{1});
    for (Index i = 0; i < n; ++i) rep[static_cast<std::size_t>(k)][i] = i;
  }

  for (std::size_t k = 0; k < identifications.size() && k <= static_cast<std::size_t>(dim); ++k) {
    for (const Identification& id : identifications[k]) {
      if (id.from < 0 || id.from >= c.count(static_cast<int>(k)) || id.to < 0 ||
          id.to >= c.count(static_cast<int>(k))) {
        throw std::out_of_range("quotient: identification index out of range");
      }
      if (id.rel_sign != 1 && id.rel_sign != -1) {
        throw std::invalid_argument("quotient: rel_sign must be +/-1");
      }
      if (id.from == id.to) continue;
      rep[k][static_cast<std::size_t>(id.from)] = id.to;
      rep_sign[k][static_cast<std::size_t>(id.from)] = id.rel_sign;
    }
  }

  ChainMap map = ChainMap::sized(c.counts());
  std::vector<Index> kept_counts(static_cast<std::size_t>(dim) + 1, 0);
  for (int k = 0; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const Index n = c.count(k);
    std::vector<Index> compact(static_cast<std::size_t>(n), invalid_index);
    Index nn = 0;
    for (Index i = 0; i < n; ++i) {
      if (rep[sk][static_cast<std::size_t>(i)] == i) compact[static_cast<std::size_t>(i)] = nn++;
    }
    kept_counts[sk] = nn;
    for (Index i = 0; i < n; ++i) {
      Index r = i;
      int s = 1;
      Index steps = 0;
      while (rep[sk][static_cast<std::size_t>(r)] != r) {
        s *= rep_sign[sk][static_cast<std::size_t>(r)];
        r = rep[sk][static_cast<std::size_t>(r)];
        if (++steps > n) throw std::runtime_error("quotient: cyclic identification chain");
      }
      map.index[sk][static_cast<std::size_t>(i)] = compact[static_cast<std::size_t>(r)];
      map.sign[sk][static_cast<std::size_t>(i)] = static_cast<Sign>(s);
    }
  }

  Complex out(dim);
  out.attach_vertices(kept_counts[0]);
  std::vector<Index> row_idx;
  std::vector<Sign> row_sg;
  for (int k = 1; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const BoundaryOperator& bnd = c.boundary(k);
    for (Index e = 0; e < c.count(k); ++e) {
      if (rep[sk][static_cast<std::size_t>(e)] != e) continue;  // identified away
      row_idx.clear();
      row_sg.clear();
      for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
        const Index b = bnd.indices[m];
        row_idx.push_back(map.index[sk - 1][static_cast<std::size_t>(b)]);
        row_sg.push_back(
            static_cast<Sign>(bnd.signs[m] * map.sign[sk - 1][static_cast<std::size_t>(b)]));
      }
      out.attach_cell(k, std::span<const Index>(row_idx), std::span<const Sign>(row_sg));
    }
  }

  return QuotientResult{std::move(out), std::move(map)};
}

// Finds parallel k-cells: cells whose boundary chains coincide as signed
// sets, i.e. candidates for identification after lower skeleta have been
// glued. Two cells match when they reference the same faces and their
// orientation patterns agree up to a uniform flip (the reported rel_sign).
//
// This is deliberately a separate query rather than something quotient does
// implicitly: in a BRep two distinct curves may legitimately span the same
// two vertices, and only the caller (with geometric knowledge) can decide.
// Degenerate boundary chains (repeated faces) are skipped.
inline std::vector<Identification> find_parallel_cells(const Complex& c, int k) {
  std::vector<Identification> out;
  if (k < 1 || k > c.dim()) {
    throw std::invalid_argument("find_parallel_cells: dimension out of range");
  }
  const BoundaryOperator& bnd = c.boundary(k);
  std::map<std::vector<Index>, std::pair<Index, std::vector<Sign>>> seen;
  std::vector<std::pair<Index, Sign>> row;
  for (Index e = 0; e < c.count(k); ++e) {
    row.clear();
    for (Index m = bnd.offsets[e]; m < bnd.offsets[e + 1]; ++m) {
      row.emplace_back(bnd.indices[m], bnd.signs[m]);
    }
    if (row.empty()) continue;
    std::sort(row.begin(), row.end(),
              [](const auto& x, const auto& y) { return x.first < y.first; });
    bool degenerate = false;
    for (std::size_t m = 1; m < row.size(); ++m) {
      if (row[m].first == row[m - 1].first) degenerate = true;
    }
    if (degenerate) continue;

    std::vector<Index> key;
    std::vector<Sign> sg;
    key.reserve(row.size());
    sg.reserve(row.size());
    for (const auto& [i, s] : row) {
      key.push_back(i);
      sg.push_back(s);
    }
    auto it = seen.find(key);
    if (it == seen.end()) {
      seen.emplace(std::move(key), std::make_pair(e, std::move(sg)));
      continue;
    }
    const auto& [rep_e, rep_sg] = it->second;
    const int r = static_cast<int>(sg[0]) * static_cast<int>(rep_sg[0]);
    bool consistent = true;
    for (std::size_t m = 0; m < sg.size(); ++m) {
      if (static_cast<int>(sg[m]) != r * static_cast<int>(rep_sg[m])) {
        consistent = false;
        break;
      }
    }
    if (consistent) out.push_back({e, rep_e, static_cast<Sign>(r)});
  }
  return out;
}

}  // namespace graphos
