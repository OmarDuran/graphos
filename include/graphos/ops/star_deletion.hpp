#pragma once

#include <stdexcept>
#include <utility>
#include <vector>

#include "graphos/core/complex.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/exec/forall.hpp"
#include "graphos/exec/memory.hpp"

namespace graphos {

struct StarDeletionResult {
  Complex complex;
  ChainMap map;
};

// Removes st(S) for the marked S: each marked cell with all its cofaces,
// cascading up the strata. The result is the largest subcomplex of c meeting
// no marked cell. Survivors are compacted; the chain map sends deleted cells
// to 0.
//
// The union of all ranks' marks is deleted, and the cascade crosses partition
// boundaries in a number of halo exchanges bounded by dim.
//
// A cell left without cofaces is kept: in a mixed-dimensional complex a bare
// vertex or edge is a legitimate maximal cell, indistinguishable from an
// orphan. Pruning is the caller's decision.
//
// The kernel form every bulk operation takes: mark → cascade → scan →
// scatter, each phase data-parallel over one stratum.
inline StarDeletionResult star_deletion(const Complex& c, const Marker& cells) {
  const int dim = c.dim();
  cells.validate_for(c);

  // mark: deleted[k][i] = 1 where S selects
  std::vector<exec::Buffer<Index>> deleted;
  deleted.reserve(static_cast<std::size_t>(dim) + 1);
  for (int k = 0; k <= dim; ++k) {
    deleted.emplace_back(static_cast<std::size_t>(c.count(k)));
    Index* g = deleted[static_cast<std::size_t>(k)].data();
    const char* mf = cells.flags(k).data();
    exec::forall(c.count(k), [=](Index i) { g[i] = mf[i] ? 1 : 0; });
  }

  // cascade: one pass per stratum suffices, since ∂ reaches exactly one
  // stratum down
  for (int k = 1; k <= dim; ++k) {
    const BoundaryOperator& bnd = c.boundary(k);
    const Index* offsets = bnd.offsets.data();
    const Index* faces = bnd.indices.data();
    const Index* dlo = deleted[static_cast<std::size_t>(k) - 1].data();
    Index* dhi = deleted[static_cast<std::size_t>(k)].data();
    exec::forall(c.count(k), [=](Index e) {
      if (dhi[e]) return;
      for (Index m = offsets[e]; m < offsets[e + 1]; ++m) {
        if (dlo[faces[m]]) {
          dhi[e] = 1;
          return;
        }
      }
    });
  }

  // scan: compact survivors into contiguous indices
  ChainMap map = ChainMap::sized(c.counts());
  std::vector<Index> kept_counts(static_cast<std::size_t>(dim) + 1, 0);
  for (int k = 0; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const Index n = c.count(k);
    exec::Buffer<Index> keep(static_cast<std::size_t>(n));
    exec::Buffer<Index> pos(static_cast<std::size_t>(n));
    const Index* g = deleted[sk].data();
    Index* kp = keep.data();
    exec::forall(n, [=](Index i) { kp[i] = g[i] ? 0 : 1; });
    kept_counts[sk] = exec::exclusive_scan(keep.data(), pos.data(), n);
    Index* mi = map.index[sk].data();
    const Index* pp = pos.data();
    exec::forall(n, [=](Index i) { mi[i] = g[i] ? invalid_index : pp[i]; });
  }

  // scatter: assemble ∂_k on the survivors
  std::vector<BoundaryOperator> strata(static_cast<std::size_t>(dim) + 1);
  for (int k = 1; k <= dim; ++k) {
    const std::size_t sk = static_cast<std::size_t>(k);
    const BoundaryOperator& bnd = c.boundary(k);
    const Index n = c.count(k);
    const Index nn = kept_counts[sk];
    BoundaryOperator& out = strata[sk];
    out.offsets.assign(static_cast<std::size_t>(nn) + 1, 0);

    const Index* g = deleted[sk].data();
    const Index* mi = map.index[sk].data();
    const Index* offs = bnd.offsets.data();
    Index* no = out.offsets.data();
    exec::forall(n, [=](Index e) {
      if (!g[e]) no[mi[e] + 1] = offs[e + 1] - offs[e];
    });
    exec::inclusive_scan_inplace(no, nn + 1);

    out.indices.resize(static_cast<std::size_t>(no[nn]));
    out.signs.resize(static_cast<std::size_t>(no[nn]));
    const Index* faces = bnd.indices.data();
    const Sign* sg = bnd.signs.data();
    const Index* mlo = map.index[sk - 1].data();
    Index* oidx = out.indices.data();
    Sign* osg = out.signs.data();
    exec::forall(n, [=](Index e) {
      if (g[e]) return;
      Index w = no[mi[e]];
      for (Index m = offs[e]; m < offs[e + 1]; ++m) {
        oidx[w] = mlo[faces[m]];
        osg[w] = sg[m];
        ++w;
      }
    });
  }

  return StarDeletionResult{Complex(std::move(kept_counts), std::move(strata)), std::move(map)};
}

}  // namespace graphos
