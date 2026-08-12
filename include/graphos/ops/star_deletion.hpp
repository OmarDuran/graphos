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

// Deletes the closed stars of the marked cells: each marked cell is removed
// together with every cell it is a face of, cascading up through the
// skeleta. The result is the largest subcomplex of c containing none of the
// marked cells. Survivors are compacted; the chain map records each cell's
// fate (invalid_index for cells sent to zero).
//
// Collective. Every rank calls this with a Marker marked over its local
// partition; the deletion is the union of all ranks' marks, and the cascade
// propagates across partition boundaries (a fixed number of halo exchanges,
// bounded by the dimension). P=1 today: the local partition is everything.
//
// Faces left without any coface are kept: in a mixed-dimensional complex a
// bare vertex or edge is a legitimate maximal cell, and graphos cannot
// distinguish an orphan from one. Pruning is a separate caller decision.
//
// This op is the exemplar of the kernel form all bulk ops converge to:
// mark -> cascade -> scan -> scatter phases through the exec seams, every
// phase data-parallel over one stratum.
inline StarDeletionResult star_deletion(const Complex& c, const Marker& cells) {
  const int dim = c.dim();
  cells.validate_for(c);

  // mark: deleted[k][i] = 1 where the marker selects
  std::vector<exec::Buffer<Index>> deleted;
  deleted.reserve(static_cast<std::size_t>(dim) + 1);
  for (int k = 0; k <= dim; ++k) {
    deleted.emplace_back(static_cast<std::size_t>(c.count(k)));
    Index* g = deleted[static_cast<std::size_t>(k)].data();
    const char* mf = cells.flags(k).data();
    exec::forall(c.count(k), [=](Index i) { g[i] = mf[i] ? 1 : 0; });
  }

  // cascade: one parallel pass per stratum suffices because ∂ reaches
  // exactly one level down
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

  // scan: compact surviving cells into contiguous indices
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

  // scatter: assemble each surviving stratum's boundary operator
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
