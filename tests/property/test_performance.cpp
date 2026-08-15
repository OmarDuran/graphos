// Constant-factor guards, which the scaling assertions cannot cover:
//
//  1. Storage: bytes per cell against the CSR model — int32 indices, int8
//     incidence numbers, both directions. Catches structural bloat exactly
//     and machine-independently.
//
//  2. Throughput: bulk operations are memory-bound, so cells/s should track
//     the machine's own memcpy bandwidth, measured in-process. Floors sit an
//     order of magnitude below the observed ratio, so only a constant-factor
//     regression of that size fails.

#include <chrono>
#include <cstring>
#include <vector>

#include "generators.hpp"
#include "graphos/graphos.hpp"
#include "graphos_test.hpp"

using Clock = std::chrono::steady_clock;
using graphos::Index;

// The floors are meaningful only in optimized, uninstrumented builds: a
// sanitizer slows the operations while memcpy stays native, breaking the
// calibration. Such builds report without asserting.
#if !defined(NDEBUG) || defined(__SANITIZE_ADDRESS__)
#define GRAPHOS_PERF_ASSERT 0
#elif defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(undefined_behavior_sanitizer) || \
    __has_feature(thread_sanitizer)
#define GRAPHOS_PERF_ASSERT 0
#else
#define GRAPHOS_PERF_ASSERT 1
#endif
#else
#define GRAPHOS_PERF_ASSERT 1
#endif

namespace {

template <typename F>
double best_of_3_ms(F&& f) {
  double best = 1e30;
  for (int r = 0; r < 3; ++r) {
    const auto t0 = Clock::now();
    f();
    const auto t1 = Clock::now();
    best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  return best;
}

volatile char g_sink = 0;  // defeats dead-code elimination of the memcpy

double memcpy_gbps() {
  const std::size_t n = std::size_t(64) << 20;  // 64 MB
  std::vector<char> a(n, 1), b(n, 0);
  const double t = best_of_3_ms([&] {
    std::memcpy(b.data(), a.data(), n);
    g_sink = b[n / 2];
  });
  return static_cast<double>(n) / (std::max(t, 0.01) * 1e6);  // GB/s
}

std::size_t stored_bytes(const graphos::Complex& c) {
  std::size_t bytes = 0;
  for (int k = 1; k <= c.dim(); ++k) {
    const graphos::BoundaryOperator& b = c.boundary(k);
    bytes += b.offsets.size() * sizeof(Index) + b.indices.size() * sizeof(Index) +
             b.signs.size() * sizeof(graphos::Sign);
  }
  return bytes;
}

std::size_t coboundary_bytes(const graphos::Complex& c) {
  std::size_t bytes = 0;
  for (int k = 0; k < c.dim(); ++k) {
    const graphos::CoboundaryOperator cob = graphos::coboundary(c, k);
    bytes += cob.offsets.size() * sizeof(Index) + cob.indices.size() * sizeof(Index) +
             cob.signs.size() * sizeof(graphos::Sign);
  }
  return bytes;
}

}  // namespace

GRAPHOS_TEST(storage_matches_the_csr_memory_model) {
  const int n = 16;
  const graphos::Complex c =
      graphos::from_simplices(3, graphos_bench::kuhn_vertex_count(n), graphos_bench::kuhn_tets(n));
  double cells = 0;
  for (int k = 0; k <= 3; ++k) cells += c.count(k);

  // model: nnz per cell ≈ 2.8 on a simplicial complex, five bytes an entry
  // plus offsets, per direction — ∂ alone ≈ 18 B/cell, ∂ with δ ≈ 36
  const double boundary_per_cell = static_cast<double>(stored_bytes(c)) / cells;
  const double both_per_cell = static_cast<double>(stored_bytes(c) + coboundary_bytes(c)) / cells;
  std::printf("  storage: boundary %.1f B/cell, with coboundary %.1f B/cell\n", boundary_per_cell,
              both_per_cell);
  CHECK(boundary_per_cell > 10.0);  // sanity: the model itself
  CHECK(boundary_per_cell < 30.0);  // ∂ storage has not bloated
  CHECK(both_per_cell < 60.0);      // frozen-complex-equivalent budget
}

GRAPHOS_TEST(bulk_op_throughput_tracks_memory_bandwidth) {
  const double gbps = memcpy_gbps();

  const int n = 32;
  const graphos::Complex c =
      graphos::from_simplices(3, graphos_bench::kuhn_vertex_count(n), graphos_bench::kuhn_tets(n));
  double cells = 0;
  for (int k = 0; k <= 3; ++k) cells += c.count(k);

  graphos::Marker one(c);
  one.mark(0, c.count(0) / 2);
  const double sd_ms = best_of_3_ms([&] { graphos::star_deletion(c, one); });
  const double fz_ms = best_of_3_ms([&] { graphos::freeze(c); });
  const auto cls = graphos::classify_facets(c);
  const double sub_ms = best_of_3_ms([&] { graphos::subcomplex(c, cls.boundary); });

  // observed: 2–3 Mcells/s per GB/s of memcpy; a floor of 0.2 leaves an
  // order of magnitude for slower memory and CI noise
  const double floor_cells_per_ms = gbps * 0.2e3;
  const auto report = [&](const char* name, double t_ms) {
    const double rate = cells / t_ms;
    std::printf("  %-14s %8.0f cells/ms  (floor %.0f, memcpy %.0f GB/s, %s)\n", name, rate,
                floor_cells_per_ms, gbps, GRAPHOS_PERF_ASSERT ? "asserted" : "report-only");
#if GRAPHOS_PERF_ASSERT
    CHECK(rate > floor_cells_per_ms);
#endif
  };
  report("star_deletion", sd_ms);
  report("freeze", fz_ms);
  report("subcomplex", sub_ms);
}

GRAPHOS_TEST_MAIN()
