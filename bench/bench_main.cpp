// Benchmark driver: structured simplicial grids, per-operation wall time and
// throughput. Usage: graphos_bench [n ...], grid sizes, default 16 32. These
// numbers are the evidence behind any performance claim.

#include <chrono>
#include <cstdio>
#include <cstdlib>

#include "generators.hpp"
#include "graphos/graphos.hpp"

using Clock = std::chrono::steady_clock;

namespace {

double ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

void run(int n) {
  using namespace graphos;
  const Index nv = graphos_bench::kuhn_vertex_count(n);
  const auto tets = graphos_bench::kuhn_tets(n);
  const double ntets = static_cast<double>(tets.size());
  std::printf("n=%d: %zu tets, %d vertices\n", n, tets.size(), nv);

  auto t0 = Clock::now();
  Complex c = from_simplices(3, nv, tets);
  auto t1 = Clock::now();
  const Index total_cells = c.count(0) + c.count(1) + c.count(2) + c.count(3);
  std::printf("  %-22s %10.1f ms  %8.2f Mtets/s\n", "from_simplices", ms(t0, t1),
              ntets / ms(t0, t1) / 1e3);

  t0 = Clock::now();
  FrozenComplex f = freeze(c);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms  %8.2f Mcells/s\n", "freeze", ms(t0, t1),
              total_cells / ms(t0, t1) / 1e3);

  t0 = Clock::now();
  const Adjacency i30 = incidence(c, 3, 0);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms  %8.2f Mtets/s\n", "incidence(3,0)", ms(t0, t1),
              ntets / ms(t0, t1) / 1e3);

  t0 = Clock::now();
  const FacetClassification cls = classify_facets(c);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms\n", "classify_facets", ms(t0, t1));

  t0 = Clock::now();
  const auto bdry = subcomplex(c, cls.boundary);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms  %8.2f Mcells/s\n", "subcomplex(boundary)", ms(t0, t1),
              total_cells / ms(t0, t1) / 1e3);

  Marker one(c);
  one.mark(0, nv / 2);
  t0 = Clock::now();
  const auto sd = star_deletion(c, one);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms  %8.2f Mcells/s\n", "star_deletion(vertex)", ms(t0, t1),
              total_cells / ms(t0, t1) / 1e3);

  t0 = Clock::now();
  const auto cc = connected_components(c);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms  %8.2f Mcells/s\n", "components(whole)", ms(t0, t1),
              total_cells / ms(t0, t1) / 1e3);

  t0 = Clock::now();
  const auto o = orient(c);
  t1 = Clock::now();
  std::printf("  %-22s %10.1f ms  %8.2f Mtets/s  (orientable=%d)\n", "orient", ms(t0, t1),
              ntets / ms(t0, t1) / 1e3, int(o.orientable));
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 1) {
    for (int a = 1; a < argc; ++a) run(std::atoi(argv[a]));
  } else {
    run(16);
    run(32);
  }
  return 0;
}
