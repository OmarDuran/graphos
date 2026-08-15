// Scaling assertions: each operation has a complexity model, and growing the
// complex 8× must not grow the time as a quadratic algorithm would, by ~64×.
// The bounds are loose by design — this catches a change of complexity class,
// not of constant factor, which bench/ measures.

#include <chrono>
#include <cstdio>

#include "generators.hpp"
#include "graphos/graphos.hpp"
#include "graphos_test.hpp"

using Clock = std::chrono::steady_clock;
using graphos::Index;

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

}  // namespace

GRAPHOS_TEST(ops_scale_subquadratically) {
  const int n_small = 8, n_large = 16;  // 8x the tets
  const auto small_tets = graphos_bench::kuhn_tets(n_small);
  const auto large_tets = graphos_bench::kuhn_tets(n_large);

  const double build_s = best_of_3_ms(
      [&] { graphos::from_simplices(3, graphos_bench::kuhn_vertex_count(n_small), small_tets); });
  const double build_l = best_of_3_ms(
      [&] { graphos::from_simplices(3, graphos_bench::kuhn_vertex_count(n_large), large_tets); });

  graphos::Complex cs =
      graphos::from_simplices(3, graphos_bench::kuhn_vertex_count(n_small), small_tets);
  graphos::Complex cl =
      graphos::from_simplices(3, graphos_bench::kuhn_vertex_count(n_large), large_tets);

  const double freeze_s = best_of_3_ms([&] { graphos::freeze(cs); });
  const double freeze_l = best_of_3_ms([&] { graphos::freeze(cl); });
  const double inc_s = best_of_3_ms([&] { graphos::incidence(cs, 3, 0); });
  const double inc_l = best_of_3_ms([&] { graphos::incidence(cl, 3, 0); });

  graphos::Marker ms(cs), ml(cl);
  ms.mark(0, cs.count(0) / 2);
  ml.mark(0, cl.count(0) / 2);
  const double sd_s = best_of_3_ms([&] { graphos::star_deletion(cs, ms); });
  const double sd_l = best_of_3_ms([&] { graphos::star_deletion(cl, ml); });

  // quadratic growth at 8× shows as ~64×; the headroom absorbs timer noise on
  // small inputs
  const auto check_ratio = [](const char* name, double s, double l, double bound) {
    const double ratio = l / std::max(s, 0.01);
    std::printf("  scaling %-16s %6.2fx (bound %.0fx)\n", name, ratio, bound);
    CHECK(ratio < bound);
  };
  check_ratio("from_simplices", build_s, build_l, 30.0);  // n log n allowed
  check_ratio("freeze", freeze_s, freeze_l, 24.0);
  check_ratio("incidence(3,0)", inc_s, inc_l, 24.0);
  check_ratio("star_deletion", sd_s, sd_l, 24.0);
}

GRAPHOS_TEST_MAIN()
