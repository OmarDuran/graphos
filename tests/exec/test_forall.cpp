#include "graphos/exec/forall.hpp"
#include "graphos/exec/memory.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(forall_visits_every_index_once) {
  graphos::exec::Buffer<Index> b(100);
  Index* p = b.data();
  graphos::exec::forall(100, [=](Index i) { p[i] = 2 * i; });
  for (Index i = 0; i < 100; ++i) CHECK(b[static_cast<std::size_t>(i)] == 2 * i);
}

GRAPHOS_TEST(forall_zero_range_is_noop) {
  bool touched = false;
  graphos::exec::forall(0, [&](Index) { touched = true; });
  CHECK(!touched);
}

GRAPHOS_TEST(exclusive_scan_returns_total) {
  const Index in[3] = {1, 2, 3};
  Index out[3] = {-1, -1, -1};
  const Index total = graphos::exec::exclusive_scan(in, out, 3);
  CHECK(total == 6);
  CHECK(out[0] == 0);
  CHECK(out[1] == 1);
  CHECK(out[2] == 3);
  CHECK(graphos::exec::exclusive_scan(in, out, 0) == 0);
}

GRAPHOS_TEST(inclusive_scan_inplace) {
  Index a[4] = {0, 1, 2, 3};  // offsets-style: leading zero then sizes
  graphos::exec::inclusive_scan_inplace(a, 4);
  CHECK(a[0] == 0);
  CHECK(a[1] == 1);
  CHECK(a[2] == 3);
  CHECK(a[3] == 6);
}

GRAPHOS_TEST_MAIN()
