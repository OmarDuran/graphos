#include "graphos/exec/memory.hpp"

#include <utility>

#include "graphos/core/types.hpp"
#include "graphos_test.hpp"

using graphos::Index;
using graphos::exec::Buffer;

GRAPHOS_TEST(buffer_allocates_and_indexes) {
  Buffer<Index> b(8);
  CHECK(b.size() == 8);
  CHECK(b.data() != nullptr);
  b[7] = 41;
  CHECK(b[7] == 41);
}

GRAPHOS_TEST(buffer_move_transfers_ownership) {
  Buffer<Index> a(4);
  a[0] = 42;
  Buffer<Index> b(std::move(a));
  CHECK(b.size() == 4);
  CHECK(b[0] == 42);
  CHECK(a.data() == nullptr);
  CHECK(a.size() == 0);

  Buffer<Index> c;
  c = std::move(b);
  CHECK(c.size() == 4);
  CHECK(c[0] == 42);
  CHECK(b.data() == nullptr);
}

GRAPHOS_TEST(buffer_zero_size_does_not_allocate) {
  const Buffer<Index> empty(0);
  CHECK(empty.size() == 0);
  CHECK(empty.data() == nullptr);
}

GRAPHOS_TEST_MAIN()
