#include "graphos/exec/array.hpp"

#include <utility>
#include <vector>

#include "graphos/core/types.hpp"
#include "graphos_test.hpp"

using graphos::Index;
using graphos::exec::Array;

GRAPHOS_TEST(constructs_from_vector_and_roundtrips) {
  const std::vector<Index> src = {5, 4, 3, 2, 1};
  const Array<Index> a(src);
  CHECK(a.size() == 5);
  for (std::size_t i = 0; i < 5; ++i) CHECK(a[i] == src[i]);
}

GRAPHOS_TEST(move_transfers_ownership) {
  Array<Index> a(std::vector<Index>{7, 8});
  Array<Index> b(std::move(a));
  CHECK(b.size() == 2);
  CHECK(b[1] == 8);
  CHECK(a.size() == 0);

  Array<Index> c;
  c = std::move(b);
  CHECK(c.size() == 2);
  CHECK(c[0] == 7);
  CHECK(b.size() == 0);
}

GRAPHOS_TEST(zero_size_has_null_data) {
  const Array<Index> empty;
  CHECK(empty.size() == 0);
  CHECK(empty.data() == nullptr);
}

GRAPHOS_TEST(view_is_kernel_indexable) {
  Array<Index> a(std::vector<Index>{1, 2, 3});
  auto v = a.view();  // ManagedArray with CHAI, raw pointer on host fallback
  CHECK(v[0] == 1);
  CHECK(v[2] == 3);
}

GRAPHOS_TEST_MAIN()
