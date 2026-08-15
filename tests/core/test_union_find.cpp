#include "graphos/core/union_find.hpp"
#include "graphos_test.hpp"

using graphos::Index;
using graphos::UnionFind;

GRAPHOS_TEST(singletons_are_their_own_roots) {
  UnionFind uf(4);
  CHECK(uf.size() == 4);
  for (Index i = 0; i < 4; ++i) CHECK(uf.find(i) == i);
}

GRAPHOS_TEST(unions_merge_and_chains_resolve) {
  UnionFind uf(6);
  uf.unite(0, 1);
  uf.unite(2, 3);
  uf.unite(1, 2);  // chain: {0,1,2,3}
  CHECK(uf.find(3) == uf.find(0));
  CHECK(uf.find(1) == uf.find(2));
  CHECK(uf.find(4) != uf.find(0));
  CHECK(uf.find(4) != uf.find(5));
  uf.unite(5, 5);  // self-union is a no-op
  CHECK(uf.find(5) == 5);
}

GRAPHOS_TEST_MAIN()
