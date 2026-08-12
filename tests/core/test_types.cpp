#include "graphos/core/types.hpp"

#include "graphos_test.hpp"

using graphos::ChainMap;
using graphos::invalid_index;
using graphos::Sign;

GRAPHOS_TEST(sized_allocates_identity_signs) {
  const ChainMap m = ChainMap::sized({2, 3});
  CHECK(m.index.size() == 2);
  CHECK(m.sign.size() == 2);
  CHECK(m.index[0].size() == 2);
  CHECK(m.index[1].size() == 3);
  for (const Sign s : m.sign[1]) CHECK(s == 1);
}

GRAPHOS_TEST(compose_follows_indices_and_multiplies_signs) {
  ChainMap first = ChainMap::sized({2});
  first.index[0] = {1, invalid_index};
  first.sign[0] = {-1, 1};
  ChainMap second = ChainMap::sized({2});
  second.index[0] = {invalid_index, 0};
  second.sign[0] = {1, -1};

  const ChainMap out = graphos::compose(first, second);
  CHECK(out.index[0][0] == 0);
  CHECK(out.sign[0][0] == 1);  // (-1) * (-1)
  CHECK(out.index[0][1] == invalid_index);
}

GRAPHOS_TEST(compose_sends_zero_to_zero) {
  ChainMap first = ChainMap::sized({1});
  first.index[0] = {0};
  ChainMap second = ChainMap::sized({1});
  second.index[0] = {invalid_index};

  const ChainMap out = graphos::compose(first, second);
  CHECK(out.index[0][0] == invalid_index);
  CHECK(out.sign[0][0] == 1);  // signs of zero cells stay neutral
}

GRAPHOS_TEST_MAIN()
