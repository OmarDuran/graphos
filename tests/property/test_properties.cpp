// Property-based fuzzing: random complexes through random operation
// sequences, asserting the laws of the calculus after every step — the layer
// that reaches the cases no fixture enumerates.

#include <random>
#include <set>
#include <vector>

#include "graphos/graphos.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

// the structural invariants, and the Euler–Poincaré identity
// χ = Σ_k (−1)^k β_k over Z₂
void assert_laws(const graphos::Complex& c) {
  c.validate();
  CHECK(graphos::d_squared_is_zero(c));
  const auto betti = graphos::betti_numbers_z2(c);
  long long chi = 0;
  for (std::size_t k = 0; k < betti.size(); ++k) {
    chi += (k % 2 == 0 ? 1LL : -1LL) * betti[k];
  }
  CHECK(chi == graphos::euler_characteristic(c));
}

void assert_chain_map_valid(const graphos::ChainMap& m, const graphos::Complex& target) {
  for (std::size_t k = 0; k < m.index.size(); ++k) {
    for (const Index i : m.index[k]) {
      CHECK(i == graphos::invalid_index || (i >= 0 && i < target.count(static_cast<int>(k))));
    }
  }
}

graphos::Complex random_simplicial(std::mt19937& rng) {
  const int dim = std::uniform_int_distribution<int>(2, 3)(rng);
  const Index nv = std::uniform_int_distribution<Index>(dim + 2, 10)(rng);
  const int ncells = std::uniform_int_distribution<int>(2, 12)(rng);
  std::set<std::vector<Index>> used;
  std::vector<std::vector<Index>> cells;
  std::uniform_int_distribution<Index> pick(0, nv - 1);
  // bounded attempts: a small vertex pool need not admit ncells distinct
  // simplices, since C(nv, dim+1) may be smaller
  for (int tries = 0; tries < 400 && static_cast<int>(cells.size()) < ncells; ++tries) {
    std::set<Index> verts;
    while (static_cast<int>(verts.size()) < dim + 1) verts.insert(pick(rng));
    std::vector<Index> cell(verts.begin(), verts.end());
    if (used.insert(cell).second) cells.push_back(std::move(cell));
  }
  return graphos::from_simplices(dim, nv, cells);
}

graphos::Marker random_marker(std::mt19937& rng, const graphos::Complex& c, double p) {
  graphos::Marker m(c);
  std::bernoulli_distribution flip(p);
  for (int k = 0; k <= c.dim(); ++k) {
    m.mark_where(k, [&](Index) { return flip(rng); });
  }
  return m;
}

}  // namespace

GRAPHOS_TEST(random_complexes_through_random_ops) {
  for (unsigned seed = 0; seed < 40; ++seed) {
    std::mt19937 rng(seed);
    graphos::Complex c = random_simplicial(rng);
    assert_laws(c);

    for (int step = 0; step < 3; ++step) {
      const int op = std::uniform_int_distribution<int>(0, 6)(rng);
      switch (op) {
        case 0: {
          const auto r = graphos::star_deletion(c, random_marker(rng, c, 0.15));
          assert_chain_map_valid(r.map, r.complex);
          c = r.complex;
          break;
        }
        case 1: {
          const auto r = graphos::subcomplex(c, random_marker(rng, c, 0.3));
          assert_chain_map_valid(r.map, r.complex);
          // the embedding inverts the map on survivors
          for (std::size_t k = 0; k < r.embedding.index.size(); ++k) {
            for (std::size_t i = 0; i < r.embedding.index[k].size(); ++i) {
              const Index parent = r.embedding.index[k][i];
              CHECK(r.map.index[k][static_cast<std::size_t>(parent)] == static_cast<Index>(i));
            }
          }
          c = r.complex;
          break;
        }
        case 2: {
          if (c.dim() < 1 || c.count(c.dim() - 1) == 0) break;
          graphos::Marker interface(c);
          std::bernoulli_distribution flip(0.2);
          interface.mark_where(c.dim() - 1, [&](Index) { return flip(rng); });
          const auto r = graphos::cut_along(c, interface);
          assert_chain_map_valid(r.ancestor, c);  // ancestors index the PARENT
          c = r.complex;
          break;
        }
        case 3: {
          const auto r = graphos::orient(c);
          c = r.complex;
          break;
        }
        case 4: {
          // agglomerate by facet-connected components; an incoherent
          // orientation is legitimately rejected, and that path is under test
          if (c.dim() < 1 || c.count(c.dim()) == 0) break;
          const auto labels = graphos::connected_components(c, c.dim(), c.dim() - 1);
          try {
            const auto r = graphos::agglomerate(c, labels.label);
            assert_chain_map_valid(r.map, r.complex);
            c = r.complex;
          } catch (const std::invalid_argument&) {
            // the documented rejection: the stratum is not coherent
          }
          break;
        }
        case 5: {
          const auto before = graphos::betti_numbers_z2(c);
          const auto r = graphos::collapse(c);
          const auto after = graphos::betti_numbers_z2(r.complex);
          CHECK(before == after);  // simple homotopy equivalence
          c = r.complex;
          break;
        }
        case 6: {
          // dimensional descent by facet-connected components; rejection on
          // orientation is a documented path
          if (c.dim() < 1 || c.count(c.dim()) == 0) break;
          const auto labels = graphos::connected_components(c, c.dim(), c.dim() - 1);
          try {
            const auto r = graphos::coarsen(c, labels.label);
            assert_chain_map_valid(r.map, r.complex);
            c = r.complex;
          } catch (const std::invalid_argument&) {
          }
          break;
        }
      }
      assert_laws(c);
    }
  }
}

GRAPHOS_TEST(disjoint_union_is_additive) {
  for (unsigned seed = 100; seed < 120; ++seed) {
    std::mt19937 rng(seed);
    const graphos::Complex a = random_simplicial(rng);
    const graphos::Complex b = random_simplicial(rng);
    const auto du = graphos::disjoint_union(a, b);
    assert_laws(du.complex);
    CHECK(graphos::euler_characteristic(du.complex) ==
          graphos::euler_characteristic(a) + graphos::euler_characteristic(b));
    CHECK(graphos::connected_components(du.complex).count ==
          graphos::connected_components(a).count + graphos::connected_components(b).count);
  }
}

GRAPHOS_TEST(subdivision_preserves_homotopy_of_random_complexes) {
  for (unsigned seed = 200; seed < 210; ++seed) {
    std::mt19937 rng(seed);
    const graphos::Complex c = random_simplicial(rng);
    const auto sd = graphos::barycentric_subdivision(c);
    assert_laws(sd.complex);
    CHECK(graphos::betti_numbers_z2(sd.complex) == graphos::betti_numbers_z2(c));
    CHECK(graphos::euler_characteristic(sd.complex) == graphos::euler_characteristic(c));

    // the transfer law on the orientable ones: the signed carrier takes a
    // coherent orientation of the parent to one of sd(C)
    const auto o = graphos::orient(c);
    if (!o.orientable) continue;
    const auto osd = graphos::barycentric_subdivision(o.complex);
    const int n = osd.complex.dim();
    const graphos::CoboundaryOperator cob = graphos::coboundary(osd.complex, n - 1);
    for (Index f = 0; f < osd.complex.count(n - 1); ++f) {
      const Index lo = cob.offsets[static_cast<std::size_t>(f)];
      const Index hi = cob.offsets[static_cast<std::size_t>(f) + 1];
      if (hi - lo != 2) continue;
      int sum = 0;
      for (Index m = lo; m < hi; ++m) {
        sum += cob.signs[static_cast<std::size_t>(m)] *
               osd.carrier_sign[static_cast<std::size_t>(n)]
                               [static_cast<std::size_t>(cob.indices[static_cast<std::size_t>(m)])];
      }
      CHECK(sum == 0);
    }
  }
}

GRAPHOS_TEST_MAIN()
