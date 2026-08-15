#include "fixtures.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/product.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos/ops/star_deletion.hpp"
#include "graphos/queries/facets.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos/queries/neighborhood.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

graphos::Complex make_circle() {
  graphos::Complex c(1);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  return c;
}

long long alternating_sum(const std::vector<Index>& betti) {
  long long chi = 0;
  for (std::size_t k = 0; k < betti.size(); ++k) {
    chi += (k % 2 == 0 ? 1LL : -1LL) * betti[k];
  }
  return chi;
}

}  // namespace

GRAPHOS_TEST(disk_is_contractible) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto b = graphos::betti_numbers_z2(c);
  CHECK(b == (std::vector<Index>{1, 0, 0}));
  CHECK(alternating_sum(b) == graphos::euler_characteristic(c));
}

GRAPHOS_TEST(circle_has_one_tunnel) {
  const auto b = graphos::betti_numbers_z2(make_circle());
  CHECK(b == (std::vector<Index>{1, 1}));
}

GRAPHOS_TEST(torus_has_two_tunnels_and_a_cavity) {
  const auto prod = graphos::product(make_circle(), make_circle());
  const auto b = graphos::betti_numbers_z2(prod.complex);
  CHECK(b == (std::vector<Index>{1, 2, 1}));
  CHECK(alternating_sum(b) == 0);
}

GRAPHOS_TEST(moebius_band_retracts_to_a_circle) {
  graphos::Complex square(2);
  square.attach_vertices(4);
  square.attach_cell(1, {0, 1}, {-1, +1});
  square.attach_cell(1, {1, 2}, {-1, +1});
  square.attach_cell(1, {3, 2}, {-1, +1});
  square.attach_cell(1, {0, 3}, {-1, +1});
  square.attach_cell(1, {0, 2}, {-1, +1});
  square.attach_cell(2, {0, 1, 4}, {+1, +1, -1});
  square.attach_cell(2, {4, 2, 3}, {+1, -1, -1});
  const auto ids = graphos::lift_identifications(square, {{0, 2, +1}, {3, 1, +1}});
  const auto moebius = graphos::quotient(square, ids);
  const auto b = graphos::betti_numbers_z2(moebius.complex);
  CHECK(b == (std::vector<Index>{1, 1, 0}));
}

GRAPHOS_TEST(relative_betti_of_disk_mod_boundary_is_a_sphere) {
  // H_k(D², S¹; Z₂) ≅ H̃_k(S²): one class in the top dimension
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto rel = graphos::betti_numbers_z2(c, graphos::classify_facets(c).boundary);
  CHECK(rel == (std::vector<Index>{0, 0, 1}));
}

GRAPHOS_TEST(relative_betti_with_empty_subcomplex_is_absolute) {
  const graphos::Complex c = graphos_test::make_fan();
  CHECK(graphos::betti_numbers_z2(c, graphos::Marker(c)) == graphos::betti_numbers_z2(c));
}

// Witnesses excision: H_k(K, cl(st S)) ≅ H_k(K ∖ st S, frontier(S)), with the
// two sides computed independently.
GRAPHOS_TEST(excision_isomorphism_holds) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker s(c);
  s.mark(2, 0);  // excise face A

  // left: relative to cl(st S) on the full complex
  const auto lhs = graphos::betti_numbers_z2(c, graphos::star_of(c, s));

  // right: excise, then relative to the frontier, transported through the
  // deletion's chain map
  const graphos::Marker frontier = graphos::frontier_of(c, s);
  const auto excised = graphos::star_deletion(c, s);
  graphos::Marker frontier_after(excised.complex);
  for (int k = 0; k <= 2; ++k) {
    for (graphos::Index i = 0; i < c.count(k); ++i) {
      if (frontier.marked(k, i)) {
        frontier_after.mark(k, excised.map.index[k][static_cast<std::size_t>(i)]);
      }
    }
  }
  const auto rhs = graphos::betti_numbers_z2(excised.complex, frontier_after);
  CHECK(lhs == rhs);
}

GRAPHOS_TEST(cut_disconnects_beta_zero) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  graphos::Marker interface(c);
  interface.mark(1, 0);
  const auto cut = graphos::cut_along(c, interface);
  const auto b = graphos::betti_numbers_z2(cut.complex);
  CHECK(b == (std::vector<Index>{3, 0, 0}));  // two triangles + the segment
  CHECK(alternating_sum(b) == graphos::euler_characteristic(cut.complex));
}

GRAPHOS_TEST_MAIN()
