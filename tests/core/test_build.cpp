#include "graphos/core/build.hpp"

#include "graphos/queries/homology.hpp"
#include "graphos/queries/facets.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos/queries/manifold.hpp"
#include "graphos/ops/orient.hpp"

#include "fixtures.hpp"
#include "graphos_test.hpp"

using graphos::Index;

GRAPHOS_TEST(edges_build_a_path) {
  const graphos::Complex c = graphos::from_edges(3, {{0, 1}, {1, 2}});
  c.validate();
  CHECK(c.count(0) == 3);
  CHECK(c.count(1) == 2);
  CHECK(graphos::euler_characteristic(c) == 1);
}

// The polytopal requirement: a quad and a triangle in ONE mesh, sharing a
// derived edge.
GRAPHOS_TEST(mixed_polygon_mesh) {
  const graphos::Complex c = graphos::from_polygons(5, {{0, 1, 2, 3}, {1, 4, 2}});
  c.validate();
  CHECK(c.count(0) == 5);
  CHECK(c.count(1) == 6);  // 4 + 3 with the shared edge derived once
  CHECK(c.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);
  CHECK(graphos::check_manifold(c).manifold_like);
}

GRAPHOS_TEST(triangle_pair_matches_hand_built_counts) {
  const graphos::Complex c = graphos::from_polygons(4, {{0, 1, 2}, {1, 0, 3}});
  c.validate();
  CHECK(c.count(0) == 4);
  CHECK(c.count(1) == 5);
  CHECK(c.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(c));
  // opposite windings across the shared edge: consistently oriented as built
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  for (std::size_t i = 0; i < o.map.sign[2].size(); ++i) CHECK(o.map.sign[2][i] == +1);
}

GRAPHOS_TEST(a_cube_from_its_face_cycles) {
  const std::vector<std::vector<Index>> cube = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  const graphos::Complex c = graphos::from_polyhedra(8, {cube});
  c.validate();
  CHECK(c.count(0) == 8);
  CHECK(c.count(1) == 12);
  CHECK(c.count(2) == 6);
  CHECK(c.count(3) == 1);
  CHECK(graphos::d_squared_is_zero(c));  // the whole sign machinery
  CHECK(graphos::euler_characteristic(c) == 1);
  CHECK(graphos::betti_numbers_z2(c) == (std::vector<Index>{1, 0, 0, 0}));
}

GRAPHOS_TEST(two_cubes_share_a_derived_face) {
  const std::vector<std::vector<Index>> lower = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                                 {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  const std::vector<std::vector<Index>> upper = {{4, 7, 6, 5}, {8, 9, 10, 11}, {4, 5, 9, 8},
                                                 {5, 6, 10, 9}, {6, 7, 11, 10}, {7, 4, 8, 11}};
  const graphos::Complex c = graphos::from_polyhedra(12, {lower, upper});
  c.validate();
  CHECK(c.count(0) == 12);
  CHECK(c.count(1) == 20);
  CHECK(c.count(2) == 11);  // 6 + 6 with the interface derived once
  CHECK(c.count(3) == 2);
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);
  // outward windings on both cells: consistent as built
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  for (std::size_t i = 0; i < o.map.sign[3].size(); ++i) CHECK(o.map.sign[3][i] == +1);
  CHECK(graphos::check_manifold(c).manifold_like);
}

GRAPHOS_TEST(tetrahedra_from_meshio_style_connectivity) {
  const graphos::Complex one = graphos::from_simplices(3, 4, {{0, 1, 2, 3}});
  one.validate();
  CHECK(one.count(0) == 4);
  CHECK(one.count(1) == 6);
  CHECK(one.count(2) == 4);
  CHECK(one.count(3) == 1);
  CHECK(graphos::d_squared_is_zero(one));

  const graphos::Complex two = graphos::from_simplices(3, 5, {{0, 1, 2, 3}, {1, 3, 2, 4}});
  two.validate();
  CHECK(two.count(0) == 5);
  CHECK(two.count(1) == 9);
  CHECK(two.count(2) == 7);  // shared face derived once
  CHECK(two.count(3) == 2);
  CHECK(graphos::d_squared_is_zero(two));
  CHECK(graphos::euler_characteristic(two) == 1);
  CHECK(graphos::orient(two).orientable);
}

// d-simplices: a single pentachoron (4-simplex), all strata derived —
// binomial counts C(5, k+1).
GRAPHOS_TEST(pentachoron_in_dimension_four) {
  const graphos::Complex c = graphos::from_simplices(4, 5, {{0, 1, 2, 3, 4}});
  c.validate();
  CHECK(c.count(0) == 5);
  CHECK(c.count(1) == 10);
  CHECK(c.count(2) == 10);
  CHECK(c.count(3) == 5);
  CHECK(c.count(4) == 1);
  CHECK(graphos::d_squared_is_zero(c));  // simplex signs at every level
  CHECK(graphos::euler_characteristic(c) == 1);
  CHECK(graphos::betti_numbers_z2(c) == (std::vector<Index>{1, 0, 0, 0, 0}));
}

// The boundary of the pentachoron is the 3-sphere: extract it with the
// generic ops and compute its homology — the whole stack in dimension 4.
GRAPHOS_TEST(pentachoron_boundary_is_the_three_sphere) {
  const graphos::Complex c = graphos::from_simplices(4, 5, {{0, 1, 2, 3, 4}});
  const auto s3 = graphos::subcomplex(c, graphos::classify_facets(c).boundary);
  s3.complex.validate();
  CHECK(s3.complex.count(3) == 5);
  CHECK(graphos::d_squared_is_zero(s3.complex));
  CHECK(graphos::euler_characteristic(s3.complex) == 0);  // χ(S³)
  CHECK(graphos::betti_numbers_z2(s3.complex) == (std::vector<Index>{1, 0, 0, 1, 0}));
}

// Two pentachora sharing a tetrahedral facet, given in orientations that
// agree — consistent as built, in dimension 4.
GRAPHOS_TEST(two_pentachora_share_a_derived_tetrahedron) {
  const graphos::Complex c =
      graphos::from_simplices(4, 6, {{0, 1, 2, 3, 4}, {1, 0, 2, 3, 5}});
  c.validate();
  CHECK(c.count(4) == 2);
  CHECK(c.count(3) == 9);  // 5 + 5 with the shared tet derived once
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  for (std::size_t i = 0; i < o.map.sign[4].size(); ++i) CHECK(o.map.sign[4][i] == +1);
}

GRAPHOS_TEST(rejects_malformed_input) {
  CHECK_THROWS(graphos::from_polygons(3, {{0, 1}}));           // 2-gon
  CHECK_THROWS(graphos::from_polygons(3, {{0, 1, 5}}));        // out of range
  CHECK_THROWS(graphos::from_polygons(4, {{0, 1, 1, 2}}));     // repeated vertex
  CHECK_THROWS(graphos::from_simplices(2, 4, {{0, 1, 2, 3}})); // wrong arity
  // same vertex set, cyclically incompatible orderings
  CHECK_THROWS(graphos::from_polyhedra(
      6, {{{0, 1, 2, 3}, {0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
          {{0, 2, 1, 3}, {0, 3, 1, 2}, {0, 1, 5}, {1, 2, 5}, {2, 3, 5}, {3, 0, 5}}}));
}

GRAPHOS_TEST_MAIN()
