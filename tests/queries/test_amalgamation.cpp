#include "fixtures.hpp"
#include "graphos/core/build.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos/queries/amalgamation.hpp"
#include "graphos_test.hpp"

using graphos::Index;

// The hypothesis of the amalgamation lemma: two 2-cells sharing a single
// edge — the common boundary is an interval (a 1-ball), the intersection
// is proper, and the union is again a cell.
GRAPHOS_TEST(cells_sharing_one_facet_amalgamate) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto shared = graphos::common_boundary(c, 2, 0, 1);
  CHECK(shared.n_facets == 1);
  CHECK(shared.components == 1);
  CHECK(shared.acyclic);
  const graphos::Marker excess = graphos::excess_intersection(c, 2, 0, 1);
  for (int k = 0; k < 2; ++k) CHECK(excess.marked_count(k) == 0);
  CHECK(graphos::amalgamates_to_cell(c, 2, 0, 1));
}

// Improper intersection: a quad and a pentagon sharing one edge, with the
// pentagon also touching the quad's far vertex. The excess intersection is
// exactly that vertex — the amalgamated cell would be pinched there.
GRAPHOS_TEST(excess_intersection_detects_the_pinch) {
  const graphos::Complex c = graphos::from_polygons(8, {{0, 1, 2, 3}, {1, 0, 6, 2, 7}});
  const auto shared = graphos::common_boundary(c, 2, 0, 1);
  CHECK(shared.n_facets == 1);  // only the edge (0,1) is common boundary
  CHECK(shared.acyclic);
  const graphos::Marker excess = graphos::excess_intersection(c, 2, 0, 1);
  CHECK(excess.marked_count(0) == 1);
  CHECK(excess.marked(0, 2));  // cl(σ) ∩ cl(τ) beyond cl(∂σ ∩ ∂τ)
  CHECK(!graphos::amalgamates_to_cell(c, 2, 0, 1));
}

// A circular common boundary: the two triangles of a cylinder share two
// edges forming a 1-cycle — connected but not acyclic (b₁ = 1), so no
// ball, no amalgamation (the union would be an annulus, not a cell).
GRAPHOS_TEST(circular_common_boundary_is_not_acyclic) {
  graphos::Complex square(2);
  square.attach_vertices(4);
  square.attach_cell(1, {0, 1}, {-1, +1});
  square.attach_cell(1, {1, 2}, {-1, +1});
  square.attach_cell(1, {3, 2}, {-1, +1});
  square.attach_cell(1, {0, 3}, {-1, +1});
  square.attach_cell(1, {0, 2}, {-1, +1});
  square.attach_cell(2, {0, 1, 4}, {+1, +1, -1});
  square.attach_cell(2, {4, 2, 3}, {+1, -1, -1});
  const auto ids = graphos::lift_identifications(square, {{0, 1, +1}, {3, 2, +1}});
  const auto cylinder = graphos::quotient(square, ids);

  const auto shared = graphos::common_boundary(cylinder.complex, 2, 0, 1);
  CHECK(shared.n_facets == 2);
  CHECK(shared.components == 1);
  CHECK(shared.betti[1] == 1);  // the common boundary is a circle
  CHECK(!shared.acyclic);
  CHECK(!graphos::amalgamates_to_cell(cylinder.complex, 2, 0, 1));
}

GRAPHOS_TEST(disjoint_cells_share_no_boundary) {
  const graphos::Complex c = graphos::from_polygons(8, {{0, 1, 2, 3}, {4, 5, 6, 7}});
  const auto shared = graphos::common_boundary(c, 2, 0, 1);
  CHECK(shared.n_facets == 0);
  CHECK(!shared.acyclic);
  CHECK(!graphos::amalgamates_to_cell(c, 2, 0, 1));
}

GRAPHOS_TEST_MAIN()
