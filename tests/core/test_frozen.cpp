#include <vector>

#include "fixtures.hpp"
#include "graphos/core/frozen.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

std::vector<Index> indices_of(const graphos::FrozenComplex::Row& r) {
  return {r.indices, r.indices + r.size};
}

std::vector<int> signs_of(const graphos::FrozenComplex::Row& r) {
  return {r.signs, r.signs + r.size};
}

}  // namespace

GRAPHOS_TEST(freeze_preserves_counts_and_boundary) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  CHECK(f.dim() == 2);
  CHECK(f.count(0) == 3);
  CHECK(f.count(1) == 3);
  CHECK(f.count(2) == 1);
  CHECK(indices_of(f.boundary_row(2, 0)) == (std::vector<Index>{0, 1, 2}));
  CHECK(signs_of(f.boundary_row(2, 0)) == (std::vector<int>{+1, +1, +1}));
  CHECK(indices_of(f.boundary_row(1, 0)) == (std::vector<Index>{0, 1}));
  CHECK(signs_of(f.boundary_row(1, 0)) == (std::vector<int>{-1, +1}));
}

GRAPHOS_TEST(freeze_derives_signed_coboundary) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  // δ_0 row of v0: tail of e0 (-1), head of e2 (+1)
  CHECK(indices_of(f.coboundary_row(0, 0)) == (std::vector<Index>{0, 2}));
  CHECK(signs_of(f.coboundary_row(0, 0)) == (std::vector<int>{-1, +1}));
  // δ_1 row of each edge: the single face, coefficient +1
  CHECK(indices_of(f.coboundary_row(1, 1)) == (std::vector<Index>{0}));
  CHECK(signs_of(f.coboundary_row(1, 1)) == (std::vector<int>{+1}));
}

GRAPHOS_TEST(views_expose_csr_arrays) {
  graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  const graphos::CsrView v = f.boundary_view(2);
  CHECK(v.offsets[1] == 3);
  CHECK(v.indices[0] == 0);
  CHECK(v.signs[0] == 1);
}

GRAPHOS_TEST(closure_of_a_face_is_the_whole_triangle) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  const auto cl = f.closure(2, 0);
  CHECK(cl[0] == (std::vector<Index>{0, 1, 2}));
  CHECK(cl[1] == (std::vector<Index>{0, 1, 2}));
  CHECK(cl[2] == (std::vector<Index>{0}));
}

GRAPHOS_TEST(star_of_a_boundary_vertex) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  const auto st = f.star(0, 0);
  CHECK(st[0] == (std::vector<Index>{0}));
  CHECK(st[1] == (std::vector<Index>{0, 2}));
  CHECK(st[2] == (std::vector<Index>{0}));
}

GRAPHOS_TEST(star_of_the_fan_center_is_the_whole_interior) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_fan());
  const auto st = f.star(0, 4);
  CHECK(st[0] == (std::vector<Index>{4}));
  CHECK(st[1] == (std::vector<Index>{0, 1, 2, 3}));  // the four spokes
  CHECK(st[2] == (std::vector<Index>{0, 1, 2, 3}));  // all four triangles
}

GRAPHOS_TEST(link_of_an_interior_vertex_is_a_cycle) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_fan());
  const auto lk = f.link(0, 4);
  CHECK(lk[0] == (std::vector<Index>{0, 1, 2, 3}));  // the corners
  CHECK(lk[1] == (std::vector<Index>{4, 5, 6, 7}));  // the boundary edges
  CHECK(lk[2].empty());
}

GRAPHOS_TEST(link_of_a_boundary_vertex_is_an_interval) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  const auto lk = f.link(0, 0);
  CHECK(lk[0] == (std::vector<Index>{1, 2}));
  CHECK(lk[1] == (std::vector<Index>{1}));  // the opposite edge
  CHECK(lk[2].empty());
}

GRAPHOS_TEST(dual_view_reverses_the_poset) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  const graphos::DualView d = graphos::dual(f);
  CHECK(d.dim() == 2);
  CHECK(d.count(0) == 1);  // the face
  CHECK(d.count(1) == 3);  // the edges
  CHECK(d.count(2) == 3);  // the vertices

  // ∂^dual of a dual 1-cell (a primal edge) reaches the dual 0-cell of its
  // primal coface
  const auto r = d.boundary_row(1, 0);
  CHECK(r.size == 1);
  CHECK(r.indices[0] == 0);
  // ∂^dual of a dual 2-cell (a primal vertex) reaches its two primal edges
  CHECK(d.boundary_row(2, 0).size == 2);
  // δ^dual of the dual 0-cell is the primal face boundary
  CHECK(d.coboundary_row(0, 0).size == 3);
}

GRAPHOS_TEST(dual_of_dual_rows_match_primal) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_two_triangle_disk());
  const graphos::DualView d = graphos::dual(f);
  // dual boundary of a dual 1-cell = primal coboundary of the mirrored edge
  const auto dr = d.boundary_row(1, 0);
  const auto pr = f.coboundary_row(1, 0);
  CHECK(dr.size == pr.size);
  for (graphos::Index m = 0; m < dr.size; ++m) {
    CHECK(dr.indices[m] == pr.indices[m]);
    CHECK(dr.signs[m] == pr.signs[m]);
  }
}

GRAPHOS_TEST(halo_depth_is_recorded_and_validated) {
  const graphos::Complex c = graphos_test::make_triangle();
  CHECK(graphos::freeze(c).halo_depth() == 1);
  CHECK(graphos::freeze(c, 3).halo_depth() == 3);
  CHECK_THROWS(graphos::freeze(c, 0));
}

GRAPHOS_TEST(queries_reject_bad_input) {
  const graphos::FrozenComplex f = graphos::freeze(graphos_test::make_triangle());
  CHECK_THROWS(f.closure(3, 0));
  CHECK_THROWS(f.star(0, 99));
  CHECK_THROWS(f.boundary_row(0, 0));    // ∂_0 does not exist
  CHECK_THROWS(f.coboundary_row(2, 0));  // δ_dim does not exist
}

GRAPHOS_TEST_MAIN()
