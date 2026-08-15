#include "fixtures.hpp"
#include "graphos/core/coboundary.hpp"
#include "graphos/ops/disjoint_union.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/orient.hpp"
#include "graphos/ops/product.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

// coherence: every interior facet is induced with opposite incidence numbers
// by its two cofaces
bool consistently_oriented(const graphos::Complex& c) {
  const int n = c.dim();
  const graphos::CoboundaryOperator cob = graphos::coboundary(c, n - 1);
  for (Index f = 0; f < c.count(n - 1); ++f) {
    if (cob.offsets[f + 1] - cob.offsets[f] != 2) continue;
    int sum = 0;
    for (Index m = cob.offsets[f]; m < cob.offsets[f + 1]; ++m) sum += cob.signs[m];
    if (sum != 0) return false;
  }
  return true;
}

graphos::Complex make_square() {
  graphos::Complex c(2);
  c.attach_vertices(4);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {3, 2}, {-1, +1});
  c.attach_cell(1, {0, 3}, {-1, +1});
  c.attach_cell(1, {0, 2}, {-1, +1});
  c.attach_cell(2, {0, 1, 4}, {+1, +1, -1});
  c.attach_cell(2, {4, 2, 3}, {+1, -1, -1});
  return c;
}

}  // namespace

// Witnesses that orient() makes a stratum coherent: the disk references its
// shared 1-cell with +1 from both 2-cells, so exactly one must be flipped.
GRAPHOS_TEST(normalizes_an_inconsistently_oriented_disk) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  CHECK(!consistently_oriented(c));

  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  o.complex.validate();
  CHECK(consistently_oriented(o.complex));
  CHECK(graphos::d_squared_is_zero(o.complex));
  // deterministic: cell 0 seeds with +1, so cell 1 carries the flip
  CHECK(o.map.sign[2][0] == +1);
  CHECK(o.map.sign[2][1] == -1);
  // the lower strata are untouched
  for (std::size_t i = 0; i < o.map.sign[1].size(); ++i) CHECK(o.map.sign[1][i] == +1);
}

GRAPHOS_TEST(consistent_input_is_untouched) {
  const graphos::Complex c = graphos_test::make_fan();
  CHECK(consistently_oriented(c));
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  for (std::size_t i = 0; i < o.map.sign[2].size(); ++i) CHECK(o.map.sign[2][i] == +1);
}

GRAPHOS_TEST(torus_is_orientable) {
  graphos::Complex circle(1);
  circle.attach_vertices(3);
  circle.attach_cell(1, {0, 1}, {-1, +1});
  circle.attach_cell(1, {1, 2}, {-1, +1});
  circle.attach_cell(1, {2, 0}, {-1, +1});
  const auto torus = graphos::product(circle, circle);
  const auto o = graphos::orient(torus.complex);
  CHECK(o.orientable);
  CHECK(consistently_oriented(o.complex));
}

GRAPHOS_TEST(cylinder_is_orientable_moebius_is_not) {
  const graphos::Complex square = make_square();

  const auto cyl_ids = graphos::lift_identifications(square, {{0, 1, +1}, {3, 2, +1}});
  const auto cylinder = graphos::quotient(square, cyl_ids);
  CHECK(graphos::orient(cylinder.complex).orientable);

  const auto mob_ids = graphos::lift_identifications(square, {{0, 2, +1}, {3, 1, +1}});
  const auto moebius = graphos::quotient(square, mob_ids);
  CHECK(!graphos::orient(moebius.complex).orientable);
}

GRAPHOS_TEST(zero_dimensional_complex_is_trivially_orientable) {
  graphos::Complex c(0);
  c.attach_vertices(5);
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  CHECK(o.complex.count(0) == 5);
}

// ---- the orientation classes, and the sign they leave free ----------------

// Witnesses that the free signs number one per class: two disjoint pieces give
// two classes, each internally coherent after orient(), with nothing relating
// them -- ∂ cannot express which way is out. Reporting the classes is what lets
// a caller with coordinates settle exactly that many signs.
GRAPHOS_TEST(disjoint_pieces_are_separate_orientation_classes) {
  const graphos::Complex one = graphos_test::make_two_triangle_disk();
  const graphos::Complex c = graphos::disjoint_union(one, one).complex;
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  CHECK(consistently_oriented(o.complex));
  CHECK(o.classes == 2);
  CHECK(o.stratum == 2);
  // the two 2-cells of a piece share a class; the pieces do not
  CHECK(o.class_of[0] == o.class_of[1]);
  CHECK(o.class_of[2] == o.class_of[3]);
  CHECK(o.class_of[0] != o.class_of[2]);
}

GRAPHOS_TEST(a_connected_mesh_is_one_class) {
  const auto o = graphos::orient(graphos_test::make_fan());
  CHECK(o.classes == 1);
  for (const Index cls : o.class_of) CHECK(cls == 0);
}

// ---- flipping a class ------------------------------------------------------

// Witnesses that reversing a whole class preserves coherence: every cell
// reverses together, so no interior facet changes its induced signs. This is
// what a caller applies once its determinant says a class is inside-out.
GRAPHOS_TEST(flipping_a_whole_class_preserves_coherence) {
  const auto o = graphos::orient(graphos_test::make_fan());
  std::vector<Index> all;
  for (Index e = 0; e < o.complex.count(2); ++e) all.push_back(e);
  const graphos::Complex f = graphos::flip_cells(o.complex, 2, all);
  f.validate();
  CHECK(consistently_oriented(f));
  CHECK(graphos::d_squared_is_zero(f));
  // an involution: flipping twice is the identity
  const graphos::Complex back = graphos::flip_cells(f, 2, all);
  for (std::size_t m = 0; m < back.boundary(2).signs.size(); ++m) {
    CHECK(back.boundary(2).signs[m] == o.complex.boundary(2).signs[m]);
  }
}

// Witnesses that the two properties are independent: reversing one cell breaks
// coherence and preserves ∂∘∂ = 0. Coherence is agreement between neighbours;
// ∂∘∂ = 0 is the complex being a chain complex, and flip_cells must preserve
// the second whatever it does to the first.
GRAPHOS_TEST(flipping_one_cell_breaks_coherence_but_not_the_complex) {
  const auto o = graphos::orient(graphos_test::make_fan());
  const graphos::Complex f = graphos::flip_cells(o.complex, 2, {0});
  f.validate();
  CHECK(graphos::d_squared_is_zero(f));
  CHECK(!consistently_oriented(f));
}

// Witnesses that reversing a non-top cell negates the column above it as well:
// the square is two 2-cells over five 1-cells, so flipping a 1-cell must touch
// both ∂₁ and ∂₂ for ∂∘∂ = 0 to survive.
GRAPHOS_TEST(flipping_a_facet_negates_the_column_above_it) {
  const graphos::Complex c = make_square();
  const graphos::Complex f = graphos::flip_cells(c, 1, {4});  // the shared diagonal
  f.validate();
  CHECK(graphos::d_squared_is_zero(f));
  // its own row of ∂₁ reversed ...
  const graphos::BoundaryOperator& b1 = c.boundary(1);
  const graphos::BoundaryOperator& g1 = f.boundary(1);
  for (Index m = b1.offsets[4]; m < b1.offsets[4 + 1]; ++m) CHECK(g1.signs[m] == -b1.signs[m]);
  // ... and every [σ : τ] referencing it reversed with it
  const graphos::BoundaryOperator& b2 = c.boundary(2);
  const graphos::BoundaryOperator& g2 = f.boundary(2);
  for (std::size_t m = 0; m < b2.indices.size(); ++m) {
    CHECK(g2.signs[m] == (b2.indices[m] == 4 ? -b2.signs[m] : b2.signs[m]));
  }
}

// ---- which stratum ---------------------------------------------------------

// Witnesses that a non-maximal stratum has nothing to orient: every 1-cell of
// the square is a facet of a 2-cell, so its orientation is the shared reference
// its cofaces agree on rather than a free choice, and no class is reported.
GRAPHOS_TEST(non_maximal_cells_are_not_reoriented) {
  const graphos::Complex c = make_square();
  const auto o = graphos::orient(c, 1);
  CHECK(o.orientable);
  CHECK(o.stratum == 1);
  CHECK(o.classes == 0);
  for (std::size_t m = 0; m < c.boundary(1).signs.size(); ++m) {
    CHECK(o.complex.boundary(1).signs[m] == c.boundary(1).signs[m]);
  }
  CHECK(graphos::d_squared_is_zero(o.complex));
}

// Witnesses the fracture case: a maximal lower stratum orients like any other.
// Two 1-cells meeting at a vertex with no coface form a path whose coherence
// is the same statement one dimension down, and orienting the top stratum of
// the same complex finds nothing to do.
GRAPHOS_TEST(a_maximal_lower_stratum_is_oriented) {
  graphos::Complex c(1);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});  // 0 -> 1
  c.attach_cell(1, {2, 1}, {-1, +1});  // 2 -> 1, so vertex 1 is induced twice with +1
  const auto o = graphos::orient(c, 1);
  CHECK(o.orientable);
  CHECK(o.classes == 1);
  CHECK(o.class_of[0] == 0);
  CHECK(o.class_of[1] == 0);
  CHECK(o.map.sign[1][0] == +1);
  CHECK(o.map.sign[1][1] == -1);  // the second edge is turned to make the path run through
  o.complex.validate();
  CHECK(graphos::d_squared_is_zero(o.complex));
}

// the default argument is the top stratum
GRAPHOS_TEST(orient_defaults_to_the_top_stratum) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  const auto a = graphos::orient(c);
  const auto b = graphos::orient(c, c.dim());
  CHECK(a.stratum == b.stratum);
  for (std::size_t m = 0; m < a.complex.boundary(2).signs.size(); ++m) {
    CHECK(a.complex.boundary(2).signs[m] == b.complex.boundary(2).signs[m]);
  }
}

GRAPHOS_TEST_MAIN()
