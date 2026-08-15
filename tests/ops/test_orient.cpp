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

// consistency: every interior facet (two distinct cofaces) is induced with
// opposite signs by its two sides
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

// The disk fixture references the shared edge with +1 from BOTH faces —
// inconsistent as imported. orient() must flip exactly one face.
GRAPHOS_TEST(normalizes_an_inconsistently_oriented_disk) {
  const graphos::Complex c = graphos_test::make_two_triangle_disk();
  CHECK(!consistently_oriented(c));

  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  o.complex.validate();
  CHECK(consistently_oriented(o.complex));
  CHECK(graphos::d_squared_is_zero(o.complex));
  // deterministic: cell 0 seeds with +1, so cell 1 takes the flip
  CHECK(o.map.sign[2][0] == +1);
  CHECK(o.map.sign[2][1] == -1);
  // lower strata untouched
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

// TWO DISJOINT PIECES ARE TWO CLASSES, and that is the number of signs
// coherence cannot fix. Each is internally consistent after orient(); nothing
// relates one to the other, because nothing in a chain complex could -- "which
// way is out" is not a statement the boundary operator can make. Reporting the
// classes is what lets a caller with coordinates settle exactly that many.
GRAPHOS_TEST(disjoint_pieces_are_separate_orientation_classes) {
  const graphos::Complex one = graphos_test::make_two_triangle_disk();
  const graphos::Complex c = graphos::disjoint_union(one, one).complex;
  const auto o = graphos::orient(c);
  CHECK(o.orientable);
  CHECK(consistently_oriented(o.complex));
  CHECK(o.classes == 2);
  CHECK(o.stratum == 2);
  // the two faces of each piece share a class, and the pieces do not
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

// TURNING A WHOLE CLASS IS STILL A COHERENT MESH. This is the operation a
// caller reaches for once its determinant says the class came out inside-out:
// every cell of it reverses together, so no interior facet changes its verdict.
GRAPHOS_TEST(flipping_a_whole_class_preserves_coherence) {
  const auto o = graphos::orient(graphos_test::make_fan());
  std::vector<Index> all;
  for (Index e = 0; e < o.complex.count(2); ++e) all.push_back(e);
  const graphos::Complex f = graphos::flip_cells(o.complex, 2, all);
  f.validate();
  CHECK(consistently_oriented(f));
  CHECK(graphos::d_squared_is_zero(f));
  // and it is an involution: turning twice is the identity
  const graphos::Complex back = graphos::flip_cells(f, 2, all);
  for (std::size_t m = 0; m < back.boundary(2).signs.size(); ++m) {
    CHECK(back.boundary(2).signs[m] == o.complex.boundary(2).signs[m]);
  }
}

// TURNING ONE CELL BREAKS COHERENCE AND KEEPS d.d = 0. The two are independent
// claims and the distinction is the whole subject: coherence is about
// neighbours agreeing, d.d = 0 is about the complex being a complex, and an
// operation that reverses a cell must respect the second whatever it does to
// the first.
GRAPHOS_TEST(flipping_one_cell_breaks_coherence_but_not_the_complex) {
  const auto o = graphos::orient(graphos_test::make_fan());
  const graphos::Complex f = graphos::flip_cells(o.complex, 2, {0});
  f.validate();
  CHECK(graphos::d_squared_is_zero(f));
  CHECK(!consistently_oriented(f));
}

// A NON-TOP CELL CARRIES A COLUMN ABOVE IT, and reversing it has to negate
// that column too or the complex stops being one. The square fixture is two
// faces over five edges, so flipping an edge touches d_1 and d_2 both.
GRAPHOS_TEST(flipping_a_facet_negates_the_column_above_it) {
  const graphos::Complex c = make_square();
  const graphos::Complex f = graphos::flip_cells(c, 1, {4});  // the shared diagonal
  f.validate();
  CHECK(graphos::d_squared_is_zero(f));
  // the edge's own row reversed ...
  const graphos::BoundaryOperator& b1 = c.boundary(1);
  const graphos::BoundaryOperator& g1 = f.boundary(1);
  for (Index m = b1.offsets[4]; m < b1.offsets[4 + 1]; ++m) CHECK(g1.signs[m] == -b1.signs[m]);
  // ... and every reference to it from a face reversed with it
  const graphos::BoundaryOperator& b2 = c.boundary(2);
  const graphos::BoundaryOperator& g2 = f.boundary(2);
  for (std::size_t m = 0; m < b2.indices.size(); ++m) {
    CHECK(g2.signs[m] == (b2.indices[m] == 4 ? -b2.signs[m] : b2.signs[m]));
  }
}

// ---- which stratum ---------------------------------------------------------

// ORIENTING A LOWER STRATUM LEAVES THE FACES OF A MESH ALONE. The 1-cells of
// the square are all facets of a face; their orientation is the shared
// reference the two sides agree on, not a free choice, so there is nothing to
// coherently orient and no class to report.
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

// BUT A MAXIMAL LOWER STRATUM IS ORIENTED LIKE ANY OTHER -- the fracture case.
// Two 1-cells meeting at a vertex, with no face above them, form a path whose
// coherence is the same statement one dimension down; orient(c, 1) fixes it,
// and orienting the top stratum of the same complex reports nothing to do.
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

// the default argument still means the top stratum
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
