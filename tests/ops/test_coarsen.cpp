// Witnesses that dimensional descent inverts tensor refinement. Being
// metric-free, graphos cannot locate a geometric corner; the tests know it,
// having built the grids, and supply it as protected markers. What is asserted
// is that the purely combinatorial descent then reproduces the coarse tensor
// complexes exactly.

#include <array>
#include <vector>

#include "fixtures.hpp"
#include "graphos/core/build.hpp"
#include "graphos/ops/coarsen.hpp"
#include "graphos/ops/product.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

graphos::Complex path(int n) {
  std::vector<std::vector<Index>> segs;
  for (Index i = 0; i < n; ++i) segs.push_back({i, i + 1});
  return graphos::from_edges(n + 1, segs);
}

graphos::Complex quad_mesh(int n) {
  const graphos::Complex p = path(n);
  return graphos::product(p, p).complex;
}

graphos::Complex hex_mesh(int n) {
  const graphos::Complex p = path(n);
  return graphos::product(quad_mesh(n), p).complex;
}

std::vector<Index> quad_parents(int n) {
  std::vector<Index> labels(static_cast<std::size_t>(n) * n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      labels[static_cast<std::size_t>(i * n + j)] = static_cast<Index>((i / 2) * (n / 2) + (j / 2));
  return labels;
}

std::vector<Index> hex_parents(int n) {
  std::vector<Index> labels(static_cast<std::size_t>(n) * n * n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        labels[static_cast<std::size_t>((i * n + j) * n + k)] =
            static_cast<Index>(((i / 2) * (n / 2) + (j / 2)) * (n / 2) + (k / 2));
  return labels;
}

// metric input: the coarse lattice corners, the vertices with all-even
// coordinates
graphos::Marker quad_corners(const graphos::Complex& c, int n) {
  graphos::Marker m(c);
  for (int i = 0; i <= n; i += 2)
    for (int j = 0; j <= n; j += 2) m.mark(0, static_cast<Index>(i * (n + 1) + j));
  return m;
}

// the coarse lattice frame: corner vertices, barriers for the 1-cell chains,
// and lattice-line 1-cells, barriers for the 2-cell patches, decoded from the
// product layout v = (i(n+1)+j)(n+1)+k
graphos::Marker hex_frame(const graphos::Complex& c, int n) {
  graphos::Marker m(c);
  const auto coords = [n](Index v) {
    const int s = n + 1;
    return std::array<int, 3>{static_cast<int>(v) / (s * s), (static_cast<int>(v) / s) % s,
                              static_cast<int>(v) % s};
  };
  for (int i = 0; i <= n; i += 2)
    for (int j = 0; j <= n; j += 2)
      for (int k = 0; k <= n; k += 2)
        m.mark(0, static_cast<Index>((i * (n + 1) + j) * (n + 1) + k));
  // a 1-cell lies on a lattice line iff its endpoints share the same two even
  // coordinates
  const graphos::BoundaryOperator& e = c.boundary(1);
  for (Index ed = 0; ed < c.count(1); ++ed) {
    const auto a = coords(e.indices[e.offsets[ed]]);
    const auto b = coords(e.indices[e.offsets[ed] + 1]);
    const bool on_x_line = a[1] % 2 == 0 && a[2] % 2 == 0 && b[1] % 2 == 0 && b[2] % 2 == 0;
    const bool on_y_line = a[0] % 2 == 0 && a[2] % 2 == 0 && b[0] % 2 == 0 && b[2] % 2 == 0;
    const bool on_z_line = a[0] % 2 == 0 && a[1] % 2 == 0 && b[0] % 2 == 0 && b[1] % 2 == 0;
    if (on_x_line || on_y_line || on_z_line) m.mark(1, ed);
  }
  return m;
}

void check_arities(const graphos::Complex& c, int k, Index expected) {
  const graphos::BoundaryOperator& b = c.boundary(k);
  for (Index e = 0; e < c.count(k); ++e) CHECK(b.offsets[e + 1] - b.offsets[e] == expected);
}

}  // namespace

// Witnesses the exact inverse of tensor refinement: with the corners
// protected, the descent reproduces the coarse complex in every stratum count
// and boundary arity.
GRAPHOS_TEST(quad_hierarchy_recovers_true_quads) {
  const graphos::Complex fine = quad_mesh(4);
  const graphos::Complex middle = quad_mesh(2);

  const auto co = graphos::coarsen(fine, quad_parents(4), quad_corners(fine, 4));
  co.complex.validate();
  CHECK(graphos::d_squared_is_zero(co.complex));
  CHECK(graphos::euler_characteristic(co.complex) == 1);
  for (int k = 0; k <= 2; ++k) CHECK(co.complex.count(k) == middle.count(k));  // 9/12/4
  check_arities(co.complex, 2, 4);  // every coarse cell is a TRUE quad
  check_arities(co.complex, 1, 2);

  // and again down to the single cell, on the independently built middle
  // complex whose layout the metric side knows
  const auto co2 = graphos::coarsen(middle, {0, 0, 0, 0}, quad_corners(middle, 2));
  co2.complex.validate();
  CHECK(graphos::d_squared_is_zero(co2.complex));
  CHECK(co2.complex.count(0) == 4);
  CHECK(co2.complex.count(1) == 4);
  CHECK(co2.complex.count(2) == 1);
  check_arities(co2.complex, 2, 4);
}

GRAPHOS_TEST(hex_hierarchy_recovers_true_hexes) {
  const graphos::Complex fine = hex_mesh(4);
  const graphos::Complex middle = hex_mesh(2);

  const auto co = graphos::coarsen(fine, hex_parents(4), hex_frame(fine, 4));
  co.complex.validate();
  CHECK(graphos::d_squared_is_zero(co.complex));
  CHECK(graphos::euler_characteristic(co.complex) == 1);
  for (int k = 0; k <= 3; ++k) CHECK(co.complex.count(k) == middle.count(k));  // 27/54/36/8
  check_arities(co.complex, 3, 6);  // every coarse cell is a TRUE hex
  check_arities(co.complex, 2, 4);
  check_arities(co.complex, 1, 2);

  const auto co2 = graphos::coarsen(middle, std::vector<Index>(8, 0), hex_frame(middle, 2));
  co2.complex.validate();
  CHECK(graphos::d_squared_is_zero(co2.complex));
  CHECK(co2.complex.count(0) == 8);
  CHECK(co2.complex.count(1) == 12);
  CHECK(co2.complex.count(2) == 6);
  CHECK(co2.complex.count(3) == 1);
  check_arities(co2.complex, 3, 6);
  check_arities(co2.complex, 2, 4);
}

// Witnesses that the chain map records the descent: cells interior to an
// aggregate are sent to 0 and merged cells share a coarse index.
GRAPHOS_TEST(chain_map_records_the_descent) {
  const graphos::Complex fine = quad_mesh(2);
  const auto co = graphos::coarsen(fine, {0, 0, 0, 0}, quad_corners(fine, 2));
  // the interior vertex (1,1) = 4 is sent to 0
  CHECK(co.map.index[0][4] == graphos::invalid_index);
  // the four corners survive
  for (const Index v : {0, 2, 6, 8}) CHECK(co.map.index[0][v] != graphos::invalid_index);
  // each surviving 1-cell pair shares a coarse cell, with coefficients ±1
  for (std::size_t e = 0; e < co.map.index[1].size(); ++e) {
    if (co.map.index[1][e] == graphos::invalid_index) continue;
    CHECK(co.map.sign[1][e] == 1 || co.map.sign[1][e] == -1);
  }
  CHECK(co.complex.count(1) == 4);
}

// Witnesses the maximal merge and its guard: unprotected, the boundary cycle
// is rejected — a closed coarse cell would have empty boundary — so the
// fallback keeps the fine cells.
GRAPHOS_TEST(unprotected_closed_cycle_falls_back_to_singletons) {
  const graphos::Complex fine = quad_mesh(2);
  const auto co = graphos::coarsen(fine, {0, 0, 0, 0});
  co.complex.validate();
  CHECK(graphos::d_squared_is_zero(co.complex));
  CHECK(graphos::euler_characteristic(co.complex) == 1);
  CHECK(co.complex.count(0) == 8);  // boundary ring survives, center gone
  CHECK(co.complex.count(1) == 8);
  CHECK(co.complex.count(2) == 1);
}

// Witnesses that topology alone cannot forbid a bigon: protecting two opposite
// corners yields a legitimate coarse cell with two faces, which is why the
// corners must be caller input.
GRAPHOS_TEST(two_protected_corners_make_a_bigon) {
  const graphos::Complex fine = quad_mesh(2);
  graphos::Marker two(fine);
  two.mark(0, 0).mark(0, 8);
  const auto co = graphos::coarsen(fine, {0, 0, 0, 0}, two);
  co.complex.validate();
  CHECK(graphos::d_squared_is_zero(co.complex));
  CHECK(graphos::euler_characteristic(co.complex) == 1);
  CHECK(co.complex.count(0) == 2);
  CHECK(co.complex.count(1) == 2);
  CHECK(co.complex.count(2) == 1);
}

// Witnesses mixed-dimensional passthrough: a maximal lower-dimensional cell
// survives the descent untouched.
GRAPHOS_TEST(detached_cells_pass_through) {
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {3, 4}, {-1, +1});  // detached segment
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});

  const auto co = graphos::coarsen(c, {0});
  co.complex.validate();
  CHECK(co.complex.count(0) == 5);
  CHECK(co.complex.count(1) == 4);  // triangle's cycle patch rejected + segment
  CHECK(co.complex.count(2) == 1);
  CHECK(co.map.index[1][3] != graphos::invalid_index);
}

GRAPHOS_TEST(rejects_bad_labels) {
  const graphos::Complex fine = quad_mesh(2);
  CHECK_THROWS(graphos::coarsen(fine, {0, 0, 0}));      // wrong size
  CHECK_THROWS(graphos::coarsen(fine, {0, 0, 2, 2}));   // id 1 unused
  CHECK_THROWS(graphos::coarsen(fine, {0, 0, -1, 1}));  // negative
}

GRAPHOS_TEST_MAIN()
