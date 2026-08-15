// Dimensional-descent coarsening. graphos is metric-free, so it cannot
// know where a geometric corner is — the TESTS know (they built the
// structured grids), and they pass that metric knowledge in as protected
// markers, then assert that graphos's purely combinatorial descent
// reproduces the exact coarse tensor meshes.

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

// METRIC knowledge, supplied as input: the coarse lattice corners of the
// structured grid (vertices with all-even coordinates)
graphos::Marker quad_corners(const graphos::Complex& c, int n) {
  graphos::Marker m(c);
  for (int i = 0; i <= n; i += 2)
    for (int j = 0; j <= n; j += 2) m.mark(0, static_cast<Index>(i * (n + 1) + j));
  return m;
}

// the full coarse lattice FRAME: corner vertices (barriers for the edge
// chains) and lattice-line edges (barriers for the face patches), decoded
// from the product layout v = (i*(n+1)+j)*(n+1)+k
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
  // an edge lies on a coarse lattice line iff both endpoints share the
  // same two even coordinates
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

// The full inverse of tensor refinement: with corners protected, the
// descent reproduces the coarse quad mesh EXACTLY — counts and arities.
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

  // and once more down to the single cell (on the independently built
  // middle mesh, whose layout the metric side knows)
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

// The chain map records the descent: sent to zero interiors go to zero,
// merged pairs share a coarse id.
GRAPHOS_TEST(chain_map_records_the_descent) {
  const graphos::Complex fine = quad_mesh(2);
  const auto co = graphos::coarsen(fine, {0, 0, 0, 0}, quad_corners(fine, 2));
  // center vertex (1,1) = index 4 is sent to zero
  CHECK(co.map.index[0][4] == graphos::invalid_index);
  // the four corners survive
  for (const Index v : {0, 2, 6, 8}) CHECK(co.map.index[0][v] != graphos::invalid_index);
  // every surviving edge pair shares a coarse edge; coefficients are ±1
  for (std::size_t e = 0; e < co.map.index[1].size(); ++e) {
    if (co.map.index[1][e] == graphos::invalid_index) continue;
    CHECK(co.map.sign[1][e] == 1 || co.map.sign[1][e] == -1);
  }
  CHECK(co.complex.count(1) == 4);
}

// Without protection the merge is maximal — and the domain-boundary cycle
// patch is REJECTED (a closed coarse edge would have empty boundary), so
// the fallback keeps the fine boundary cells.
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

// Combinatorial freedom, demonstrated: protecting only two opposite
// corners yields a legitimate two-edge "bigon" cell — topology alone
// cannot forbid it, which is exactly why corners are caller input.
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

// Mixed-dimensional passthrough: maximal lower-dimensional cells (a
// detached segment) survive the descent untouched.
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
