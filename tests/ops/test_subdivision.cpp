#include "fixtures.hpp"
#include "graphos/core/build.hpp"
#include "graphos/core/coboundary.hpp"
#include "graphos/ops/orient.hpp"
#include "graphos/ops/subdivision.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

// the transfer law: for a consistently oriented parent, multiplying each
// top sd cell by its carrier_sign yields a consistent orientation of the
// subdivision — every interior sd facet is induced with opposite signs
bool orientation_transfers(const graphos::Complex& oriented) {
  const auto sd = graphos::barycentric_subdivision(oriented);
  const int n = sd.complex.dim();
  const graphos::CoboundaryOperator cob = graphos::coboundary(sd.complex, n - 1);
  for (Index f = 0; f < sd.complex.count(n - 1); ++f) {
    const Index lo = cob.offsets[static_cast<std::size_t>(f)];
    const Index hi = cob.offsets[static_cast<std::size_t>(f) + 1];
    if (hi - lo != 2) continue;
    int sum = 0;
    for (Index m = lo; m < hi; ++m) {
      sum += cob.signs[static_cast<std::size_t>(m)] *
             sd.carrier_sign[static_cast<std::size_t>(n)]
                            [static_cast<std::size_t>(cob.indices[static_cast<std::size_t>(m)])];
    }
    if (sum != 0) return false;
  }
  return true;
}

}  // namespace

GRAPHOS_TEST(segment_subdivides_into_two_edges) {
  const auto sd = graphos::barycentric_subdivision(graphos_test::make_segment());
  sd.complex.validate();
  CHECK(sd.complex.count(0) == 3);  // 2 vertices + 1 edge barycenter
  CHECK(sd.complex.count(1) == 2);
  CHECK(graphos::d_squared_is_zero(sd.complex));
  CHECK(graphos::euler_characteristic(sd.complex) == 1);
}

GRAPHOS_TEST(triangle_subdivides_into_six_triangles) {
  const auto sd = graphos::barycentric_subdivision(graphos_test::make_triangle());
  sd.complex.validate();
  CHECK(sd.complex.count(0) == 7);                // 3 + 3 + 1 barycenters
  CHECK(sd.complex.count(1) == 12);               // v<e: 6, v<f: 3, e<f: 3
  CHECK(sd.complex.count(2) == 6);                // full flags v<e<f
  CHECK(graphos::d_squared_is_zero(sd.complex));  // the (−1)^m chain signs
  CHECK(graphos::euler_characteristic(sd.complex) == 1);

  // sd-vertex layout: vertices, then edges, then faces
  CHECK(sd.vertex_offset[0] == 0);
  CHECK(sd.vertex_offset[1] == 3);
  CHECK(sd.vertex_offset[2] == 6);

  // every 2-simplex is carried by the original face
  for (std::size_t i = 0; i < sd.carrier_dim[2].size(); ++i) {
    CHECK(sd.carrier_dim[2][i] == 2);
    CHECK(sd.carrier_index[2][i] == 0);
  }
}

GRAPHOS_TEST(subdivision_preserves_the_disk) {
  const auto sd = graphos::barycentric_subdivision(graphos_test::make_two_triangle_disk());
  sd.complex.validate();
  CHECK(sd.complex.count(0) == 11);  // 4 + 5 + 2
  CHECK(sd.complex.count(1) == 22);  // v<e: 10, v<f: 6, e<f: 6
  CHECK(sd.complex.count(2) == 12);  // 6 per triangle
  CHECK(graphos::d_squared_is_zero(sd.complex));
  CHECK(graphos::euler_characteristic(sd.complex) == 1);
}

// Mixed-dimensional input: a triangle with a dangling maximal edge still
// subdivides — chains through the bare edge stop at dimension 1.
GRAPHOS_TEST(signed_carrier_is_the_incidence_product) {
  const auto sd = graphos::barycentric_subdivision(graphos_test::make_triangle());
  // sd_0: original vertices carry themselves (+1); barycenters of
  // positive-dimensional cells lie outside the chain-map image (0)
  for (Index i = 0; i < 3; ++i) CHECK(sd.carrier_sign[0][i] == +1);
  for (Index i = 3; i < 7; ++i) CHECK(sd.carrier_sign[0][i] == 0);
  // sd_1: full flags (v < e) carry the incidence sign; dimension jumps
  // (v < F) and chains starting above dimension 0 (e < F) are 0
  CHECK(sd.carrier_sign[1][0] == -1);  // (v0 < e0): [v0 : e0] = −1
  CHECK(sd.carrier_sign[1][1] == +1);  // (v0 < e2): [v0 : e2] = +1
  CHECK(sd.carrier_sign[1][2] == 0);   // (v0 < F): jump
  Index zeros = 0, units = 0;
  for (const graphos::Sign s : sd.carrier_sign[1]) {
    if (s == 0) ++zeros;
    if (s == 1 || s == -1) ++units;
  }
  CHECK(zeros == 6);  // 3 jumps (v<F) + 3 chains (e<F)
  CHECK(units == 6);  // the full flags (v<e)
  // every top sd cell is a full flag: coefficient ±1
  for (const graphos::Sign s : sd.carrier_sign[2]) CHECK(s == 1 || s == -1);
  CHECK(sd.carrier_sign[2][0] == -1);  // (v0<e0<F) = (−1)·(+1)
}

GRAPHOS_TEST(orientation_transfers_through_subdivision) {
  // the fan is consistently oriented as built
  CHECK(orientation_transfers(graphos_test::make_fan()));
  // the disk is not — orient it first; the check then covers the
  // cross-parent seam through the shared edge
  const auto o2 = graphos::orient(graphos_test::make_two_triangle_disk());
  CHECK(o2.orientable);
  CHECK(orientation_transfers(o2.complex));
  // and in 3D: two tetrahedra sharing a face
  const auto o3 = graphos::orient(graphos::from_simplices(3, 5, {{0, 1, 2, 3}, {1, 3, 2, 4}}));
  CHECK(o3.orientable);
  CHECK(orientation_transfers(o3.complex));
}

GRAPHOS_TEST(subdivides_mixed_dimensional_complexes) {
  graphos::Complex c(2);
  c.attach_vertices(4);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {0, 3}, {-1, +1});  // whisker: maximal 1-cell
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  const auto sd = graphos::barycentric_subdivision(c);
  sd.complex.validate();
  CHECK(graphos::d_squared_is_zero(sd.complex));
  // χ preserved: disk with a whisker is contractible
  CHECK(graphos::euler_characteristic(sd.complex) == 1);
  // the whisker contributes exactly two sd edges
  Index whisker_faces = 0;
  for (std::size_t i = 0; i < sd.carrier_dim[1].size(); ++i) {
    if (sd.carrier_dim[1][i] == 1 && sd.carrier_index[1][i] == 3) ++whisker_faces;
  }
  CHECK(whisker_faces == 2);
}

GRAPHOS_TEST_MAIN()
