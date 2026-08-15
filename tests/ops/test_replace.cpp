// Incremental imprinting with replace(), simulated as a 2D cut-cell
// workflow. ALL metric content lives in this test — the line inclination,
// the intersection points, the trapezoid areas, and the bad-cut rule
// (area < 1/3 of the quad) — and enters graphos only as markers, patches,
// and glue. graphos does the topology: excision + amalgamation, side
// classification, healing by agglomeration, and the invariant checks.
//
// Scenario: 3x3 unit-quad mesh, vertices (x, y) = (i, j), cut by the line
//   y = 1.2 + 0.2 x
// which crosses the middle row, splitting each of its quads into a lower
// and an upper trapezoid with areas 0.3/0.5/0.7 and 0.7/0.5/0.3. The two
// area-0.3 pieces are bad cuts and are healed into their same-side
// neighbors below/above.

#include <vector>

#include "fixtures.hpp"
#include "graphos/core/build.hpp"
#include "graphos/ops/agglomerate.hpp"
#include "graphos/ops/orient.hpp"
#include "graphos/ops/product.hpp"
#include "graphos/ops/replace.hpp"
#include "graphos/queries/components.hpp"
#include "graphos_test.hpp"

using graphos::Index;

namespace {

graphos::Complex path(int n) {
  std::vector<std::vector<Index>> segs;
  for (Index i = 0; i < n; ++i) segs.push_back({i, i + 1});
  return graphos::from_edges(n + 1, segs);
}

// product layout: vertex (i,j) = i*4+j; vertical edge at x=i, y∈[j,j+1] is
// i*3+j; horizontal edge y=j, x∈[i,i+1] is 12+i*4+j; quad (i,j) = i*3+j
graphos::Complex quad_mesh_3x3() {
  const graphos::Complex p = path(3);
  return graphos::product(p, p).complex;
}

}  // namespace

// Sanity: replace one face of the two-triangle disk with a 3-triangle fan
// over a new center vertex; the frontier edges are reused via lifted glue.
GRAPHOS_TEST(replace_face_with_fan) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::Marker region(disk);
  region.mark(2, 1);  // face B only; its edges survive as the frontier

  // patch: fan over B's triangle (vertices 0,1,3 -> patch 0,1,2; center 3)
  const graphos::Complex fan = graphos::from_simplices(2, 4, {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}});
  const auto r = graphos::replace(disk, region, fan, {{0, 0, +1}, {1, 1, +1}, {2, 3, +1}});
  r.complex.validate();
  CHECK(graphos::d_squared_is_zero(r.complex));
  CHECK(r.complex.count(0) == 5);  // 4 host + 1 new center
  CHECK(r.complex.count(1) == 8);  // 5 host + 3 spokes (outer 3 lifted onto host)
  CHECK(r.complex.count(2) == 4);  // A + 3 fan triangles
  CHECK(graphos::euler_characteristic(r.complex) == 1);
  CHECK(r.map.index[2][1] == graphos::invalid_index);  // B is gone
  CHECK(r.map.index[2][0] != graphos::invalid_index);  // A survives
  CHECK(r.patch_map.index[0][3] != graphos::invalid_index);
}

GRAPHOS_TEST(glue_into_excised_region_throws) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::Marker region(disk);
  region.mark(0, 3);  // B's apex: its star (edges + face B) is excised
  const graphos::Complex patch = graphos_test::make_triangle();
  CHECK_THROWS(graphos::replace(disk, region, patch, {{0, 3, +1}}));  // target gone
  CHECK_THROWS(graphos::replace(disk, region, patch, {{0, 99, +1}}));
}

// The cut-cell scenario, end to end.
GRAPHOS_TEST(inclined_cut_imprint_and_bad_cut_healing) {
  const graphos::Complex host = quad_mesh_3x3();

  // --- metric side: the line y = 1.2 + 0.2x crosses the middle row,
  // cutting the vertical edges at x = 0..3; no horizontal edge is crossed
  // --- topology side: the region is the middle-row quads and the crossed
  // vertical edges (their stars are exactly those quads)
  graphos::Marker region(host);
  for (Index i = 0; i < 3; ++i) region.mark(2, i * 3 + 1);  // quads (i,1)
  for (Index i = 0; i < 4; ++i) region.mark(1, i * 3 + 1);  // vertical edges at row 1

  // the imprinted patch for the swath y∈[1,2]: bottom vertices b_i (0..3),
  // top t_i (4..7), cut points p_i (8..11); three lower and three upper
  // trapezoids, wound CCW
  const graphos::Complex patch = graphos::from_polygons(12, {{0, 1, 9, 8},
                                                             {1, 2, 10, 9},
                                                             {2, 3, 11, 10},  // lower pieces
                                                             {8, 9, 5, 4},
                                                             {9, 10, 6, 5},
                                                             {10, 11, 7, 6}});  // upper pieces

  // glue: patch bottom/top rows onto the surviving host lattice rows
  std::vector<graphos::Identification> glue;
  for (Index i = 0; i <= 3; ++i) {
    glue.push_back({i, i * 4 + 1, +1});      // b_i -> (i, 1)
    glue.push_back({4 + i, i * 4 + 2, +1});  // t_i -> (i, 2)
  }

  const auto imprint = graphos::replace(host, region, patch, glue);
  imprint.complex.validate();
  CHECK(graphos::d_squared_is_zero(imprint.complex));
  CHECK(graphos::euler_characteristic(imprint.complex) == 1);
  CHECK(imprint.complex.count(0) == 20);  // 16 host + 4 cut points
  CHECK(imprint.complex.count(1) == 31);  // 24 - 4 split + 8 halves + 3 cut
  CHECK(imprint.complex.count(2) == 12);  // 6 host quads + 6 pieces
  for (Index i = 0; i < 3; ++i) CHECK(imprint.map.index[2][i * 3 + 1] == graphos::invalid_index);
  for (Index i = 0; i < 4; ++i) CHECK(imprint.map.index[1][i * 3 + 1] == graphos::invalid_index);

  // the fault topology: the three cut edges (both endpoints are cut
  // points), located combinatorially in the patch and transported
  graphos::Marker interface(imprint.complex);
  const graphos::BoundaryOperator& pe = patch.boundary(1);
  for (Index e = 0; e < patch.count(1); ++e) {
    const bool cut_edge = pe.indices[pe.offsets[e]] >= 8 && pe.indices[pe.offsets[e] + 1] >= 8;
    if (cut_edge) interface.mark(1, imprint.patch_map.index[1][e]);
  }
  CHECK(interface.marked_count(1) == 3);
  // conformity: every cut edge separates exactly two cells
  const graphos::CoboundaryOperator cob = graphos::coboundary(imprint.complex, 1);
  for (Index e = 0; e < imprint.complex.count(1); ++e) {
    if (interface.marked(1, e)) CHECK(cob.offsets[e + 1] - cob.offsets[e] == 2);
  }

  // --- metric side: piece areas from the trapezoid formula
  // lower: 0.3, 0.5, 0.7   upper: 0.7, 0.5, 0.3  -> bad = area < 1/3
  const double lower_area[3] = {0.3, 0.5, 0.7};
  const double upper_area[3] = {0.7, 0.5, 0.3};
  // result-cell areas, indexed by result id
  std::vector<double> area(12, 0.0);
  std::vector<Index> host_res(9, graphos::invalid_index);
  for (Index q = 0; q < 9; ++q) {
    host_res[q] = imprint.map.index[2][q];
    if (host_res[q] != graphos::invalid_index) area[host_res[q]] = 1.0;
  }
  std::vector<Index> lower_res(3), upper_res(3);
  for (Index i = 0; i < 3; ++i) {
    lower_res[i] = imprint.patch_map.index[2][i];
    upper_res[i] = imprint.patch_map.index[2][3 + i];
    area[lower_res[i]] = lower_area[i];
    area[upper_res[i]] = upper_area[i];
  }

  // sides: excluding the cut edges as connectors separates below from above
  const auto sides = graphos::connected_components(imprint.complex, 2, 1, interface);
  CHECK(sides.count == 2);
  // healing policy (metric): each bad piece merges into the host quad
  // across its unsplit horizontal edge — same side, verified topologically
  CHECK(sides.label[lower_res[0]] == sides.label[host_res[0]]);  // (0,0) below
  CHECK(sides.label[upper_res[2]] == sides.label[host_res[8]]);  // (2,2) above
  CHECK(sides.label[lower_res[0]] != sides.label[upper_res[2]]);

  std::vector<Index> labels(12);
  for (Index r = 0; r < 12; ++r) labels[r] = r;
  labels[lower_res[0]] = labels[host_res[0]];
  labels[upper_res[2]] = labels[host_res[8]];
  // compact label ids
  {
    std::vector<Index> remap(12, graphos::invalid_index);
    Index next = 0;
    for (Index r = 0; r < 12; ++r) {
      if (remap[labels[r]] == graphos::invalid_index) remap[labels[r]] = next++;
      labels[r] = remap[labels[r]];
    }
  }

  const auto oriented = graphos::orient(imprint.complex);
  CHECK(oriented.orientable);
  const auto healed = graphos::agglomerate(oriented.complex, labels);
  healed.complex.validate();
  CHECK(graphos::d_squared_is_zero(healed.complex));
  CHECK(graphos::euler_characteristic(healed.complex) == 1);
  CHECK(healed.complex.count(2) == 10);  // 12 cells - 2 merges

  // metric verdict: every healed cell now has area >= 1/3
  std::vector<double> healed_area(10, 0.0);
  for (Index r = 0; r < 12; ++r) healed_area[labels[r]] += area[r];
  for (const double a : healed_area) CHECK(a > 1.0 / 3.0 - 1e-12);

  // the merged cells are hexagon-shaped (4 + 4 - 2 shared); the fault
  // edges survive the healing untouched
  const graphos::BoundaryOperator& hb = healed.complex.boundary(2);
  const Index m1 = healed.map.index[2][lower_res[0]];
  const Index m2 = healed.map.index[2][upper_res[2]];
  CHECK(hb.offsets[m1 + 1] - hb.offsets[m1] == 6);
  CHECK(hb.offsets[m2 + 1] - hb.offsets[m2] == 6);
  CHECK(healed.map.index[2][lower_res[0]] == healed.map.index[2][host_res[0]]);
  for (Index e = 0; e < imprint.complex.count(1); ++e) {
    if (interface.marked(1, e)) CHECK(healed.map.index[1][e] != graphos::invalid_index);
  }
}

GRAPHOS_TEST_MAIN()
