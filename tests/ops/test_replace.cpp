// Witnesses that excision and amalgamation compose into imprinting: replace()
// followed by agglomerate() preserves ∂∘∂ = 0, conformity of the interface,
// and the chain maps that transport cochains across both.
//
// Every metric quantity — the line, the intersection points, the areas, the
// rule area < 1/3 — lives in the test and reaches graphos only as markers,
// patches and glue. The complex is a 3×3 quad mesh cut by y = 1.2 + 0.2x,
// which splits each middle-row 2-cell into trapezoids of area 0.3/0.5/0.7 and
// 0.7/0.5/0.3; the two smallest are healed into their same-side neighbours.

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

// product layout: vertex (i,j) = 4i+j; the 1-cell at x = i over y ∈ [j, j+1]
// is 3i+j; the one at y = j over x ∈ [i, i+1] is 12+4i+j; the 2-cell is 3i+j
graphos::Complex quad_mesh_3x3() {
  const graphos::Complex p = path(3);
  return graphos::product(p, p).complex;
}

}  // namespace

// Witnesses that the glue lifts: replacing one 2-cell of the disk by a fan of
// three over a new vertex reuses the frontier 1-cells through the lifted
// vertex correspondence, rather than duplicating them.
GRAPHOS_TEST(replace_face_with_fan) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::Marker region(disk);
  region.mark(2, 1);  // face B only; its edges survive as the frontier

  // patch: a fan over B's 2-cell (vertices 0,1,3 → 0,1,2; apex 3)
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

// The imprint, end to end.
GRAPHOS_TEST(inclined_cut_imprint_and_bad_cut_healing) {
  const graphos::Complex host = quad_mesh_3x3();

  // metric: y = 1.2 + 0.2x crosses the middle row, cutting the 1-cells at
  // x = 0..3 and no others
  // topology: the region is the middle-row 2-cells with the cut 1-cells,
  // whose stars are exactly those 2-cells
  graphos::Marker region(host);
  for (Index i = 0; i < 3; ++i) region.mark(2, i * 3 + 1);  // quads (i,1)
  for (Index i = 0; i < 4; ++i) region.mark(1, i * 3 + 1);  // vertical edges at row 1

  // the patch over y ∈ [1,2]: bottom vertices 0..3, top 4..7, cut points
  // 8..11; three lower and three upper 2-cells, coherently wound
  const graphos::Complex patch = graphos::from_polygons(12, {{0, 1, 9, 8},
                                                             {1, 2, 10, 9},
                                                             {2, 3, 11, 10},  // lower pieces
                                                             {8, 9, 5, 4},
                                                             {9, 10, 6, 5},
                                                             {10, 11, 7, 6}});  // upper pieces

  // glue: the patch's bottom and top rows onto the surviving host rows
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

  // the interface: the three 1-cells with both endpoints cut points, located
  // combinatorially in the patch and transported through the chain map
  graphos::Marker interface(imprint.complex);
  const graphos::BoundaryOperator& pe = patch.boundary(1);
  for (Index e = 0; e < patch.count(1); ++e) {
    const bool cut_edge = pe.indices[pe.offsets[e]] >= 8 && pe.indices[pe.offsets[e] + 1] >= 8;
    if (cut_edge) interface.mark(1, imprint.patch_map.index[1][e]);
  }
  CHECK(interface.marked_count(1) == 3);
  // conformity: each interface 1-cell has exactly two cofaces
  const graphos::CoboundaryOperator cob = graphos::coboundary(imprint.complex, 1);
  for (Index e = 0; e < imprint.complex.count(1); ++e) {
    if (interface.marked(1, e)) CHECK(cob.offsets[e + 1] - cob.offsets[e] == 2);
  }

  // metric: areas 0.3/0.5/0.7 below and 0.7/0.5/0.3 above; bad is < 1/3
  const double lower_area[3] = {0.3, 0.5, 0.7};
  const double upper_area[3] = {0.7, 0.5, 0.3};
  // areas indexed by result cell
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

  // sides: removing the interface as connectors separates below from above
  const auto sides = graphos::connected_components(imprint.complex, 2, 1, interface);
  CHECK(sides.count == 2);
  // healing (metric): each bad cell merges across its unsplit 1-cell into the
  // neighbour on the same side, which is checked topologically
  CHECK(sides.label[lower_res[0]] == sides.label[host_res[0]]);  // (0,0) below
  CHECK(sides.label[upper_res[2]] == sides.label[host_res[8]]);  // (2,2) above
  CHECK(sides.label[lower_res[0]] != sides.label[upper_res[2]]);

  std::vector<Index> labels(12);
  for (Index r = 0; r < 12; ++r) labels[r] = r;
  labels[lower_res[0]] = labels[host_res[0]];
  labels[upper_res[2]] = labels[host_res[8]];
  // compact the labels
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

  // metric: every healed cell has area ≥ 1/3
  std::vector<double> healed_area(10, 0.0);
  for (Index r = 0; r < 12; ++r) healed_area[labels[r]] += area[r];
  for (const double a : healed_area) CHECK(a > 1.0 / 3.0 - 1e-12);

  // a merged cell has six faces (4 + 4 − 2 cancelled) and the interface
  // survives the agglomeration
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
