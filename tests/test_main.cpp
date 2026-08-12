#include <cstdio>

#include "graphos/graphos.hpp"

namespace {

int failures = 0;

#define CHECK(cond)                                                     \
  do {                                                                  \
    if (!(cond)) {                                                      \
      ++failures;                                                       \
      std::printf("FAILED %s:%d: %s\n", __FILE__, __LINE__, #cond);     \
    }                                                                   \
  } while (0)

// Oriented triangle: edges e0=[0,1], e1=[1,2], e2=[2,0]; the 2-cell is
// attached along the cycle e0 + e1 + e2.
graphos::Complex make_triangle() {
  graphos::Complex c(2);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  return c;
}

// Oriented segment: two vertices, one edge.
graphos::Complex make_segment() {
  graphos::Complex c(1);
  c.attach_vertices(2);
  c.attach_cell(1, {0, 1}, {-1, +1});
  return c;
}

void test_triangle_invariants() {
  const graphos::Complex c = make_triangle();
  c.validate();
  CHECK(c.count(0) == 3);
  CHECK(c.count(1) == 3);
  CHECK(c.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);  // disk
}

void test_disjoint_union_mixed_dimensional() {
  const graphos::Complex a = make_triangle();
  const graphos::Complex b = make_segment();
  const auto du = graphos::disjoint_union(a, b);
  du.complex.validate();
  CHECK(du.complex.dim() == 2);
  CHECK(du.complex.count(0) == 5);
  CHECK(du.complex.count(1) == 4);
  CHECK(du.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(du.complex));
  CHECK(graphos::euler_characteristic(du.complex) == 2);  // disk + segment
  // B's lone edge lands after A's edges, its vertices after A's vertices
  CHECK(du.b_map.index[0][0] == 3);
  CHECK(du.b_map.index[1][0] == 3);
  const graphos::BoundaryOperator& e = du.complex.boundary(1);
  CHECK(e.indices[e.offsets[3]] == 3);
  CHECK(e.indices[e.offsets[3] + 1] == 4);
}

// Pushout of two triangles along a shared edge, with B attached through
// REVERSED vertex identifications (B0'~A1, B1'~A0) so the shared edge is
// discovered with relative orientation -1. The orientation stress test.
void test_pushout_with_orientation_flip() {
  const graphos::Complex a = make_triangle();
  const graphos::Complex b = make_triangle();
  const std::vector<graphos::Identification> vertex_ids = {{0, 1, +1}, {1, 0, +1}};
  const auto po = graphos::pushout(a, b, vertex_ids, /*deduplicate=*/true);
  po.complex.validate();

  CHECK(po.complex.count(0) == 4);
  CHECK(po.complex.count(1) == 5);
  CHECK(po.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(po.complex));
  CHECK(graphos::euler_characteristic(po.complex) == 1);  // still a disk

  // every cell of A survives with orientation intact
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < po.a_map.index[k].size(); ++i) {
      CHECK(po.a_map.index[k][i] != graphos::invalid_index);
      CHECK(po.a_map.sign[k][i] == 1);
    }
  }
  // B's e0 identified with A's e0 through a flip; B's other cells survive
  CHECK(po.b_map.index[1][0] == po.a_map.index[1][0]);
  CHECK(po.b_map.sign[1][0] == -1);
  CHECK(po.b_map.index[0][0] == po.a_map.index[0][1]);
  CHECK(po.b_map.index[0][1] == po.a_map.index[0][0]);
  CHECK(po.b_map.index[0][2] != graphos::invalid_index);
  CHECK(po.b_map.index[2][0] != graphos::invalid_index);
  CHECK(po.b_map.index[2][0] != po.a_map.index[2][0]);

  // B's 2-cell must reference the shared edge with the flipped coefficient
  const graphos::Index shared_edge = po.a_map.index[1][0];
  const graphos::Index b_face = po.b_map.index[2][0];
  const graphos::BoundaryOperator& f = po.complex.boundary(2);
  bool found = false;
  for (graphos::Index m = f.offsets[b_face]; m < f.offsets[b_face + 1]; ++m) {
    if (f.indices[m] == shared_edge) {
      found = true;
      CHECK(f.signs[m] == -1);
    }
  }
  CHECK(found);
}

// Without deduplication both copies of the shared edge survive as parallel
// cells — the BRep "two curves over the same vertex pair" scenario.
void test_pushout_without_dedupe() {
  const graphos::Complex a = make_triangle();
  const graphos::Complex b = make_triangle();
  const std::vector<graphos::Identification> vertex_ids = {{0, 1, +1}, {1, 0, +1}};
  const auto po = graphos::pushout(a, b, vertex_ids, /*deduplicate=*/false);
  po.complex.validate();
  CHECK(po.complex.count(0) == 4);
  CHECK(po.complex.count(1) == 6);
  CHECK(po.complex.count(2) == 2);
  CHECK(graphos::d_squared_is_zero(po.complex));
}

// Deleting the star of the glued-on triangle's apex must sweep out its two
// private edges and its 2-cell, restoring a single triangle.
void test_star_deletion_cascade() {
  const graphos::Complex a = make_triangle();
  const graphos::Complex b = make_triangle();
  const std::vector<graphos::Identification> vertex_ids = {{0, 1, +1}, {1, 0, +1}};
  const auto po = graphos::pushout(a, b, vertex_ids, /*deduplicate=*/true);

  std::vector<std::vector<graphos::Index>> cells(3);
  cells[0].push_back(po.b_map.index[0][2]);  // B's apex in the pushout
  const auto sd = graphos::star_deletion(po.complex, cells);
  sd.complex.validate();

  CHECK(sd.complex.count(0) == 3);
  CHECK(sd.complex.count(1) == 3);
  CHECK(sd.complex.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(sd.complex));
  CHECK(graphos::euler_characteristic(sd.complex) == 1);

  // all of A survives the deletion; B's private cells are sent to zero
  const auto a_final = graphos::compose(po.a_map, sd.map);
  for (int k = 0; k <= 2; ++k) {
    for (std::size_t i = 0; i < a_final.index[k].size(); ++i) {
      CHECK(a_final.index[k][i] != graphos::invalid_index);
    }
  }
  const auto b_final = graphos::compose(po.b_map, sd.map);
  CHECK(b_final.index[0][2] == graphos::invalid_index);
  CHECK(b_final.index[1][1] == graphos::invalid_index);
  CHECK(b_final.index[1][2] == graphos::invalid_index);
  CHECK(b_final.index[2][0] == graphos::invalid_index);
  // ...but B's share of the interface persists, with the flip preserved
  CHECK(b_final.index[0][0] != graphos::invalid_index);
  CHECK(b_final.index[1][0] != graphos::invalid_index);
  CHECK(b_final.sign[1][0] == -1);
}

void test_find_parallel_cells_reports_orientation() {
  graphos::Complex c(1);
  c.attach_vertices(2);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {0, 1}, {+1, -1});  // same faces, reversed
  c.attach_cell(1, {0, 1}, {-1, +1});  // same faces, same orientation
  const auto dups = graphos::find_parallel_cells(c, 1);
  CHECK(dups.size() == 2);
  CHECK(dups[0].from == 1 && dups[0].to == 0 && dups[0].rel_sign == -1);
  CHECK(dups[1].from == 2 && dups[1].to == 0 && dups[1].rel_sign == +1);
}

void test_quotient_chain_resolution() {
  graphos::Complex c(0);
  c.attach_vertices(3);
  // 2 ~ 1 ~ 0: identification chains must resolve to the root
  std::vector<std::vector<graphos::Identification>> ids(1);
  ids[0] = {{2, 1, +1}, {1, 0, +1}};
  const auto q = graphos::quotient(c, ids);
  CHECK(q.complex.count(0) == 1);
  CHECK(q.map.index[0][0] == 0);
  CHECK(q.map.index[0][1] == 0);
  CHECK(q.map.index[0][2] == 0);
}

}  // namespace

int main() {
  test_triangle_invariants();
  test_disjoint_union_mixed_dimensional();
  test_pushout_with_orientation_flip();
  test_pushout_without_dedupe();
  test_star_deletion_cascade();
  test_find_parallel_cells_reports_orientation();
  test_quotient_chain_resolution();

  if (failures == 0) {
    std::printf("all tests passed\n");
    return 0;
  }
  std::printf("%d check(s) failed\n", failures);
  return 1;
}
