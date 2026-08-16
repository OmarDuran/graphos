#include "fixtures.hpp"
#include "graphos/core/complex.hpp"
#include "graphos_test.hpp"

GRAPHOS_TEST(triangle_invariants) {
  const graphos::Complex c = graphos_test::make_triangle();
  c.validate();
  CHECK(c.count(0) == 3);
  CHECK(c.count(1) == 3);
  CHECK(c.count(2) == 1);
  CHECK(graphos::d_squared_is_zero(c));
  CHECK(graphos::euler_characteristic(c) == 1);  // disk
}

GRAPHOS_TEST(segment_invariants) {
  const graphos::Complex c = graphos_test::make_segment();
  c.validate();
  CHECK(graphos::euler_characteristic(c) == 1);
}

GRAPHOS_TEST(attach_cell_rejects_bad_input) {
  graphos::Complex c(2);
  c.attach_vertices(2);
  CHECK_THROWS(c.attach_cell(1, {0, 7}, {-1, +1}));  // face out of range
  CHECK_THROWS(c.attach_cell(1, {0, 1}, {-1, +2}));  // sign not +/-1
  CHECK_THROWS(c.attach_cell(3, {0, 1}, {-1, +1}));  // dimension out of range
  CHECK(c.count(1) == 0);
}

GRAPHOS_TEST(d_squared_detects_broken_orientation) {
  graphos::Complex c(2);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  // e0 + e1 − e2 does not close: ∂∂ = 2v2 − 2v0 ≠ 0
  c.attach_cell(2, {0, 1, 2}, {+1, +1, -1});
  c.validate();                           // structurally fine...
  CHECK(!graphos::d_squared_is_zero(c));  // ...but not a chain complex
}

GRAPHOS_TEST(bulk_constructor_checks_stratification) {
  std::vector<graphos::BoundaryOperator> strata(2);
  strata[1].append_row(std::vector<graphos::Index>{0, 1}, std::vector<graphos::Sign>{-1, +1});
  // N₁ claims two 1-cells but ∂₁ has one row
  CHECK_THROWS(graphos::Complex({2, 2}, std::move(strata)));
}

// -- commutes_with_boundary: ∂π = π∂ ---------------------------------------

namespace {

// π = id, the chain map every complex carries
graphos::ChainMap identity_map(const graphos::Complex& c) {
  graphos::ChainMap m = graphos::ChainMap::sized(c.counts());
  for (int k = 0; k <= c.dim(); ++k) {
    for (graphos::Index i = 0; i < c.count(k); ++i) {
      m.index[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = i;
    }
  }
  return m;
}

// π = 0: every generator to zero
graphos::ChainMap zero_map(const graphos::Complex& c) {
  graphos::ChainMap m = graphos::ChainMap::sized(c.counts());
  for (int k = 0; k <= c.dim(); ++k) {
    for (graphos::Index i = 0; i < c.count(k); ++i) {
      m.index[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = graphos::invalid_index;
    }
  }
  return m;
}

}  // namespace

GRAPHOS_TEST(identity_and_zero_are_chain_maps) {
  const graphos::Complex tri = graphos_test::make_triangle();
  CHECK(graphos::commutes_with_boundary(tri, tri, identity_map(tri)));
  // the zero map is a chain map: ∂0 = 0 = 0∂
  CHECK(graphos::commutes_with_boundary(tri, tri, zero_map(tri)));

  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  CHECK(graphos::commutes_with_boundary(disk, disk, identity_map(disk)));
  CHECK(graphos::commutes_with_boundary(disk, disk, zero_map(disk)));
}

// −id is a chain map (∂ is linear); −id on ONE stratum is not, since the sign
// then fails to pass through ∂ between that stratum and its neighbour
GRAPHOS_TEST(a_global_sign_flip_commutes_but_a_stratum_wise_one_does_not) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  graphos::ChainMap neg = identity_map(disk);
  for (int k = 0; k <= disk.dim(); ++k) {
    for (auto& sg : neg.sign[static_cast<std::size_t>(k)]) sg = graphos::Sign{-1};
  }
  CHECK(graphos::commutes_with_boundary(disk, disk, neg));

  graphos::ChainMap one_stratum = identity_map(disk);
  for (auto& sg : one_stratum.sign[1]) sg = graphos::Sign{-1};  // edges only
  CHECK(!graphos::commutes_with_boundary(disk, disk, one_stratum));
}

GRAPHOS_TEST(a_perturbed_map_is_rejected) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();

  // one edge sent to the wrong edge: ∂ of the image no longer matches
  graphos::ChainMap wrong_target = identity_map(disk);
  wrong_target.index[1][0] = (disk.count(1) > 1) ? 1 : 0;
  CHECK(!graphos::commutes_with_boundary(disk, disk, wrong_target));

  // one edge sent to zero while its cofaces are kept: π∂ loses a term
  graphos::ChainMap partial = identity_map(disk);
  partial.index[1][0] = graphos::invalid_index;
  CHECK(!graphos::commutes_with_boundary(disk, disk, partial));

  // one vertex flipped: ∂ of an incident edge changes sign on one endpoint only
  graphos::ChainMap flip_vertex = identity_map(disk);
  flip_vertex.sign[0][0] = graphos::Sign{-1};
  CHECK(!graphos::commutes_with_boundary(disk, disk, flip_vertex));
}

GRAPHOS_TEST(a_malformed_map_is_rejected_rather_than_read) {
  const graphos::Complex disk = graphos_test::make_two_triangle_disk();
  const graphos::Complex tri = graphos_test::make_triangle();

  graphos::ChainMap short_rows = identity_map(disk);
  short_rows.index[1].pop_back();  // one row short of the source
  CHECK(!graphos::commutes_with_boundary(disk, disk, short_rows));

  graphos::ChainMap ragged = identity_map(disk);
  ragged.sign[1].pop_back();  // index and sign disagree
  CHECK(!graphos::commutes_with_boundary(disk, disk, ragged));

  // the disk's identity indexes cells the triangle does not have
  CHECK(!graphos::commutes_with_boundary(disk, tri, identity_map(disk)));

  graphos::ChainMap too_few_strata = identity_map(disk);
  too_few_strata.index.pop_back();
  too_few_strata.sign.pop_back();
  CHECK(!graphos::commutes_with_boundary(disk, disk, too_few_strata));
}

// A degenerate source: with no 1-cells there is nothing for ∂ to test, so any
// shape-valid map passes. Pinned so the predicate is not read as stronger.
GRAPHOS_TEST(a_zero_dimensional_source_imposes_no_condition) {
  graphos::Complex pts(0);
  pts.attach_vertices(3);
  graphos::ChainMap m = graphos::ChainMap::sized(pts.counts());
  m.index[0] = {2, 0, 1};
  CHECK(graphos::commutes_with_boundary(pts, pts, m));
}

GRAPHOS_TEST_MAIN()
