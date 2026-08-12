#pragma once

#include "graphos/core/complex.hpp"

namespace graphos_test {

// Oriented triangle: edges e0=[0,1], e1=[1,2], e2=[2,0]; the 2-cell is
// attached along the cycle e0 + e1 + e2.
inline graphos::Complex make_triangle() {
  graphos::Complex c(2);
  c.attach_vertices(3);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  return c;
}

// Oriented segment: two vertices, one edge.
inline graphos::Complex make_segment() {
  graphos::Complex c(1);
  c.attach_vertices(2);
  c.attach_cell(1, {0, 1}, {-1, +1});
  return c;
}

// Disk made of two triangles sharing edge e0=[0,1]:
// A over vertex 2 (edges e0, e1=[1,2], e2=[2,0]),
// B over vertex 3 (edges e0, e3=[1,3], e4=[3,0]).
inline graphos::Complex make_two_triangle_disk() {
  graphos::Complex c(2);
  c.attach_vertices(4);
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 0}, {-1, +1});
  c.attach_cell(1, {1, 3}, {-1, +1});
  c.attach_cell(1, {3, 0}, {-1, +1});
  c.attach_cell(2, {0, 1, 2}, {+1, +1, +1});
  c.attach_cell(2, {0, 3, 4}, {+1, +1, +1});
  return c;
}

// Square fan: corners 0..3, center 4, four triangles around the center.
// Edges: spokes s0=[0,4](0) .. s3=[3,4](3), boundary b01(4), b12(5),
// b23(6), b30(7).
inline graphos::Complex make_fan() {
  graphos::Complex c(2);
  c.attach_vertices(5);
  c.attach_cell(1, {0, 4}, {-1, +1});
  c.attach_cell(1, {1, 4}, {-1, +1});
  c.attach_cell(1, {2, 4}, {-1, +1});
  c.attach_cell(1, {3, 4}, {-1, +1});
  c.attach_cell(1, {0, 1}, {-1, +1});
  c.attach_cell(1, {1, 2}, {-1, +1});
  c.attach_cell(1, {2, 3}, {-1, +1});
  c.attach_cell(1, {3, 0}, {-1, +1});
  c.attach_cell(2, {4, 1, 0}, {+1, +1, -1});
  c.attach_cell(2, {5, 2, 1}, {+1, +1, -1});
  c.attach_cell(2, {6, 3, 2}, {+1, +1, -1});
  c.attach_cell(2, {7, 0, 3}, {+1, +1, -1});
  return c;
}

}  // namespace graphos_test
