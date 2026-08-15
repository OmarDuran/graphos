#pragma once

// core: the complex, its operators ∂_k and δ_k, and its constructors
#include "graphos/core/build.hpp"
#include "graphos/core/coboundary.hpp"
#include "graphos/core/complex.hpp"
#include "graphos/core/frozen.hpp"
#include "graphos/core/incidence.hpp"
#include "graphos/core/marker.hpp"
#include "graphos/core/types.hpp"
#include "graphos/core/union_find.hpp"

// exec: hardware-portability seams
#include "graphos/exec/array.hpp"
#include "graphos/exec/forall.hpp"
#include "graphos/exec/memory.hpp"

// ops: the operation calculus — complex in, complex and chain maps out
#include "graphos/ops/agglomerate.hpp"
#include "graphos/ops/coarsen.hpp"
#include "graphos/ops/collapse.hpp"
#include "graphos/ops/cut.hpp"
#include "graphos/ops/disjoint_union.hpp"
#include "graphos/ops/join.hpp"
#include "graphos/ops/lift_identifications.hpp"
#include "graphos/ops/orient.hpp"
#include "graphos/ops/product.hpp"
#include "graphos/ops/pushout.hpp"
#include "graphos/ops/quotient.hpp"
#include "graphos/ops/replace.hpp"
#include "graphos/ops/star_deletion.hpp"
#include "graphos/ops/subcomplex.hpp"
#include "graphos/ops/subdivision.hpp"

// queries: homology, connectivity, manifoldness, and the structural invariants
#include "graphos/queries/adjacency.hpp"
#include "graphos/queries/amalgamation.hpp"
#include "graphos/queries/components.hpp"
#include "graphos/queries/facets.hpp"
#include "graphos/queries/homology.hpp"
#include "graphos/queries/manifold.hpp"
#include "graphos/queries/neighborhood.hpp"
