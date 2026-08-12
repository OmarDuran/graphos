# graphos

A metric-free computational topology engine. `graphos` manages finite cell
complexes — cells stratified by dimension, boundary operators ∂_k, chain maps
— and, eventually, distributed communication maps, completely decoupled from
physical coordinates, spatial metrics, and mesh geometry. Geometry, metrics,
Hodge stars, and PDE solvers belong to the companion project `exokalk`.

## Vocabulary

graphos speaks computational topology:

| Term | Meaning in graphos |
|---|---|
| **k-cell** | An abstract topological entity of dimension k; nothing but an index |
| **complex** | Cell counts per dimension + boundary operators; may be mixed-dimensional (maximal cells in any skeleton) |
| **boundary operator ∂_k** | Signed CSR matrix from k-cells to (k-1)-cells, entries ±1; `d_squared_is_zero` checks ∂∘∂ = 0 |
| **attaching a cell** | Adding a k-cell along its boundary chain (`attach_cell`) |
| **chain map** | The output of every operation: where each cell went, with orientation coefficient; cells can be sent to zero |
| **coproduct** | `disjoint_union(a, b)` |
| **quotient** | `quotient(c, identifications)` — glue cells together, orientation flips propagate through stars |
| **parallel cells** | `find_parallel_cells(c, k)` — cells with equal boundary chains up to a uniform flip |
| **pushout A ⊔_C B** | `pushout(a, b, vertex_identifications)` — union glued along the shared subcomplex; the combinatorial core of BRep/domain union and mesh conformity |
| **star deletion** | `star_deletion(c, cells)` — remove cells and their closed stars (upward cascade); the combinatorial core of domain difference |

Decisions that require geometry ("these two vertices are the same point")
enter as explicit `Identification` inputs; everything downstream is
combinatorial.

## Architecture

- **Metric-free isolation.** No coordinates, lengths, areas, or normals ever
  enter graphos. Topology is integer index spaces plus signed CSR incidence.
- **Frozen complexes, bulk edits.** Operations are out-of-place
  transformations: complex in, complex + chain map out. Fields (cochains)
  attached to cells are transported through the chain maps.
- **Flat arrays only.** The poset is CSR strata, never a pointer graph — the
  precondition for GPU portability and zero-copy NetworkX views.
- **Portability seams.** The RAJA-CHAI-Umpire triplet connects at three
  points:
  1. `exec::forall` / `exec::*_scan` ([exec/forall.hpp](include/graphos/exec/forall.hpp)) —
     kernel phases dispatch through RAJA policies (`GRAPHOS_ENABLE_RAJA`),
     with an identical-semantics serial fallback. `star_deletion` is the
     exemplar op in full kernel form (mark → cascade → scan → scatter).
  2. `exec::Buffer<T>` ([exec/memory.hpp](include/graphos/exec/memory.hpp)) —
     kernel scratch drawn from Umpire pools (`GRAPHOS_ENABLE_UMPIRE`), plain
     heap otherwise.
  3. Persistent complex storage as CHAI-managed arrays — planned; lands with
     the frozen/builder split, since CHAI's copy-on-context semantics
     presuppose frozen storage. Device execution policies follow from it.

## Build and test

```bash
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

Requires CMake ≥ 3.20 and a C++20 compiler; the default build has no other
dependencies (header-only, `graphos::graphos` CMake target). Portability
options:

```bash
cmake -S . -B build -DGRAPHOS_ENABLE_RAJA=ON -DGRAPHOS_FETCH_TPL=ON
```

`GRAPHOS_FETCH_TPL` downloads and builds the third-party libs via
FetchContent; without it they are located with `find_package`.

## Roadmap

1. Interface split / conformity cut (vertex duplication along a cut set) with
   topological side classification.
2. Derived operators: coboundary (transposes), multi-level closures, star and
   link queries, fixed-arity strided storage for single-cell-type complexes.
3. Frozen/builder split; persistent storage on CHAI/Umpire; device execution
   policies (CUDA/HIP/SYCL) for the bulk kernels.
4. Python bindings (pybind11) with opaque handles; NetworkX 3.x backend
   exposing zero-copy structural views.
5. Distributed layer: global IDs, METIS/ParMETIS partitions, MPI halo
   exchange.
