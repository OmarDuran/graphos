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
| **star deletion** | `star_deletion(c, marker)` — remove marked cells and their closed stars (upward cascade); the combinatorial core of domain difference |
| **cutting along a subcomplex** | `cut_along(c, interface_marker)` — split the complex along a marked interface: sides are connected components of each closure cell's cut star, each side gets a copy, originals survive as the detached interface (fracture) domain. Tips/rims (one side) are not copied; junctions (3+ sides) get one copy per side. Purely topological — no geometric side test needed |
| **coboundary δ_k** | `coboundary(c, k)` — the signed transpose of ∂_{k+1}; applying it to a k-cochain is the discrete differential |
| **frozen complex** | `freeze(c)` → `FrozenComplex` — the immutable query object: ∂ and δ in device-capable storage (`exec::Array`, the CHAI seam), host row access, kernel views |
| **star / closure / link** | `FrozenComplex::star/closure/link(k, cell)` — st(σ), cl(σ), lk(σ) = cl(st(σ)) \ st(cl(σ)); the queries NetworkX views will sit on |
| **marker** | `Marker` — locally-evaluated, collectively-meaningful cell selection (`mark`, `mark_where(k, pred)`); the argument form of every cell-selecting op |
| **closed subcomplex** | `subcomplex(c, marker)` — cl(marked): the marked cells with their full closure, extracted as its own complex with both chain maps (parent→sub and the embedding sub→parent). Fracture domains, material regions, boundary complexes, k-skeleta |
| **facet classification** | `classify_facets(c)` — every (n−1)-cell into exactly one of: free (0 cofaces — detached interface domains), boundary (1 — ∂Ω), interior (2), nonmanifold (3+ — DFN junctions). Compose with `subcomplex` to extract ∂Ω |

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
  3. `exec::Array<T>` ([exec/array.hpp](include/graphos/exec/array.hpp)) —
     persistent storage for frozen complexes, `chai::ManagedArray`-backed
     under `GRAPHOS_ENABLE_CHAI` (host fallback otherwise). `Complex` is the
     mutable builder; `freeze()` produces the immutable `FrozenComplex`
     whose ∂/δ arrays and `CsrView`s are what device kernels will capture.
     Device execution policies (CUDA/HIP/SYCL) are the remaining step.

## Collective semantics (the SPMD contract)

graphos is written so the same program runs unchanged on a laptop and on a
cluster: every public operation is *specified* as collective, and the
current serial implementation is the P = 1 special case. Distribution
(global IDs, ParMETIS partitioning, halo exchange) will land inside
`freeze()` and the ops as an implementation detail — not as an API change.

Every public operation carries one of three contracts (PETSc vocabulary):

| Contract | Meaning | Operations |
|---|---|---|
| **Collective** | All ranks call it, same order; result is globally consistent | `disjoint_union`, `cut_along`, `star_deletion`, `subcomplex`, `classify_facets`, `freeze`, `count`, `validate`, `d_squared_is_zero`, `euler_characteristic` |
| **Logically collective** | All ranks participate; arguments are supplied per-rank for locally owned cells | `quotient`, `pushout` (identifications) |
| **Local** | Per-rank, no communication; correct within the ghost ring | `star`, `closure`, `link`, row access, views |

Two disciplines follow:

1. **The one law**: a distributed program must call the collective
   operations in the same order on every rank — never branch a
   topology-changing call on rank-local data. This is the only way
   distribution is visible in user code.
2. **Selection is by marking, not by index lists**: ops take a `Marker`,
   marked locally (`mark_where(k, predicate)` is the canonical form). A
   predicate evaluates on each rank over its own cells, so the same program
   text is meaningful at any rank count.

`freeze(c, halo_depth)` records the ghost-ring depth (default 1) that Local
queries are guaranteed correct within.

## Build, test, install

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build
cmake --install build --prefix ~/opt/graphos
```

Consume the installed package with `find_package(graphos CONFIG REQUIRED)`
and link `graphos::graphos`. **Full instructions — including the
RAJA/Umpire/CHAI stack via the bundled Spack environment
([spack.yaml](spack.yaml)), GPU variants, and troubleshooting — are in
[INSTALL.md](INSTALL.md).**

Requires CMake ≥ 3.20 and a C++20 compiler; the default build has no other
dependencies (header-only, `graphos::graphos` CMake target).

Unit tests mirror the include tree — one test file per header, one ctest
entry per file (`tests/core/test_complex.cpp` ↔
`include/graphos/core/complex.hpp`, registered as `core.test_complex`), with
a zero-dependency harness in `tests/graphos_test.hpp` and shared complexes
in `tests/fixtures.hpp`. Run a single suite with e.g.
`ctest --test-dir build -R ops.test_cut`. Portability
options:

```bash
cmake -S . -B build -DGRAPHOS_ENABLE_RAJA=ON -DGRAPHOS_FETCH_TPL=ON
```

`GRAPHOS_FETCH_TPL` downloads and builds the third-party libs via
FetchContent; without it they are located with `find_package`.

## Roadmap

1. Kernel-form port of `cut_along` (side labeling as parallel label
   propagation) and `quotient`; 3D cut coverage (branching fracture tests).
2. Device execution policies (CUDA/HIP/SYCL) for the bulk kernels, on the
   frozen storage; CHAI-enabled CI build. Fixed-arity strided storage for
   single-cell-type complexes.
4. Python bindings (pybind11) with opaque handles; NetworkX 3.x backend
   exposing zero-copy structural views.
5. Distributed layer: global IDs, METIS/ParMETIS partitions, MPI halo
   exchange.
