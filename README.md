# graphos

A metric-free computational topology engine. The formal specification —
objects, morphisms, the operation calculus, laws with their witnessing
tests, and design principles — is in [THEORY.md](THEORY.md). `graphos` manages finite cell
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
| **mesh ingestion** | `from_polygons` (vertex cycles), `from_polyhedra` (polygonal face lists — the general polymesh form), `from_edges`, and `from_simplices(d, …)` for d-simplices in ANY dimension (all strata derived, top cells oriented by input-order parity). Deterministic orientation conventions, ∂∘∂ = 0 by construction. Polytopal first: any polygon, any polyhedron, any mix |
| **chain map** | The output of every operation: where each cell went, with orientation coefficient; cells can be sent to zero |
| **coproduct** | `disjoint_union(a, b)` |
| **quotient** | `quotient(c, identifications)` — glue cells together, orientation flips propagate through stars |
| **parallel cells** | `find_parallel_cells(c, k)` — cells with equal boundary chains up to a uniform flip |
| **pushout A ⊔_C B** | `pushout(a, b, vertex_identifications)` — union glued along the shared subcomplex; the combinatorial core of BRep/domain union and mesh conformity |
| **star deletion** | `star_deletion(c, marker)` — remove marked cells and their closed stars (upward cascade); the combinatorial core of domain difference |
| **cutting along a subcomplex** | `cut_along(c, interface_marker)` — split the complex along a marked interface: sides are connected components of each closure cell's cut star, each side gets a copy, originals survive as the detached interface (fracture) domain. Tips/rims (one side) are not copied; junctions (3+ sides) get one copy per side. Purely topological — no geometric side test needed |
| **coboundary δ_k** | `coboundary(c, k)` — the signed transpose of ∂_{k+1}; applying it to a k-cochain is the discrete differential |
| **incidence I(k, j)** | `incidence(c, k, j)` — unsigned transitive incidence: per k-cell, the j-cells of its closure (j < k), star (j > k), or itself. The substrate of DoF gathers and NetworkX incidence views; deliberately unsigned (multi-level sign compositions telescope — that is ∂∘∂ = 0) |
| **connected components** | `connected_components(c, k, via[, exclude])` — components of k-cells through shared via-cells (excluding marked connectors = the sides-of-a-cut computation); `connected_components(c)` — β₀ of the whole mixed-dimensional complex |
| **product A × B** | `product(a, b)` — Cartesian product with Leibniz-rule boundary, ∂(α×β) = ∂α×β + (−1)^{dim α} α×∂β; `product(mesh, segment)` is extrusion. χ multiplies; factor structure returned as blocks |
| **lifting identifications** | `lift_identifications(c, vertex_pairs)` — extend a vertex pairing upward through the strata by boundary-chain matching (orientation flips reported); with `quotient`, periodic boundary conditions in two calls |
| **barycentric subdivision** | `barycentric_subdivision(c)` — the order complex of the face poset: one vertex per cell, one k-simplex per strict chain of k+1 cells. Simplicial output from ANY complex; the carrier map (maximal chain element) is the refinement/prolongation relation, and the SIGNED carrier (incidence product along the flag) is the subdivision chain map — refinement transfers orientation rather than leaving it to the caller |
| **dual complex** | `dual(frozen)` → `DualView` — the poset order-reversed, zero-copy: dual k-cells are primal (n−k)-cells, ∂^dual_k = δ_{n−k}. The combinatorial half of the DEC dual mesh; exokalk adds geometry via the subdivision |
| **orientation** | `orient(c)` — propagate a consistent global orientation across top cells (interior facets induced with opposite signs), or report non-orientability; the chain map records the flips for cochain transport |
| **Betti numbers (Z₂)** | `betti_numbers_z2(c)` — β_k by boundary-matrix reduction mod 2: components, tunnels, cavities; the dimension counts behind well-posedness of mixed formulations on multiply-connected domains |
| **manifoldness** | `check_manifold(c)` — necessary conditions with an offending-cell marker: purity, facet coface counts, link connectivity (pinch detection) |
| **agglomeration** | `agglomerate(c, labels)` — merge top cells into polytopal aggregates (interior facets cancel — inconsistent orientation is detected); the inverse of refinement, for multigrid/multiscale coarse spaces |
| **replacement (incremental imprinting)** | `replace(c, region, patch, vertex_glue)` — excision + amalgamation fused: excise the region's open star, amalgamate a patch along the frontier via a vertex glue lifted through the strata. The primitive of cut-cell workflows: imprint a fault through a cell without rebuilding the complex |
| **coarsening (dimensional descent)** | `coarsen(c, labels[, protected])` — agglomeration applied down the dimension ladder: cells, then interface patches, then edge chains (vertices only relabel — merging them is `quotient`'s job). Patches merge only when orientable, coefficient-clean, open, and coface-consistent — otherwise graceful singleton fallback. `protected` cells are caller-declared barriers (metric-free graphos cannot know corners); protecting the coarse lattice frame recovers refined tensor meshes EXACTLY |
| **relative homology** | `betti_numbers_z2(c, rel)` — β_k(K, A; Z₂) with A = cl(marked); makes the excision isomorphism H_k(K, cl st S) ≅ H_k(K∖st S, frontier S) computable on both sides |
| **join A ∗ B** | `join(a, b)` — ∂(α∗β) = ∂α∗β + (−1)^{dim α+1} α∗∂β, factors embedded as subcomplexes; join(point, X) is the cone, join(S⁰, X) the suspension |
| **neighborhood markers** | `star_of`, `closure_of`, `link_of`, `frontier_of` — set-level st(S), cl(S), lk(S) = cl(st S)∖st(cl S), and the excision frontier cl(st S)∖st(S) (the void boundary a pushout patch glues to) |
| **elementary collapse** | `free_faces(c)` (unique coface, single occurrence) and `collapse(c)` — greedy Whitehead collapsing to a free-face-free core; every step a simple homotopy equivalence |
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

## Queries

Read-only questions live in `include/graphos/queries/`, separated from the
transformations in `ops/`. The taxonomy:

| Question | Where |
|---|---|
| Boundary ∂_k(σ) / coboundary δ_k(σ) rows | `Complex::boundary`, `coboundary()`, `FrozenComplex::boundary_row/coboundary_row` (core) |
| Star, closure, link, excision frontier (per cell and per set) | `FrozenComplex::star/closure/link`, `queries/neighborhood.hpp` |
| Generalized adjacency (facet-adjacent, node-adjacent, any via-stratum) | `queries/adjacency.hpp` |
| Facet co-degree: boundary / interior / nonmanifold / free | `queries/facets.hpp` (`classify_facets`) |
| Closedness ∂K = ∅ (and ∂∂K = ∅ by composition) | `queries/facets.hpp` (`is_closed`) |
| Purity, facet condition, link connectivity | `queries/manifold.hpp` (`check_manifold`) |
| Link sphere/ball property per vertex (exact for dim ≤ 3) | `queries/manifold.hpp` (`classify_vertex_links`) |
| χ, Betti numbers (absolute and relative, Z₂) | `euler_characteristic` (core), `queries/homology.hpp` |
| b₀ / component labels (with excluded connectors) | `queries/components.hpp` |
| Global orientability | `ops/orient.hpp` (`orient().orientable` — an op, since it also constructs) |
| Cellular amalgamation criteria: common boundary ∂σ ∩ ∂τ and its acyclicity, excess intersection cl σ ∩ cl τ ∖ cl(∂σ ∩ ∂τ), and the union-is-a-cell verdict | `queries/amalgamation.hpp` (`common_boundary`, `excess_intersection`, `amalgamates_to_cell`) |

## Collective semantics (the SPMD contract)

graphos is written so the same program runs unchanged on a laptop and on a
cluster: every public operation is *specified* as collective, and the
current serial implementation is the P = 1 special case. Distribution
(global IDs, ParMETIS partitioning, halo exchange) will land inside
`freeze()` and the ops as an implementation detail — not as an API change.

Every public operation carries one of three contracts (PETSc vocabulary):

| Contract | Meaning | Operations |
|---|---|---|
| **Collective** | All ranks call it, same order; result is globally consistent | `from_edges`/`from_polygons`/`from_polyhedra`/`from_simplices`, `disjoint_union`, `product`, `join`, `cut_along`, `star_deletion`, `subcomplex`, `classify_facets`, `connected_components`, `barycentric_subdivision`, `orient`, `betti_numbers_z2`, `check_manifold`, `star_of`/`closure_of`/`link_of`/`frontier_of`, `free_faces`, `collapse`, `freeze`, `count`, `validate`, `d_squared_is_zero`, `euler_characteristic` |
| **Logically collective** | All ranks participate; arguments are supplied per-rank for locally owned cells | `quotient`, `pushout`, `lift_identifications` (identifications), `agglomerate`, `coarsen` (labels, protections), `replace` (region, glue) |
| **Local** | Per-rank, no communication; correct within the ghost ring | `star`, `closure`, `link`, `incidence`, `dual`, row access, views |

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

## Python bindings

`-DGRAPHOS_BUILD_PYTHON=ON` builds `graphos._core` (pybind11, fetched
automatically) into `python/graphos/`; run with `PYTHONPATH=python`. The
whole calculus and query surface is exposed with the same names; markers
take Python callables as predicates (`marker.mark_where(1, lambda e: ...)`),
identifications and glue are lists of `(from, to[, sign])` tuples, and
chain maps come back as NumPy arrays. Lifetime follows the builder/frozen
split: accessors on `Complex` and on operation results return array
copies; `FrozenComplex.boundary/coboundary` return zero-copy views whose
base keeps the frozen complex alive. `python/tests/test_graphos.py`
re-asserts the mathematical certificates through the bindings and runs as
the `python.test_graphos` ctest entry.

### NetworkX integration (`graphos.nx`)

Read-only graph **views** duck-type the NetworkX protocol directly over
the CSR arrays — `hasse_graph(frozen)` (nodes `(dim, index)`, arcs cell →
face with the orientation sign as an edge attribute), `incidence_graph(c,
k, j)` (the nfempy `build_graph(dim, codim)` shape), and
`adjacency_graph(c, k, via[, exclude])` (the mesh dual graph and its
generalizations). Each has `.to_networkx()` for guaranteed interop.

The views carry `__networkx_backend__ = "graphos"`, and the package
registers the **"graphos" NetworkX backend** (via the packaging entry
point; the source tree ships discovery metadata so `PYTHONPATH=python`
works too): plain `nx.connected_components(adjacency_graph(...))`
dispatches to the native C++ engine — including the sides-of-a-cut
computation via the `exclude` marker — and algorithms without a native
implementation fall back through conversion when
`nx.config.fallback_to_nx` is enabled. `python/tests/test_networkx.py`
covers protocol, dispatch, and fallback.

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
`ctest --test-dir build -R ops.test_cut`.

Beyond the unit suites:

- **Property-based fuzzing** (`property.test_properties`): random complexes
  through random op sequences, asserting the algebra's laws (validate,
  ∂∘∂ = 0, Euler–Poincaré, chain-map validity, collapse homotopy
  invariance) after every step.
- **Scaling assertions** (`property.test_scaling`): each op's runtime is
  checked against its complexity model on 8×-grown meshes — a complexity
  regression turns the suite red.
- **Performance guards** (`property.test_performance`): stored bytes/cell
  asserted against the CSR memory model (structural bloat fails
  deterministically), and bulk-op throughput asserted against the
  machine's own measured memcpy bandwidth — a constant-factor regression
  fails without hardcoding machine-specific numbers. Floors assert only in
  optimized, uninstrumented builds (sanitized/debug builds report only).
- **Multilevel hierarchies** (`property.test_multilevel`): three levels of
  quad (2D) and hex (3D) meshes built by tensor refinement (`product` of
  path complexes); the cross-level topology is asserted through
  `agglomerate` — parent-label coarsening reproduces the next level's
  structure, every seam facet connects exactly the parents of its fine
  cofaces, and level-by-level coarsening equals one-step coarsening with
  identical chain maps.
- **Benchmarks** (`bench/graphos_bench [n ...]`): structured tet grids,
  per-op wall time and throughput; the numbers behind any performance claim.
- **Sanitizers**: `-DGRAPHOS_SANITIZE=address,undefined` applies ASan/UBSan
  to everything in the tree.
- **CI** (`.github/workflows/ci.yml`): Linux gcc Release, Linux clang
  Debug+sanitizers, macOS Release, and a RAJA-enabled build, on every push
  and pull request. Portability
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
