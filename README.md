# graphos

A metric-free computational topology engine.

graphos represents a finite stratified signed incidence structure
C = (N₀, …, N_n ; ∂₁, …, ∂_n) — cell counts per stratum, boundary operators
∂_k with entries in {−1, 0, +1}, and the chain maps its operations induce. A
k-cell is an index: no realization, coordinate or metric is stored, and none
is needed to state or check anything below. Complexes may be
mixed-dimensional, polytopal, non-regular, non-manifold and non-orientable;
validity is queried rather than presupposed.

The formal specification — objects, morphisms, the operation calculus, the
laws and the tests that witness them — is [THEORY.md](THEORY.md). Geometry,
metrics, Hodge stars and PDE solvers belong to the companion project
`exokal`.

## Vocabulary

The surface, in the vocabulary of THEORY.md:

| Term | Meaning in graphos |
|---|---|
| **k-cell** | A cell of stratum k; an index, carrying no realization |
| **complex** | C = (N₀, …, N_n ; ∂₁, …, ∂_n); a chain complex exactly when ∂∘∂ = 0, which `d_squared_is_zero` decides |
| **boundary operator ∂_k** | ∂_k ∈ {−1, 0, +1}^{N_{k−1} × N_k} in signed CSR; row σ holds the (k−1)-faces τ of σ with [σ : τ] |
| **attaching a cell** | `attach_cell` — a k-cell along the (k−1)-chain ∂σ, the CW paradigm |
| **mesh ingestion** | `from_edges`, `from_polygons` (vertex cycles), `from_polyhedra` (face lists), `from_simplices(d, …)` in any dimension. Intermediate strata are derived, orientation conventions are deterministic, and ∂∘∂ = 0 holds by construction |
| **chain map** | f_* : C_k(C) → C_k(C′), induced on generators by every operation; a generator may be sent to 0 |
| **coproduct** | `disjoint_union(a, b)` — A ⊔ B, strata concatenated; the factors may differ in dimension |
| **quotient** | `quotient(c, identifications)` — C/~; an orientation flip propagates into st(σ) |
| **parallel cells** | `find_parallel_cells(c, k)` — k-cells whose boundary chains agree as signed sets up to a uniform flip |
| **pushout A ⊔_C B** | `pushout(a, b, vertex_identifications)` — glued along the subcomplex C the identified vertices generate; deduplication extends C upward and makes the interface conforming |
| **star deletion** | `star_deletion(c, marker)` — removes st(S), leaving the largest subcomplex meeting no marked cell |
| **cutting along a subcomplex** | `cut_along(c, interface_marker)` — for x ∈ cl(S), the sides are the incidence-connected components of st(x) ∖ cl(S), each receiving a copy of x. cl(S) survives as a detached subcomplex; a one-sided rim is not copied, giving the crack front, and a junction gets one copy per side. No geometric side test enters |
| **coboundary δ_k** | `coboundary(c, k)` — δ_k = ∂_{k+1}^T; applied to a k-cochain it is the discrete exterior derivative |
| **incidence I(k, j)** | `incidence(c, k, j)` — per k-cell σ, the j-cells of cl(σ) (j < k), of st(σ) (j > k), or σ itself. Unsigned by necessity: a multi-level composite of incidence numbers telescopes to 0 by ∂∘∂ = 0 |
| **connected components** | `connected_components(c, k, via[, exclude])` — components of k-cells adjacent through a common via-cell; excluding marked connectors gives the sides of a cut. `connected_components(c)` is β₀ of the whole complex |
| **product A × B** | `product(a, b)` — ∂(α × β) = ∂α × β + (−1)^{dim α} α × ∂β. `product(mesh, segment)` is extrusion; χ(A × B) = χ(A) · χ(B); the factor structure survives in `blocks` |
| **lifting identifications** | `lift_identifications(c, vertex_pairs)` — extends a vertex pairing upward by boundary-chain matching, reporting rel_sign. With `quotient`, a periodic identification in two calls: a square becomes a cylinder, a twisted pairing a Möbius band |
| **barycentric subdivision** | `barycentric_subdivision(c)` — sd(C), the order complex of the face poset: one vertex per cell, one k-simplex per strict chain of k+1 cells. Simplicial from any complex. The carrier is the refinement relation; the signed carrier ∏ᵢ [cᵢ : cᵢ₊₁] is the subdivision chain map, so refinement transfers orientation |
| **dual complex** | `dual(frozen)` → `DualView` — the face poset order-reversed, zero-copy: the dual k-cell is the primal (n−k)-cell and ∂^dual_k = δ_{n−k}. The combinatorial half of the DEC dual mesh |
| **orientation** | `orient(c, k)` — coherently orients the maximal cells of stratum k, so every interior facet is induced with opposite incidence numbers, or reports non-orientability. The orientation classes are returned: coherence leaves exactly one free sign per class, which only a geometric notion of outward can settle |
| **Betti numbers (Z₂)** | `betti_numbers_z2(c)` — β_k = N_k − rank ∂_k − rank ∂_{k+1} by reduction mod 2. These fix the dimension of the harmonic space, hence well-posedness of mixed formulations on multiply-connected domains |
| **manifoldness** | `check_manifold(c)` — the decidable necessary conditions, with an offending-cell marker: purity, facet coface counts, connectivity of st(σ) for dim σ ≤ n−2 |
| **agglomeration** | `agglomerate(c, labels)` — merges top cells into aggregates with boundary Σ ∂(members); interior facets cancel exactly when the stratum is coherently oriented, and an incoherent one is detected |
| **replacement** | `replace(c, region, patch, vertex_glue)` — excision and amalgamation fused: excise st(S), glue a patch along the frontier through a vertex correspondence lifted by boundary-chain matching. Imprints an interface through a cell without rebuilding the complex |
| **coarsening** | `coarsen(c, labels[, protected])` — agglomeration by dimensional descent: top cells, then interface patches, then 1-chains; vertices only relabel. A patch merges only if coherently orientable, coefficient-clean, open and coface-consistent, else its members stay singletons. `protected` cells are declared barriers, since a metric-free engine cannot locate a corner; protecting the coarse lattice frame inverts tensor refinement exactly |
| **relative homology** | `betti_numbers_z2(c, rel)` — β_k(K, A; Z₂) with A = cl(marked), making the excision isomorphism H_k(K, cl st S) ≅ H_k(K ∖ st S, frontier S) computable on both sides |
| **join A ∗ B** | `join(a, b)` — ∂(α ∗ β) = ∂α ∗ β + (−1)^{dim α+1} α ∗ ∂β, both factors embedded. join(point, X) is the cone, join(S⁰, X) the suspension, and cl(st σ) ≅ σ ∗ lk(σ) becomes computable on both sides |
| **neighbourhood markers** | `star_of`, `closure_of`, `link_of`, `frontier_of` — st(S), cl(S), lk(S) = cl(st S) ∖ st(cl S), and the frontier cl(st S) ∖ st(S), which bounds the void a pushout patch glues to |
| **elementary collapse** | `free_faces(c)` — τ with one proper coface, occurring once in ∂σ — and `collapse(c)`, which removes free pairs until none remains. Every step is a simple homotopy equivalence, so the homotopy type is preserved |
| **frozen complex** | `freeze(c)` → `FrozenComplex` — the immutable object queries sit on: ∂ and δ in device-capable storage (`exec::Array`, the CHAI seam), with host rows and kernel views |
| **star / closure / link** | `FrozenComplex::star/closure/link(k, cell)` — st(σ), cl(σ), lk(σ) = cl(st σ) ∖ st(cl σ), per cell |
| **marker** | `Marker` — a selection of cells, evaluated locally and meaningful collectively (`mark`, `mark_where(k, pred)`); the argument form of every operation on a subcomplex |
| **closed subcomplex** | `subcomplex(c, marker)` — cl(S) as its own complex, with both chain maps. Interface domains, material regions, ∂K, and the k-skeleton are all instances |
| **facet classification** | `classify_facets(c)` — partitions the facets by top-coface count: `maximal` (0, a detached interface — no proper coface, *not* free in the Whitehead sense, which is `free_faces`), `boundary` (1, ∂K), `interior` (2), `nonmanifold` (3+, a junction). With `subcomplex` it extracts ∂K |

A decision that requires geometry — that two vertices are the same point —
enters as an explicit `Identification`; everything downstream of it is
combinatorial.

## Architecture

- **Metric-free.** No coordinate, length, area or normal enters graphos. A
  complex is integer index spaces with signed CSR incidence.
- **Out-of-place operations.** Complex in, complex and chain map out; cochains
  are transported by gathering through the chain map.
- **Flat arrays.** The face poset is CSR strata, never a pointer graph — the
  precondition for device portability and zero-copy views.
- **Portability seams.** The RAJA-CHAI-Umpire triplet connects at three
  points:
  1. `exec::forall` / `exec::*_scan` ([exec/forall.hpp](include/graphos/exec/forall.hpp)) —
     kernel phases dispatch through RAJA policies (`GRAPHOS_ENABLE_RAJA`),
     with a serial fallback of identical semantics. `star_deletion` and
     `subcomplex` are in full kernel form — mark → cascade → scan → scatter —
     and are the shape the remaining operations are being moved to.
  2. `exec::Buffer<T>` ([exec/memory.hpp](include/graphos/exec/memory.hpp)) —
     kernel scratch drawn from Umpire pools (`GRAPHOS_ENABLE_UMPIRE`), plain
     heap otherwise.
  3. `exec::Array<T>` ([exec/array.hpp](include/graphos/exec/array.hpp)) —
     persistent storage for frozen complexes, `chai::ManagedArray`-backed
     under `GRAPHOS_ENABLE_CHAI` (host fallback otherwise). `Complex` is the
     mutable builder; `freeze()` produces the immutable `FrozenComplex`
     whose ∂ and δ arrays and `CsrView`s a device kernel captures. Device
     policies (CUDA/HIP/SYCL) do not exist yet — see the Roadmap.

## Queries

Read-only questions live in `queries/`, apart from the transformations in
`ops/`:

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
| Orientability of a stratum, and its orientation classes | `ops/orient.hpp` (`orient().orientable`, `.class_of`) — an operation, since it also constructs |
| Amalgamation: ∂σ ∩ ∂τ and its Z₂-acyclicity, the excess cl σ ∩ cl τ ∖ cl(∂σ ∩ ∂τ), and whether σ ∪ τ is a cell | `queries/amalgamation.hpp` |

## Collective semantics (the SPMD contract)

Every public operation is *specified* as collective; the present serial
implementation is the P = 1 case. Distribution — global IDs, ParMETIS
partitioning, halo exchange — is intended to land inside `freeze()` and the
operations without an API change. None of it exists yet (see the Roadmap);
what follows is the contract it will honour.

Three contracts, in PETSc vocabulary:

| Contract | Meaning | Operations |
|---|---|---|
| **Collective** | All ranks call it, same order; result is globally consistent | `from_edges`/`from_polygons`/`from_polyhedra`/`from_simplices`, `disjoint_union`, `product`, `join`, `cut_along`, `star_deletion`, `subcomplex`, `classify_facets`, `connected_components`, `barycentric_subdivision`, `orient`, `betti_numbers_z2`, `check_manifold`, `star_of`/`closure_of`/`link_of`/`frontier_of`, `free_faces`, `collapse`, `freeze`, `count`, `validate`, `d_squared_is_zero`, `euler_characteristic` |
| **Logically collective** | All ranks participate; arguments are supplied per-rank for locally owned cells | `quotient`, `pushout`, `lift_identifications` (identifications), `agglomerate`, `coarsen` (labels, protections), `replace` (region, glue) |
| **Local** | Per-rank, no communication; correct within the ghost ring | `star`, `closure`, `link`, `incidence`, `dual`, row access, views |

Two disciplines follow:

1. **Same order on every rank.** A distributed program must call the
   collective operations in one order — never branch a topology-changing call
   on rank-local data. This is the only way distribution appears in user code.
2. **Selection by marking, not by index list.** Operations take a `Marker`,
   marked locally; `mark_where(k, predicate)` is the canonical form. The
   predicate runs on each rank over its own cells, so the same program text is
   meaningful at any rank count.

`freeze(c)` takes a complex and nothing else. The ghost-ring depth that bounds
the Local queries is a property of the distributed layer, and arrives as an
argument when that layer does.

## Python bindings

`-DGRAPHOS_BUILD_PYTHON=ON` builds `graphos._core` (pybind11, pinned) into
`python/graphos/`; run with `PYTHONPATH=python`. The whole calculus and query
surface carries the same names: a marker takes a Python predicate
(`marker.mark_where(1, lambda e: ...)`), identifications and glue are
`(from, to[, sign])` tuples, and chain maps return as NumPy arrays.

Lifetime follows the builder/frozen split — accessors on `Complex` and on
operation results return copies, while `FrozenComplex.boundary/coboundary`
return zero-copy views whose base keeps the complex alive.
`python/tests/test_graphos.py` re-asserts the laws through the bindings as the
`python.test_graphos` ctest entry.

### NetworkX integration (`graphos.nx`)

Read-only views duck-type the NetworkX protocol directly over the CSR
arrays: `hasse_graph(frozen)` — the face poset, nodes `(dim, index)` and arcs
σ → τ carrying [σ : τ] as an edge attribute — `incidence_graph(c, k, j)` for
I(k, j), and `adjacency_graph(c, k, via[, exclude])` for adjacency through a
via-stratum. Each has `.to_networkx()` for guaranteed interop.

The views carry `__networkx_backend__ = "graphos"` and the package registers
the backend through the packaging entry point, with discovery metadata in the
source tree so `PYTHONPATH=python` works too. A plain
`nx.connected_components(adjacency_graph(...))` then dispatches to the C++
engine — the sides-of-a-cut computation included, through the `exclude`
marker — and an algorithm with no native implementation falls back through
conversion when `nx.config.fallback_to_nx` is set.
`python/tests/test_networkx.py` covers protocol, dispatch and fallback.

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

Unit tests mirror the include tree: one file per header, one ctest entry per
file (`tests/core/test_complex.cpp` ↔ `include/graphos/core/complex.hpp`,
registered as `core.test_complex`), on a zero-dependency harness in
`tests/graphos_test.hpp` with shared complexes in `tests/fixtures.hpp`. Run
one suite with `ctest --test-dir build -R ops.test_cut`.

Beyond the unit suites:

- **Property-based fuzzing** (`property.test_properties`): random complexes
  through random operation sequences, asserting the laws after every step —
  the structural invariants, ∂∘∂ = 0, Euler–Poincaré, chain-map validity, and
  homotopy invariance under collapse.
- **Scaling assertions** (`property.test_scaling`): each operation's runtime
  against its complexity model on 8×-grown complexes, so a change of
  complexity class fails.
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
- **Benchmarks** (`bench/graphos_bench [n ...]`): structured simplicial
  grids, per-operation wall time and throughput — the evidence behind any
  performance claim.
- **Sanitizers**: `-DGRAPHOS_SANITIZE=address,undefined` applies ASan/UBSan
  to everything in the tree.
- **CI** (`.github/workflows/ci.yml`): a `format` check
  (`clang-format --dry-run --Werror`, version pinned), a dependency-free
  `baseline` matrix (Linux gcc Release, Linux clang Debug + sanitizers, macOS
  Release), and four staged jobs that answer whether the *installation* is
  reliable — `1. tpls` builds camp/RAJA/Umpire/CHAI from source
  ([scripts/build_tpls.sh](scripts/build_tpls.sh), cached on the script's own
  hash), `2. library` builds and installs against them, `3. tests` runs ctest
  and then compiles a fresh consumer against the install alone, and
  `4. container` runs the same four stages in [Dockerfile](Dockerfile) from a
  bare image. Each phase is a separate job, so a red run names its cause.

Portability options:

```bash
cmake -S . -B build -DGRAPHOS_ENABLE_RAJA=ON -DGRAPHOS_FETCH_TPL=ON
```

`GRAPHOS_FETCH_TPL` downloads and builds the third-party libs via
FetchContent; without it they are located with `find_package`. CHAI is
reached through `find_package` only, so it comes from the Spack environment
or from `scripts/build_tpls.sh`.

## Roadmap

Shipped since this list was last written: the Python bindings and the
NetworkX 3.x backend (§ Python bindings), and the staged installation CI with
its container (§ Build, test, install).

1. **Kernel form for the operation calculus.** Two of the fifteen operations
   go through the exec seam today — `star_deletion` and `subcomplex`, in the
   mark → cascade → scan → scatter shape. The other thirteen are serial host
   code: `cut_along` (side labelling as parallel label propagation) and
   `quotient` are the load-bearing ones, then `agglomerate`, `coarsen`,
   `orient`, `subdivision` and the rest. Also 3-dimensional cut coverage,
   branching interfaces included.
2. **Device execution.** No CUDA/HIP/SYCL policy exists yet; `exec::forall`
   dispatches to RAJA host policies or to plain loops. Device policies over
   the frozen storage, and a CHAI-enabled CI build, follow (1) — there is
   little to move to a device until the operations are in kernel form.
3. **Fixed-arity strided storage.** A complex of a single cell type needs no
   offset array, only a stride — the specialization `BoundaryOperator`
   anticipates.
4. **Distributed layer.** Global IDs, METIS/ParMETIS partitioning and MPI halo
   exchange inside `freeze()` and the operations. The contracts in
   § Collective semantics are written for it; today every one is the P = 1
   case, and no MPI or ParMETIS call exists in the tree.
