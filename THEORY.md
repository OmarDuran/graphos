# graphos — Formal Specification

This document states the mathematics graphos implements: its objects, its
morphisms, the operation calculus, the laws the library guarantees, and the
principles that constrain its design. Every law cites the test that
witnesses it — a claim without a citation is not made.

## 1. Objects

A **complex** is a finite stratified signed incidence structure

  C = (N₀, …, N_n ; ∂₁, …, ∂_n)

where N_k ∈ ℕ counts the *k-cells* (a k-cell is nothing but an index — no
geometric realization is stored) and each **boundary operator**
∂_k ∈ {−1, 0, +1}^{N_{k−1} × N_k} is a sparse signed matrix, stored CSR
with entries strictly ±1 (absence encodes 0). The structural invariants —
row counts, index ranges, sign values — are `validate()`.

C is a **chain complex** when ∂_{k−1} ∘ ∂_k = 0 for all k
(`d_squared_is_zero`). The data structure does not force this; the
constructors establish it, and every operation preserves it (Law L1).

Deliberate generalities: complexes may be **mixed-dimensional** (maximal
cells in any stratum), **polytopal** (no cell-type assumption anywhere —
a cell is its boundary chain), non-regular (loop edges, repeated faces are
representable; queries that require regularity skip degenerate cells and
say so), non-manifold, and non-orientable. Nothing in the core assumes
otherwise; validity is *queried*, not presupposed.

The **dual** of a complex is its face poset order-reversed: the dual k-cell
is the primal (n−k)-cell and ∂^dual_k = δ_{n−k} (`DualView`, zero-copy).

## 2. Morphisms

A **chain map** (`ChainMap`) records, per stratum, where each cell of a
source complex went in a target complex: to a cell with an orientation
coefficient ±1, or to zero (`invalid_index`). These are the morphisms
induced on generators by every operation; cochains (fields, DoFs) are
transported by gathering through them.

Laws of morphisms:

- **L0 (composition).** `compose(f, g)` is associative, sends
  zero-to-zero, and multiplies orientation coefficients.
  *Witness:* `core.test_types`, and end-to-end in
  `ops.test_star_deletion` (chain_maps_compose_across_pushout_and_deletion).

## 3. Constructors

Complexes enter graphos three ways, all establishing ∂∘∂ = 0 by
construction:

- **Cell attachment** (`attach_cell`): a k-cell is attached along an
  explicit boundary chain — the CW paradigm.
- **Mesh ingestion** (`from_edges`, `from_polygons`, `from_polyhedra`,
  `from_simplices(d, …)` for any d ≥ 1): connectivity in the forms meshes
  arrive in; intermediate strata derived; deterministic orientation
  conventions (edges low→high vertex; faces by canonical cycle; simplices
  by sorted tuples with permutation parity). Inconsistent input (a face
  set described with cyclically incompatible orderings) is an error, not a
  guess. *Witness:* `core.test_build` (cube, mixed polygon mesh,
  pentachoron; ∂Δ⁴ ≅ S³).
- **Freezing** (`freeze`): the immutable query object; derives δ_k and
  moves storage into device-capable arrays. `Complex` is the builder,
  `FrozenComplex` the queryable epoch — topology-changing operations are
  out-of-place, producing new complexes plus morphisms (Principle P3).

## 4. The operation calculus

Each operation is a construction with a formal characterization; each
returns the induced chain map(s).

| Operation | Characterization | Law / witness |
|---|---|---|
| `disjoint_union(A,B)` | coproduct A ⊔ B | **L2**: χ and b₀ additive — `property.test_properties` (disjoint_union_is_additive) |
| `quotient(C, ~)` | quotient by cell identifications with orientation | sign propagation through stars — `ops.test_pushout` |
| `pushout(A,B,φ)` | A ⊔_C B along a vertex-generated shared subcomplex | orientation-reversed gluing keeps ∂∘∂ = 0 — `ops.test_pushout` |
| `lift_identifications` | extension of a vertex map through the strata by boundary-chain matching | square → cylinder / Möbius — `ops.test_lift_identifications` |
| `product(A,B)` | A × B, ∂(a×b) = ∂a×b + (−1)^{dim a} a×∂b | **L3**: χ(A×B) = χ(A)·χ(B); S¹×S¹ = T² — `ops.test_product` |
| `join(A,B)` | A ∗ B, ∂(a∗b) = ∂a∗b + (−1)^{dim a+1} a∗∂b | cone(S¹) = D², ΣS¹ = S² — `ops.test_join` |
| `star_deletion(C,S)` | excision of the open star st(S) | **L4 (excision)**: H_k(K, cl st S) ≅ H_k(K∖st S, frontier S), both sides computed — `queries.test_homology` |
| `subcomplex(C,S)` | the closed subcomplex cl(S), with embedding | embedding inverts the map on survivors — `ops.test_subcomplex` |
| `cut_along(C,S)` | cutting along a subcomplex; sides = components of st(x)∖cl(S) under full poset comparability | crack-tip topology; nonmanifold-sound — `ops.test_cut`, fuzzer |
| `replace(C,S,P,φ)` | excision + amalgamation fused (local surgery) | cut-cell imprint — `ops.test_replace` |
| `barycentric_subdivision` | the order complex of the face poset, with the SIGNED carrier: the incidence product ∏[cᵢ:cᵢ₊₁] along each flag (the subdivision chain-map coefficient; nonzero exactly on full flags) | **L5 (homotopy invariance)**: Betti preserved — `property.test_properties`. **L8 (orientation transfer)**: a consistent parent orientation, pushed through the signed carrier, is a consistent orientation of the subdivision (interior sd facets cancel — within cells by ∂∘∂ = 0 and the diamond property, across cells by the parent's consistency) — `ops.test_subdivision`, fuzzer |
| `collapse` / `free_faces` | Whitehead collapses at free faces (unique coface, single occurrence) | **L5**: simple homotopy equivalence; Möbius ↘ S¹ — `ops.test_collapse` |
| `orient` | construction of a fundamental-chain orientation, or refutation | interior facets induced with opposite signs — `ops.test_orient` |
| `agglomerate(C,ℓ)` | cellular amalgamation of top cells (numerical name: agglomeration); ∂(Σcᵢ) = Σ∂cᵢ, interior cancellation | ±2 coefficient ⇒ inconsistent orientation, rejected — `ops.test_agglomerate` |
| `coarsen(C,ℓ,π)` | dimensional descent of amalgamation; barriers π declared (P1) | **L6 (functoriality)**: level-by-level = one-step, morphisms equal element-wise — `property.test_multilevel`, `ops.test_coarsen` |

**L1 (chain condition preservation).** Every operation maps chain
complexes to chain complexes. *Witness:* the property fuzzer asserts
∂∘∂ = 0 after every step of random operation sequences on random
complexes (`property.test_properties`), in addition to every unit suite.

## 5. Queries and decidability

Read-only questions live in `queries/`. Where a question is undecidable or
only partially decidable, the API says which proxy it computes:

- **Homology** is over Z₂ (`betti_numbers_z2`, absolute and relative).
  **L7 (Euler–Poincaré):** Σ(−1)^k b_k = χ, asserted for every random
  complex in the fuzzer. Caveat: b_k(Z₂) ≥ b_k(ℝ), equal absent 2-torsion.
- **Link classification** (`classify_vertex_links`) is *exact* for
  complexes of dimension ≤ 3 (links are curves and surfaces, classified
  completely by closedness, connectedness, χ); for higher-dimensional
  links it degrades to the homology profile — necessary, not sufficient,
  since sphere recognition is undecidable in high dimension.
- **Amalgamation criteria** (`queries/amalgamation.hpp`): `acyclic` names
  what Z₂ homology certifies — acyclicity, not contractibility (homology
  cannot see π₁); exact ball detection for 1-dimensional common
  boundaries. `excess_intersection` = cl σ ∩ cl τ ∖ cl(∂σ ∩ ∂τ), the
  witness of improper intersection.
- **Manifoldness** (`check_manifold`) states necessary conditions
  (purity, facet condition, link connectivity), never claims sufficiency.

Terminological discipline: *free face* means the Whitehead condition
(exactly one coface, single occurrence — `free_faces`); a cell with *no*
proper coface is **maximal** (`FacetClassification::maximal`). The two are
never conflated.

## 6. Principles

- **P1 — metric-freeness.** No coordinates, measures, or geometric
  predicates exist in graphos. Every geometric decision enters as declared
  data: `Identification` lists (which cells are the same),
  `Marker`s (which cells are selected/protected), labels (how cells
  aggregate). Corners, features, and cut qualities are caller knowledge.
- **P2 — collective semantics.** Every public operation is specified
  collectively (Collective / Logically collective / Local); the serial
  implementation is the P = 1 case, and distribution is an implementation
  detail of `freeze()` and the ops, never an API change. See README,
  “Collective semantics”.
- **P3 — out-of-place epochs.** Complexes are never mutated by operations;
  each operation yields a new complex and the induced morphisms. History
  is carried by chain maps, and the efficient inverse of an operation is
  the remembered morphism, not a recomputation.
- **P4 — no silent topology change.** An operation either satisfies its
  stated hypotheses, rejects loudly (inconsistent orientation, ill-posed
  glue, malformed input), or falls back detectably (coarsening patches
  degrade to singletons; every fallback is observable in the morphism).
- **P5 — stated proxies.** Where the exact question is undecidable, the
  implemented proxy is named for what it is (acyclic, homology profile,
  necessary condition), in the API and its documentation.
