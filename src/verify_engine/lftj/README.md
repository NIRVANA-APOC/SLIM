# Vendored LFTJ engine (SubgraphMatching)

Enumeration backend of the online stage. The code in
this directory is a trimmed, namespace-isolated copy of the **LFTJ** engine of

* upstream: `https://github.com/RapidsAtHKUST/SubgraphMatching`
  (SIGMOD 2020, "In-Memory Subgraph Matching: An In-Depth Study"); the file
  identity is pinned by the sha256 sums below
* license: **MIT**, Copyright (c) 2021 Rapids @ HKUST — see `LICENSE` in this
  directory, copied verbatim from the upstream tree. This is the only vendored
  engine in `src/verify_engine/` that carries an explicit license.
* paper of the algorithm: SIGMOD 2020, Section 3.3 — LFTJ, the worst-case
  optimal join over the edge-relation tables `C(u, u')`. The SubgraphMatching
  implementation is the one studied in that paper; SLIM drives it with the
  GQL plan of He & Singh, SIGMOD 2008 (see `generate_query_plan.h`), which is
  what the upstream study uses for the LFTJ engine as well.

## Upstream files and their hashes

| upstream file | sha256 (upstream) |
| ------------- | ----------------- |
| `LICENSE` | `d68a769a91f332e5bf6acaeb78f0bffc3a2e936103c7151dacfce1c34c5ae122` |
| `configuration/types.h` | `f6cbbaf48091fd31fe4997bd9a73679a7a9f959c2e2625b219e8a4a0cdbc8761` |
| `graph/graph.h` | `ea3252c3ccc7e80d5b8427e6a26514a36974e09cb117de1cfeb420bae57394c0` |
| `graph/graph.cpp` | `a9b40bd3bce873f0967e8e20c82c86705429ebab1d1e871624cf86f615b6f033` |
| `matching/BuildTable.h` | `53fbd07a9126ed883ef58f12765e21384690f054e08a4d388574e4feacaf823c` |
| `matching/BuildTable.cpp` | `a075ed47b452907e8a267c0cfc9389e28f0f3210f3021e2dabc9aff8f1cec21a` |
| `matching/GenerateQueryPlan.h` | `09909a00d2a014255dafc4a65583ec1aef96cd2b8c826e97ee5b7205b4abfb68` |
| `matching/GenerateQueryPlan.cpp` | `147864837e63d675008154fd88ba1b9b029fd3623f414117a78e355609a1b2a9` |
| `matching/EvaluateQuery.h` | `74a0f592508939501a9a5c573960658affb11dd1f01f14d2d7a7e747c4a720a5` |
| `matching/EvaluateQuery.cpp` | `66e038865a98e6bf578dbcab9ab337e64946dcf1bfa1d867384d2d28c849ca29` |
| `utility/computesetintersection.h` | `0cd24414535fcae49841a98fc9f34ecfc4f0b70dd435e4495ed5cb5696455a15` |
| `utility/computesetintersection.cpp` | `210ccb3ea0393a669552794b6b1fcf04fb150655d27a816b9681e9488c470626` |

## File map

| file here | upstream origin | content |
| --------- | --------------- | ------- |
| `graph.h` | `graph/graph.h` + `graph/graph.cpp`, `configuration/types.h` | `Graph` (CSR), `ui` / `VertexID` / `LabelID`, `Edges` |
| `build_table.h` | `matching/BuildTable.{h,cpp}` | `BuildTable::buildTables` |
| `generate_query_plan.h` | `matching/GenerateQueryPlan.{h,cpp}` | `GenerateQueryPlan::generateGQLQueryPlan` |
| `evaluate_query.h` | `matching/EvaluateQuery.{h,cpp}` | `EvaluateQuery::LFTJ` and its four helpers |
| `computesetintersection.h` | `utility/computesetintersection.{h,cpp}` | `ComputeSetIntersection::ComputeCandidates` and its AVX2 helpers |
| `LICENSE` | `LICENSE` | upstream MIT license, verbatim |
| `../lftj_verifier.h` | — (SLIM code) | the `EnumerationBackend` wrapper that injects C(u) |

## SLIM-side changes

Whole-directory:

1. **Namespace isolation.** Everything is wrapped in `namespace slim_lftj`.
   Upstream defines its symbols in the global namespace (`Graph`, `Edges`, `ui`,
   `VertexID`, ...); all three backends are compiled into the single translation
   unit `src/backtracking_verifier.cpp`, where IVE's symbols are global.
2. **Header-only.** Every function is defined in these headers (upstream's
   `.cpp` bodies moved in) and marked `inline`. The engines are included by
   exactly one translation unit, and `CMakeLists.txt` stays untouched.
3. **Deleted:** `main`, the command-line/configuration machinery, all file IO
   (graph loaders, compressed loaders, store, the k-core table that pulls in
   `GraphOperations`), all printing, the filter and the other engines (`GQL`,
   `QSI`, `DPiso`, `CECI`, `Spectrum`), the `DISTRIBUTION` /
   `ENABLE_FAILING_SET` / `SPECTRUM` / `ENABLE_QFLITER` variants, and the
   `configuration/config.h` macros they select on (`HYBRID`, `SI`: the vendored
   set intersection hardcodes the upstream defaults `HYBRID == 0`, `SI == 0`).
   **Minimal closure:** the LFTJ path needs neither `utility/graphoperations` nor
   `utility/bitsetoperation`, and no `immintrin`-free fallback — the pilot
   recorded this (equivalence report, sec. 1).
4. **Preprocessor directives dropped inside the vendored bodies:** the
   conditional blocks listed above are resolved at vendoring time, so the
   remaining lines are exactly the branch upstream compiles by default.

`graph.h`:

5. Trimmed to the CSR and the accessors the LFTJ path calls:
   `getVerticesCount`, `getEdgesCount`, `getVertexDegree`,
   `getVertexNeighbors`, `getVertexLabel`, `checkEdgeExistence` (upstream's
   label-free variant, kept verbatim including its swap). Dropped the reverse
   index, the k-core table, the NLF / label-offset structures
   (`OPTIMIZED_LABELED_GRAPH`) and the `enable_label_offset_` flag. The
   accessor set is identical to the previously reviewed GQL copy, so the two
   files differ only in the namespace, the header guard and the `Edges` class.
6. Added `build_from_csr(...)` (marked `SLIM integration`) so the engine can be
   fed SLIM's in-memory CSR; upstream only builds a `Graph` by reading a file.
   Vertex identifiers are preserved one to one, which is what makes the injected
   C(u) meaningful to the engine.
7. `Edges` (upstream `configuration/types.h`) is vendored here: `offset_`,
   `edge_`, `vertex_count_`, `edge_count_`, `max_degree_` and the destructor
   that frees both arrays. Upstream's other `types.h` entities (`TreeNode`, the
   `MatchingIndexType` enum) belong to engines that are not vendored.

`build_table.h`:

8. `BuildTable::buildTables` is upstream's body with three documented edits (the
   `ENABLE_QFLITER` block dropped, as everywhere):
   * the two per-call `flag` / `updated_flag` arrays are freed before returning.
     Upstream never frees them (it is written for one query per process), so an
     online process answering many queries would leak `2 * |V| ui` per query
     (ASan/LeakSanitizer on the hand-written cases; ~30 MB per query on
     us_patents).
   * every `edge_` array is allocated with `kAvx2Slack` extra elements (see
     `computesetintersection.h`).
   The printing / memory-estimation helpers, the `QFliter` include and the
   `Tool`/`print` code are not vendored.

`generate_query_plan.h`:

9. Carried over from the GQL backend that this engine replaces (same upstream
   function, same patch), with the namespace renamed and the unused `pivot`
   output removed. **Multi-component fix (SLIM integration).** Upstream walks
   the greedy rule "least candidates among the vertices adjacent to the selected
   set" and reads `next_vertex` uninitialized when no unvisited vertex is
   adjacent to the selected set — i.e. for every disconnected query graph, and
   for an isolated query vertex that is not the start vertex. The pilot measured
   a garbage order (`0,1,1` for a 3-vertex query) and a SIGSEGV
   (`triangle_plus_isolated_vertex`); the vendored copy initializes
   `next_vertex = 0`, records whether the greedy found anything (`found`) and
   falls back to `selectGQLComponentRoot`, which applies the same key (fewest
   candidates, then higher degree, then smaller id) restricted to unvisited
   vertices.
10. `ui` → `VertexID` and redundant `const` on value returns removed; the body
    is otherwise upstream's, unindented to this repository's style.

`evaluate_query.h`:

11. `LFTJ`, `generateBN`, `allocateBuffer` and `releaseBuffer` are upstream's
    bodies apart from the namespace, the dropped macro blocks,
    `size_t` → `std::size_t` in the signature, and one added `kAvx2Slack` on the
    buffers `allocateBuffer` allocates (the AVX2 intersection reads up to
    `kAvx2Slack` elements past the logical end of an input). `generateBN` is the
    3-argument overload (earlier neighbours only), because the plan no longer
    produces `pivot`.
12. **Multi-component fix (SLIM integration).**
    `generateValidCandidateIndex` takes `candidates_count` and
    starts a vertex without an earlier neighbour (a component root, or an
    isolated query vertex) from its full candidate list, exactly like depth 0.
    Upstream reads `bn[depth][0]` unconditionally there, which is out of range:
    the pilot measured a SIGSEGV as soon as the plan patch made such an order
    reachable (equivalence report, sec. 1.2).
13. `output_limit_num == 0` stops after the first match upstream, so the wrapper
    passes `SIZE_MAX` for "unlimited"; the search checks the limit before
    extending, so it stops exactly at `limit`.

`computesetintersection.h`:

14. `ComputeCandidates` is upstream's default-path body with the `HYBRID` / `SI`
    dispatch resolved to the two AVX2 implementations (`ComputeCNGallopingAVX2`
    for very skewed inputs, `ComputeCNMergeBasedAVX2` otherwise) plus the
    `GallopingSearchAVX2` / `BinarySearchForGallopingSearchAVX2` helpers. The
    AVX512 and scalar variants and the count-only overloads are not vendored
    (LFTJ calls only the AVX2 pair). Requires `-mavx2`-capable hardware, which
    the release flags (`-march=native`) and the rest of the engine already
    assume.
15. **`kAvx2Slack` (SLIM integration).** The AVX2 kernels load full 8-element
    vectors and read up to 7 elements past the logical end of an input
    (`ComputeCNMergeBasedAVX2` merges 8 at a time; `GallopingSearchAVX2` loads 8
    elements even when fewer remain). Upstream allocates exactly sized arrays
    and gets away with it only because the allocator rounds up; ASan caught the
    read past the end of an edge table on dblp. The vendored copy defines
    `constexpr ui kAvx2Slack = 8` and `build_table.h` / `evaluate_query.h`
    allocate every array the intersection can read with that much slack.
16. **Two UBSan/`-Wsign-compare` cleanups (SLIM integration).** `int real_mask`
    is cast to `unsigned int` before `<< 16` (UBSan: "left shift of ... cannot be
    represented in type 'int'"), and the `mask != 0xffffffff` tests are written
    `mask != -1` (same bits, no signed/unsigned comparison warning). The
    arithmetic is unchanged.

## Restrictions

None for query shape: connected, disconnected, isolated query vertices,
isolated data vertices and empty candidate sets are all supported (see
the engine equivalence report for the hand-written cases).

The engine allocates its per-query buffers with `new`/`delete`: the edge tables
(`buildTables`) and the search scratch (`allocateBuffer`, including one
`bool[|V|]`). The wrapper releases everything before returning, so nothing is
kept alive across queries, but the allocation cost is part of the per-query
timing (upstream behaviour, kept as-is).
