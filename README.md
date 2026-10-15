# SLIM: Spectral Laplacian Interlacing Matcher

A C++17 implementation of *SLIM: Exact Subgraph Matching via Spectral Dominance*:
`slim_index` builds the spectral signatures and the multi-rank dominance index
offline, `slim_query` runs spectral filtering followed by exact enumeration
online.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Artifacts: `build/slim_index` (offline) and `build/slim_query` (online);
`build/slim_offline` and `build/slim_online` are compatibility symlinks to them.

## Index

`slim_index` writes one index file per data graph: the local spectral
signatures, the label histograms, the `(label, layer, dimension)` partitions and
the adaptive multi-rank (AMSD) retrieval section the online stage runs on.

That section indexes, per label partition, the threshold families of the top
spectral coordinates of each signature: for rank `j` the family is a ladder of
decreasing thresholds, and rung `b` labels the vertices whose `j`-th largest
stored eigenvalue is at least that threshold. A query threshold that falls
between two rungs uses the next rung down, which can only keep extra vertices
that the exact test then removes. The `j`-th coordinate inequality is a
necessary condition of spectral dominance for every `j`, so a query vertex can
intersect any subset of the ranks, in any order, and stop at any time; the exact
dominance test decides the final candidate set, which therefore does not depend
on how much of the retrieval ran before it. Besides the rank families the
section carries each partition's degrees (a contiguous read instead of one
random access into the data graph) and the dimension condition
`|sigma(v)| >= |sigma(u)|`, materialized from the partition metadata on demand
and therefore free in storage.

For every query vertex the online stage picks the mechanism its cost model
expects to pay for: the threshold bitmaps, a per-candidate rank test inside the
partition scan, or the plain partition scan. The index is built with up to 6
ranks and up to 32 rungs per rank; `slim_index --verbose` prints the ranks, the
rungs and the bitmap bytes it produced.

Index files carry format version 2. A file written by an older build has no
multi-rank section and is refused at load time (`the index has no multi-rank
index section`): rebuild it with this `slim_index`.

## Contents

```text
CMakeLists.txt                    build script for slim_index and slim_query
include/slim/*.h                  library headers
src/*.cpp                         library sources and the two stage entry points
src/verify_engine/                enumeration engines (IVE default; see below)
datasets/real/*.graph             9 real-world labeled graphs
datasets/synthetic/*.graph        3 synthetic graphs with 100K vertices
queries/<dataset>/s<size>/        query graphs, 1710 files in total
```

Both graph and query files use the CHIP text format: a `t <name> <n> <m>` header,
one `v <id> L<label>` line per vertex, one `e <u> <v>` line per edge.

## Datasets

Nine real-world labeled graphs:

| Dataset | \|V\| | \|E\| | \|Σ\| | d_avg | File |
|---|---:|---:|---:|---:|---:|
| yeast | 3,112 | 12,519 | 71 | 8.05 | 169 KB |
| human | 4,674 | 86,282 | 44 | 36.92 | 972 KB |
| hprd | 9,460 | 34,998 | 307 | 7.40 | 486 KB |
| wordnet | 82,670 | 120,399 | 5 | 2.91 | 2.5 MB |
| dblp | 317,080 | 1,049,866 | 15 | 6.62 | 19 MB |
| youtube | 1,134,890 | 2,987,624 | 25 | 5.27 | 56 MB |
| us_patents | 3,774,768 | 16,518,947 | 20 | 8.75 | 324 MB |
| reddit | 232,965 | 11,606,919 | 41 | 99.65 | 170 MB |
| ogbn_products | 2,449,029 | 61,859,012 | 47 | 50.52 | 1,010 MB |

`wordnet` keeps 5,817 isolated vertices and `reddit` keeps 1,522, so their `|V|`
includes them. The first seven graphs are the labeled graphs of the in-memory
subgraph matching benchmark of Sun and Luo (SIGMOD 2020); `reddit` is the
GraphSAINT/SGC version of the Reddit graph; `ogbn_products` is the OGB graph.

Three synthetic graphs, one per label distribution:

| Dataset | Label distribution | File |
|---|---|---|
| Syn-Uniform | Uniform over `[1, \|Σ\|]` | `datasets/synthetic/e4l_n100k_d20_nws_l20_uni_s0.graph` |
| Syn-Gaussian | Gaussian over `[1, \|Σ\|]` | `datasets/synthetic/e4l_n100k_d20_nws_l20_gau_s0.graph` |
| Syn-Zipf | Zipf over `[1, \|Σ\|]` | `datasets/synthetic/e4l_n100k_d20_nws_l20_zipf_s0.graph` |

The three share one Newman-Watts-Strogatz topology (ring degree 8 plus random
shortcuts, structure seed 20260815): `|V| = 100,000`, `|E| = 1,000,000`
(`d_avg = 20`), `|Σ| = 20`; only the label assignment differs, so the label
skew is varied while the structure stays fixed. File names are the identifiers
of the experiment registry (`100k` vertices, `d20` average degree, `nws` model,
`l20` labels, `s0` replica), which keeps every number traceable to the runs
reported in the paper.

## Queries

All queries are connected induced subgraphs sampled with a random walk with
restart (restart probability 0.15) and deduplicated by vertex set and
isomorphism:

- The nine real graphs have 180 queries each: query size
  `|V(q)| ∈ {8, 16, 24, 32, 48, 64}` (6 sizes) × three classes (sparse / medium /
  dense, 10 queries each). The class of a query is decided by the edge count of
  the sampling pool of its `(dataset, size)` pair (bottom 30% / middle 40% / top
  30%), so `dense` queries are the structurally hardest ones. That is 1,620
  files under `queries/<dataset>/s<size>/`.
- The three synthetic graphs have 30 queries each, `|V(q)| = 32`, again 10
  sparse, 10 medium and 10 dense, for 90 files in total.

Query files are named `<dataset>_s<size>_<class>_<NN>.query`, e.g.
`queries/human/s16/human_s16_dense_03.query`.

## Usage examples

```bash
# Offline: build the dominance-aware index (paper defaults: alpha=1.0, k=1, tau_scale=1.0)
./build/slim_index -i datasets/real/yeast.graph -o yeast.index

# Online: a single query
./build/slim_query -i yeast.index -d datasets/real/yeast.graph \
    -q queries/yeast/s8/yeast_s8_sparse_00.query

# Online: batch, metrics written to CSV (the 30 yeast s8 queries take a few
# seconds in total; the slowest query alone takes seconds)
find queries/yeast/s8 -name '*.query' | sort > queries.txt
./build/slim_query -i yeast.index -d datasets/real/yeast.graph \
    --query-list queries.txt --output-csv metrics.csv

# Online: stop after the first match (full enumeration of s64 dense is slow)
./build/slim_query -i yeast.index -d datasets/real/yeast.graph \
    -q queries/yeast/s64/yeast_s64_dense_00.query \
    -l 1

# Synthetic graph: same two stages
./build/slim_index -i datasets/synthetic/e4l_n100k_d20_nws_l20_uni_s0.graph \
    -o syn_uniform.index
./build/slim_query -i syn_uniform.index \
    -d datasets/synthetic/e4l_n100k_d20_nws_l20_uni_s0.graph \
    -q queries/e4l_n100k_d20_nws_l20_uni_s0/s32/e4l_n100k_d20_nws_l20_uni_s0_s32_sparse_00.query
```

Query sizes range from 8 to 64; the larger and denser ones take longer — a size
64 dense query can enumerate for minutes on the larger graphs. Use `--limit` to
cap enumeration, or `--filter-only` when only the filtering stage is of
interest. Other datasets work the same way: replace `-d` with
`datasets/real/<dataset>.graph` (or the synthetic file) and point `-q` /
`--query-list` at the queries under `queries/<dataset>/`.

## Metrics CSV

`--output-csv <file>` (and `--metrics` for stdout) writes one row per query:

| Column | Meaning |
|---|---|
| `query_path` | query file of the row |
| `end_to_end_us` | wall clock of the query, from reading the query file to the end of the run |
| `match_core_us` | online time of the matcher: covers `query_spectral_us`, `filter_us` and `verify_us` and the glue between them; this is the per-query online time the paper reports |
| `index_load_us`, `data_load_us`, `query_load_us` | loading of the index, the data graph and the query graph |
| `query_spectral_us` | query-side signatures (`Tsig`) |
| `filter_us` | spectral filtering and candidate propagation (`Tfilter`) |
| `verify_us` | exact enumeration (`Tenum`) |
| `verify_propagation_us`, `verify_adjacency_us`, `verify_enumeration_us` | the three stages of `verify_us` |
| `candidates_after_filter` | `Σ_u \|C(u)\|` after filtering |
| `candidates_after_propagation` | `Σ_u \|C(u)\|` after propagation |
| `match_count` | embeddings found (capped by `--limit`) |
| `early_stopped` | 1 when enumeration stopped at the match limit |
| `coarse_dominance_checks`, `full_dominance_checks` | pairs that entered the index-level retrieval and the full dominance test |
| `data_vertices`, `query_vertices`, `query_edges` | sizes of the data graph and the query graph |

## Command-line reference

### `slim_index`

| Option | Meaning | Default |
|---|---|---|
| `-i, --input` | input graph file | required |
| `-o, --output` | output index file | required |
| `-a, --alpha` | label weight, must be `> 0` | `1.0` |
| `-k, --k` | neighborhood radius k | `1` |
| `--tau-scale c` | `tau = c * lambda_max / sqrt(\|V\|)` | `1.0` |
| `--neighborhood-size-limit N` | vertices whose neighborhood is larger than N skip the eigen-decomposition and are kept conservatively; `0` disables the guard | `2048` |
| `--threads N` | worker threads for the signature pass | `1` |
| `-v, --verbose` | verbose output | off |
| `-h, --help` | print usage | — |

### `slim_query`

| Option | Meaning |
|---|---|
| `-i, --index` | SLIM index file (required) |
| `-d, --data` | data graph file (required) |
| `-q, --query` | single query file (exactly one of `-q` / `--query-list`) |
| `--query-list` | file with one query path per line (exactly one of `-q` / `--query-list`) |
| `-l, --limit` | stop after N matches; `0` = enumerate all (default `0`) |
| `--output-csv` | write a machine-readable metrics CSV |
| `--metrics` | print the metrics CSV to stdout |
| `--filter-only` | stop after spectral filtering (for measurement) |
| `--trace-dominance` | also print coarse/full dominance test counts to stderr |
| `-v, --verbose` | verbose output |
| `-h, --help` | print usage |

`alpha`, `k` and `tau` are always taken from the index; the online stage has no
algorithm parameters. `--threads` only affects offline build speed: an index
written with any N is bit-identical to the single-threaded one.

## Third-party code

`src/verify_engine/` holds the enumeration engines of the online stage. The
command-line interface drives the default one (IVE); the other two are
library-level alternatives (`slim::VerifierBackend`) that the experiments used
for cross-checking.

| Engine | Directory | Provenance |
|---|---|---|
| IVE | `src/verify_engine/ive/` | vendored from IVE (ICDE 2024, "IVE: Accelerating Enumeration-based Subgraph Matching via Exploring Isolated Vertices"); the upstream project ships no license file, so this copy is for local research use only |
| VEQ | `src/verify_engine/veq/` | vendored from VEQ; the upstream repository publishes prebuilt binaries only and ships no license file, so this copy is for local research use only |
| LFTJ | `src/verify_engine/lftj/` | trimmed copy of the LFTJ engine of SubgraphMatching (SIGMOD 2020), MIT license, Copyright (c) 2021 Rapids @ HKUST; the license text is in `src/verify_engine/lftj/LICENSE` and the provenance of every file is in `src/verify_engine/lftj/README.md` |
