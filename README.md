# SLIM: Spectral-based Subgraph Matching

## Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Run

### Build Index
```bash
./offline/slim_offline \
    -i ../dataset/yeast.graph \
    -o yeast.index \
    -a 1.0
```

### Single Query
```bash
./online/slim_online \
    -i yeast.index \
    -d ../dataset/yeast.graph \
    -q ../dataset/queries/query_5_000.graph
```

### Batch Test (All 150 Queries)
```bash
for q in ../dataset/queries/*.graph; do
    ./online/slim_online -i yeast.index -d ../dataset/yeast.graph -q "$q"
done
```

## Output Format

Per-query output:
```
Query time: 2453 us (2.453 ms)
Matches found: 296114
```

## System Requirements
- C++17 compiler (GCC 9+, Clang 10+)
- CMake 3.15+
- OpenMP (optional, for offline parallelization)
