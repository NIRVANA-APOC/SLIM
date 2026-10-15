#ifndef SLIM_VERIFY_ENGINE_IVE_VERIFIER_H
#define SLIM_VERIFY_ENGINE_IVE_VERIFIER_H

#include <slim/backtracking_verifier.h>
#include <slim/labeled_graph.h>
#include <slim/types.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

#include "ive/BuildTable.h"
#include "ive/enumeration.h"
#include "ive/graph.h"
#include "ive/reordering.h"

namespace slim {
class IveVerifier : public EnumerationBackend {
public:
    explicit IveVerifier(const LabeledGraph& data_graph);
    ~IveVerifier() override;

    IveVerifier(const IveVerifier&) = delete;
    IveVerifier& operator=(const IveVerifier&) = delete;

    Result verify(const LabeledGraph& query_graph,
                  std::vector<std::vector<VertexId>>& candidates,
                  std::size_t limit) override;

private:
    struct Impl {
        ::Graph data_graph;
        ::Graph query_graph;
        ui previous_query_size = 0;
        std::vector<ui*> candidate_ptrs;
        std::vector<ui> candidate_counts;
        std::unique_ptr<bool[]> component_roots;
    };

    std::unique_ptr<Impl> impl_;
};

namespace ive_detail {
static_assert(sizeof(ui) == sizeof(VertexId), "IVE ui must match slim::VertexId");

inline std::uint64_t micros_since(std::chrono::high_resolution_clock::time_point start) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start)
            .count());
}

inline void fill_ive_graph(::Graph& out, const LabeledGraph& graph) {
    const std::size_t n = graph.num_vertices();
    std::vector<ui> offsets(n + 1, 0);
    for (std::size_t v = 0; v < n; ++v) {
        offsets[v + 1] = offsets[v] + static_cast<ui>(graph.degree(static_cast<VertexId>(v)));
    }
    std::vector<ui> neighbors(offsets[n]);
    std::vector<ui> labels(n);
    for (std::size_t v = 0; v < n; ++v) {
        std::size_t count = 0;
        const VertexId* row = graph.neighbors(static_cast<VertexId>(v), count);
        std::memcpy(neighbors.data() + offsets[v], row, count * sizeof(ui));
        labels[v] = static_cast<ui>(graph.label(static_cast<VertexId>(v)));
    }
    out.buildFromArrays(static_cast<ui>(n), static_cast<ui>(graph.num_edges()),
                        offsets.data(), neighbors.data(), labels.data());
}
}

inline IveVerifier::IveVerifier(const LabeledGraph& data_graph) : impl_(new Impl) {
    ive_detail::fill_ive_graph(impl_->data_graph, data_graph);
}

inline IveVerifier::~IveVerifier() {
    if (impl_->previous_query_size != 0) {
        freeBsrGraph(bsr_graph_, impl_->previous_query_size);
    }
}

inline IveVerifier::Result IveVerifier::verify(const LabeledGraph& query_graph,
                                               std::vector<std::vector<VertexId>>& candidates,
                                               std::size_t limit) {
    Result result;
    const std::size_t query_size = query_graph.num_vertices();
    if (query_size == 0 || candidates.size() != query_size) {
        return result;
    }
    for (const auto& set : candidates) {
        if (set.empty()) {
            return result;
        }
    }

    Impl& impl = *impl_;
    ive_detail::fill_ive_graph(impl.query_graph, query_graph);
    const ui qn = static_cast<ui>(query_size);

    impl.candidate_ptrs.resize(query_size);
    impl.candidate_counts.resize(query_size);
    for (std::size_t u = 0; u < query_size; ++u) {
        impl.candidate_ptrs[u] = candidates[u].data();
        impl.candidate_counts[u] = static_cast<ui>(candidates[u].size());
    }

    if (impl.previous_query_size != 0) {
        freeBsrGraph(bsr_graph_, impl.previous_query_size);
    }
    impl.previous_query_size = qn;

    ive_reset_query_state(qn);

    ui* order = nullptr;
    impl.component_roots = std::make_unique<bool[]>(query_size);
    MDE(&impl.data_graph, &impl.query_graph, order, impl.candidate_counts.data(),
        impl.component_roots.get());

    const auto t_tables = std::chrono::high_resolution_clock::now();
    Edges*** edge_matrix = new Edges**[query_size];
    for (std::size_t i = 0; i < query_size; ++i) {
        edge_matrix[i] = new Edges*[query_size];
        for (std::size_t j = 0; j < query_size; ++j) {
            edge_matrix[i][j] = nullptr;
        }
    }
    buildTables(&impl.data_graph, &impl.query_graph, impl.candidate_ptrs.data(),
                impl.candidate_counts.data(), edge_matrix, order);
    freeEdgeMatrix(edge_matrix, qn);
    for (std::size_t i = 0; i < query_size; ++i) {
        delete[] edge_matrix[i];
    }
    delete[] edge_matrix;
    qfliter_bsr_graph_ = bsr_graph_;
    result.adjacency_us = ive_detail::micros_since(t_tables);

    const auto t_search = std::chrono::high_resolution_clock::now();
    std::size_t call_count = 0;
    const std::size_t engine_limit =
        (limit == 0) ? std::numeric_limits<std::size_t>::max() : limit;
    result.count = explore(&impl.data_graph, &impl.query_graph, impl.candidate_ptrs.data(),
                           impl.candidate_counts.data(), engine_limit, call_count, order,
                           impl.component_roots.get());
    result.enumeration_us = ive_detail::micros_since(t_search);

    result.early_stopped = (limit > 0) && (result.count >= limit);
    if (limit > 0 && result.count > limit) {
        result.count = limit;
    }

    delete[] order;
    return result;
}
}

#endif
