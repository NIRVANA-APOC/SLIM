#ifndef SLIM_VERIFY_ENGINE_LFTJ_VERIFIER_H
#define SLIM_VERIFY_ENGINE_LFTJ_VERIFIER_H

#include <slim/backtracking_verifier.h>
#include <slim/labeled_graph.h>
#include <slim/types.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "lftj/build_table.h"
#include "lftj/evaluate_query.h"
#include "lftj/generate_query_plan.h"
#include "lftj/graph.h"

namespace slim {
class LftjVerifier : public EnumerationBackend {
public:
    explicit LftjVerifier(const LabeledGraph& data_graph);
    ~LftjVerifier() override = default;

    LftjVerifier(const LftjVerifier&) = delete;
    LftjVerifier& operator=(const LftjVerifier&) = delete;

    Result verify(const LabeledGraph& query_graph,
                  std::vector<std::vector<VertexId>>& candidates,
                  std::size_t limit) override;

private:
    struct Impl {
        slim_lftj::Graph data;
        std::vector<slim_lftj::ui> offsets;
        std::vector<slim_lftj::VertexID> neighbors;
        std::vector<slim_lftj::LabelID> labels;
        std::vector<std::vector<slim_lftj::ui>> candidate_storage;
        std::vector<slim_lftj::ui*> candidate_ptrs;
        std::vector<slim_lftj::ui> candidate_counts;
    };

    std::unique_ptr<Impl> impl_;
};

namespace lftj_detail {
using Clock = std::chrono::high_resolution_clock;

inline std::uint64_t micros_since(Clock::time_point start, Clock::time_point end) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
}

inline void build_csr(const LabeledGraph& source, slim_lftj::Graph& target,
                      std::vector<slim_lftj::ui>& offsets,
                      std::vector<slim_lftj::VertexID>& neighbors,
                      std::vector<slim_lftj::LabelID>& labels) {
    const std::size_t num_vertices = source.num_vertices();
    offsets.resize(num_vertices + 1);
    labels.resize(num_vertices);
    neighbors.clear();
    neighbors.reserve(2 * source.num_edges());

    slim_lftj::ui max_degree = 0;
    std::size_t cursor = 0;
    for (VertexId v = 0; v < num_vertices; ++v) {
        offsets[v] = static_cast<slim_lftj::ui>(cursor);
        std::size_t degree = 0;
        const VertexId* const row = source.neighbors(v, degree);
        for (std::size_t i = 0; i < degree; ++i) {
            neighbors.push_back(row[i]);
        }
        cursor += degree;
        if (degree > max_degree) {
            max_degree = static_cast<slim_lftj::ui>(degree);
        }
        labels[v] = source.label(v);
    }
    offsets[num_vertices] = static_cast<slim_lftj::ui>(cursor);

    slim_lftj::ui labels_count = 0;
    for (const slim_lftj::LabelID label : labels) {
        if (label + 1 > labels_count) {
            labels_count = label + 1;
        }
    }

    target.build_from_csr(static_cast<slim_lftj::ui>(num_vertices),
                          static_cast<slim_lftj::ui>(neighbors.size() / 2), labels_count,
                          max_degree, offsets, neighbors, labels);
}

inline slim_lftj::Edges*** allocate_edge_matrix(std::size_t query_size) {
    slim_lftj::Edges*** edge_matrix = new slim_lftj::Edges**[query_size];
    for (std::size_t i = 0; i < query_size; ++i) {
        edge_matrix[i] = new slim_lftj::Edges*[query_size];
    }
    return edge_matrix;
}

inline void release_edge_matrix(slim_lftj::Edges*** edge_matrix, std::size_t query_size) {
    for (std::size_t i = 0; i < query_size; ++i) {
        for (std::size_t j = 0; j < query_size; ++j) {
            delete edge_matrix[i][j];
        }
        delete[] edge_matrix[i];
    }
    delete[] edge_matrix;
}
}

inline LftjVerifier::LftjVerifier(const LabeledGraph& data_graph) : impl_(new Impl) {
    lftj_detail::build_csr(data_graph, impl_->data, impl_->offsets, impl_->neighbors,
                           impl_->labels);
}

inline EnumerationBackend::Result LftjVerifier::verify(
    const LabeledGraph& query_graph, std::vector<std::vector<VertexId>>& candidates,
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

    const auto t_adjacency_start = lftj_detail::Clock::now();

    slim_lftj::Graph query;
    {
        std::vector<slim_lftj::ui> query_offsets;
        std::vector<slim_lftj::VertexID> query_neighbors;
        std::vector<slim_lftj::LabelID> query_labels;
        lftj_detail::build_csr(query_graph, query, query_offsets, query_neighbors, query_labels);
    }

    impl.candidate_storage.resize(query_size);
    impl.candidate_ptrs.resize(query_size);
    impl.candidate_counts.resize(query_size);
    for (std::size_t u = 0; u < query_size; ++u) {
        std::vector<slim_lftj::ui>& set = impl.candidate_storage[u];
        set.assign(candidates[u].begin(), candidates[u].end());
        impl.candidate_ptrs[u] = set.data();
        impl.candidate_counts[u] = static_cast<slim_lftj::ui>(set.size());
    }

    slim_lftj::Edges*** edge_matrix = lftj_detail::allocate_edge_matrix(query_size);
    slim_lftj::BuildTable::buildTables(&impl.data, &query, impl.candidate_ptrs.data(),
                                       impl.candidate_counts.data(), edge_matrix);

    slim_lftj::ui* order = nullptr;
    slim_lftj::GenerateQueryPlan::generateGQLQueryPlan(&impl.data, &query,
                                                       impl.candidate_counts.data(), order);

    result.adjacency_us =
        lftj_detail::micros_since(t_adjacency_start, lftj_detail::Clock::now());

    const std::size_t engine_limit =
        (limit == 0) ? std::numeric_limits<std::size_t>::max() : limit;
    std::size_t call_count = 0;
    const auto t_enum_start = lftj_detail::Clock::now();
    const std::size_t engine_count = slim_lftj::EvaluateQuery::LFTJ(
        &impl.data, &query, edge_matrix, impl.candidate_ptrs.data(),
        impl.candidate_counts.data(), order, engine_limit, call_count);
    result.enumeration_us =
        lftj_detail::micros_since(t_enum_start, lftj_detail::Clock::now());

    delete[] order;
    lftj_detail::release_edge_matrix(edge_matrix, query_size);

    result.count = engine_count;
    result.early_stopped = (limit > 0) && (result.count >= limit);
    return result;
}
}

#endif
