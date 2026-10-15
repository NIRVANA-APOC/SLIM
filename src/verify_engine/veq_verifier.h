#ifndef SLIM_VERIFY_ENGINE_VEQ_VERIFIER_H
#define SLIM_VERIFY_ENGINE_VEQ_VERIFIER_H

#include <slim/backtracking_verifier.h>
#include <slim/labeled_graph.h>
#include <slim/types.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "veq/veq_entry.h"

namespace slim {
class VeqVerifier : public EnumerationBackend {
public:
    explicit VeqVerifier(const LabeledGraph& data_graph);
    ~VeqVerifier() override;

    VeqVerifier(const VeqVerifier&) = delete;
    VeqVerifier& operator=(const VeqVerifier&) = delete;

    Result verify(const LabeledGraph& query_graph,
                  std::vector<std::vector<VertexId>>& candidates,
                  std::size_t limit) override;

private:
    struct Impl {
        std::unordered_map<int, int> label_ids;
        std::vector<std::vector<int>> candidate_sets;
        std::vector<const int*> candidate_ptrs;
        std::vector<int> candidate_sizes;
        VEQ::Graph* query = nullptr;
    };

    std::unique_ptr<Impl> impl_;
};

namespace veq_detail {
inline void check_query_supported(const LabeledGraph& query_graph) {
    const std::size_t query_size = query_graph.num_vertices();
    if (query_size > static_cast<std::size_t>(VEQ::MAX_NUM_VERTEX)) {
        throw std::runtime_error(
            "VEQ backend supports at most " + std::to_string(VEQ::MAX_NUM_VERTEX) +
            " query vertices; use the IVE backend instead");
    }
    std::size_t max_degree = 0;
    for (VertexId u = 0; u < query_size; ++u) {
        const std::size_t degree = query_graph.degree(u);
        if (degree > max_degree) {
            max_degree = degree;
        }
    }
    if (max_degree >= static_cast<std::size_t>(VEQ::MAX_QUERY_DEGREE)) {
        throw std::runtime_error(
            "VEQ backend supports query degree < " +
            std::to_string(VEQ::MAX_QUERY_DEGREE) + " only; use the IVE backend instead");
    }
    std::vector<char> seen(query_size, 0);
    std::vector<VertexId> stack;
    stack.push_back(0);
    seen[0] = 1;
    std::size_t reached = 1;
    while (!stack.empty()) {
        const VertexId u = stack.back();
        stack.pop_back();
        std::size_t degree = 0;
        const VertexId* neighbors = query_graph.neighbors(u, degree);
        if (degree == 0) {
            throw std::runtime_error(
                "VEQ backend does not support isolated query vertices; use the IVE backend instead");
        }
        for (std::size_t i = 0; i < degree; ++i) {
            const VertexId v = neighbors[i];
            if (seen[v] == 0) {
                seen[v] = 1;
                ++reached;
                stack.push_back(v);
            }
        }
    }
    if (reached != query_size) {
        throw std::runtime_error(
            "VEQ backend does not support disconnected query graphs (upstream DAG "
            "covers the root's component only); use the IVE backend instead");
    }
}
}

inline VeqVerifier::VeqVerifier(const LabeledGraph& data_graph) : impl_(new Impl) {
    impl_->label_ids = VEQ::BuildLabelIndex(data_graph);
    VEQ::init_veq_data(data_graph, impl_->label_ids);
}

inline VeqVerifier::~VeqVerifier() {
    VEQ::dispose_veq_query(impl_->query);
    impl_->query = nullptr;
    VEQ::dispose_veq_data();
}

inline EnumerationBackend::Result VeqVerifier::verify(
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
    veq_detail::check_query_supported(query_graph);

    Impl& impl = *impl_;

    VEQ::dispose_veq_query(impl.query);
    impl.query = nullptr;

    impl.candidate_sets.resize(query_size);
    impl.candidate_ptrs.resize(query_size);
    impl.candidate_sizes.resize(query_size);
    for (std::size_t u = 0; u < query_size; ++u) {
        const int query_label =
            static_cast<int>(query_graph.label(static_cast<VertexId>(u)));
        const auto label = impl.label_ids.find(query_label);
        if (label == impl.label_ids.end()) {
            return result;
        }
        std::vector<int>& set = impl.candidate_sets[u];
        set.assign(candidates[u].begin(), candidates[u].end());
        impl.candidate_ptrs[u] = set.data();
        impl.candidate_sizes[u] = static_cast<int>(set.size());
    }

    impl.query = new VEQ::Graph();
    VEQ::FillVeqGraph(*impl.query, query_graph, impl.label_ids,
                      static_cast<int>(impl.label_ids.size()));

    long long engine_count = 0;
    int64_t time_total_ns = 0;
    const bool exceeded =
        VEQ::run_veq_query(impl.query, impl.candidate_ptrs.data(),
                           impl.candidate_sizes.data(), limit,
                           false, engine_count, time_total_ns);
    if (exceeded) {
        throw std::runtime_error("VEQ backend hit its internal time limit");
    }

    result.adjacency_us = static_cast<std::uint64_t>(VEQ::filteringTime / 1000);
    result.enumeration_us = static_cast<std::uint64_t>(VEQ::verificationTime / 1000);

    result.count = (engine_count < 0) ? 0 : static_cast<std::size_t>(engine_count);
    result.early_stopped = (limit > 0) && (result.count >= limit);
    if (limit > 0 && result.count > limit) {
        result.count = limit;
    }
    return result;
}
}

#endif
