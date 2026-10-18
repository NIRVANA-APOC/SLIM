#include <slim/backtracking_verifier.h>

#include "verify_engine/ive_verifier.h"
#include "verify_engine/veq_verifier.h"
#include "verify_engine/lftj_verifier.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <utility>

namespace slim {
namespace {
std::uint64_t micros_since(std::chrono::high_resolution_clock::time_point start,
                           std::chrono::high_resolution_clock::time_point end) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
}
}

BacktrackingVerifier::BacktrackingVerifier(const LabeledGraph& data_graph,
                                           VerifierBackend backend)
    : data_graph_(data_graph), engine_(nullptr) {
    switch (backend) {
        case VerifierBackend::Veq:
            engine_ = std::make_unique<VeqVerifier>(data_graph);
            break;
        case VerifierBackend::Ive:
            engine_ = std::make_unique<IveVerifier>(data_graph);
            break;
        case VerifierBackend::Lftj:
            engine_ = std::make_unique<LftjVerifier>(data_graph);
            break;
    }
}

BacktrackingVerifier::~BacktrackingVerifier() = default;

void BacktrackingVerifier::propagate_edge_constraints(
    const LabeledGraph& query_graph,
    std::vector<std::vector<VertexId>>& candidates) {
    const std::size_t query_size = query_graph.num_vertices();
    const std::size_t words = (data_graph_.num_vertices() + 63) / 64;

    if (support_.size() < query_size * words) {
        support_.resize(query_size * words);
    }
    if (support_fresh_.size() < query_size) {
        support_fresh_.resize(query_size);
    }
    std::uint64_t* const support = support_.data();
    char* const fresh = support_fresh_.data();

    auto rebuild_support = [&](VertexId u) {
        std::uint64_t* row = support + static_cast<std::size_t>(u) * words;
        std::fill(row, row + words, 0);
        for (VertexId w : candidates[u]) {
            std::size_t degree = 0;
            const VertexId* neighbors = data_graph_.neighbors(w, degree);
            for (std::size_t j = 0; j < degree; ++j) {
                const VertexId x = neighbors[j];
                row[x >> 6] |= (std::uint64_t{1} << (x & 63));
            }
        }
    };

    constexpr int kPropagationRounds = 3;
    const int max_rounds = kPropagationRounds;

    for (int round = 0; round < max_rounds; ++round) {
        bool changed = false;
        std::fill(fresh, fresh + query_size, 0);
        for (VertexId u = 0; u < query_size; ++u) {
            std::vector<VertexId>& current = candidates[u];
            if (current.empty()) {
                continue;
            }

            std::size_t neighbor_count = 0;
            const VertexId* query_neighbors = query_graph.neighbors(u, neighbor_count);
            if (neighbor_count == 0) {
                continue;
            }

            for (std::size_t n = 0; n < neighbor_count; ++n) {
                const VertexId neighbor_u = query_neighbors[n];
                if (fresh[neighbor_u] == 0) {
                    rebuild_support(neighbor_u);
                    fresh[neighbor_u] = 1;
                }
            }

            std::vector<VertexId> next_candidates;
            next_candidates.reserve(current.size());

            for (VertexId v : current) {
                const std::size_t word = v >> 6;
                const std::uint64_t bit = std::uint64_t{1} << (v & 63);
                bool valid = true;
                for (std::size_t n = 0; n < neighbor_count; ++n) {
                    const std::uint64_t* row =
                        support + static_cast<std::size_t>(query_neighbors[n]) * words;
                    if ((row[word] & bit) == 0) {
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    next_candidates.push_back(v);
                } else {
                    changed = true;
                }
            }

            const bool u_changed = next_candidates.size() != current.size();
            current = std::move(next_candidates);
            if (u_changed) {
                rebuild_support(u);
                fresh[u] = 1;
            }
        }

        if (!changed) {
            break;
        }
    }
}

BacktrackingVerifier::Result BacktrackingVerifier::verify(const LabeledGraph& query_graph,
                                  CandidateSets& candidates,
                                  std::size_t limit) {
    Result result;
    const std::size_t query_size = query_graph.num_vertices();
    if (query_size == 0 || query_size > MAX_QUERY_VERTICES) {
        return result;
    }

    for (const auto& set : candidates) {
        if (set.empty()) {
            return result;
        }
    }

    const auto t_prop_start = std::chrono::high_resolution_clock::now();
    propagate_edge_constraints(query_graph, candidates);
    result.candidates_after_propagation = 0;
    for (const auto& set : candidates) {
        if (set.empty()) {
            result.candidates_after_propagation = 0;
            result.propagation_us = micros_since(t_prop_start,
                std::chrono::high_resolution_clock::now());
            return result;
        }
        result.candidates_after_propagation += set.size();
    }
    const auto t_prop_end = std::chrono::high_resolution_clock::now();
    result.propagation_us = micros_since(t_prop_start, t_prop_end);

    EnumerationBackend::Result engine_result = engine_->verify(query_graph, candidates, limit);
    result.count = engine_result.count;
    result.early_stopped = engine_result.early_stopped;
    result.adjacency_us = engine_result.adjacency_us;
    result.enumeration_us = engine_result.enumeration_us;
    return result;
}
}
