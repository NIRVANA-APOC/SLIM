#include <slim/slim_matcher.h>

#include <chrono>
#include <stdexcept>
#include <utility>
#include <vector>

namespace slim {
namespace {
std::uint64_t micros_since(std::chrono::high_resolution_clock::time_point start,
                           std::chrono::high_resolution_clock::time_point end) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
}

std::size_t count_candidates(const CandidateSets& candidates) {
    std::size_t total = 0;
    for (const auto& set : candidates) {
        total += set.size();
    }
    return total;
}
}

SlimMatcher::SlimMatcher(const std::string& index_path) {
    index_ = SpectralDominanceIndex::load(index_path);
}

void SlimMatcher::set_data_graph(LabeledGraph data_graph, VerifierBackend backend) {
    data_graph_ = std::move(data_graph);
    has_data_graph_ = true;
    encoder_ = LabelEncoder::from_graph(data_graph_);
    filter_ = std::make_unique<SpectralFilter>(*index_, data_graph_, encoder_);
    verifier_ = std::make_unique<BacktrackingVerifier>(data_graph_, backend);
}

const LabeledGraph& SlimMatcher::data_graph() const {
    if (!has_data_graph_) {
        throw std::runtime_error("data graph is not loaded");
    }
    return data_graph_;
}

MatchResult SlimMatcher::match(const LabeledGraph& query_graph,
                               std::size_t limit,
                               bool filter_only,
                               MatchStats* stats) {
    MatchStats local_stats;
    if (stats == nullptr) {
        stats = &local_stats;
    }
    *stats = MatchStats();

    if (!has_data_graph_) {
        throw std::runtime_error("data graph is not loaded");
    }
    if (query_graph.num_vertices() > MAX_QUERY_VERTICES) {
        throw std::runtime_error("query graph is too large");
    }

    const auto start_total = std::chrono::high_resolution_clock::now();

    const auto start_spectral = std::chrono::high_resolution_clock::now();
    SpectralFilter& filter = *filter_;
    std::vector<LocalSpectralSignature> query_signatures(query_graph.num_vertices());
    KHopNeighborhood query_neighborhood(query_graph);
    for (VertexId u = 0; u < query_graph.num_vertices(); ++u) {
        query_signatures[u] = filter.compute_query_signature(query_graph, u,
                                                             query_neighborhood);
    }
    const auto end_spectral = std::chrono::high_resolution_clock::now();
    stats->query_spectral_us = micros_since(start_spectral, end_spectral);

    const auto start_filter = std::chrono::high_resolution_clock::now();
    filter.filter_into(query_graph, query_signatures, candidate_buffer_);
    const auto end_filter = std::chrono::high_resolution_clock::now();
    stats->filter_us = micros_since(start_filter, end_filter);
    stats->coarse_dominance_checks = filter.coarse_dominance_checks();
    stats->full_dominance_checks = filter.full_dominance_checks();
    stats->candidates_after_filter = count_candidates(candidate_buffer_);
    stats->label_degree_candidates = filter.count_label_degree_candidates(query_graph);

    if (filter_only) {
        stats->total_us = micros_since(start_total,
                                       std::chrono::high_resolution_clock::now());
        return {0, false};
    }

    const auto start_verify = std::chrono::high_resolution_clock::now();
    BacktrackingVerifier::Result result = verifier_->verify(query_graph, candidate_buffer_, limit);
    const auto end_verify = std::chrono::high_resolution_clock::now();
    stats->candidates_after_propagation = result.candidates_after_propagation;
    stats->verify_propagation_us = result.propagation_us;
    stats->verify_adjacency_us = result.adjacency_us;
    stats->verify_enumeration_us = result.enumeration_us;
    stats->verify_us = micros_since(start_verify, end_verify);
    stats->match_count = result.count;
    stats->early_stopped = result.early_stopped;
    stats->total_us = micros_since(start_total, end_verify);
    return {result.count, result.early_stopped};
}
}
