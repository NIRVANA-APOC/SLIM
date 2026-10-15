#include <slim/slim_matcher.h>

#include <slim/candidate_dump.h>

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

std::size_t count_candidates(const Candidates& candidates) {
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

SlimMatcher::SlimMatcher(const std::vector<std::string>& index_paths, int max_radius) {
    pending_indexes_ = load_cumulative_index_set(index_paths, max_radius);
}

void SlimMatcher::set_data_graph(LabeledGraph data_graph, VerifierBackend backend) {
    data_graph_ = std::move(data_graph);
    has_data_graph_ = true;
    encoder_ = LabelEncoder::from_graph(data_graph_);
    if (index_ != nullptr) {
        filter_ = std::make_unique<SpectralFilter>(*index_, data_graph_, encoder_);
        filter_->set_channel_stats(channel_stats_);
        filter_->set_label_filter(label_filter_);
        filter_->set_spectral_filter(spectral_filter_);
        filter_->set_query_radius(query_radius_);
    } else {
        validate_cumulative_index_set(pending_indexes_, data_graph_);
        cumulative_filter_ = std::make_unique<CumulativeSpectralFilter>(
            std::move(pending_indexes_), data_graph_);
        cumulative_filter_->set_channel_stats(channel_stats_);
        cumulative_filter_->set_label_filter(label_filter_);
        cumulative_filter_->set_spectral_filter(spectral_filter_);
    }
    verifier_ = std::make_unique<BacktrackingVerifier>(data_graph_, backend);
}

int SlimMatcher::max_radius() const {
    if (cumulative_filter_ != nullptr) {
        return cumulative_filter_->max_radius();
    }
    if (!pending_indexes_.empty()) {
        return static_cast<int>(pending_indexes_.size());
    }
    return index_ != nullptr ? index_->k() : 0;
}

const SpectralDominanceIndex& SlimMatcher::index_for_radius(int radius) const {
    const int layers = max_radius();
    if (radius < 1 || radius > layers) {
        throw std::runtime_error("radius " + std::to_string(radius) +
                                 " is outside the loaded index set 1.." +
                                 std::to_string(layers));
    }
    if (index_ != nullptr) {
        return *index_;
    }
    const std::vector<std::unique_ptr<SpectralDominanceIndex>>& indexes =
        cumulative_filter_ != nullptr ? cumulative_filter_->indexes() : pending_indexes_;
    return *indexes[static_cast<std::size_t>(radius) - 1];
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
    stats->index_tau = index_ != nullptr ? index_->tau() : index_for_radius(1).tau();
    if (index_ != nullptr && index_->has_identity()) {
        stats->index_neighborhood_size_limit =
            index_->identity().neighborhood_size_limit;
    }
    if (query_graph.num_vertices() > MAX_QUERY_VERTICES) {
        throw std::runtime_error("query graph is too large");
    }

    const auto start_total = std::chrono::high_resolution_clock::now();

    if (cumulative_filter_ != nullptr) {
        CumulativeSpectralFilter& cumulative = *cumulative_filter_;
        cumulative.set_channel_stats(channel_stats_);
        cumulative.set_label_filter(label_filter_);
        cumulative.set_spectral_filter(spectral_filter_);
        cumulative.set_layer_dump_dir(layer_dump_dir_, layer_dump_stem_);
        CumulativeFilterStats cumulative_stats;
        cumulative.filter_into(query_graph, candidate_buffer_, &cumulative_stats);
        stats->radius_mode = RadiusMode::Cumulative;
        stats->max_radius = cumulative.max_radius();
        stats->layers = cumulative_stats.layers;
        stats->cumulative_released = cumulative_stats.released;
        stats->cumulative_missing_signature_keeps = cumulative_stats.missing_signature_keeps;
        stats->query_spectral_us = cumulative_stats.spectral_us;
        stats->filter_us = cumulative_stats.filter_us;
        stats->base_ld_candidates = cumulative_stats.base_candidates;
        stats->candidates_after_filter = count_candidates(candidate_buffer_);
        stats->coarse_dominance_checks = cumulative_stats.coarse_dominance_checks;
        stats->full_dominance_checks = cumulative_stats.full_dominance_checks;
        stats->channels_enabled = channel_stats_;
        stats->ld_candidates = cumulative_stats.ld_candidates;
        stats->label_rejects = cumulative_stats.label_rejects;
        stats->dim_rejects = cumulative_stats.dim_rejects;
        stats->key_rejects = cumulative_stats.key_rejects;
        stats->win_rejects = cumulative_stats.win_rejects;
        stats->empty_accepts = cumulative_stats.empty_accepts;
        stats->label_filter_active = label_filter_;
        if (!candidate_dump_path_.empty()) {
            write_candidate_dump(candidate_dump_path_, candidate_buffer_);
        }
        if (filter_only) {
            stats->total_us = micros_since(start_total,
                                           std::chrono::high_resolution_clock::now());
            return {0, false};
        }
        const auto start_verify = std::chrono::high_resolution_clock::now();
        BacktrackingVerifier::Result result =
            verifier_->verify(query_graph, candidate_buffer_, limit);
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

    filter.set_channel_stats(channel_stats_);
    filter.set_label_filter(label_filter_);
    filter.set_spectral_filter(spectral_filter_);
    const auto start_filter = std::chrono::high_resolution_clock::now();
    filter.filter_into(query_graph, query_signatures, candidate_buffer_);
    const auto end_filter = std::chrono::high_resolution_clock::now();
    stats->filter_us = micros_since(start_filter, end_filter);
    stats->coarse_dominance_checks = filter.coarse_dominance_checks();
    stats->full_dominance_checks = filter.full_dominance_checks();
    stats->candidates_after_filter = count_candidates(candidate_buffer_);
    stats->base_ld_candidates = filter.count_label_degree_candidates(query_graph);
    stats->channels_enabled = filter.channel_stats();
    stats->ld_candidates = filter.ld_candidates();
    stats->label_rejects = filter.label_rejects();
    stats->dim_rejects = filter.dim_rejects();
    stats->label_filter_active = filter.label_filter_enabled() && index_->has_ball_labels();
    stats->key_rejects = filter.key_rejects();
    stats->win_rejects = filter.win_rejects();
    stats->empty_accepts = filter.empty_accepts();

    if (!candidate_dump_path_.empty()) {
        write_candidate_dump(candidate_dump_path_, candidate_buffer_);
    }

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
