#ifndef SLIM_SLIM_MATCHER_H
#define SLIM_SLIM_MATCHER_H

#include <slim/backtracking_verifier.h>
#include <slim/config.h>
#include <slim/cumulative_spectral_filter.h>
#include <slim/k_hop_neighborhood.h>
#include <slim/label_encoder.h>
#include <slim/labeled_graph.h>
#include <slim/spectral_dominance_index.h>
#include <slim/spectral_filter.h>
#include <slim/spectral_operator.h>
#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace slim {
struct MatchStats {
    std::uint64_t query_spectral_us = 0;
    std::uint64_t filter_us = 0;
    std::uint64_t verify_us = 0;
    std::uint64_t verify_propagation_us = 0;
    std::uint64_t verify_adjacency_us = 0;
    std::uint64_t verify_enumeration_us = 0;
    std::uint64_t total_us = 0;
    std::size_t candidates_after_filter = 0;
    std::size_t candidates_after_propagation = 0;
    std::size_t match_count = 0;
    std::uint64_t coarse_dominance_checks = 0;
    std::uint64_t full_dominance_checks = 0;
    bool early_stopped = false;

    std::size_t base_ld_candidates = 0;

    RadiusMode radius_mode = RadiusMode::Single;
    int max_radius = 0;
    std::vector<CumulativeLayerStats> layers;
    std::uint64_t cumulative_released = 0;
    std::uint64_t cumulative_missing_signature_keeps = 0;

    double index_tau = 0.0;

    std::uint64_t index_neighborhood_size_limit = 0;

    bool channels_enabled = false;
    std::uint64_t ld_candidates = 0;
    std::uint64_t label_rejects = 0;
    std::uint64_t dim_rejects = 0;
    std::uint64_t key_rejects = 0;
    std::uint64_t win_rejects = 0;
    std::uint64_t empty_accepts = 0;

    bool amsd_adaptive = true;
    std::uint64_t amsd_active_vertices = 0;
    std::uint64_t amsd_rank_ands = 0;
    std::uint64_t amsd_rank_ands_skipped = 0;
    std::uint64_t amsd_ranks_unavailable = 0;
    std::uint64_t amsd_direct_scans = 0;
    std::uint64_t amsd_region = 0;
    std::uint64_t amsd_enumerated = 0;
    std::uint64_t amsd_preempted = 0;
    std::uint64_t amsd_bucket_scan_vertices = 0;
    std::uint64_t amsd_bit_test_vertices = 0;
    std::uint64_t amsd_bitmap_vertices = 0;
    std::uint64_t amsd_chain_tests = 0;
    std::uint64_t amsd_chain_rejections = 0;
    std::uint64_t amsd_chain_enabled_scans = 0;
    std::uint64_t dispatch_ns = 0;

    bool label_filter_active = false;
};

struct MatchResult {
    std::size_t match_count = 0;
    bool early_stopped = false;
};

class SlimMatcher {
public:
    explicit SlimMatcher(const std::string& index_path);

    SlimMatcher(const std::vector<std::string>& index_paths, int max_radius);

    void set_data_graph(LabeledGraph data_graph,
                        VerifierBackend backend = VerifierBackend::Ive);

    const LabeledGraph& data_graph() const;

    const SpectralDominanceIndex& index() const { return *index_; }

    bool cumulative() const { return index_ == nullptr; }

    int max_radius() const;

    const SpectralDominanceIndex& index_for_radius(int radius) const;

    MatchResult match(const LabeledGraph& query_graph,
                      std::size_t limit,
                      bool filter_only,
                      MatchStats* stats = nullptr);

    void set_channel_stats(bool enabled) { channel_stats_ = enabled; }

    void set_stage_timings(bool enabled) {
        stage_timings_ = enabled;
        if (filter_ != nullptr) {
            filter_->set_stage_timings(enabled);
        }
    }

    void set_label_filter(bool enabled) { label_filter_ = enabled; }
    void set_spectral_filter(bool enabled) { spectral_filter_ = enabled; }

    void set_query_radius(int radius) {
        query_radius_ = radius;
        if (filter_ != nullptr) {
            filter_->set_query_radius(radius);
        }
    }

    void set_amsd_adaptive(bool adaptive) {
        amsd_adaptive_ = adaptive;
        if (filter_ != nullptr) {
            filter_->set_amsd_adaptive(adaptive);
        }
        if (cumulative_filter_ != nullptr) {
            cumulative_filter_->set_amsd_adaptive(adaptive);
        }
    }

    void set_candidate_dump_path(std::string path) { candidate_dump_path_ = std::move(path); }

    void set_layer_dump_dir(std::string dir, std::string stem) {
        layer_dump_dir_ = std::move(dir);
        layer_dump_stem_ = std::move(stem);
    }

private:
    std::unique_ptr<SpectralDominanceIndex> index_;
    std::vector<std::unique_ptr<SpectralDominanceIndex>> pending_indexes_;
    std::unique_ptr<CumulativeSpectralFilter> cumulative_filter_;
    LabeledGraph data_graph_;
    bool has_data_graph_ = false;
    LabelEncoder encoder_;
    std::unique_ptr<SpectralFilter> filter_;
    std::unique_ptr<BacktrackingVerifier> verifier_;
    CandidateSets candidate_buffer_;
    bool channel_stats_ = false;
    int query_radius_ = 0;
    bool label_filter_ = true;
    bool spectral_filter_ = true;
    bool amsd_adaptive_ = true;
    bool stage_timings_ = false;
    std::string candidate_dump_path_;
    std::string layer_dump_dir_;
    std::string layer_dump_stem_;
};
}

#endif
