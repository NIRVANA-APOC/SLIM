#ifndef SLIM_CUMULATIVE_SPECTRAL_FILTER_H
#define SLIM_CUMULATIVE_SPECTRAL_FILTER_H

#include <slim/k_hop_neighborhood.h>
#include <slim/label_encoder.h>
#include <slim/labeled_graph.h>
#include <slim/spectral_dominance_index.h>
#include <slim/spectral_filter.h>
#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace slim {
struct CumulativeLayerStats {
    int radius = 0;
    int index_k = 0;
    double index_tau = 0.0;
    std::uint64_t index_neighborhood_size_limit = 0;
    std::uint64_t entry = 0;
    std::uint64_t exit = 0;
    std::uint64_t missing_signature_keeps = 0;
    std::uint64_t released = 0;
    bool ran = false;
    std::uint64_t spectral_us = 0;
    std::uint64_t filter_us = 0;
};

struct CumulativeFilterStats {
    std::vector<CumulativeLayerStats> layers;
    std::uint64_t base_candidates = 0;
    std::uint64_t candidates_after_filter = 0;
    std::uint64_t spectral_us = 0;
    std::uint64_t filter_us = 0;
    std::uint64_t released = 0;
    std::uint64_t missing_signature_keeps = 0;

    std::uint64_t coarse_dominance_checks = 0;
    std::uint64_t full_dominance_checks = 0;
    std::uint64_t ld_candidates = 0;
    std::uint64_t label_rejects = 0;
    std::uint64_t dim_rejects = 0;
    std::uint64_t key_rejects = 0;
    std::uint64_t win_rejects = 0;
    std::uint64_t empty_accepts = 0;
};

std::vector<std::unique_ptr<SpectralDominanceIndex>> load_cumulative_index_set(
    const std::vector<std::string>& paths, int max_radius);

void validate_cumulative_index_set(
    const std::vector<std::unique_ptr<SpectralDominanceIndex>>& indexes,
    const LabeledGraph& data_graph);

class CumulativeSpectralFilter {
public:
    CumulativeSpectralFilter(std::vector<std::unique_ptr<SpectralDominanceIndex>> indexes,
                             const LabeledGraph& data_graph);

    int max_radius() const { return static_cast<int>(layers_.size()); }

    const std::vector<std::unique_ptr<SpectralDominanceIndex>>& indexes() const {
        return indexes_;
    }

    void set_label_filter(bool enabled) { label_filter_ = enabled; }
    void set_spectral_filter(bool enabled) { spectral_filter_ = enabled; }
    void set_channel_stats(bool enabled) { channel_stats_ = enabled; }

    void set_amsd_adaptive(bool adaptive) { amsd_adaptive_ = adaptive; }

    void set_layer_dump_dir(std::string dir, std::string stem) {
        layer_dump_dir_ = std::move(dir);
        layer_dump_stem_ = std::move(stem);
    }

    void filter_into(const LabeledGraph& query_graph, CandidateSets& out,
                     CumulativeFilterStats* stats);

private:
    struct Layer {
        SpectralDominanceIndex* index = nullptr;
        std::unique_ptr<SpectralFilter> filter;
        std::vector<LocalSpectralSignature> query_signatures;
    };

    std::vector<std::unique_ptr<SpectralDominanceIndex>> indexes_;
    LabelEncoder encoder_;
    std::vector<Layer> layers_;
    const LabeledGraph& data_graph_;
    bool label_filter_ = true;
    bool spectral_filter_ = true;
    bool channel_stats_ = false;
    bool amsd_adaptive_ = true;
    std::string layer_dump_dir_;
    std::string layer_dump_stem_;

    CandidateSets accumulate_;
    CandidateSets layer_out_;
    CandidateSets scratch_;
};
}

#endif
