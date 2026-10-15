#ifndef SLIM_SPECTRAL_FILTER_H
#define SLIM_SPECTRAL_FILTER_H

#include <slim/k_hop_neighborhood.h>
#include <slim/label_encoder.h>
#include <slim/labeled_graph.h>
#include <slim/spectral_dominance_index.h>
#include <slim/spectral_operator.h>
#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstdint>
#include <utility>
#include <vector>

namespace slim {
using CandidateSets = std::vector<std::vector<VertexId>>;

using Candidates = CandidateSets;

class SpectralFilter {
public:
    SpectralFilter(const SpectralDominanceIndex& index,
                   const LabeledGraph& data_graph,
                   const LabelEncoder& encoder);

    CandidateSets filter(const LabeledGraph& query_graph,
                         const std::vector<LocalSpectralSignature>& query_signatures) const;

    std::size_t count_label_degree_candidates(const LabeledGraph& query_graph) const;

    void filter_into(const LabeledGraph& query_graph,
                     const std::vector<LocalSpectralSignature>& query_signatures,
                     CandidateSets& out) const;

    LocalSpectralSignature compute_query_signature(const LabeledGraph& query_graph,
                                                   VertexId u,
                                                   KHopNeighborhood& neighborhood) const;

    std::uint64_t coarse_dominance_checks() const { return coarse_dominance_checks_; }
    std::uint64_t full_dominance_checks() const { return full_dominance_checks_; }

    std::uint64_t scanned_buckets() const { return scanned_buckets_; }
    std::uint64_t scanned_cells() const { return scanned_cells_; }

    void set_label_filter(bool enabled) { label_filter_enabled_ = enabled; }
    bool label_filter_enabled() const { return label_filter_enabled_; }

    void set_spectral_filter(bool enabled) { spectral_test_enabled_ = enabled; }

    void set_query_radius(int radius) { query_radius_ = radius; }
    int query_radius() const {
        return query_radius_ > 0 ? query_radius_ : static_cast<int>(index_.k());
    }
    bool spectral_filter_enabled() const { return spectral_test_enabled_; }

    void set_channel_stats(bool enabled) { channel_stats_ = enabled; }
    bool channel_stats() const { return channel_stats_; }
    std::uint64_t ld_candidates() const { return ld_candidates_; }
    std::uint64_t label_rejects() const { return label_rejects_; }
    std::uint64_t dim_rejects() const { return dim_rejects_; }
    std::uint64_t key_rejects() const { return key_rejects_; }
    std::uint64_t win_rejects() const { return win_rejects_; }
    std::uint64_t empty_accepts() const { return empty_accepts_; }

private:
    using QueryLabelCount = std::pair<std::uint32_t, std::uint8_t>;

    void classify_into(const LabeledGraph& query_graph,
                       const LocalSpectralSignature& query_sig,
                       VertexId u,
                       const QueryLabelCount* query_labels,
                       const QueryLabelCount* query_labels_end,
                       std::vector<VertexId>& accepted) const;

    const SpectralDominanceIndex& index_;
    const LabeledGraph& data_graph_;
    const LabelEncoder& encoder_;

    mutable std::uint64_t coarse_dominance_checks_ = 0;
    mutable std::uint64_t full_dominance_checks_ = 0;
    mutable std::uint64_t scanned_buckets_ = 0;
    mutable std::uint64_t scanned_cells_ = 0;

    bool channel_stats_ = false;
    int query_radius_ = 0;
    bool label_filter_enabled_ = true;
    bool spectral_test_enabled_ = true;
    mutable std::uint64_t ld_candidates_ = 0;
    mutable std::uint64_t label_rejects_ = 0;
    mutable std::uint64_t dim_rejects_ = 0;
    mutable std::uint64_t key_rejects_ = 0;
    mutable std::uint64_t win_rejects_ = 0;
    mutable std::uint64_t empty_accepts_ = 0;

    mutable std::vector<VertexId> radix_scratch_;

    mutable std::vector<QueryLabelCount> query_ball_labels_;
    mutable std::vector<std::uint64_t> query_ball_begin_;
};
}

#endif
