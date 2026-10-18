#ifndef SLIM_SPECTRAL_FILTER_H
#define SLIM_SPECTRAL_FILTER_H

#include <slim/amsd_index.h>
#include <slim/amsd_retrieval.h>
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

    void set_amsd_adaptive(bool adaptive) {
        amsd_adaptive_ = adaptive;
        amsd_provider_adaptive_ = adaptive;
    }
    bool amsd_adaptive() const { return amsd_adaptive_; }

    void set_stage_timings(bool enabled) { stage_timings_ = enabled; }
    bool stage_timings() const { return stage_timings_; }

    std::uint64_t amsd_rank_ands() const { return amsd_stats_.rank_ands; }
    std::uint64_t amsd_rank_ands_skipped() const { return amsd_stats_.rank_ands_skipped; }
    std::uint64_t amsd_ranks_unavailable() const { return amsd_stats_.ranks_unavailable; }
    std::uint64_t amsd_direct_scans() const { return amsd_stats_.direct_scans; }
    std::uint64_t amsd_region() const { return amsd_stats_.region; }
    std::uint64_t amsd_enumerated() const { return amsd_stats_.enumerated; }
    std::uint64_t amsd_preempted() const {
        return amsd_stats_.region >= amsd_stats_.enumerated
                   ? amsd_stats_.region - amsd_stats_.enumerated
                   : 0;
    }
    std::uint64_t amsd_active_vertices() const { return amsd_active_vertices_; }
    std::uint64_t amsd_bucket_scan_vertices() const { return amsd_bucket_vertices_; }
    std::uint64_t amsd_bit_test_vertices() const { return amsd_bit_test_vertices_; }
    std::uint64_t amsd_bitmap_vertices() const { return amsd_bitmap_vertices_; }
    std::uint64_t amsd_chain_tests() const { return amsd_stats_.chain_tests; }
    std::uint64_t amsd_chain_rejections() const { return amsd_stats_.chain_rejections; }
    std::uint64_t amsd_chain_enabled_scans() const { return amsd_stats_.chain_enabled_scans; }
    std::uint64_t amsd_dimension_builds() const { return amsd_stats_.dimension_bitmap_builds; }
    std::uint64_t amsd_dimension_hits() const { return amsd_stats_.dimension_bitmap_hits; }
    std::uint64_t dispatch_ns() const { return dispatch_ns_; }

private:
    using QueryLabelCount = std::pair<std::uint32_t, std::uint8_t>;

    void classify_into(const LabeledGraph& query_graph,
                       const LocalSpectralSignature& query_sig,
                       VertexId u,
                       const QueryLabelCount* query_labels,
                       const QueryLabelCount* query_labels_end,
                       std::vector<VertexId>& accepted) const;

    void retrieve_amsd_into(const SpectralDominanceIndex::LabelGroup& group,
                            const AmsdGroup& amsd_group,
                            const LocalSpectralSignature& query_sig,
                            std::size_t query_degree,
                            double query_lambda_max,
                            const QueryLabelCount* query_labels,
                            const QueryLabelCount* query_labels_end,
                            std::vector<VertexId>& accepted) const;

    bool verify_amsd_degrees() const;

    enum class Provider { BucketScan, BitTest, Bitmap };

    struct RankChannel {
        const BitVector* suffix = nullptr;
        std::uint32_t rank = 0;
        std::size_t suffix_size = 0;
    };

    Provider choose_provider(const SpectralDominanceIndex::LabelGroup& group,
                             const AmsdGroup& amsd_group,
                             const LocalSpectralSignature& query_sig,
                             std::size_t* region_out) const;

    void collect_rank_channels(const AmsdGroup& amsd_group,
                               const LocalSpectralSignature& query_sig) const;

    void scan_buckets_into(const SpectralDominanceIndex::LabelGroup& group,
                           const LocalSpectralSignature& query_sig,
                           std::size_t query_degree,
                           double query_lambda_max,
                           const QueryLabelCount* query_labels,
                           const QueryLabelCount* query_labels_end,
                           std::size_t region_hint,
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

    bool amsd_adaptive_ = true;
    bool amsd_provider_adaptive_ = true;
    bool stage_timings_ = false;
    bool use_index_degrees_ = false;

    mutable AmsdRetriever amsd_retriever_;
    mutable BitVector amsd_result_;
    mutable AmsdQueryPlan amsd_plan_;
    mutable AmsdRetrievalStats amsd_stats_;
    mutable std::vector<VertexId> amsd_survivors_;
    mutable std::vector<RankChannel> amsd_channels_;
    mutable std::uint64_t amsd_active_vertices_ = 0;
    mutable std::uint64_t amsd_bucket_vertices_ = 0;
    mutable std::uint64_t amsd_bit_test_vertices_ = 0;
    mutable std::uint64_t amsd_bitmap_vertices_ = 0;
    mutable std::uint64_t dispatch_ns_ = 0;

    mutable std::vector<VertexId> radix_scratch_;

    mutable std::vector<QueryLabelCount> query_ball_labels_;
    mutable std::vector<std::uint64_t> query_ball_begin_;
};
}

#endif
