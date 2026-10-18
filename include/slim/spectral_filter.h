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

private:
    using QueryLabelCount = std::pair<std::uint32_t, std::uint8_t>;

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

    bool use_index_degrees_ = false;

    mutable AmsdRetriever amsd_retriever_;
    mutable BitVector amsd_result_;
    mutable AmsdQueryPlan amsd_plan_;
    mutable std::vector<VertexId> amsd_survivors_;
    mutable std::vector<RankChannel> amsd_channels_;

    mutable std::vector<VertexId> radix_scratch_;
    mutable std::vector<QueryLabelCount> query_ball_labels_;
    mutable std::vector<std::uint64_t> query_ball_begin_;
};
}

#endif
