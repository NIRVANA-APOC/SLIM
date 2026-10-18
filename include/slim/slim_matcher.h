#ifndef SLIM_SLIM_MATCHER_H
#define SLIM_SLIM_MATCHER_H

#include <slim/backtracking_verifier.h>
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
    std::size_t label_degree_candidates = 0;
    std::size_t candidates_after_filter = 0;
    std::size_t candidates_after_propagation = 0;
    std::size_t match_count = 0;
    std::uint64_t coarse_dominance_checks = 0;
    std::uint64_t full_dominance_checks = 0;
    bool early_stopped = false;
};

struct MatchResult {
    std::size_t match_count = 0;
    bool early_stopped = false;
};

class SlimMatcher {
public:
    explicit SlimMatcher(const std::string& index_path);

    void set_data_graph(LabeledGraph data_graph,
                        VerifierBackend backend = VerifierBackend::Ive);

    const LabeledGraph& data_graph() const;

    const SpectralDominanceIndex& index() const { return *index_; }

    MatchResult match(const LabeledGraph& query_graph,
                      std::size_t limit,
                      bool filter_only,
                      MatchStats* stats = nullptr);

private:
    std::unique_ptr<SpectralDominanceIndex> index_;
    LabeledGraph data_graph_;
    bool has_data_graph_ = false;
    LabelEncoder encoder_;
    std::unique_ptr<SpectralFilter> filter_;
    std::unique_ptr<BacktrackingVerifier> verifier_;
    CandidateSets candidate_buffer_;
};
}

#endif
