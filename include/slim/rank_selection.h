#ifndef SLIM_RANK_SELECTION_H
#define SLIM_RANK_SELECTION_H

#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstdint>
#include <vector>

namespace slim {
constexpr double kRankBudgetRatio = 4.0;
constexpr std::uint32_t kMaxRankCandidate = 32;

struct RankSelectionOptions {
    double tau = 0.0;
    double budget_entries = 0.0;
    std::uint32_t max_rank = kMaxRankCandidate;
};

struct RankSelectionResult {
    std::vector<std::uint32_t> selected_ranks;
    std::uint64_t space_entries = 0;
};

RankSelectionResult select_rank_set(const std::vector<LocalSpectralSignature>& signatures,
                                    const RankSelectionOptions& options);
}

#endif
