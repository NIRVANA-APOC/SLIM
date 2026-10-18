#ifndef SLIM_RANK_SELECTION_H
#define SLIM_RANK_SELECTION_H

#include <slim/config.h>
#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstdint>
#include <vector>

namespace slim {
constexpr double kMinRankBudgetRatio = 1.0;

constexpr double kDefaultRankBudgetRatio = 4.0;

constexpr std::uint32_t kDefaultMaxRankCandidate = 32;

struct RankSelectionOptions {
    RankSelectMode mode = RankSelectMode::Adaptive;
    double tau = 0.0;
    std::size_t num_vertices = 0;
    double budget_entries = 0.0;
    std::uint32_t max_rank = kDefaultMaxRankCandidate;
    std::vector<std::uint32_t> fallback_ranks;
};

struct RankNmiRecord {
    std::uint32_t rank_a = 0;
    std::uint32_t rank_b = 0;
    double value = 0.0;
};

struct RankSelectionResult {
    RankSelectMode mode = RankSelectMode::Adaptive;
    double tau = 0.0;
    std::size_t num_vertices = 0;
    double budget_entries = 0.0;
    std::uint32_t max_rank = 0;

    std::vector<std::uint64_t> c_j;
    std::vector<double> h_j;
    std::vector<double> d_j;

    std::vector<std::uint32_t> selected_ranks;
    std::vector<std::uint32_t> selection_order;
    std::vector<std::uint64_t> selected_entries;
    std::uint64_t space_entries = 0;

    std::vector<RankNmiRecord> nmi_calls;
};

RankSelectionResult select_rank_set(const std::vector<LocalSpectralSignature>& signatures,
                                    const RankSelectionOptions& options);

double normalized_conditional_entropy(const std::vector<LocalSpectralSignature>& signatures,
                                      std::size_t num_vertices, double tau, std::uint32_t rank,
                                      std::uint32_t* count_out = nullptr);

double normalized_conditional_mutual_information(
    const std::vector<LocalSpectralSignature>& signatures, std::size_t num_vertices, double tau,
    std::uint32_t rank_i, std::uint32_t rank_j);
}

#endif
