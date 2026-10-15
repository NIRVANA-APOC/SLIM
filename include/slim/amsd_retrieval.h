#ifndef SLIM_AMSD_RETRIEVAL_H
#define SLIM_AMSD_RETRIEVAL_H

#include <slim/amsd_index.h>
#include <slim/spectral_dominance_index.h>
#include <slim/types.h>

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace slim {
struct AmsdQueryPlan {
    std::uint32_t dimension = 0;

    std::vector<std::pair<std::uint32_t, double>> ranks;
};

struct AmsdRetrievalStats {
    std::uint64_t dimension_bitmap_builds = 0;
    std::uint64_t dimension_bitmap_hits = 0;
    std::uint64_t rank_ands = 0;
    std::uint64_t rank_ands_skipped = 0;
    std::uint64_t ranks_unavailable = 0;
    std::uint64_t direct_scans = 0;
    std::uint64_t region = 0;
    std::uint64_t enumerated = 0;
    std::uint64_t chain_tests = 0;
    std::uint64_t chain_rejections = 0;
    std::uint64_t chain_enabled_scans = 0;
};

constexpr double kScanNsPerCandidate = 12.0;
constexpr double kBitmapNsPerWord = 1.0;

constexpr double kAmsdSetupNs = 300.0;

constexpr double kChannelCorrelationPenalty = 4.0;

constexpr std::size_t kChainSampleStride = 16;

constexpr double kSurvivorExtractNs = 1.5;

inline double corrected_keep_fraction(double keep_fraction) {
    const double channel_keep = std::max(0.0, std::min(1.0, keep_fraction));
    return 1.0 - (1.0 - channel_keep) / kChannelCorrelationPenalty;
}

inline bool amsd_pays_off(std::size_t region, double keep_fraction, std::size_t words) {
    const double keep = corrected_keep_fraction(keep_fraction);
    const double preempted = static_cast<double>(region) * (1.0 - keep);
    const double survivors = static_cast<double>(region) * keep;
    return preempted * kScanNsPerCandidate >
           kAmsdSetupNs + static_cast<double>(words) * kBitmapNsPerWord +
               survivors * kSurvivorExtractNs;
}

constexpr double kRankBitTestNs = 3.0;
constexpr double kPreemptedCandidateNs = 12.0;

inline bool amsd_bit_test_pays_off(double keep_fraction, std::size_t rank_count) {
    if (rank_count == 0) {
        return false;
    }
    return (1.0 - keep_fraction) * kPreemptedCandidateNs >
           static_cast<double>(rank_count) * kRankBitTestNs;
}

class AmsdRetriever {
public:
    const BitVector& dimension_bitmap(const SpectralDominanceIndex::LabelGroup& label_group, std::uint32_t dimension,
                                      AmsdRetrievalStats* stats);

    void retrieve(const SpectralDominanceIndex::LabelGroup& label_group, const AmsdGroup& group,
                  const AmsdQueryPlan& plan, bool adaptive, BitVector& result,
                  AmsdRetrievalStats* stats);

private:
    struct DimensionCacheEntry {
        const SpectralDominanceIndex::LabelGroup* group = nullptr;
        std::uint32_t dimension = 0;
        BitVector bits;
    };

    struct Channel {
        const BitVector* suffix = nullptr;
        std::size_t suffix_size = 0;
        std::uint32_t rank = 0;
    };

    static constexpr std::size_t kDimensionCacheEntries = 512;
    std::array<DimensionCacheEntry, kDimensionCacheEntries> dimension_cache_{};
    std::size_t cache_next_ = 0;

    std::vector<Channel> channels_;
};
}

#endif
