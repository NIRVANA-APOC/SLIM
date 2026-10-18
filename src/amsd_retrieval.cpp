#include <slim/amsd_retrieval.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace slim {
namespace {
double expected_reduction(std::size_t current, std::size_t suffix_size,
                          std::size_t partition_size) {
    if (partition_size == 0) {
        return 0.0;
    }
    const double fraction = static_cast<double>(suffix_size) /
                            static_cast<double>(partition_size);
    const double survivors = static_cast<double>(current) * std::min(1.0, fraction);
    return std::max(0.0, static_cast<double>(current) - survivors);
}

bool worthwhile(double reduction, std::size_t words) {
    return reduction * kScanNsPerCandidate > static_cast<double>(words) * kBitmapNsPerWord;
}
}

const BitVector& AmsdRetriever::dimension_bitmap(
    const SpectralDominanceIndex::LabelGroup& label_group, std::uint32_t dimension) {
    for (const DimensionCacheEntry& entry : dimension_cache_) {
        if (entry.group == &label_group && entry.dimension == dimension) {
            return entry.bits;
        }
    }

    DimensionCacheEntry& slot = dimension_cache_[cache_next_];
    cache_next_ = (cache_next_ + 1) % kDimensionCacheEntries;
    slot.group = &label_group;
    slot.dimension = dimension;
    slot.bits.resize_bits(label_group.vertices.size());
    slot.bits.clear();
    for (const SpectralDominanceIndex::Bucket& bucket : label_group.buckets) {
        if (bucket.dimension < dimension) {
            continue;
        }
        slot.bits.set_range(static_cast<std::size_t>(bucket.begin),
                            static_cast<std::size_t>(bucket.begin + bucket.count));
    }
    return slot.bits;
}

void AmsdRetriever::retrieve(const SpectralDominanceIndex::LabelGroup& label_group,
                             const AmsdGroup& group, const AmsdQueryPlan& plan,
                             BitVector& result) {
    const std::size_t partition_size = label_group.vertices.size();
    const BitVector& base = dimension_bitmap(label_group, plan.dimension);
    const std::size_t words = base.words();
    const std::size_t region = base.count();

    result.resize_bits(base.bits());
    result.clear();
    for (std::size_t word = 0; word < words; ++word) {
        result.word_data()[word] = base.word_data()[word];
    }

    if (!worthwhile(static_cast<double>(region), words)) {
        return;
    }

    channels_.clear();
    channels_.reserve(plan.ranks.size());
    for (const std::pair<std::uint32_t, double>& entry : plan.ranks) {
        const AmsdRankIndex* rank = group.find_rank(entry.first);
        if (rank == nullptr) {
            continue;
        }
        const int rung = rank->rung_for(entry.second);
        if (rung < 0) {
            continue;
        }
        const std::size_t index = static_cast<std::size_t>(rung);
        channels_.push_back({&rank->suffix[index], rank->suffix_sets[index], entry.first});
    }

    std::sort(channels_.begin(), channels_.end(),
              [](const Channel& a, const Channel& b) {
                  return a.suffix_size < b.suffix_size;
              });

    std::size_t count = region;
    for (const Channel& channel : channels_) {
        if (count <= 1) {
            break;
        }
        const double reduction = expected_reduction(count, channel.suffix_size, partition_size);
        if (!worthwhile(reduction, words)) {
            continue;
        }
        result.and_with(*channel.suffix);
        count = result.count();
    }
}
}
