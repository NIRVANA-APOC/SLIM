#include <slim/amsd_index.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <vector>

namespace slim {
namespace {
inline int popcount64(std::uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(value);
#else
    int count = 0;
    while (value != 0) {
        value &= value - 1;
        ++count;
    }
    return count;
#endif
}

inline int ctz64(std::uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(value);
#else
    int count = 0;
    while ((value & 1u) == 0) {
        value >>= 1;
        ++count;
    }
    return count;
#endif
}
}

void BitVector::resize_bits(std::size_t bits) {
    bits_ = bits;
    words_.assign((bits + 63) / 64, 0);
}

void BitVector::clear() {
    std::fill(words_.begin(), words_.end(), std::uint64_t{0});
}

bool BitVector::none() const {
    for (std::uint64_t word : words_) {
        if (word != 0) {
            return false;
        }
    }
    return true;
}

std::size_t BitVector::count() const {
    std::size_t total = 0;
    for (std::uint64_t word : words_) {
        total += static_cast<std::size_t>(popcount64(word));
    }
    return total;
}

void BitVector::set_range(std::size_t begin, std::size_t end) {
    if (begin >= end || begin >= bits_) {
        return;
    }
    if (end > bits_) {
        end = bits_;
    }
    const std::size_t first_word = begin >> 6;
    const std::size_t last_word = (end - 1) >> 6;
    const std::uint64_t head = ~std::uint64_t{0} << (begin & 63u);
    const std::uint64_t tail = (end & 63u) == 0 ? ~std::uint64_t{0}
                                                : ((std::uint64_t{1} << (end & 63u)) - 1);
    if (first_word == last_word) {
        words_[first_word] |= (head & tail);
        return;
    }
    words_[first_word] |= head;
    for (std::size_t word = first_word + 1; word < last_word; ++word) {
        words_[word] = ~std::uint64_t{0};
    }
    words_[last_word] |= tail;
}

void BitVector::and_with(const BitVector& other) {
    const std::size_t n = std::min(words_.size(), other.words_.size());
    for (std::size_t i = 0; i < n; ++i) {
        words_[i] &= other.words_[i];
    }
}

std::size_t BitVector::next_set_bit(std::size_t from) const {
    if (from >= bits_) {
        return bits_;
    }
    std::size_t word_index = from >> 6;
    std::uint64_t word = words_[word_index] & (~std::uint64_t{0} << (from & 63u));
    for (;;) {
        if (word != 0) {
            const std::size_t index = (word_index << 6) + static_cast<std::size_t>(ctz64(word));
            return index < bits_ ? index : bits_;
        }
        if (++word_index >= words_.size()) {
            return bits_;
        }
        word = words_[word_index];
    }
}

int AmsdRankIndex::rung_for(double threshold) const {
    const auto it = std::lower_bound(rungs.begin(), rungs.end(), threshold,
                                     std::greater<double>());
    if (it == rungs.end()) {
        return -1;
    }
    return static_cast<int>(it - rungs.begin());
}

const AmsdRankIndex* AmsdGroup::find_rank(std::uint32_t rank) const {
    const auto it = std::lower_bound(
        ranks.begin(), ranks.end(), rank,
        [](const AmsdRankIndex& entry, std::uint32_t value) { return entry.rank < value; });
    if (it == ranks.end() || it->rank != rank) {
        return nullptr;
    }
    return &*it;
}

const AmsdGroup* AmsdIndex::find(Label label) const {
    const auto it = std::lower_bound(
        groups.begin(), groups.end(), label,
        [](const AmsdGroup& group, Label value) { return group.label < value; });
    if (it == groups.end() || it->label != label) {
        return nullptr;
    }
    return &*it;
}

AmsdIndex build_amsd_index(const std::vector<AmsdBuildGroup>& groups,
                           const std::vector<LocalSpectralSignature>& signatures,
                           const std::vector<std::uint32_t>& degrees,
                           std::uint32_t max_ranks, std::uint32_t max_bins) {
    AmsdIndex amsd;
    amsd.present = true;
    amsd.configured_ranks = max_ranks;
    amsd.configured_bins = max_bins;
    amsd.groups.resize(groups.size());
    amsd.degrees = degrees;

    std::vector<VertexId> sorted_local;
    std::vector<double> sorted_value;
    std::vector<double> rungs;
    BitVector accumulated;

    for (std::size_t g = 0; g < groups.size(); ++g) {
        const AmsdBuildGroup& input = groups[g];
        AmsdGroup& group = amsd.groups[g];
        group.label = input.label;
        group.size = static_cast<std::uint32_t>(input.size);

        std::uint32_t max_dimension = 0;
        for (std::size_t i = 0; i < input.size; ++i) {
            const VertexId v = input.vertices[i];
            max_dimension = std::max(max_dimension,
                                     static_cast<std::uint32_t>(signatures[v].dimension()));
        }
        if (input.size == 0) {
            continue;
        }

        const std::uint32_t ranks_to_build = std::min(max_ranks, max_dimension);
        for (std::uint32_t rank = 1; rank <= ranks_to_build; ++rank) {
            sorted_local.clear();
            sorted_value.clear();
            std::uint32_t present = 0;
            for (std::size_t i = 0; i < input.size; ++i) {
                const LocalSpectralSignature& signature = signatures[input.vertices[i]];
                const std::size_t dimension = signature.dimension();
                if (dimension < rank) {
                    continue;
                }
                sorted_local.push_back(static_cast<VertexId>(i));
                sorted_value.push_back(
                    stored_spectral_value(signature.eigenvalues[dimension - rank]));
                ++present;
            }
            if (present == 0) {
                continue;
            }

            std::vector<std::size_t> order(present);
            std::iota(order.begin(), order.end(), std::size_t{0});
            std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
                if (sorted_value[a] != sorted_value[b]) {
                    return sorted_value[a] < sorted_value[b];
                }
                return sorted_local[a] < sorted_local[b];
            });
            std::vector<VertexId> by_value(present);
            std::vector<double> values(present);
            for (std::size_t i = 0; i < present; ++i) {
                by_value[i] = sorted_local[order[i]];
                values[i] = sorted_value[order[i]];
            }

            const std::size_t rung_budget =
                std::min<std::size_t>(static_cast<std::size_t>(max_bins),
                                      static_cast<std::size_t>(present));
            rungs.clear();
            for (std::size_t k = 0; k < rung_budget; ++k) {
                std::size_t step = 1;
                if (rung_budget > 1) {
                    const double exponent =
                        static_cast<double>(k) / static_cast<double>(rung_budget - 1);
                    step = static_cast<std::size_t>(
                        std::ceil(std::pow(static_cast<double>(present), exponent)));
                }
                step = std::min(step, static_cast<std::size_t>(present));
                const double value = values[present - step];
                if (!rungs.empty() && value >= rungs.back()) {
                    continue;
                }
                const std::size_t first =
                    static_cast<std::size_t>(std::lower_bound(values.begin(), values.end(),
                                                              value) - values.begin());
                if (first == 0) {
                    break;
                }
                rungs.push_back(value);
            }
            if (rungs.empty()) {
                continue;
            }

            AmsdRankIndex entry;
            entry.rank = rank;
            entry.rungs = rungs;
            entry.suffix.resize(rungs.size());
            entry.suffix_sets.resize(rungs.size());

            accumulated.resize_bits(input.size);
            accumulated.clear();
            std::size_t next = present;
            for (std::size_t b = 0; b < rungs.size(); ++b) {
                const std::size_t first =
                    static_cast<std::size_t>(std::lower_bound(values.begin(), values.end(),
                                                              rungs[b]) - values.begin());
                for (std::size_t i = first; i < next; ++i) {
                    accumulated.set(by_value[i]);
                }
                next = std::min(next, first);
                entry.suffix[b] = accumulated;
                entry.suffix_sets[b] = static_cast<std::uint32_t>(accumulated.count());
            }
            amsd.bitmap_bytes += static_cast<std::uint64_t>(rungs.size()) *
                                 static_cast<std::uint64_t>(accumulated.words()) * sizeof(std::uint64_t);
            group.ranks.push_back(std::move(entry));
        }
    }
    return amsd;
}
}
