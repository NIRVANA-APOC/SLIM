#ifndef SLIM_AMSD_INDEX_H
#define SLIM_AMSD_INDEX_H

#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstdint>
#include <vector>

namespace slim {
class BitVector {
public:
    void resize_bits(std::size_t bits);
    void clear();
    std::size_t bits() const { return bits_; }
    std::size_t words() const { return words_.size(); }
    bool empty() const { return bits_ == 0; }
    bool none() const;
    std::size_t count() const;

    bool test(std::size_t index) const {
        return (words_[index >> 6] >> (index & 63u)) & 1u;
    }
    void set(std::size_t index) {
        words_[index >> 6] |= (std::uint64_t{1} << (index & 63u));
    }
    void reset(std::size_t index) {
        words_[index >> 6] &= ~(std::uint64_t{1} << (index & 63u));
    }
    void set_range(std::size_t begin, std::size_t end);

    void and_with(const BitVector& other);

    const std::uint64_t* word_data() const { return words_.data(); }
    std::uint64_t* word_data() { return words_.data(); }

    std::size_t next_set_bit(std::size_t from) const;

private:
    std::size_t bits_ = 0;
    std::vector<std::uint64_t> words_;
};

struct AmsdRankIndex {
    std::uint32_t rank = 0;

    std::vector<double> rungs;

    std::vector<std::uint32_t> suffix_sets;
    std::vector<BitVector> suffix;

    int rung_for(double threshold) const;
};

struct AmsdGroup {
    Label label = 0;
    std::uint32_t size = 0;
    std::vector<AmsdRankIndex> ranks;

    const AmsdRankIndex* find_rank(std::uint32_t rank) const;
};

struct AmsdIndex {
    bool present = false;
    std::uint32_t configured_ranks = 0;
    std::uint32_t configured_bins = 0;
    std::uint64_t bitmap_bytes = 0;
    std::vector<AmsdGroup> groups;

    std::vector<std::uint32_t> degrees;

    const AmsdGroup* find(Label label) const;
};

struct AmsdBuildGroup {
    Label label = 0;
    const VertexId* vertices = nullptr;
    std::size_t size = 0;
};

AmsdIndex build_amsd_index(const std::vector<AmsdBuildGroup>& groups,
                           const std::vector<LocalSpectralSignature>& signatures,
                           const std::vector<std::uint32_t>& degrees,
                           std::uint32_t max_ranks, std::uint32_t max_bins);
}

#endif
