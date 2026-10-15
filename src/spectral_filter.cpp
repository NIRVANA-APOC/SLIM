#include <slim/spectral_filter.h>

#include <chrono>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#define SLIM_PREFETCH(addr) __builtin_prefetch(static_cast<const void*>(addr))
#else
#define SLIM_PREFETCH(addr) ((void)0)
#endif

namespace slim {
namespace {
constexpr std::uint64_t kHeaderLookahead = 16;
constexpr std::uint64_t kDataLookahead = 8;

void sort_vertex_ids(std::vector<VertexId>& ids, std::vector<VertexId>& scratch,
                     std::uint32_t max_id) {
    const std::size_t count = ids.size();
    if (count < 2) {
        return;
    }
    if (count < 256) {
        std::sort(ids.begin(), ids.end());
        return;
    }

    unsigned passes = 1;
    for (std::uint32_t range = max_id; range > 0xFFu; range >>= 8) {
        ++passes;
    }

    scratch.resize(count);
    std::array<std::uint32_t, 256> offset;
    std::vector<VertexId>* src = &ids;
    std::vector<VertexId>* dst = &scratch;
    for (unsigned pass = 0; pass < passes; ++pass) {
        const unsigned shift = 8u * pass;
        offset.fill(0);
        for (std::size_t i = 0; i < count; ++i) {
            ++offset[((*src)[i] >> shift) & 0xFFu];
        }
        std::uint32_t running = 0;
        for (std::uint32_t& bucket : offset) {
            const std::uint32_t size = bucket;
            bucket = running;
            running += size;
        }
        for (std::size_t i = 0; i < count; ++i) {
            const VertexId v = (*src)[i];
            (*dst)[offset[(v >> shift) & 0xFFu]++] = v;
        }
        std::swap(src, dst);
    }
    if (src != &ids) {
        ids.swap(scratch);
    }
}

using QueryLabelCount = std::pair<std::uint32_t, std::uint8_t>;

bool label_covered(const std::uint8_t* data_counts, const QueryLabelCount* query_begin,
                   const QueryLabelCount* query_end) {
    if (query_begin == query_end || data_counts == nullptr) {
        return true;
    }
    for (const QueryLabelCount* q = query_begin; q != query_end; ++q) {
        const std::uint8_t have = data_counts[q->first];
        if (have != 255 && have < q->second) {
            return false;
        }
    }
    return true;
}
}

SpectralFilter::SpectralFilter(const SpectralDominanceIndex& index,
                               const LabeledGraph& data_graph,
                               const LabelEncoder& encoder)
    : index_(index), data_graph_(data_graph), encoder_(encoder) {
    if (!index_.has_amsd()) {
        throw std::invalid_argument(
            "the index has no multi-rank index section; rebuild it with this slim_index");
    }
    use_index_degrees_ = verify_amsd_degrees();
}

bool SpectralFilter::verify_amsd_degrees() const {
    const AmsdIndex& amsd = index_.amsd();
    if (amsd.degrees.size() != index_.num_vertices()) {
        return false;
    }
    for (std::size_t position = 0; position < amsd.groups.size(); ++position) {
        const AmsdGroup& group = amsd.groups[position];
        const SpectralDominanceIndex::LabelGroup* label_group = index_.find_group(group.label);
        if (label_group == nullptr ||
            index_.group_position(group.label) != position ||
            group.size != label_group->vertices.size()) {
            return false;
        }
    }
    for (VertexId v = 0; v < amsd.degrees.size(); ++v) {
        if (amsd.degrees[v] != static_cast<std::uint32_t>(data_graph_.degree(v))) {
            return false;
        }
    }
    return !amsd.groups.empty();
}

CandidateSets SpectralFilter::filter(
    const LabeledGraph& query_graph,
    const std::vector<LocalSpectralSignature>& query_signatures) const {
    CandidateSets candidates;
    filter_into(query_graph, query_signatures, candidates);
    return candidates;
}

void SpectralFilter::filter_into(
    const LabeledGraph& query_graph,
    const std::vector<LocalSpectralSignature>& query_signatures,
    CandidateSets& candidates) const {
    const std::size_t query_size = query_graph.num_vertices();
    if (query_signatures.size() != query_size) {
        throw std::runtime_error("query signature count mismatch");
    }
    if (candidates.size() != query_size) {
        candidates.resize(query_size);
    }

    coarse_dominance_checks_ = 0;
    full_dominance_checks_ = 0;
    scanned_buckets_ = 0;
    scanned_cells_ = 0;
    ld_candidates_ = 0;
    label_rejects_ = 0;
    dim_rejects_ = 0;
    key_rejects_ = 0;
    win_rejects_ = 0;
    empty_accepts_ = 0;
    amsd_stats_ = AmsdRetrievalStats();
    amsd_active_vertices_ = 0;
    amsd_bucket_vertices_ = 0;
    amsd_bit_test_vertices_ = 0;
    amsd_bitmap_vertices_ = 0;
    dispatch_ns_ = 0;

    const bool label_filter = label_filter_enabled_ && index_.has_ball_labels();
    const bool spectral_test = spectral_test_enabled_;
    query_ball_labels_.clear();
    query_ball_begin_.assign(query_size + 1, 0);
    if (label_filter) {
        const std::vector<Label>& alphabet = index_.alphabet();
        KHopNeighborhood neighborhood(query_graph);
        std::vector<std::uint32_t> ids;
        for (VertexId u = 0; u < query_size; ++u) {
            const std::vector<VertexId>& ball = neighborhood.extract(u, query_radius());
            ids.clear();
            ids.reserve(ball.size());
            for (VertexId w : ball) {
                const Label label = query_graph.label(w);
                const auto it = std::lower_bound(alphabet.begin(), alphabet.end(), label);
                if (it == alphabet.end() || *it != label) {
                    throw std::runtime_error("query label outside the data alphabet");
                }
                ids.push_back(static_cast<std::uint32_t>(it - alphabet.begin()));
            }
            std::sort(ids.begin(), ids.end());
            std::size_t i = 0;
            while (i < ids.size()) {
                std::size_t j = i;
                while (j < ids.size() && ids[j] == ids[i]) {
                    ++j;
                }
                const std::size_t count = j - i;
                query_ball_labels_.emplace_back(
                    ids[i], static_cast<std::uint8_t>(count < 255 ? count : 255));
                i = j;
            }
            query_ball_begin_[u + 1] = query_ball_labels_.size();
        }
    }

    const double tau = index_.tau();
    const std::uint32_t max_id = static_cast<std::uint32_t>(index_.num_vertices());

    for (VertexId u = 0; u < query_size; ++u) {
        std::vector<VertexId>& accepted = candidates[u];
        accepted.clear();

        const LocalSpectralSignature& query_sig = query_signatures[u];
        const Label query_label = query_graph.label(u);
        const std::size_t query_degree = query_graph.degree(u);
        const std::size_t query_dim = query_sig.dimension();
        const double query_lambda_max = query_sig.lambda_max();
        const QueryLabelCount* query_labels =
            label_filter ? query_ball_labels_.data() + query_ball_begin_[u] : nullptr;
        const QueryLabelCount* query_labels_end =
            label_filter ? query_ball_labels_.data() + query_ball_begin_[u + 1] : nullptr;

        if (channel_stats_ && spectral_test) {
            classify_into(query_graph, query_sig, u, query_labels, query_labels_end, accepted);
            sort_vertex_ids(accepted, radix_scratch_, max_id);
            continue;
        }

        if (const SpectralDominanceIndex::EmptyGroup* empty =
                index_.find_empty_group(query_label)) {
            for (VertexId v : empty->vertices) {
                if (data_graph_.degree(v) < query_degree) {
                    continue;
                }
                if (label_filter &&
                    !label_covered(index_.ball_counts(v), query_labels, query_labels_end)) {
                    continue;
                }
                accepted.push_back(v);
            }
        }

        if (spectral_test && (query_dim == 0 || tau <= 0.0)) {
            continue;
        }

        const SpectralDominanceIndex::LabelGroup* group = index_.find_group(query_label);
        if (group == nullptr) {
            continue;
        }

        const AmsdGroup* amsd_group = nullptr;
        Provider provider = Provider::BucketScan;
        std::size_t region_hint = 0;
        const std::size_t position = index_.group_position(query_label);
        const AmsdGroup* candidate_group = index_.amsd_group_at(position);
        if (candidate_group != nullptr && candidate_group->size == group->vertices.size()) {
            amsd_group = candidate_group;
            if (spectral_test) {
                provider = choose_provider(*group, *candidate_group, query_sig, &region_hint);
            }
        }
        if (provider != Provider::BitTest) {
            amsd_channels_.clear();
        }

        if (amsd_group != nullptr) {
            ++amsd_active_vertices_;
            if (provider == Provider::Bitmap) {
                ++amsd_bitmap_vertices_;
            } else if (provider == Provider::BitTest) {
                ++amsd_bit_test_vertices_;
            } else {
                ++amsd_bucket_vertices_;
            }
        }

        const auto dispatch_start =
            stage_timings_ ? std::chrono::steady_clock::now()
                           : std::chrono::steady_clock::time_point();
        if (provider == Provider::Bitmap) {
            retrieve_amsd_into(*group, *amsd_group, query_sig, query_degree,
                               query_lambda_max, query_labels, query_labels_end, accepted);
        } else {
            scan_buckets_into(*group, query_sig, query_degree, query_lambda_max, query_labels,
                              query_labels_end, region_hint, accepted);
        }
        if (stage_timings_) {
            const auto dispatch_end = std::chrono::steady_clock::now();
            dispatch_ns_ += static_cast<std::uint64_t>(std::chrono::duration_cast<
                std::chrono::nanoseconds>(dispatch_end - dispatch_start).count());
        }

        sort_vertex_ids(accepted, radix_scratch_, max_id);
    }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
void SpectralFilter::collect_rank_channels(const AmsdGroup& amsd_group,
                                           const LocalSpectralSignature& query_sig) const {
    amsd_channels_.clear();
    const std::size_t query_dim = query_sig.dimension();
    const std::size_t group_size = amsd_group.size;
    for (const AmsdRankIndex& rank : amsd_group.ranks) {
        if (rank.rank > query_dim) {
            continue;
        }
        const double threshold =
            query_sig.eigenvalues[query_dim - rank.rank] - DOMINANCE_EPSILON;
        const int rung = rank.rung_for(threshold);
        if (rung < 0) {
            continue;
        }
        const std::size_t index = static_cast<std::size_t>(rung);
        amsd_channels_.push_back({&rank.suffix[index], rank.rank});
        amsd_channels_.back().suffix_size = rank.suffix_sets[index];
        (void)group_size;
    }
    if (amsd_channels_.size() > 1) {
        std::size_t best = 0;
        for (std::size_t i = 1; i < amsd_channels_.size(); ++i) {
            if (amsd_channels_[i].suffix_size < amsd_channels_[best].suffix_size) {
                best = i;
            }
        }
        std::swap(amsd_channels_[0], amsd_channels_[best]);
    }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
SpectralFilter::Provider SpectralFilter::choose_provider(
    const SpectralDominanceIndex::LabelGroup& group, const AmsdGroup& amsd_group,
    const LocalSpectralSignature& query_sig, std::size_t* region_out) const {
    const std::size_t group_size = group.vertices.size();
    if (region_out != nullptr) {
        *region_out = 0;
    }
    collect_rank_channels(amsd_group, query_sig);
    if (amsd_channels_.empty()) {
        return Provider::BucketScan;
    }
    const double keep_fraction = group_size == 0
                                     ? 1.0
                                     : static_cast<double>(amsd_channels_.front().suffix_size) /
                                           static_cast<double>(group_size);

    const double words = static_cast<double>((group_size + 63) / 64);
    const double group_potential =
        static_cast<double>(group_size) * (1.0 - corrected_keep_fraction(keep_fraction)) /
        kChannelCorrelationPenalty;
    if (group_potential * kScanNsPerCandidate > kAmsdSetupNs + words) {
        std::size_t region = 0;
        for (const SpectralDominanceIndex::Bucket& bucket : group.buckets) {
            if (bucket.dimension >= static_cast<std::uint32_t>(query_sig.dimension())) {
                region += static_cast<std::size_t>(bucket.count);
            }
        }
        if (region_out != nullptr) {
            *region_out = region;
        }
        if (amsd_pays_off(region, keep_fraction, (group_size + 63) / 64)) {
            return Provider::Bitmap;
        }
    }
    if (amsd_bit_test_pays_off(keep_fraction, amsd_channels_.size())) {
        return Provider::BitTest;
    }
    return Provider::BucketScan;
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
void SpectralFilter::retrieve_amsd_into(
    const SpectralDominanceIndex::LabelGroup& group, const AmsdGroup& amsd_group,
    const LocalSpectralSignature& query_sig, std::size_t query_degree,
    double query_lambda_max, const QueryLabelCount* query_labels,
    const QueryLabelCount* query_labels_end, std::vector<VertexId>& accepted) const {
    const std::size_t group_size = group.vertices.size();
    const std::size_t query_dim = query_sig.dimension();
    const bool spectral_test = spectral_test_enabled_;
    const bool label_filter = query_labels != nullptr;
    ++amsd_active_vertices_;

    if (!spectral_test) {
        accepted.reserve(accepted.size() + group_size);
        for (std::size_t i = 0; i < group_size; ++i) {
            const VertexId v = group.vertices[i];
            if (data_graph_.degree(v) < query_degree) {
                continue;
            }
            if (label_filter &&
                !label_covered(index_.ball_counts(v), query_labels, query_labels_end)) {
                continue;
            }
            accepted.push_back(v);
        }
        return;
    }

    amsd_plan_.dimension = static_cast<std::uint32_t>(query_dim);
    amsd_plan_.ranks.clear();
    for (const AmsdRankIndex& rank : amsd_group.ranks) {
        if (rank.rank > query_dim) {
            continue;
        }
        const double threshold =
            query_sig.eigenvalues[query_dim - rank.rank] - DOMINANCE_EPSILON;
        amsd_plan_.ranks.emplace_back(rank.rank, threshold);
    }

    amsd_retriever_.retrieve(group, amsd_group, amsd_plan_, amsd_adaptive_, amsd_result_,
                             &amsd_stats_);

    amsd_survivors_.clear();
    const std::size_t survivor_count = amsd_result_.count();
    amsd_survivors_.reserve(survivor_count);
    for (std::size_t i = amsd_result_.next_set_bit(0); i < group_size;
         i = amsd_result_.next_set_bit(i + 1)) {
        amsd_survivors_.push_back(group.vertices[i]);
    }
    accepted.reserve(accepted.size() + survivor_count);

    const bool use_index_degrees = use_index_degrees_;
    const std::uint32_t* degrees = use_index_degrees ? index_.amsd().degrees.data() : nullptr;

    for (std::size_t pos = 0; pos < survivor_count; ++pos) {
        if (pos + kHeaderLookahead < survivor_count) {
            SLIM_PREFETCH(&index_.signature(amsd_survivors_[pos + kHeaderLookahead]));
        }
        if (pos + kDataLookahead < survivor_count) {
            const LocalSpectralSignature& ahead =
                index_.signature(amsd_survivors_[pos + kDataLookahead]);
            if (!ahead.eigenvalues.empty()) {
                SLIM_PREFETCH(ahead.eigenvalues.data() + ahead.eigenvalues.size() - 1);
            }
        }
        if (label_filter && pos + kHeaderLookahead < survivor_count) {
            SLIM_PREFETCH(index_.ball_counts(amsd_survivors_[pos + kHeaderLookahead]));
        }

        const VertexId v = amsd_survivors_[pos];
        const std::size_t degree = degrees != nullptr ? degrees[v] : data_graph_.degree(v);
        if (degree < query_degree) {
            continue;
        }
        ++coarse_dominance_checks_;
        if (label_filter &&
            !label_covered(index_.ball_counts(v), query_labels, query_labels_end)) {
            continue;
        }
        const LocalSpectralSignature& data_sig = index_.signature(v);
        if (data_sig.lambda_max() + DOMINANCE_EPSILON < query_lambda_max) {
            continue;
        }
        ++full_dominance_checks_;
        if (data_sig.dominates(query_sig)) {
            accepted.push_back(v);
        }
    }
}

void SpectralFilter::scan_buckets_into(
    const SpectralDominanceIndex::LabelGroup& group, const LocalSpectralSignature& query_sig,
    std::size_t query_degree, double query_lambda_max, const QueryLabelCount* query_labels,
    const QueryLabelCount* query_labels_end, std::size_t region_hint,
    std::vector<VertexId>& accepted) const {
    const std::vector<RankChannel>& rank_channels = amsd_channels_;
    const bool label_filter = query_labels != nullptr;
    const bool spectral_test = spectral_test_enabled_;
    const std::size_t query_dim = query_sig.dimension();
    const double tau = index_.tau();
    const bool has_ball_labels = label_filter;
    const bool chain_available = !rank_channels.empty() && spectral_test;
    bool chain_enabled = false;
    const double chain_cost_ns =
        static_cast<double>(rank_channels.size()) * kRankBitTestNs;
    std::size_t chain_sampled = 0;
    std::size_t chain_rejected = 0;
    std::size_t chain_until_sample = 0;

    const double shifted_lambda_max =
        std::max(0.0, query_lambda_max - DOMINANCE_EPSILON);
    const std::uint32_t query_layer =
        static_cast<std::uint32_t>(std::floor(shifted_lambda_max / tau));
    const std::uint32_t query_dimension = static_cast<std::uint32_t>(query_dim);

    auto bucket_it =
        spectral_test
            ? std::lower_bound(
                  group.buckets.begin(), group.buckets.end(), query_layer,
                  [](const SpectralDominanceIndex::Bucket& bucket, std::uint32_t layer) {
                      return bucket.layer < layer;
                  })
            : group.buckets.begin();

    std::size_t capacity = accepted.size() + region_hint;
    if (region_hint == 0) {
        for (auto it = bucket_it; it != group.buckets.end(); ++it) {
            if (it->dimension >= query_dimension) {
                capacity += it->count;
            }
        }
    }
    accepted.reserve(capacity);

    if (chain_enabled) {
        ++amsd_stats_.chain_enabled_scans;
    }
    for (; bucket_it != group.buckets.end(); ++bucket_it) {
        const SpectralDominanceIndex::Bucket& bucket = *bucket_it;
        ++scanned_buckets_;
        scanned_cells_ += bucket.count;
        if (spectral_test && bucket.dimension < query_dimension) {
            continue;
        }
        const VertexId* cells = group.vertices.data() + bucket.begin;
        const std::size_t bucket_begin = static_cast<std::size_t>(bucket.begin);
        for (std::uint64_t pos = 0; pos < bucket.count; ++pos) {
            if (pos + kHeaderLookahead < bucket.count) {
                SLIM_PREFETCH(&index_.signature(cells[pos + kHeaderLookahead]));
            }
            if (pos + kDataLookahead < bucket.count) {
                const LocalSpectralSignature& ahead = index_.signature(cells[pos + kDataLookahead]);
                if (!ahead.eigenvalues.empty()) {
                    SLIM_PREFETCH(ahead.eigenvalues.data() + ahead.eigenvalues.size() - 1);
                }
            }
            if (has_ball_labels && pos + kHeaderLookahead < bucket.count) {
                SLIM_PREFETCH(index_.ball_counts(cells[pos + kHeaderLookahead]));
            }
            const VertexId v = cells[pos];
            if (chain_available) {
                bool test_now = chain_enabled;
                if (!test_now) {
                    if (chain_until_sample == 0) {
                        test_now = true;
                        chain_until_sample = kChainSampleStride;
                    }
                    --chain_until_sample;
                }
                if (test_now) {
                    bool passed = true;
                    const std::size_t local = bucket_begin + static_cast<std::size_t>(pos);
                    for (const RankChannel& channel : rank_channels) {
                        if (!channel.suffix->test(local)) {
                            passed = false;
                            break;
                        }
                    }
                    ++chain_sampled;
                    ++amsd_stats_.chain_tests;
                    if (!passed) {
                        ++chain_rejected;
                    }
                    chain_enabled = chain_rejected * kPreemptedCandidateNs >
                                    chain_cost_ns * static_cast<double>(chain_sampled);
                    if (!passed) {
                        ++amsd_stats_.chain_rejections;
                        continue;
                    }
                }
            }
            if (data_graph_.degree(v) < query_degree) {
                continue;
            }
            if (spectral_test) {
                ++coarse_dominance_checks_;
            }
            if (label_filter &&
                !label_covered(index_.ball_counts(v), query_labels, query_labels_end)) {
                continue;
            }
            if (!spectral_test) {
                accepted.push_back(v);
                continue;
            }
            const LocalSpectralSignature& data_sig = index_.signature(v);
            if (data_sig.lambda_max() + DOMINANCE_EPSILON < query_lambda_max) {
                continue;
            }
            ++full_dominance_checks_;
            if (data_sig.dominates(query_sig)) {
                accepted.push_back(v);
            }
        }
    }
}

std::size_t SpectralFilter::count_label_degree_candidates(
    const LabeledGraph& query_graph) const {
    std::size_t total = 0;
    for (VertexId u = 0; u < query_graph.num_vertices(); ++u) {
        const Label query_label = query_graph.label(u);
        const std::size_t query_degree = query_graph.degree(u);
        if (const SpectralDominanceIndex::LabelGroup* group = index_.find_group(query_label)) {
            for (VertexId v : group->vertices) {
                if (data_graph_.degree(v) >= query_degree) {
                    ++total;
                }
            }
        }
        if (const SpectralDominanceIndex::EmptyGroup* empty =
                index_.find_empty_group(query_label)) {
            for (VertexId v : empty->vertices) {
                if (data_graph_.degree(v) >= query_degree) {
                    ++total;
                }
            }
        }
    }
    return total;
}

void SpectralFilter::classify_into(const LabeledGraph& query_graph,
                                   const LocalSpectralSignature& query_sig,
                                   VertexId u,
                                   const QueryLabelCount* query_labels,
                                   const QueryLabelCount* query_labels_end,
                                   std::vector<VertexId>& accepted) const {
    const bool label_filter = query_labels != nullptr;
    const Label query_label = query_graph.label(u);
    const std::size_t query_degree = query_graph.degree(u);
    const std::size_t query_dim = query_sig.dimension();
    const double query_lambda_max = query_sig.lambda_max();
    const double tau = index_.tau();

    if (const SpectralDominanceIndex::EmptyGroup* empty = index_.find_empty_group(query_label)) {
        for (VertexId v : empty->vertices) {
            if (data_graph_.degree(v) < query_degree) {
                continue;
            }
            ++ld_candidates_;
            if (label_filter &&
                !label_covered(index_.ball_counts(v), query_labels, query_labels_end)) {
                ++label_rejects_;
                continue;
            }
            accepted.push_back(v);
            ++empty_accepts_;
        }
    }

    if (query_dim == 0 || tau <= 0.0) {
        return;
    }
    const SpectralDominanceIndex::LabelGroup* group = index_.find_group(query_label);
    if (group == nullptr) {
        return;
    }

    const std::uint32_t query_dimension = static_cast<std::uint32_t>(query_dim);
    const std::uint32_t query_layer = static_cast<std::uint32_t>(std::floor(
        std::max(0.0, query_lambda_max - DOMINANCE_EPSILON) / tau));

    for (const SpectralDominanceIndex::Bucket& bucket : group->buckets) {
        const bool reachable = bucket.layer >= query_layer;
        if (reachable) {
            ++scanned_buckets_;
            scanned_cells_ += bucket.count;
        }
        const VertexId* cells = group->vertices.data() + bucket.begin;
        for (std::uint64_t pos = 0; pos < bucket.count; ++pos) {
            const VertexId v = cells[pos];
            if (data_graph_.degree(v) < query_degree) {
                continue;
            }
            ++ld_candidates_;
            if (bucket.dimension < query_dimension) {
                ++dim_rejects_;
                continue;
            }
            if (reachable) {
                ++coarse_dominance_checks_;
            }
            if (label_filter &&
                !label_covered(index_.ball_counts(v), query_labels, query_labels_end)) {
                ++label_rejects_;
                continue;
            }
            const LocalSpectralSignature& data_sig = index_.signature(v);
            if (data_sig.lambda_max() + DOMINANCE_EPSILON < query_lambda_max) {
                ++key_rejects_;
                continue;
            }
            if (reachable) {
                ++full_dominance_checks_;
            }
            if (data_sig.dominates(query_sig)) {
                accepted.push_back(v);
            } else {
                ++win_rejects_;
            }
        }
    }
}

LocalSpectralSignature SpectralFilter::compute_query_signature(
    const LabeledGraph& query_graph,
    VertexId u,
    KHopNeighborhood& neighborhood) const {
    const LabelSpectralOperator op(index_.alpha(), encoder_);
    return compute_local_spectral_signature(query_graph, op, query_radius(), u, neighborhood);
}
}
