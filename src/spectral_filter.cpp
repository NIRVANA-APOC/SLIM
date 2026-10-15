#include <slim/spectral_filter.h>

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
    : index_(index), data_graph_(data_graph), encoder_(encoder) {}

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
    const bool has_ball_labels = label_filter;

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

        const double shifted_lambda_max =
            std::max(0.0, query_lambda_max - DOMINANCE_EPSILON);
        const std::uint32_t query_layer =
            static_cast<std::uint32_t>(std::floor(shifted_lambda_max / tau));
        const std::uint32_t query_dimension = static_cast<std::uint32_t>(query_dim);

        auto bucket_it =
            spectral_test
                ? std::lower_bound(
                      group->buckets.begin(), group->buckets.end(), query_layer,
                      [](const SpectralDominanceIndex::Bucket& bucket, std::uint32_t layer) {
                          return bucket.layer < layer;
                      })
                : group->buckets.begin();

        std::size_t capacity = accepted.size();
        for (auto it = bucket_it; it != group->buckets.end(); ++it) {
            if (it->dimension >= query_dimension) {
                capacity += it->count;
            }
        }
        accepted.reserve(capacity);

        for (; bucket_it != group->buckets.end(); ++bucket_it) {
            const SpectralDominanceIndex::Bucket& bucket = *bucket_it;
            ++scanned_buckets_;
            scanned_cells_ += bucket.count;
            if (spectral_test && bucket.dimension < query_dimension) {
                continue;
            }
            const VertexId* cells = group->vertices.data() + bucket.begin;
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

        sort_vertex_ids(accepted, radix_scratch_, max_id);
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
