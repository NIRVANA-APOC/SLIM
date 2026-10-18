#include <slim/rank_selection.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

namespace slim {
namespace {
constexpr double kTieEps = 1e-12;

std::uint32_t stratum_of(std::size_t dimension) {
    std::uint32_t z = 0;
    std::size_t value = dimension;
    while (value > 1) {
        value >>= 1;
        ++z;
    }
    return z;
}

std::uint32_t bucket_of(double value, double tau, std::uint32_t cols) {
    if (tau <= 0.0) {
        return 0;
    }
    const double raw = std::floor(value / tau);
    if (raw <= 0.0) {
        return 0;
    }
    if (raw >= static_cast<double>(cols)) {
        return cols == 0 ? 0 : cols - 1;
    }
    return static_cast<std::uint32_t>(raw);
}

struct RankLayout {
    const std::vector<LocalSpectralSignature>* signatures = nullptr;
    std::vector<std::uint32_t> dims;
    std::vector<VertexId> order;
    std::vector<std::uint64_t> c;
    double tau = 0.0;
    std::size_t num_vertices = 0;
    std::uint32_t max_rank = 0;
    std::uint32_t max_dim = 0;
    std::uint32_t num_strata = 0;
    std::uint32_t cols = 0;
    double max_lambda = 0.0;
};

RankLayout build_layout(const std::vector<LocalSpectralSignature>& signatures, double tau,
                        std::uint32_t max_rank_cap) {
    RankLayout layout;
    layout.signatures = &signatures;
    layout.num_vertices = signatures.size();
    layout.tau = tau;

    layout.dims.resize(layout.num_vertices);
    layout.max_dim = 0;
    layout.max_lambda = 0.0;
    for (std::size_t v = 0; v < layout.num_vertices; ++v) {
        const std::uint32_t dimension = static_cast<std::uint32_t>(signatures[v].dimension());
        layout.dims[v] = dimension;
        layout.max_dim = std::max(layout.max_dim, dimension);
        if (dimension > 0) {
            layout.max_lambda =
                std::max(layout.max_lambda, stored_spectral_value(signatures[v].lambda_max()));
        }
    }
    layout.max_rank = std::min(max_rank_cap, layout.max_dim);

    layout.order.resize(layout.num_vertices);
    std::iota(layout.order.begin(), layout.order.end(), VertexId{0});
    std::sort(layout.order.begin(), layout.order.end(), [&layout](VertexId a, VertexId b) {
        if (layout.dims[a] != layout.dims[b]) {
            return layout.dims[a] > layout.dims[b];
        }
        return a < b;
    });

    layout.c.assign(static_cast<std::size_t>(layout.max_rank) + 1, 0);
    std::vector<std::uint64_t> histogram(static_cast<std::size_t>(layout.max_dim) + 1, 0);
    for (std::size_t v = 0; v < layout.num_vertices; ++v) {
        ++histogram[layout.dims[v]];
    }
    std::uint64_t running = 0;
    for (std::uint32_t dimension = layout.max_dim; dimension >= 1; --dimension) {
        running += histogram[dimension];
        if (dimension <= layout.max_rank) {
            layout.c[dimension] = running;
        }
        if (dimension == 1) {
            break;
        }
    }

    layout.num_strata = stratum_of(layout.max_dim) + 1;
    layout.cols = 0;
    if (tau > 0.0 && layout.max_lambda > 0.0) {
        const double raw_cols = std::floor(layout.max_lambda / tau) + 1.0;
        constexpr double kMaxCols = 1u << 21;
        layout.cols = static_cast<std::uint32_t>(std::min(raw_cols, kMaxCols));
        if (layout.cols == 0) {
            layout.cols = 1;
        }
    }
    if (layout.cols == 0) {
        layout.cols = 1;
    }
    return layout;
}

double conditional_entropy(const RankLayout& layout, std::uint32_t rank, std::size_t limit,
                           std::vector<std::uint32_t>& cells,
                           std::vector<std::uint64_t>& row_sums) {
    if (limit < 2) {
        return 0.0;
    }
    const std::size_t cell_count =
        static_cast<std::size_t>(layout.num_strata) * static_cast<std::size_t>(layout.cols);
    cells.assign(cell_count, 0);
    row_sums.assign(layout.num_strata, 0);
    const auto& signatures = *layout.signatures;
    for (std::size_t i = 0; i < limit; ++i) {
        const VertexId v = layout.order[i];
        const std::uint32_t dimension = layout.dims[v];
        const std::uint32_t z = stratum_of(dimension);
        const double value =
            stored_spectral_value(signatures[v].eigenvalues[dimension - rank]);
        const std::uint32_t bucket = bucket_of(value, layout.tau, layout.cols);
        ++cells[static_cast<std::size_t>(z) * layout.cols + bucket];
    }
    for (std::uint32_t z = 0; z < layout.num_strata; ++z) {
        std::uint64_t row = 0;
        for (std::uint32_t b = 0; b < layout.cols; ++b) {
            row += cells[static_cast<std::size_t>(z) * layout.cols + b];
        }
        row_sums[z] = row;
    }
    double sum = 0.0;
    for (std::uint32_t z = 0; z < layout.num_strata; ++z) {
        const std::uint64_t row = row_sums[z];
        if (row == 0) {
            continue;
        }
        for (std::uint32_t b = 0; b < layout.cols; ++b) {
            const std::uint64_t count = cells[static_cast<std::size_t>(z) * layout.cols + b];
            if (count != 0) {
                sum += static_cast<double>(count) *
                       std::log(static_cast<double>(row) / static_cast<double>(count));
            }
        }
    }
    return sum / static_cast<double>(limit);
}

void compute_rank_stats(const RankLayout& layout, RankSelectionResult& result) {
    const std::size_t rank_count = static_cast<std::size_t>(layout.max_rank);
    result.c_j.assign(rank_count, 0);
    result.h_j.assign(rank_count, 0.0);
    result.d_j.assign(rank_count, 0.0);
    for (std::uint32_t rank = 1; rank <= layout.max_rank; ++rank) {
        result.c_j[rank - 1] = layout.c[rank];
    }
    std::vector<std::uint32_t> cells;
    std::vector<std::uint64_t> row_sums;
    for (std::uint32_t rank = 1; rank <= layout.max_rank; ++rank) {
        const std::uint64_t count = layout.c[rank];
        if (count < 2) {
            continue;
        }
        const double entropy = conditional_entropy(layout, rank, static_cast<std::size_t>(count),
                                                   cells, row_sums);
        result.h_j[rank - 1] = entropy;
        result.d_j[rank - 1] = entropy > 0.0 ? entropy / std::log(static_cast<double>(count)) : 0.0;
    }
}

double conditional_nmi(const RankLayout& layout, std::uint32_t rank_i, std::uint32_t rank_j) {
    const std::uint32_t max_rank = std::max(rank_i, rank_j);
    const std::size_t limit = layout.c[max_rank];
    if (limit < 2) {
        return 0.0;
    }
    const auto& signatures = *layout.signatures;
    const std::uint32_t cols = layout.cols;
    std::vector<std::uint64_t> keys;
    keys.reserve(limit);
    for (std::size_t i = 0; i < limit; ++i) {
        const VertexId v = layout.order[i];
        const std::uint32_t dimension = layout.dims[v];
        const std::uint32_t z = stratum_of(dimension);
        const std::uint32_t bucket_i = bucket_of(
            stored_spectral_value(signatures[v].eigenvalues[dimension - rank_i]), layout.tau, cols);
        const std::uint32_t bucket_j = bucket_of(
            stored_spectral_value(signatures[v].eigenvalues[dimension - rank_j]), layout.tau, cols);
        const std::uint64_t key =
            (static_cast<std::uint64_t>(z) * cols + bucket_i) * cols + bucket_j;
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());

    std::vector<std::uint32_t> hist_i(static_cast<std::size_t>(layout.num_strata) * cols, 0);
    std::vector<std::uint32_t> hist_j(hist_i.size(), 0);
    std::vector<std::uint64_t> row_sums(layout.num_strata, 0);
    struct Run {
        std::uint64_t key;
        std::uint32_t count;
    };
    std::vector<Run> runs;
    runs.reserve(std::min<std::size_t>(keys.size(), 1u << 16));
    for (std::size_t begin = 0; begin < keys.size();) {
        std::size_t end = begin + 1;
        while (end < keys.size() && keys[end] == keys[begin]) {
            ++end;
        }
        const std::uint64_t key = keys[begin];
        const std::uint32_t bucket_j = static_cast<std::uint32_t>(key % cols);
        const std::uint32_t bucket_i = static_cast<std::uint32_t>((key / cols) % cols);
        const std::uint32_t z = static_cast<std::uint32_t>(key / (static_cast<std::uint64_t>(cols) * cols));
        const std::uint32_t count = static_cast<std::uint32_t>(end - begin);
        hist_i[static_cast<std::size_t>(z) * cols + bucket_i] += count;
        hist_j[static_cast<std::size_t>(z) * cols + bucket_j] += count;
        row_sums[z] += count;
        runs.push_back({key, count});
        begin = end;
    }

    const double n = static_cast<double>(limit);
    double entropy_i = 0.0;
    double entropy_j = 0.0;
    for (std::uint32_t z = 0; z < layout.num_strata; ++z) {
        const std::uint64_t row = row_sums[z];
        if (row == 0) {
            continue;
        }
        for (std::uint32_t b = 0; b < cols; ++b) {
            const std::uint64_t ci = hist_i[static_cast<std::size_t>(z) * cols + b];
            if (ci != 0) {
                entropy_i += static_cast<double>(ci) *
                             std::log(static_cast<double>(row) / static_cast<double>(ci));
            }
            const std::uint64_t cj = hist_j[static_cast<std::size_t>(z) * cols + b];
            if (cj != 0) {
                entropy_j += static_cast<double>(cj) *
                             std::log(static_cast<double>(row) / static_cast<double>(cj));
            }
        }
    }
    entropy_i /= n;
    entropy_j /= n;
    double joint_entropy = 0.0;
    for (const Run& run : runs) {
        const std::uint32_t z =
            static_cast<std::uint32_t>(run.key / (static_cast<std::uint64_t>(cols) * cols));
        const std::uint64_t row = row_sums[z];
        joint_entropy += static_cast<double>(run.count) *
                         std::log(static_cast<double>(row) / static_cast<double>(run.count));
    }
    joint_entropy /= n;
    const double mutual = entropy_i + entropy_j - joint_entropy;
    const double denominator = std::sqrt(std::max(0.0, entropy_i) * std::max(0.0, entropy_j));
    if (denominator <= kTieEps) {
        return 0.0;
    }
    const double nmi = mutual / denominator;
    if (!(nmi > 0.0)) {
        return 0.0;
    }
    return std::min(1.0, nmi);
}

std::vector<std::uint32_t> uniform_ranks(const RankLayout& layout, double budget) {
    std::vector<std::uint32_t> best;
    std::uint64_t best_space = 0;
    for (std::uint32_t k = 1; k <= layout.max_rank; ++k) {
        std::vector<std::uint32_t> ranks;
        ranks.reserve(k);
        if (k == 1) {
            ranks.push_back(1);
        } else {
            const double step = static_cast<double>(layout.max_rank - 1) / static_cast<double>(k - 1);
            for (std::uint32_t i = 0; i < k; ++i) {
                const auto rank = static_cast<std::uint32_t>(
                    1.0 + std::llround(static_cast<double>(i) * step));
                if (rank >= 1 && rank <= layout.max_rank) {
                    ranks.push_back(rank);
                }
            }
            std::sort(ranks.begin(), ranks.end());
            ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
        }
        std::uint64_t space = 0;
        for (std::uint32_t rank : ranks) {
            space += layout.c[rank];
        }
        const bool fits = static_cast<double>(space) <= budget + 0.5;
        if (fits && (ranks.size() > best.size() ||
                     (ranks.size() == best.size() && space < best_space))) {
            best = std::move(ranks);
            best_space = space;
        }
    }
    return best;
}
}

RankSelectionResult select_rank_set(const std::vector<LocalSpectralSignature>& signatures,
                                    const RankSelectionOptions& options) {
    RankSelectionResult result;
    result.mode = options.mode;
    result.tau = options.tau;
    result.num_vertices = options.num_vertices;
    result.budget_entries = options.budget_entries;
    result.max_rank = std::min<std::uint32_t>(options.max_rank, 1u << 20);

    if (options.mode == RankSelectMode::All) {
        result.selected_ranks = options.fallback_ranks;
        std::sort(result.selected_ranks.begin(), result.selected_ranks.end());
        result.selected_ranks.erase(
            std::unique(result.selected_ranks.begin(), result.selected_ranks.end()),
            result.selected_ranks.end());
        result.selection_order = result.selected_ranks;
        const RankLayout layout = build_layout(signatures, options.tau, result.max_rank);
        result.max_rank = layout.max_rank;
        compute_rank_stats(layout, result);
        std::uint64_t space = 0;
        for (std::uint32_t rank : result.selected_ranks) {
            if (rank <= layout.max_rank) {
                space += layout.c[rank];
                result.selected_entries.push_back(layout.c[rank]);
            }
        }
        result.space_entries = space;
        return result;
    }

    const RankLayout layout = build_layout(signatures, options.tau, result.max_rank);
    result.max_rank = layout.max_rank;
    compute_rank_stats(layout, result);

    if (layout.max_rank == 0 || layout.num_vertices == 0) {
        return result;
    }

    if (options.mode == RankSelectMode::Rank1) {
        if (layout.c[1] > 0) {
            result.selected_ranks = {1};
            result.selected_entries = {layout.c[1]};
            result.space_entries = layout.c[1];
            result.selection_order = {1};
        }
        return result;
    }

    if (options.mode == RankSelectMode::Adaptive) {
        std::map<std::uint64_t, double> nmi_cache;
        const auto nmi_of = [&](std::uint32_t rank_a, std::uint32_t rank_b) {
            const std::uint32_t lo = std::min(rank_a, rank_b);
            const std::uint32_t hi = std::max(rank_a, rank_b);
            const std::uint64_t key =
                (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
            const auto it = nmi_cache.find(key);
            if (it != nmi_cache.end()) {
                return it->second;
            }
            const double value = conditional_nmi(layout, lo, hi);
            nmi_cache.emplace(key, value);
            result.nmi_calls.push_back({lo, hi, value});
            return value;
        };

        struct Candidate {
            std::uint32_t rank;
            double bound;
        };
        std::vector<Candidate> candidates;
        candidates.reserve(layout.max_rank);
        for (std::uint32_t rank = 1; rank <= layout.max_rank; ++rank) {
            const std::uint64_t count = layout.c[rank];
            if (count < 2 || result.d_j[rank - 1] <= 0.0) {
                continue;
            }
            if (static_cast<double>(count) > options.budget_entries + 0.5) {
                continue;
            }
            const double bound = static_cast<double>(count) /
                                 static_cast<double>(layout.num_vertices) *
                                 result.d_j[rank - 1];
            candidates.push_back({rank, bound});
        }
        std::sort(candidates.begin(), candidates.end(), [&result](const Candidate& a,
                                                                  const Candidate& b) {
            if (a.bound != b.bound) {
                return a.bound > b.bound;
            }
            const std::uint64_t ca = result.c_j[a.rank - 1];
            const std::uint64_t cb = result.c_j[b.rank - 1];
            if (ca != cb) {
                return ca > cb;
            }
            return a.rank < b.rank;
        });

        std::vector<bool> selected(static_cast<std::size_t>(layout.max_rank) + 1, false);
        std::uint64_t space = 0;
        while (true) {
            std::uint32_t best_rank = 0;
            double best_delta = 0.0;
            for (const Candidate& candidate : candidates) {
                if (selected[candidate.rank]) {
                    continue;
                }
                if (static_cast<double>(space) + static_cast<double>(layout.c[candidate.rank]) >
                    options.budget_entries + 0.5) {
                    continue;
                }
                if (best_rank != 0 && candidate.bound < best_delta - kTieEps) {
                    break;
                }
                double max_nmi = 0.0;
                for (std::uint32_t chosen : result.selected_ranks) {
                    max_nmi = std::max(max_nmi, nmi_of(chosen, candidate.rank));
                }
                const double delta = candidate.bound * (1.0 - max_nmi);
                if (delta > best_delta + kTieEps) {
                    best_delta = delta;
                    best_rank = candidate.rank;
                }
            }
            if (best_rank == 0) {
                break;
            }
            selected[best_rank] = true;
            result.selected_ranks.push_back(best_rank);
            result.selection_order.push_back(best_rank);
            result.selected_entries.push_back(layout.c[best_rank]);
            space += layout.c[best_rank];
            if (static_cast<double>(space) + 1.0 > options.budget_entries + 0.5) {
                break;
            }
        }
        std::sort(result.selected_ranks.begin(), result.selected_ranks.end());
        result.space_entries = space;
        return result;
    }

    result.selected_ranks = uniform_ranks(layout, options.budget_entries);
    result.selection_order = result.selected_ranks;
    std::uint64_t space = 0;
    for (std::uint32_t rank : result.selected_ranks) {
        space += layout.c[rank];
        result.selected_entries.push_back(layout.c[rank]);
    }
    result.space_entries = space;
    return result;
}

double normalized_conditional_entropy(const std::vector<LocalSpectralSignature>& signatures,
                                      std::size_t num_vertices, double tau, std::uint32_t rank,
                                      std::uint32_t* count_out) {
    RankLayout layout = build_layout(signatures, tau, rank);
    if (rank == 0 || rank > layout.max_rank || layout.num_vertices != num_vertices) {
        if (count_out != nullptr) {
            *count_out = 0;
        }
        return 0.0;
    }
    if (count_out != nullptr) {
        *count_out = static_cast<std::uint32_t>(layout.c[rank]);
    }
    std::vector<std::uint32_t> cells;
    std::vector<std::uint64_t> row_sums;
    return conditional_entropy(layout, rank, static_cast<std::size_t>(layout.c[rank]), cells,
                               row_sums);
}

double normalized_conditional_mutual_information(
    const std::vector<LocalSpectralSignature>& signatures, std::size_t num_vertices, double tau,
    std::uint32_t rank_i, std::uint32_t rank_j) {
    if (rank_i == 0 || rank_j == 0) {
        return 0.0;
    }
    RankLayout layout = build_layout(signatures, tau, std::max(rank_i, rank_j));
    if (layout.num_vertices != num_vertices || layout.max_rank < std::max(rank_i, rank_j)) {
        return 0.0;
    }
    return conditional_nmi(layout, rank_i, rank_j);
}
}
