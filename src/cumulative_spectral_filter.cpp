#include <slim/cumulative_spectral_filter.h>

#include <slim/candidate_dump.h>
#include <slim/graph_fingerprint.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace slim {
namespace {
std::uint64_t micros_since(std::chrono::high_resolution_clock::time_point start,
                           std::chrono::high_resolution_clock::time_point end) {
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
}

std::size_t count_candidates(const CandidateSets& candidates) {
    std::size_t total = 0;
    for (const auto& set : candidates) {
        total += set.size();
    }
    return total;
}

void intersect_into(const CandidateSets& left, const CandidateSets& right, CandidateSets& out) {
    if (left.size() != right.size()) {
        throw std::runtime_error("candidate set size mismatch between layers");
    }
    out.resize(left.size());
    for (std::size_t u = 0; u < left.size(); ++u) {
        std::vector<VertexId>& dst = out[u];
        dst.clear();
        const std::vector<VertexId>& a = left[u];
        const std::vector<VertexId>& b = right[u];
        std::size_t i = 0;
        std::size_t j = 0;
        while (i < a.size() && j < b.size()) {
            if (a[i] < b[j]) {
                ++i;
            } else if (a[i] > b[j]) {
                ++j;
            } else {
                dst.push_back(a[i]);
                ++i;
                ++j;
            }
        }
    }
}

std::size_t count_not_contained(const CandidateSets& subset, const CandidateSets& superset) {
    std::size_t missing = 0;
    const std::size_t sets = std::min(subset.size(), superset.size());
    for (std::size_t u = 0; u < sets; ++u) {
        const std::vector<VertexId>& sub = subset[u];
        const std::vector<VertexId>& sup = superset[u];
        std::size_t i = 0;
        std::size_t j = 0;
        while (i < sub.size()) {
            while (j < sup.size() && sup[j] < sub[i]) {
                ++j;
            }
            if (j == sup.size() || sup[j] != sub[i]) {
                ++missing;
            }
            ++i;
        }
    }
    return missing;
}

std::size_t count_missing_signature_keeps(const CandidateSets& sets,
                                          const SpectralDominanceIndex& index) {
    std::size_t keeps = 0;
    for (const std::vector<VertexId>& set : sets) {
        for (VertexId v : set) {
            if (index.signature(v).eigenvalues.empty()) {
                ++keeps;
            }
        }
    }
    return keeps;
}

std::string radius_list(const std::vector<int>& radii) {
    std::string text;
    for (std::size_t i = 0; i < radii.size(); ++i) {
        if (i != 0) {
            text += ",";
        }
        text += std::to_string(radii[i]);
    }
    return text;
}
}

std::vector<std::unique_ptr<SpectralDominanceIndex>> load_cumulative_index_set(
    const std::vector<std::string>& paths, int max_radius) {
    if (paths.empty()) {
        throw std::runtime_error("cumulative mode needs one --index per radius (1..K)");
    }
    if (max_radius < 0) {
        throw std::runtime_error("--max-radius must be >= 1");
    }

    std::vector<std::pair<int, std::unique_ptr<SpectralDominanceIndex>>> loaded;
    loaded.reserve(paths.size());
    for (const std::string& path : paths) {
        std::unique_ptr<SpectralDominanceIndex> index = SpectralDominanceIndex::load(path);
        if (!index->has_identity()) {
            throw std::runtime_error(
                "index " + path + " has no identity section (built by an older slim_index): its radius, "
                "data graph, alpha and label encoding cannot be verified. Rebuild it with the current "
                "slim_index, or load a single index instead.");
        }
        const SpectralDominanceIndex::Identity& identity = index->identity();
        if (identity.format_version != INDEX_FORMAT_VERSION) {
            throw std::runtime_error(
                "index " + path + " has format version " + std::to_string(identity.format_version) +
                ", this binary reads version " + std::to_string(INDEX_FORMAT_VERSION));
        }
        loaded.emplace_back(index->k(), std::move(index));
    }

    std::sort(loaded.begin(), loaded.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<int> radii;
    radii.reserve(loaded.size());
    for (const auto& entry : loaded) {
        radii.push_back(entry.first);
    }
    for (std::size_t i = 1; i < loaded.size(); ++i) {
        if (loaded[i].first == loaded[i - 1].first) {
            throw std::runtime_error(
                "cumulative index set lists radius " + std::to_string(loaded[i].first) +
                " twice (radii " + radius_list(radii) + "): one layer per radius");
        }
    }
    if (loaded.front().first != 1) {
        throw std::runtime_error(
            "cumulative index set must include the r=1 layer, got radii " + radius_list(radii));
    }

    if (max_radius > 0) {
        const bool present = std::find(radii.begin(), radii.end(), max_radius) != radii.end();
        if (!present) {
            throw std::runtime_error(
                "--max-radius " + std::to_string(max_radius) +
                " is not one of the index set's radii (" + radius_list(radii) + ")");
        }
        loaded.erase(std::remove_if(loaded.begin(), loaded.end(),
                                    [max_radius](const auto& entry) {
                                        return entry.first > max_radius;
                                    }),
                     loaded.end());
        radii.erase(std::remove_if(radii.begin(), radii.end(),
                                   [max_radius](int radius) { return radius > max_radius; }),
                    radii.end());
    }

    const auto* first = loaded.front().second.get();
    const SpectralDominanceIndex::Identity& first_identity = first->identity();
    for (const auto& entry : loaded) {
        const SpectralDominanceIndex* index = entry.second.get();
        const SpectralDominanceIndex::Identity& identity = index->identity();
        if (identity.data_fingerprint != first_identity.data_fingerprint) {
            throw std::runtime_error(
                "index set mixes data graphs: k=" + std::to_string(index->k()) + " has fingerprint " +
                fingerprint_hex(identity.data_fingerprint) + " while k=" + std::to_string(first->k()) +
                " has " + fingerprint_hex(first_identity.data_fingerprint));
        }
        if (identity.label_fingerprint != first_identity.label_fingerprint) {
            throw std::runtime_error(
                "index set mixes label encodings: k=" + std::to_string(index->k()) + " has alphabet " +
                fingerprint_hex(identity.label_fingerprint) + " while k=" + std::to_string(first->k()) +
                " has " + fingerprint_hex(first_identity.label_fingerprint));
        }
        if (index->alpha() != first->alpha()) {
            throw std::runtime_error(
                "index set mixes alpha: k=" + std::to_string(index->k()) + " was built with alpha=" +
                std::to_string(index->alpha()) + " while k=" + std::to_string(first->k()) +
                " used alpha=" + std::to_string(first->alpha()));
        }
        if (identity.tau_scale != first_identity.tau_scale) {
            throw std::runtime_error(
                "index set mixes tau-scale: k=" + std::to_string(index->k()) + " used " +
                std::to_string(identity.tau_scale) + " while k=" + std::to_string(first->k()) + " used " +
                std::to_string(first_identity.tau_scale));
        }
    }

    std::vector<std::unique_ptr<SpectralDominanceIndex>> indexes;
    indexes.reserve(loaded.size());
    for (auto& entry : loaded) {
        indexes.push_back(std::move(entry.second));
    }
    return indexes;
}

void validate_cumulative_index_set(
    const std::vector<std::unique_ptr<SpectralDominanceIndex>>& indexes,
    const LabeledGraph& data_graph) {
    const std::uint64_t graph_fingerprint = data_graph_fingerprint(data_graph);
    const std::uint64_t label_fingerprint = label_alphabet_fingerprint(data_graph);
    for (const auto& index : indexes) {
        const SpectralDominanceIndex::Identity& identity = index->identity();
        if (index->num_vertices() != data_graph.num_vertices() ||
            identity.data_vertices != data_graph.num_vertices() ||
            identity.data_edges != data_graph.num_edges() ||
            identity.data_fingerprint != graph_fingerprint) {
            throw std::runtime_error(
                "index set does not describe the loaded data graph: its k=" + std::to_string(index->k()) +
                " layer was built from " + std::to_string(identity.data_vertices) + " vertices / " +
                std::to_string(identity.data_edges) + " edges (fingerprint " +
                fingerprint_hex(identity.data_fingerprint) + "), the data graph has " +
                std::to_string(data_graph.num_vertices()) + " vertices / " +
                std::to_string(data_graph.num_edges()) + " edges (fingerprint " +
                fingerprint_hex(graph_fingerprint) + ")");
        }
        if (identity.label_fingerprint != label_fingerprint) {
            throw std::runtime_error(
                "index set was built with a different label encoding than the loaded data graph: index k=" +
                std::to_string(index->k()) + " has alphabet " + fingerprint_hex(identity.label_fingerprint) +
                ", the data graph has " + fingerprint_hex(label_fingerprint));
        }
    }
}

CumulativeSpectralFilter::CumulativeSpectralFilter(
    std::vector<std::unique_ptr<SpectralDominanceIndex>> indexes,
    const LabeledGraph& data_graph)
    : indexes_(std::move(indexes)), data_graph_(data_graph) {
    if (indexes_.empty()) {
        throw std::runtime_error("cumulative filter needs at least one index");
    }
    encoder_ = LabelEncoder::from_graph(data_graph_);
    layers_.reserve(indexes_.size());
    for (const auto& index : indexes_) {
        Layer layer;
        layer.index = index.get();
        layer.filter = std::make_unique<SpectralFilter>(*index, data_graph_, encoder_);
        layers_.push_back(std::move(layer));
    }
}

void CumulativeSpectralFilter::filter_into(const LabeledGraph& query_graph, CandidateSets& out,
                                           CumulativeFilterStats* stats) {
    CumulativeFilterStats local;
    if (stats == nullptr) {
        stats = &local;
    }
    *stats = CumulativeFilterStats();

    const std::size_t query_size = query_graph.num_vertices();
    stats->layers.resize(layers_.size());
    for (std::size_t i = 0; i < layers_.size(); ++i) {
        CumulativeLayerStats& layer_stats = stats->layers[i];
        layer_stats.radius = layers_[i].index->k();
        layer_stats.index_k = layers_[i].index->k();
        layer_stats.index_tau = layers_[i].index->tau();
        layer_stats.index_neighborhood_size_limit =
            layers_[i].index->identity().neighborhood_size_limit;
    }

    KHopNeighborhood neighborhood(query_graph);
    for (std::size_t i = 0; i < layers_.size(); ++i) {
        Layer& layer = layers_[i];
        const auto start = std::chrono::high_resolution_clock::now();
        layer.query_signatures.resize(query_size);
        for (VertexId u = 0; u < query_size; ++u) {
            layer.query_signatures[u] =
                layer.filter->compute_query_signature(query_graph, u, neighborhood);
        }
        const std::uint64_t elapsed =
            micros_since(start, std::chrono::high_resolution_clock::now());
        stats->layers[i].spectral_us = elapsed;
        stats->spectral_us += elapsed;
    }

    stats->base_candidates = layers_.front().filter->count_label_degree_candidates(query_graph);

    bool emptied = false;
    for (std::size_t i = 0; i < layers_.size(); ++i) {
        Layer& layer = layers_[i];
        CumulativeLayerStats& layer_stats = stats->layers[i];
        if (emptied) {
            layer_stats.entry = 0;
            layer_stats.exit = 0;
            layer_stats.ran = false;
            continue;
        }
        layer_stats.ran = true;
        layer_stats.entry = i == 0 ? stats->base_candidates : count_candidates(accumulate_);

        const auto start = std::chrono::high_resolution_clock::now();
        layer.filter->set_label_filter(label_filter_);
        layer.filter->set_spectral_filter(spectral_filter_);
        layer.filter->set_channel_stats(channel_stats_);
        layer.filter->filter_into(query_graph, layer.query_signatures, layer_out_);

        if (i == 0) {
            accumulate_.swap(layer_out_);
        } else {
            intersect_into(accumulate_, layer_out_, scratch_);
            layer_stats.released = count_not_contained(scratch_, accumulate_);
            accumulate_.swap(scratch_);
        }

        layer_stats.exit = count_candidates(accumulate_);
        layer_stats.missing_signature_keeps =
            count_missing_signature_keeps(accumulate_, *layer.index);
        layer_stats.filter_us = micros_since(start, std::chrono::high_resolution_clock::now());
        stats->filter_us += layer_stats.filter_us;
        stats->released += layer_stats.released;
        stats->missing_signature_keeps += layer_stats.missing_signature_keeps;
        stats->coarse_dominance_checks += layer.filter->coarse_dominance_checks();
        stats->full_dominance_checks += layer.filter->full_dominance_checks();
        stats->ld_candidates += layer.filter->ld_candidates();
        stats->label_rejects += layer.filter->label_rejects();
        stats->dim_rejects += layer.filter->dim_rejects();
        stats->key_rejects += layer.filter->key_rejects();
        stats->win_rejects += layer.filter->win_rejects();
        stats->empty_accepts += layer.filter->empty_accepts();

        if (!layer_dump_dir_.empty()) {
            write_candidate_dump(layer_dump_dir_ + "/" + layer_dump_stem_ + ".r" +
                                     std::to_string(layer_stats.radius) + ".cand",
                                 accumulate_);
        }

        if (layer_stats.exit == 0) {
            emptied = true;
        }
    }

    stats->candidates_after_filter = count_candidates(accumulate_);
    out.swap(accumulate_);
}
}
