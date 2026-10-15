#ifndef SLIM_SPECTRAL_DOMINANCE_INDEX_H
#define SLIM_SPECTRAL_DOMINANCE_INDEX_H

#include <slim/amsd_index.h>
#include <slim/config.h>
#include <slim/label_encoder.h>
#include <slim/labeled_graph.h>
#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace slim {
constexpr std::uint32_t INDEX_FORMAT_VERSION = 2;

class SpectralDominanceIndex {
public:
    struct DominanceKey {
        std::uint32_t layer = 0;
        std::uint32_t dimension = 0;
    };

    struct Bucket {
        std::uint32_t layer = 0;
        std::uint32_t dimension = 0;
        std::uint64_t begin = 0;
        std::uint64_t count = 0;
    };

    struct LabelGroup {
        Label label = 0;
        std::vector<Bucket> buckets;
        std::vector<VertexId> vertices;
    };

    struct EmptyGroup {
        Label label = 0;
        std::vector<VertexId> vertices;
    };

    const std::uint8_t* ball_counts(VertexId v) const {
        return has_ball_labels_
                   ? ball_counts_.data() + static_cast<std::size_t>(v) * num_labels_
                   : nullptr;
    }
    bool has_ball_labels() const { return has_ball_labels_; }
    std::uint32_t num_labels() const { return num_labels_; }
    const std::vector<Label>& alphabet() const { return alphabet_; }

    struct Identity {
        std::uint32_t format_version = 0;
        std::uint64_t data_fingerprint = 0;
        std::uint64_t label_fingerprint = 0;
        std::uint64_t neighborhood_size_limit = 0;
        double tau_scale = 0.0;
        std::uint64_t data_vertices = 0;
        std::uint64_t data_edges = 0;
    };

    bool has_identity() const { return has_identity_; }
    const Identity& identity() const { return identity_; }

    double alpha() const { return alpha_; }
    int k() const { return k_; }

    double tau() const { return tau_; }

    std::size_t num_vertices() const { return signatures_.size(); }

    const LocalSpectralSignature& signature(VertexId v) const { return signatures_[v]; }

    const LabelGroup* find_group(Label label) const;
    const EmptyGroup* find_empty_group(Label label) const;

    std::size_t group_position(Label label) const;

    bool has_amsd() const { return amsd_.present; }
    const AmsdIndex& amsd() const { return amsd_; }
    const AmsdGroup* amsd_group_at(std::size_t position) const {
        return position < amsd_.groups.size() ? &amsd_.groups[position] : nullptr;
    }

    void save(const std::string& path) const;
    static std::unique_ptr<SpectralDominanceIndex> load(const std::string& path);

private:
    SpectralDominanceIndex() = default;

    double alpha_ = 1.0;
    int k_ = 1;
    double tau_ = 0.0;
    bool has_identity_ = false;
    Identity identity_;
    std::vector<LocalSpectralSignature> signatures_;
    std::vector<LabelGroup> groups_;
    std::vector<EmptyGroup> empty_groups_;
    AmsdIndex amsd_;

    bool has_ball_labels_ = false;
    std::uint32_t num_labels_ = 0;
    std::vector<Label> alphabet_;
    std::vector<std::uint8_t> ball_counts_;

    friend class OfflineIndexBuilder;
};

class OfflineIndexBuilder {
public:
    struct BuildTimings {
        std::uint64_t preprocess_ms = 0;
        std::uint64_t assemble_ms = 0;
        std::uint64_t amsd_ms = 0;
    };

    explicit OfflineIndexBuilder(const LabeledGraph& graph, SlimConfig config);

    std::unique_ptr<SpectralDominanceIndex> build();

    std::size_t skipped_vertices() const { return skipped_vertices_; }

    std::size_t label_skipped_vertices() const { return label_skipped_vertices_; }

    const BuildTimings& timings() const { return timings_; }

private:
    const LabeledGraph& graph_;
    SlimConfig config_;
    std::size_t skipped_vertices_ = 0;
    std::size_t label_skipped_vertices_ = 0;
    BuildTimings timings_;
};
}

#endif
