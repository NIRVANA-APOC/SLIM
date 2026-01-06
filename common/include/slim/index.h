

#ifndef SLIM_INDEX_H
#define SLIM_INDEX_H

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <string>
#include <fstream>

namespace slim {

class QueryTree;
class Graph;

using VertexId = uint32_t;

using Label = uint32_t;

constexpr VertexId INVALID_VERTEX = static_cast<VertexId>(-1);

constexpr uint32_t INDEX_VERSION = 1;

constexpr double EPSILON = 1e-6;

class QueryTree {
public:

    VertexId root() const { return root_; }

    VertexId order(int depth) const {
        if (depth < 0 || static_cast<size_t>(depth) >= order_.size()) {
            return INVALID_VERTEX;
        }
        return order_[depth];
    }

    VertexId parent(VertexId u) const {
        auto it = parent_map_.find(u);
        return (it == parent_map_.end()) ? INVALID_VERTEX : it->second;
    }

    size_t numVertices() const { return order_.size(); }

    static QueryTree build(const Graph& query_graph);

    static QueryTree buildWithCandidates(const Graph& query_graph,
                                         const std::vector<size_t>& candidates_count);

private:
    VertexId root_;
    std::vector<VertexId> order_;
    std::unordered_map<VertexId, VertexId> parent_map_;
};

struct LabelEncoder {
    std::unordered_map<Label, double> encoding;
    size_t num_labels;

    double encode(Label label) const {
        auto it = encoding.find(label);
        if (it == encoding.end()) {
            throw std::runtime_error("Label not found in encoder: " + std::to_string(label));
        }
        return it->second;
    }

    bool contains(Label label) const {
        return encoding.find(label) != encoding.end();
    }
};

struct SpectralSignature {
    std::vector<double> eigenvalues;

    size_t dimension() const { return eigenvalues.size(); }

    double lambdaMax() const {
        return eigenvalues.empty() ? 0.0 : eigenvalues.back();
    }

    bool dominatedBy(const SpectralSignature& other) const {
        size_t m = eigenvalues.size();
        size_t n = other.eigenvalues.size();

        if (m > n) {
            return false;
        }

        size_t delta = n - m;

        for (size_t k = 0; k < m; ++k) {
            if (eigenvalues[k] > other.eigenvalues[k + delta] + EPSILON) {
                return false;
            }
        }

        return true;
    }
};

struct HierarchicalKey {
    Label label;
    uint32_t layer;
    uint32_t dimension;

    bool operator==(const HierarchicalKey& other) const {
        return label == other.label &&
               layer == other.layer &&
               dimension == other.dimension;
    }
};

struct HierarchicalKeyHash {
    size_t operator()(const HierarchicalKey& key) const {
        return std::hash<Label>()(key.label) ^
               (std::hash<uint32_t>()(key.layer) << 1) ^
               (std::hash<uint32_t>()(key.dimension) << 2);
    }
};

class HierarchicalIndex {
public:

    static std::unique_ptr<HierarchicalIndex> deserialize(std::ifstream& file);

    void serialize(std::ofstream& file) const;

    double tau() const { return tau_; }

    uint32_t maxLayer() const { return max_layer_; }

    uint32_t maxDimension() const { return max_dimension_; }

    uint32_t layer(VertexId v) const {
        if (v >= layers_.size()) {
            return 0;
        }
        return layers_[v];
    }

    const std::vector<VertexId>& getByKey(
        Label label,
        uint32_t layer,
        uint32_t dimension
    ) const {
        HierarchicalKey key{label, layer, dimension};
        auto it = hierarchical_index_.find(key);
        if (it != hierarchical_index_.end()) {
            return it->second;
        }
        static const std::vector<VertexId> empty;
        return empty;
    }

    uint32_t maxLayerForLabel(Label label) const {
        auto it = label_max_layer_.find(label);
        if (it != label_max_layer_.end()) {
            return it->second;
        }
        return 0;
    }

    uint32_t maxDimensionForLabelLayer(Label label, uint32_t layer) const {
        auto it = label_layer_max_dim_.find({label, layer});
        if (it != label_layer_max_dim_.end()) {
            return it->second;
        }
        return 0;
    }

    friend std::unique_ptr<HierarchicalIndex> std::make_unique<HierarchicalIndex>();

    HierarchicalIndex(double tau,
                     uint32_t max_layer,
                     uint32_t max_dimension,
                     std::vector<uint32_t> layers,
                     std::unordered_map<HierarchicalKey, std::vector<VertexId>, HierarchicalKeyHash> hierarchical_index,
                     std::unordered_map<Label, uint32_t> label_max_layer,
                     std::map<std::pair<Label, uint32_t>, uint32_t> label_layer_max_dim)
        : tau_(tau), max_layer_(max_layer), max_dimension_(max_dimension),
          layers_(std::move(layers)),
          hierarchical_index_(std::move(hierarchical_index)),
          label_max_layer_(std::move(label_max_layer)),
          label_layer_max_dim_(std::move(label_layer_max_dim)) {}

private:
    HierarchicalIndex() = default;

    double tau_;
    uint32_t max_layer_;
    uint32_t max_dimension_;
    std::vector<uint32_t> layers_;
    std::unordered_map<HierarchicalKey, std::vector<VertexId>, HierarchicalKeyHash> hierarchical_index_;
    std::unordered_map<Label, uint32_t> label_max_layer_;
    std::map<std::pair<Label, uint32_t>, uint32_t> label_layer_max_dim_;
};

class Index {
public:

    static std::unique_ptr<Index> load(const std::string& path);

    void save(const std::string& path) const;

    static std::unique_ptr<Index> create(
        const LabelEncoder& label_encoder,
        double alpha,
        std::vector<SpectralSignature> signatures,
        std::vector<Label> labels,
        std::unordered_map<Label, std::vector<VertexId>> label_index,
        std::unique_ptr<HierarchicalIndex> hierarchical_index = nullptr);

    uint32_t version() const { return version_; }

    const LabelEncoder& labelEncoder() const { return label_encoder_; }

    double alpha() const { return alpha_; }

    size_t numVertices() const { return signatures_.size(); }

    const SpectralSignature& signature(VertexId v) const {
        if (v >= signatures_.size()) {
            throw std::runtime_error("Vertex ID out of range: " + std::to_string(v));
        }
        return signatures_[v];
    }

    Label label(VertexId v) const {
        if (v >= labels_.size()) {
            throw std::runtime_error("Vertex ID out of range: " + std::to_string(v));
        }
        return labels_[v];
    }

    const std::vector<VertexId>& getByLabel(Label label) const {
        auto it = label_index_.find(label);
        if (it == label_index_.end()) {
            static const std::vector<VertexId> empty;
            return empty;
        }
        return it->second;
    }

    bool useHierarchical() const { return use_hierarchical_; }

    const HierarchicalIndex* hierarchicalIndex() const {
        return hierarchical_index_.get();
    }

private:

    Index() = default;

    friend struct IndexLoader;

    uint32_t version_;
    LabelEncoder label_encoder_;
    double alpha_;
    std::vector<SpectralSignature> signatures_;
    std::vector<Label> labels_;
    std::unordered_map<Label, std::vector<VertexId>> label_index_;
    std::unique_ptr<HierarchicalIndex> hierarchical_index_;
    bool use_hierarchical_;
};

}

#endif
