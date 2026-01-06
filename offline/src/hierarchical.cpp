#include <slim/offline/hierarchical.h>
#include <algorithm>
#include <map>

namespace slim {
namespace offline {

std::unique_ptr<HierarchicalIndex> HierarchicalIndexBuilder::build(
    const std::vector<SpectralSignature>& signatures,
    const std::vector<Label>& labels,
    const std::unordered_map<Label, std::vector<VertexId>>& label_index) {

    if (signatures.empty() || labels.empty()) {
        return nullptr;
    }

    double tau = computeTau(signatures);
    if (tau == 0.0) {
        return nullptr;
    }

    std::vector<uint32_t> layers = computeLayers(signatures, tau);

    std::unordered_map<HierarchicalKey, std::vector<VertexId>, HierarchicalKeyHash> hierarchical_index;
    std::unordered_map<Label, uint32_t> label_max_layer;
    std::map<std::pair<Label, uint32_t>, uint32_t> label_layer_max_dim;

    uint32_t global_max_layer = 0;
    uint32_t global_max_dimension = 0;

    for (size_t v = 0; v < signatures.size(); ++v) {
        Label label = labels[v];
        uint32_t layer = layers[v];
        uint32_t dimension = static_cast<uint32_t>(signatures[v].dimension());

        HierarchicalKey key{label, layer, dimension};
        hierarchical_index[key].push_back(static_cast<VertexId>(v));

        if (label_max_layer.find(label) == label_max_layer.end()) {
            label_max_layer[label] = layer;
        } else {
            label_max_layer[label] = std::max(label_max_layer[label], layer);
        }

        auto key_pair = std::make_pair(label, layer);
        if (label_layer_max_dim.find(key_pair) == label_layer_max_dim.end()) {
            label_layer_max_dim[key_pair] = dimension;
        } else {
            label_layer_max_dim[key_pair] = std::max(label_layer_max_dim[key_pair], dimension);
        }

        global_max_layer = std::max(global_max_layer, layer);
        global_max_dimension = std::max(global_max_dimension, dimension);
    }

    return std::unique_ptr<HierarchicalIndex>(new HierarchicalIndex(
        tau,
        global_max_layer,
        global_max_dimension,
        std::move(layers),
        std::move(hierarchical_index),
        std::move(label_max_layer),
        std::move(label_layer_max_dim)
    ));
}

double HierarchicalIndexBuilder::computeTau(const std::vector<SpectralSignature>& signatures) {
    if (signatures.empty()) {
        return 0.0;
    }

    double lambda_max_global = 0.0;
    for (const auto& sig : signatures) {
        lambda_max_global = std::max(lambda_max_global, sig.lambdaMax());
    }

    double tau = lambda_max_global / std::sqrt(static_cast<double>(signatures.size()));
    return tau;
}

std::vector<uint32_t> HierarchicalIndexBuilder::computeLayers(
    const std::vector<SpectralSignature>& signatures, double tau) {

    std::vector<uint32_t> layers(signatures.size());

    for (size_t v = 0; v < signatures.size(); ++v) {
        double lambda_max = signatures[v].lambdaMax();
        layers[v] = static_cast<uint32_t>(std::floor(lambda_max / tau));
    }

    return layers;
}

}
}
