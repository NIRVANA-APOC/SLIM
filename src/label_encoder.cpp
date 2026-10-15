#include <slim/label_encoder.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace slim {
LabelEncoder LabelEncoder::from_graph(const LabeledGraph& graph) {
    std::vector<Label> unique;
    unique.reserve(graph.num_vertices());
    for (VertexId v = 0; v < graph.num_vertices(); ++v) {
        unique.push_back(graph.label(v));
    }
    std::sort(unique.begin(), unique.end());
    unique.erase(std::unique(unique.begin(), unique.end()), unique.end());

    LabelEncoder encoder;
    encoder.sorted_labels_ = std::move(unique);
    return encoder;
}

double LabelEncoder::encode(Label label) const {
    const auto it = std::lower_bound(sorted_labels_.begin(), sorted_labels_.end(), label);
    if (it == sorted_labels_.end() || *it != label) {
        throw std::runtime_error("label is not present in the label alphabet: " +
                                 std::to_string(label));
    }
    const std::size_t index = static_cast<std::size_t>(it - sorted_labels_.begin());
    return static_cast<double>(index) / static_cast<double>(sorted_labels_.size());
}
}
