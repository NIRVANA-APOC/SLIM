#include <slim/offline/encoder.h>

namespace slim {
namespace offline {

LabelEncoder LabelEncoderBuilder::fromGraph(const Graph& graph) {
    std::vector<Label> unique_labels = collectUniqueLabels(graph);

    std::sort(unique_labels.begin(), unique_labels.end());

    LabelEncoder encoder;
    encoder.num_labels = unique_labels.size();

    if (unique_labels.size() == 1) {

        encoder.encoding[unique_labels[0]] = 0.0;
    } else {
        for (size_t i = 0; i < unique_labels.size(); ++i) {
            encoder.encoding[unique_labels[i]] = static_cast<double>(i) / unique_labels.size();
        }
    }

    return encoder;
}

std::vector<Label> LabelEncoderBuilder::collectUniqueLabels(const Graph& graph) {
    std::unordered_set<Label> label_set;

    for (size_t v = 0; v < graph.numVertices(); ++v) {
        label_set.insert(graph.label(v));
    }

    return std::vector<Label>(label_set.begin(), label_set.end());
}

}
}
