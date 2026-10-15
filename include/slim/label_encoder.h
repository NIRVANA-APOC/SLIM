#ifndef SLIM_LABEL_ENCODER_H
#define SLIM_LABEL_ENCODER_H

#include <slim/labeled_graph.h>
#include <slim/types.h>

#include <cstddef>
#include <vector>

namespace slim {
class LabelEncoder {
public:
    static LabelEncoder from_graph(const LabeledGraph& graph);

    double encode(Label label) const;

    std::size_t num_labels() const { return sorted_labels_.size(); }

private:
    std::vector<Label> sorted_labels_;
};
}

#endif
