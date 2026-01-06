#ifndef SLIM_OFFLINE_ENCODER_H
#define SLIM_OFFLINE_ENCODER_H

#include <slim/index.h>
#include <slim/graph.h>
#include <vector>
#include <unordered_set>
#include <algorithm>

namespace slim {
namespace offline {

class LabelEncoderBuilder {
public:
    static LabelEncoder fromGraph(const Graph& graph);

private:
    static std::vector<Label> collectUniqueLabels(const Graph& graph);
};

}
}

#endif
