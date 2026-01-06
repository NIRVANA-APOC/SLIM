#ifndef SLIM_OFFLINE_HIERARCHICAL_H
#define SLIM_OFFLINE_HIERARCHICAL_H

#include <slim/index.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>

namespace slim {
namespace offline {

class HierarchicalIndexBuilder {
public:
    static std::unique_ptr<HierarchicalIndex> build(
        const std::vector<SpectralSignature>& signatures,
        const std::vector<Label>& labels,
        const std::unordered_map<Label, std::vector<VertexId>>& label_index);

private:
    static double computeTau(const std::vector<SpectralSignature>& signatures);
    static std::vector<uint32_t> computeLayers(
        const std::vector<SpectralSignature>& signatures, double tau);
};

}
}

#endif
