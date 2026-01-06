#ifndef SLIM_OFFLINE_BUILDER_H
#define SLIM_OFFLINE_BUILDER_H

#include <slim/index.h>
#include <slim/graph.h>
#include <slim/spectral/laplacian.h>
#include <slim/spectral/eigen.h>
#include <slim/offline/encoder.h>
#include <slim/offline/hierarchical.h>
#include <memory>
#include <vector>
#include <unordered_map>

namespace slim {
namespace offline {

constexpr double HIERARCHICAL_INDEX_THRESHOLD = 100.0;

class IndexBuilder {
public:
    IndexBuilder(const Graph& graph, double alpha);

    void setNumThreads(int threads);

    std::unique_ptr<Index> build();

private:
    const Graph& graph_;
    double alpha_;
    int num_threads_;

    std::vector<SpectralSignature> computeAllSignatures(
        const LaplacianBuilder& builder);

    std::unordered_map<Label, std::vector<VertexId>> buildLabelIndex() const;

    bool shouldUseHierarchical(
        const std::unordered_map<Label, std::vector<VertexId>>& label_index) const;
};

}
}

#endif
