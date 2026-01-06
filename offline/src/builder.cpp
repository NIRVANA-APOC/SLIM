#include <slim/offline/builder.h>
#include <iostream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace slim {
namespace offline {

IndexBuilder::IndexBuilder(const Graph& graph, double alpha)
    : graph_(graph), alpha_(alpha), num_threads_(0) {

    if (graph_.numVertices() == 0) {
        throw std::runtime_error("Cannot build index for empty graph");
    }

    if (alpha_ <= 0.0) {
        throw std::runtime_error("Alpha must be positive");
    }
}

void IndexBuilder::setNumThreads(int threads) {
    num_threads_ = threads;
}

std::unique_ptr<Index> IndexBuilder::build() {
    std::cout << "Building index...\n";

    std::cout << "Building label encoder...\n";
    LabelEncoder label_encoder = LabelEncoderBuilder::fromGraph(graph_);

    LaplacianBuilder laplacian_builder(alpha_, label_encoder);

    std::cout << "Computing spectral signatures for " << graph_.numVertices() << " vertices...\n";
    std::vector<SpectralSignature> signatures = computeAllSignatures(laplacian_builder);

    std::cout << "Building label index...\n";
    std::unordered_map<Label, std::vector<VertexId>> label_index = buildLabelIndex();

    std::vector<Label> labels(graph_.numVertices());
    for (size_t v = 0; v < graph_.numVertices(); ++v) {
        labels[v] = graph_.label(v);
    }

    bool use_hierarchical = shouldUseHierarchical(label_index);
    std::unique_ptr<HierarchicalIndex> hierarchical_index = nullptr;

    if (use_hierarchical) {
        std::cout << "Building hierarchical index...\n";
        hierarchical_index = HierarchicalIndexBuilder::build(signatures, labels, label_index);
        if (hierarchical_index) {
            std::cout << "  Hierarchical index built (tau=" << hierarchical_index->tau()
                      << ", max_layer=" << hierarchical_index->maxLayer()
                      << ", max_dimension=" << hierarchical_index->maxDimension() << ")\n";
        } else {
            std::cout << "  Warning: Hierarchical index build failed, proceeding without it\n";
        }
    } else {
        std::cout << "Hierarchical index not needed (avg vertices per label <= "
                  << HIERARCHICAL_INDEX_THRESHOLD << ")\n";
    }

    return Index::create(
        label_encoder,
        alpha_,
        std::move(signatures),
        std::move(labels),
        std::move(label_index),
        std::move(hierarchical_index)
    );
}

std::vector<SpectralSignature> IndexBuilder::computeAllSignatures(
    const LaplacianBuilder& builder) {

    size_t num_vertices = graph_.numVertices();
    std::vector<SpectralSignature> signatures(num_vertices);

#ifdef _OPENMP
    int threads = num_threads_;
    if (threads <= 0) {
        threads = omp_get_max_threads();
    }
    omp_set_num_threads(threads);
#endif

#pragma omp parallel for schedule(dynamic) if(num_threads_ != 1)
    for (size_t v = 0; v < num_vertices; ++v) {
        signatures[v] = computeNeighborhoodSpectrum(
            graph_, static_cast<VertexId>(v), builder);
    }

    std::cout << "    Processed " << num_vertices << "/" << num_vertices << " vertices\n";

    return signatures;
}

std::unordered_map<Label, std::vector<VertexId>> IndexBuilder::buildLabelIndex() const {
    std::unordered_map<Label, std::vector<VertexId>> label_index;

    for (size_t v = 0; v < graph_.numVertices(); ++v) {
        Label label = graph_.label(v);
        label_index[label].push_back(static_cast<VertexId>(v));
    }

    return label_index;
}

bool IndexBuilder::shouldUseHierarchical(
    const std::unordered_map<Label, std::vector<VertexId>>& label_index) const {

    if (label_index.empty()) {
        return false;
    }

    double avg_vertices_per_label = static_cast<double>(graph_.numVertices()) / label_index.size();

    return avg_vertices_per_label > HIERARCHICAL_INDEX_THRESHOLD;
}

}
}
