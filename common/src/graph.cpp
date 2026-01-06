

#include <slim/graph.h>
#include <algorithm>
#include <unordered_set>
#include <stdexcept>

namespace slim {

Graph::Graph(size_t num_vertices,
             const std::vector<std::pair<VertexId, VertexId>>& edges,
             const std::vector<Label>& labels)
    : num_vertices_(num_vertices) {

    if (labels.size() != num_vertices) {
        throw std::runtime_error("Label count mismatch: expected " +
                                std::to_string(num_vertices) +
                                ", got " + std::to_string(labels.size()));
    }

    offsets_.resize(num_vertices + 1, 0);

    std::vector<std::vector<VertexId>> temp_edges(num_vertices);

    for (const auto& edge : edges) {
        VertexId u = edge.first;
        VertexId v = edge.second;

        if (u == v) {
            throw std::runtime_error("Self-loop detected: edge (" +
                                    std::to_string(u) + ", " +
                                    std::to_string(v) + ") violates GT2 constraint");
        }

        if (u >= num_vertices || v >= num_vertices) {
            throw std::runtime_error("Vertex ID out of range: u=" +
                                    std::to_string(u) + ", v=" +
                                    std::to_string(v));
        }

        temp_edges[u].push_back(v);
        temp_edges[v].push_back(u);
    }

    for (size_t i = 0; i < num_vertices; ++i) {
        auto& neighbors = temp_edges[i];
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }

    for (size_t i = 0; i < num_vertices; ++i) {
        offsets_[i + 1] = offsets_[i] + temp_edges[i].size();
    }

    edges_.resize(offsets_[num_vertices]);
    size_t pos = 0;

    for (size_t i = 0; i < num_vertices; ++i) {
        for (VertexId v : temp_edges[i]) {
            edges_[pos++] = v;
        }
    }

    labels_ = labels;

    num_edges_ = edges_.size() / 2;
}

}
