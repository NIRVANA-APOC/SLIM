

#ifndef SLIM_GRAPH_H
#define SLIM_GRAPH_H

#include <cstdint>
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <algorithm>

namespace slim {

using VertexId = uint32_t;

using Label = uint32_t;

class Graph {
public:

    Graph(size_t num_vertices,
          const std::vector<std::pair<VertexId, VertexId>>& edges,
          const std::vector<Label>& labels);

    size_t numVertices() const { return num_vertices_; }

    size_t numEdges() const { return num_edges_; }

    const VertexId* neighbors(VertexId v, size_t& count) const {
        if (v >= num_vertices_) {
            throw std::runtime_error("Vertex ID out of range: " + std::to_string(v));
        }
        size_t start = offsets_[v];
        size_t end = offsets_[v + 1];
        count = end - start;
        return edges_.data() + start;
    }

    size_t degree(VertexId v) const {
        if (v >= num_vertices_) {
            throw std::runtime_error("Vertex ID out of range: " + std::to_string(v));
        }
        return offsets_[v + 1] - offsets_[v];
    }

    Label label(VertexId v) const {
        if (v >= num_vertices_) {
            throw std::runtime_error("Vertex ID out of range: " + std::to_string(v));
        }
        return labels_[v];
    }

    inline bool hasEdge(VertexId u, VertexId v) const {
        if (u >= num_vertices_ || v >= num_vertices_) {
            return false;
        }
        const size_t start = offsets_[u];
        const size_t end = offsets_[u + 1];

        const VertexId* data = edges_.data();
        return std::binary_search(data + start, data + end, v);
    }

private:
    size_t num_vertices_;
    size_t num_edges_;
    std::vector<size_t> offsets_;
    std::vector<VertexId> edges_;
    std::vector<Label> labels_;
};

}

#endif
