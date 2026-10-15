#ifndef SLIM_LABELED_GRAPH_H
#define SLIM_LABELED_GRAPH_H

#include <slim/types.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace slim {
class LabeledGraph {
public:
    LabeledGraph() = default;

    LabeledGraph(std::size_t num_vertices,
                 const std::vector<std::pair<VertexId, VertexId>>& edges,
                 const std::vector<Label>& labels);

    std::size_t num_vertices() const { return num_vertices_; }
    std::size_t num_edges() const { return edges_.size() / 2; }

    const VertexId* neighbors(VertexId v, std::size_t& count) const {
        if (v >= num_vertices_) {
            throw std::runtime_error("vertex id out of range: " + std::to_string(v));
        }
        const std::size_t begin = offsets_[v];
        const std::size_t end = offsets_[v + 1];
        count = end - begin;
        return edges_.data() + begin;
    }

    std::size_t degree(VertexId v) const {
        if (v >= num_vertices_) {
            throw std::runtime_error("vertex id out of range: " + std::to_string(v));
        }
        return offsets_[v + 1] - offsets_[v];
    }

    Label label(VertexId v) const {
        if (v >= num_vertices_) {
            throw std::runtime_error("vertex id out of range: " + std::to_string(v));
        }
        return labels_[v];
    }

private:
    std::size_t num_vertices_ = 0;
    std::vector<std::uint32_t> offsets_;
    std::vector<VertexId> edges_;
    std::vector<Label> labels_;
};
}

#endif
