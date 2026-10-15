#include <slim/labeled_graph.h>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace slim {
LabeledGraph::LabeledGraph(
    std::size_t num_vertices,
    const std::vector<std::pair<VertexId, VertexId>>& edges,
    const std::vector<Label>& labels)
    : num_vertices_(num_vertices), labels_(labels) {
    if (labels_.size() != num_vertices_) {
        throw std::runtime_error("label count does not match vertex count");
    }

    std::vector<std::uint32_t> end(num_vertices_ + 1, 0);

    for (const auto& [u, v] : edges) {
        if (u == v) {
            throw std::runtime_error("self-loop is not allowed: " + std::to_string(u));
        }
        if (u >= num_vertices_ || v >= num_vertices_) {
            throw std::runtime_error("edge endpoint is out of range");
        }
        ++end[u + 1];
        ++end[v + 1];
    }

    for (std::size_t i = 0; i < num_vertices_; ++i) {
        end[i + 1] += end[i];
    }
    offsets_ = end;

    std::vector<VertexId> adjacency(edges.size() * 2);
    std::vector<std::uint32_t> cursor = offsets_;
    for (const auto& [u, v] : edges) {
        adjacency[cursor[u]++] = v;
        adjacency[cursor[v]++] = u;
    }

    std::vector<VertexId> compact;
    compact.reserve(adjacency.size());
    for (VertexId u = 0; u < num_vertices_; ++u) {
        const std::size_t begin = offsets_[u];
        const std::size_t finish = end[u + 1];
        std::sort(adjacency.begin() + begin, adjacency.begin() + finish);
        const auto last = std::unique(adjacency.begin() + begin,
                                      adjacency.begin() + finish);
        offsets_[u] = static_cast<std::uint32_t>(compact.size());
        compact.insert(compact.end(), adjacency.begin() + begin, last);
    }
    offsets_[num_vertices_] = static_cast<std::uint32_t>(compact.size());
    edges_ = std::move(compact);
}
}
