#include <slim/k_hop_neighborhood.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace slim {
KHopNeighborhood::KHopNeighborhood(const LabeledGraph& graph)
    : graph_(graph) {
    mark_.assign(graph.num_vertices(), 0);
}

const std::vector<VertexId>& KHopNeighborhood::extract(VertexId v, int k) {
    if (++generation_ >> 32 != 0) {
        std::fill(mark_.begin(), mark_.end(), 0);
        generation_ = 1;
    }
    const std::uint64_t tag = generation_ << 32;

    result_.clear();
    result_.push_back(v);
    mark_[v] = tag | 1;

    if (k <= 1) {
        std::size_t count = 0;
        const VertexId* neighbors = graph_.neighbors(v, count);
        for (std::size_t i = 0; i < count; ++i) {
            const VertexId w = neighbors[i];
            mark_[w] = tag | static_cast<std::uint64_t>(result_.size() + 1);
            result_.push_back(w);
        }
        return result_;
    }

    frontier_.clear();
    frontier_.push_back(v);
    for (int hop = 0; hop < k; ++hop) {
        next_.clear();
        for (VertexId u : frontier_) {
            std::size_t count = 0;
            const VertexId* neighbors = graph_.neighbors(u, count);
            for (std::size_t i = 0; i < count; ++i) {
                const VertexId w = neighbors[i];
                if ((mark_[w] >> 32) != generation_) {
                    mark_[w] = tag | static_cast<std::uint64_t>(result_.size() + 1);
                    result_.push_back(w);
                    next_.push_back(w);
                }
            }
        }
        if (next_.empty()) {
            break;
        }
        frontier_.swap(next_);
    }
    return result_;
}
}
