#ifndef SLIM_K_HOP_NEIGHBORHOOD_H
#define SLIM_K_HOP_NEIGHBORHOOD_H

#include <slim/labeled_graph.h>
#include <slim/types.h>

#include <cstdint>
#include <vector>

namespace slim {
class KHopNeighborhood {
public:
    explicit KHopNeighborhood(const LabeledGraph& graph);

    const std::vector<VertexId>& extract(VertexId v, int k);

    std::uint32_t local_index(VertexId w) const {
        return w < mark_.size() ? static_cast<std::uint32_t>(mark_[w]) : 0;
    }

    bool contains(VertexId w) const {
        return w < mark_.size() && (mark_[w] >> 32) == generation_;
    }

private:
    const LabeledGraph& graph_;
    std::vector<std::uint64_t> mark_;
    std::uint64_t generation_ = 0;
    std::vector<VertexId> result_;
    std::vector<VertexId> frontier_;
    std::vector<VertexId> next_;
};
}

#endif
