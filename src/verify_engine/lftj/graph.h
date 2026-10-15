#ifndef SLIM_VERIFY_ENGINE_LFTJ_GRAPH_H
#define SLIM_VERIFY_ENGINE_LFTJ_GRAPH_H

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace slim_lftj {
using ui = unsigned int;
using VertexID = std::uint32_t;
using LabelID = ui;

class Graph {
public:
    Graph() = default;

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;
    }

    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;

    void build_from_csr(ui vertices_count, ui edges_count, ui labels_count, ui max_degree,
                        const std::vector<ui>& offsets,
                        const std::vector<VertexID>& neighbors,
                        const std::vector<LabelID>& labels) {
        vertices_count_ = vertices_count;
        edges_count_ = edges_count;
        labels_count_ = labels_count;
        max_degree_ = max_degree;

        offsets_ = new ui[offsets.size()];
        for (std::size_t i = 0; i < offsets.size(); ++i) {
            offsets_[i] = offsets[i];
        }
        neighbors_ = new VertexID[neighbors.size()];
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
            neighbors_[i] = neighbors[i];
        }
        labels_ = new LabelID[labels.size()];
        for (std::size_t i = 0; i < labels.size(); ++i) {
            labels_[i] = labels[i];
        }
    }

    ui getVerticesCount() const { return vertices_count_; }

    ui getEdgesCount() const { return edges_count_; }

    ui getLabelsCount() const { return labels_count_; }

    ui getGraphMaxDegree() const { return max_degree_; }

    ui getVertexDegree(const VertexID id) const { return offsets_[id + 1] - offsets_[id]; }

    LabelID getVertexLabel(const VertexID id) const { return labels_[id]; }

    const VertexID* getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors = getVertexNeighbors(v, count);

        int begin = 0;
        int end = static_cast<int>(count) - 1;
        while (begin <= end) {
            const int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            } else if (neighbors[mid] > u) {
                end = mid - 1;
            } else {
                begin = mid + 1;
            }
        }
        return false;
    }

private:
    ui vertices_count_ = 0;
    ui edges_count_ = 0;
    ui labels_count_ = 0;
    ui max_degree_ = 0;

    ui* offsets_ = nullptr;
    VertexID* neighbors_ = nullptr;
    LabelID* labels_ = nullptr;
};

class Edges {
public:
    Edges() = default;

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }

    Edges(const Edges&) = delete;
    Edges& operator=(const Edges&) = delete;

    ui* offset_ = nullptr;
    ui* edge_ = nullptr;
    ui vertex_count_ = 0;
    ui edge_count_ = 0;
    ui max_degree_ = 0;
};
}

#endif
