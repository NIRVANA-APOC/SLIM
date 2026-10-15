#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <vector>

typedef unsigned int ui;

typedef uint32_t VertexID;
typedef ui LabelID;

class Edges {
public:
    ui* offset_;
    ui* edge_;
    ui vertex_count_;
    ui edge_count_;
    ui max_degree_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};

class Graph {
private:
    ui vertices_count_;
    ui edges_count_;
    ui labels_count_;
    ui max_degree_;
    ui max_label_frequency_;

    ui* offsets_;
    VertexID * neighbors_;
    LabelID* labels_;

    std::unordered_map<LabelID, ui> labels_frequency_;

public:
    Graph() {
        vertices_count_ = 0;
        edges_count_ = 0;
        labels_count_ = 0;
        max_degree_ = 0;
        max_label_frequency_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        labels_ = NULL;
        labels_frequency_.clear();
    }

    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;
    }

    void buildFromArrays(ui vertices_count, ui edges_count, const ui* offsets,
                         const VertexID* neighbors, const LabelID* labels) {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;

        vertices_count_ = vertices_count;
        edges_count_ = edges_count;

        offsets_ = new ui[vertices_count_ + 1];
        std::memcpy(offsets_, offsets, sizeof(ui) * (vertices_count_ + 1));

        neighbors_ = new VertexID[static_cast<size_t>(edges_count_) * 2];
        std::memcpy(neighbors_, neighbors, sizeof(VertexID) * static_cast<size_t>(edges_count_) * 2);

        labels_ = new LabelID[vertices_count_];
        std::memcpy(labels_, labels, sizeof(LabelID) * vertices_count_);

        max_degree_ = 0;
        LabelID max_label_id = 0;
        labels_frequency_.clear();
        for (ui v = 0; v < vertices_count_; ++v) {
            const ui degree = offsets_[v + 1] - offsets_[v];
            if (degree > max_degree_) max_degree_ = degree;
            labels_frequency_[labels_[v]] += 1;
            if (labels_[v] > max_label_id) max_label_id = labels_[v];
        }
        labels_count_ = static_cast<ui>(labels_frequency_.size()) > (max_label_id + 1)
                            ? static_cast<ui>(labels_frequency_.size())
                            : max_label_id + 1;
        max_label_frequency_ = 0;
        for (auto element : labels_frequency_) {
            if (element.second > max_label_frequency_) max_label_frequency_ = element.second;
        }
    }
    const ui getLabelsCount() const {
        return labels_count_;
    }

    const ui getVerticesCount() const {
        return vertices_count_;
    }

    const ui getEdgesCount() const {
        return edges_count_;
    }

    const ui getGraphMaxDegree() const {
        return max_degree_;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency_;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets_[id + 1] - offsets_[id];
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency_.find(label) == labels_frequency_.end() ? 0 : labels_frequency_.at(label);
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }
};
