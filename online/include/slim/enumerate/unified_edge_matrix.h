#ifndef SLIM_UNIFIED_EDGE_MATRIX_H
#define SLIM_UNIFIED_EDGE_MATRIX_H

#include <slim/graph.h>
#include <slim/filter.h>
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>

namespace slim
{

    struct UnifiedEdge
    {
        uint32_t *offset_ = nullptr;
        uint32_t *edge_ = nullptr;
        uint32_t vertex_count_ = 0;
        uint32_t edge_count_ = 0;
        uint32_t max_degree_ = 0;

        UnifiedEdge() = default;
        ~UnifiedEdge()
        {
            delete[] offset_;
            delete[] edge_;
        }

        UnifiedEdge(const UnifiedEdge &) = delete;
        UnifiedEdge &operator=(const UnifiedEdge &) = delete;

        UnifiedEdge(UnifiedEdge &&other) noexcept
            : offset_(other.offset_), edge_(other.edge_),
              vertex_count_(other.vertex_count_), edge_count_(other.edge_count_),
              max_degree_(other.max_degree_)
        {
            other.offset_ = nullptr;
            other.edge_ = nullptr;
        }

        UnifiedEdge &operator=(UnifiedEdge &&other) noexcept
        {
            if (this != &other)
            {
                delete[] offset_;
                delete[] edge_;
                offset_ = other.offset_;
                edge_ = other.edge_;
                vertex_count_ = other.vertex_count_;
                edge_count_ = other.edge_count_;
                max_degree_ = other.max_degree_;
                other.offset_ = nullptr;
                other.edge_ = nullptr;
            }
            return *this;
        }
    };

    class UnifiedEdgeMatrix
    {
    public:
        void build(const Graph &data_graph,
                   const Graph &query_graph,
                   const CandidateSet &candidates);

        inline const UnifiedEdge *getEdge(VertexId u, VertexId v) const
        {
            if (u >= edge_matrix_.size() || v >= edge_matrix_[u].size())
            {
                return nullptr;
            }
            return edge_matrix_[u][v].get();
        }

        inline const uint32_t *getCandidatesByIdx(VertexId u_bn, VertexId u,
                                                  uint32_t bn_idx, uint32_t &count) const
        {
            const auto *edge = getEdge(u_bn, u);
            if (edge == nullptr || edge->offset_ == nullptr)
            {
                count = 0;
                return nullptr;
            }
            count = edge->offset_[bn_idx + 1] - edge->offset_[bn_idx];
            return edge->edge_ + edge->offset_[bn_idx];
        }

        const std::vector<VertexId> &getCandidates(VertexId u) const
        {
            return candidates_[u];
        }

        uint32_t getCandidatesCount(VertexId u) const
        {
            return static_cast<uint32_t>(candidates_[u].size());
        }

        const std::vector<std::vector<VertexId>> &getAllCandidates() const
        {
            return candidates_;
        }

    private:
        std::vector<std::vector<std::unique_ptr<UnifiedEdge>>> edge_matrix_;
        std::vector<std::vector<VertexId>> candidates_;
    };

    size_t unified_enumerate(
        const Graph &query_graph,
        const Graph &data_graph,
        const UnifiedEdgeMatrix &edge_matrix,
        const VertexId *order,
        size_t order_size,
        size_t limit,
        bool *early_stopped);

}

#endif
