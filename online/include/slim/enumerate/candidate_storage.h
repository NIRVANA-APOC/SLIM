#ifndef SLIM_CANDIDATE_STORAGE_H
#define SLIM_CANDIDATE_STORAGE_H

#include <slim/index.h>
#include <slim/graph.h>
#include <slim/filter.h>
#include <vector>
#include <cstdint>

namespace slim
{

    struct EdgeTable
    {
        uint32_t *offset_ = nullptr;
        uint32_t *edge_ = nullptr;
        uint32_t vertex_count_ = 0;
        uint32_t edge_count_ = 0;

        EdgeTable() = default;
        ~EdgeTable()
        {
            delete[] offset_;
            delete[] edge_;
        }

        EdgeTable(const EdgeTable &) = delete;
        EdgeTable &operator=(const EdgeTable &) = delete;

        EdgeTable(EdgeTable &&other) noexcept
            : offset_(other.offset_), edge_(other.edge_),
              vertex_count_(other.vertex_count_), edge_count_(other.edge_count_)
        {
            other.offset_ = nullptr;
            other.edge_ = nullptr;
        }

        EdgeTable &operator=(EdgeTable &&other) noexcept
        {
            if (this != &other)
            {
                delete[] offset_;
                delete[] edge_;
                offset_ = other.offset_;
                edge_ = other.edge_;
                vertex_count_ = other.vertex_count_;
                edge_count_ = other.edge_count_;
                other.offset_ = nullptr;
                other.edge_ = nullptr;
            }
            return *this;
        }
    };

    class TreeEdgeCandidates
    {
    public:
        void build(const Graph &data_graph,
                   const QueryTree &tree,
                   const CandidateSet &candidates);

        void buildOptimized(const Graph &data_graph,
                            const Graph &query_graph,
                            const QueryTree &tree,
                            const CandidateSet &candidates);

        inline const uint32_t *getCandidatesByIdx(VertexId u, VertexId u_bn, uint32_t bn_idx,
                                                  uint32_t &count) const
        {
            const auto &edge = edge_matrix_[u_bn][u];
            if (edge.offset_ == nullptr)
            {
                count = 0;
                return nullptr;
            }
            count = edge.offset_[bn_idx + 1] - edge.offset_[bn_idx];
            return edge.edge_ + edge.offset_[bn_idx];
        }

        inline const EdgeTable &getEdgeTable(VertexId u_bn, VertexId u) const
        {
            return edge_matrix_[u_bn][u];
        }

        VertexId getParent(VertexId u) const
        {
            if (u >= parent_.size())
                return INVALID_VERTEX;
            return parent_[u];
        }

        const std::vector<std::vector<VertexId>> &getAllCandidates() const
        {
            return candidates_;
        }

    private:
        std::vector<std::vector<EdgeTable>> edge_matrix_;
        std::vector<VertexId> parent_;
        std::vector<std::vector<VertexId>> candidates_;
    };

}

#endif
