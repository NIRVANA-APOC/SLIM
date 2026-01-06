#include <slim/enumerate/unified_edge_matrix.h>
#include <slim/enumerate/set_intersection.h>
#include <slim/filter.h>
#include <cstring>
#include <algorithm>

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

namespace slim
{

    alignas(64) static thread_local uint32_t g_buffer_pool[64][65536];
    alignas(64) static thread_local uint32_t g_temp_pool[65536];
    static thread_local uint32_t *g_valid_candidates[64];
    static thread_local uint32_t *g_temp_buffer = nullptr;
    static thread_local bool g_init = false;

    static inline void initBuffers()
    {
        if (UNLIKELY(!g_init))
        {
            for (int i = 0; i < 64; ++i)
            {
                g_valid_candidates[i] = g_buffer_pool[i];
            }
            g_temp_buffer = g_temp_pool;
            g_init = true;
        }
    }

    class UnifiedEnumerator
    {
    public:
        UnifiedEnumerator(const Graph &query_graph, const Graph &data_graph,
                          const UnifiedEdgeMatrix &edge_matrix)
            : query_graph_(query_graph), data_graph_(data_graph),
              edge_matrix_(edge_matrix), count_(0), early_stopped_(false)
        {
            data_size_ = data_graph.numVertices();
            visited_ = new bool[data_size_];
            embedding_ = new VertexId[64];
        }

        ~UnifiedEnumerator()
        {
            delete[] visited_;
            delete[] embedding_;
        }

        size_t enumerate(const VertexId *order, size_t order_size, size_t limit);

        bool wasEarlyStopped() const { return early_stopped_; }
        size_t count() const { return count_; }

    private:
        const Graph &query_graph_;
        const Graph &data_graph_;
        const UnifiedEdgeMatrix &edge_matrix_;

        bool *visited_;
        VertexId *embedding_;
        size_t data_size_;
        size_t count_;
        bool early_stopped_;
    };

    size_t UnifiedEnumerator::enumerate(const VertexId *order, size_t order_size, size_t limit)
    {
        count_ = 0;
        early_stopped_ = false;

        const int max_depth = static_cast<int>(order_size);
        if (UNLIKELY(max_depth == 0))
            return 0;

        initBuffers();
        std::memset(visited_, false, data_size_);

        VertexId bn[64][64];
        uint32_t bn_count[64];
        std::memset(bn_count, 0, sizeof(bn_count));

        const size_t query_size = query_graph_.numVertices();
        std::vector<bool> visited_query(query_size, false);
        visited_query[order[0]] = true;

        for (int depth = 1; depth < max_depth; ++depth)
        {
            VertexId u = order[depth];
            size_t u_nbrs_count;
            const VertexId *u_nbrs = query_graph_.neighbors(u, u_nbrs_count);

            for (size_t j = 0; j < u_nbrs_count; ++j)
            {
                VertexId nbr = u_nbrs[j];
                if (visited_query[nbr])
                {
                    bn[depth][bn_count[depth]++] = nbr;
                }
            }
            visited_query[u] = true;
        }

        const auto &all_candidates = edge_matrix_.getAllCandidates();
        const VertexId *candidates[64];
        uint32_t candidates_count[64];

        for (int d = 0; d < max_depth; ++d)
        {
            VertexId u = order[d];
            candidates[u] = all_candidates[u].data();
            candidates_count[u] = static_cast<uint32_t>(all_candidates[u].size());
            if (UNLIKELY(candidates_count[u] == 0))
                return 0;
        }

        uint32_t idx[64];
        uint32_t idx_count[64];
        uint32_t idx_embedding[64];
        std::memset(idx, 0, sizeof(idx));

        const VertexId u0 = order[0];
        idx_count[0] = candidates_count[u0];
        for (uint32_t i = 0; i < idx_count[0]; ++i)
        {
            g_valid_candidates[0][i] = i;
        }

        int cur_depth = 0;
        size_t embedding_cnt = 0;
        const int last_depth = max_depth - 1;

        while (true)
        {
            while (idx[cur_depth] < idx_count[cur_depth])
            {
                const uint32_t valid_idx = g_valid_candidates[cur_depth][idx[cur_depth]];
                const VertexId u = order[cur_depth];
                const VertexId v = candidates[u][valid_idx];

                if (visited_[v])
                {
                    idx[cur_depth] += 1;
                    continue;
                }

                embedding_[u] = v;
                idx_embedding[u] = valid_idx;
                visited_[v] = true;
                idx[cur_depth] += 1;

                if (cur_depth == last_depth)
                {
                    embedding_cnt += 1;
                    visited_[v] = false;

                    if (limit > 0 && embedding_cnt >= limit)
                    {
                        early_stopped_ = true;
                        count_ = embedding_cnt;
                        return count_;
                    }
                }
                else
                {
                    const int next_depth = cur_depth + 1;
                    const VertexId u_next = order[next_depth];
                    const uint32_t num_bn = bn_count[next_depth];

                    if (UNLIKELY(num_bn == 0))
                    {

                        idx[next_depth] = 0;
                        idx_count[next_depth] = candidates_count[u_next];
                        for (uint32_t i = 0; i < idx_count[next_depth]; ++i)
                        {
                            g_valid_candidates[next_depth][i] = i;
                        }
                        cur_depth = next_depth;
                    }
                    else
                    {

                        const VertexId first_bn = bn[next_depth][0];
                        const uint32_t first_bn_idx = idx_embedding[first_bn];

                        uint32_t first_count = 0;
                        const uint32_t *first_cands = edge_matrix_.getCandidatesByIdx(
                            first_bn, u_next, first_bn_idx, first_count);

                        if (first_cands == nullptr || first_count == 0)
                        {
                            visited_[v] = false;
                            continue;
                        }

                        uint32_t valid_count = first_count;
                        uint32_t *cur_valid = g_valid_candidates[next_depth];
                        std::memcpy(cur_valid, first_cands, first_count * sizeof(uint32_t));

                        for (uint32_t i = 1; i < num_bn && valid_count > 0; ++i)
                        {
                            const VertexId cur_bn = bn[next_depth][i];
                            const uint32_t cur_bn_idx = idx_embedding[cur_bn];

                            uint32_t cur_count = 0;
                            const uint32_t *cur_cands = edge_matrix_.getCandidatesByIdx(
                                cur_bn, u_next, cur_bn_idx, cur_count);

                            if (cur_cands == nullptr || cur_count == 0)
                            {
                                valid_count = 0;
                                break;
                            }

                            valid_count = static_cast<uint32_t>(
                                intersect_adaptive(cur_cands, cur_count,
                                                   cur_valid, valid_count, g_temp_buffer));

                            std::swap(g_temp_buffer, cur_valid);
                            g_valid_candidates[next_depth] = cur_valid;
                        }

                        if (valid_count > 0)
                        {
                            idx[next_depth] = 0;
                            idx_count[next_depth] = valid_count;
                            cur_depth = next_depth;
                        }
                        else
                        {
                            visited_[v] = false;
                        }
                    }
                }
            }

            cur_depth -= 1;
            if (cur_depth < 0)
                break;

            visited_[embedding_[order[cur_depth]]] = false;
        }

        count_ = embedding_cnt;
        return count_;
    }

    size_t unified_enumerate(
        const Graph &query_graph,
        const Graph &data_graph,
        const UnifiedEdgeMatrix &edge_matrix,
        const VertexId *order,
        size_t order_size,
        size_t limit,
        bool *early_stopped)
    {
        UnifiedEnumerator enumerator(query_graph, data_graph, edge_matrix);
        size_t count = enumerator.enumerate(order, order_size, limit);
        if (early_stopped)
        {
            *early_stopped = enumerator.wasEarlyStopped();
        }
        return count;
    }

}
