#include <slim/enumerate/enumerator.h>
#include <slim/enumerate/set_intersection.h>
#include <cstring>
#include <algorithm>

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

namespace slim
{

    alignas(64) static thread_local uint32_t g_valid_buffer[64][8192];
    alignas(64) static thread_local uint32_t g_temp_buffer[8192];

    size_t Enumerator::enumerate(const QueryTree &tree,
                                 const CandidateSet &initial_candidates,
                                 size_t limit)
    {
        count_ = 0;
        early_stopped_ = false;

        const int max_depth = static_cast<int>(tree.numVertices());
        if (UNLIKELY(max_depth == 0))
            return 0;

        VertexId order[MAX_QUERY_SIZE];
        for (int d = 0; d < max_depth; ++d)
        {
            order[d] = tree.order(d);
        }

        VertexId bn[MAX_QUERY_SIZE][MAX_QUERY_SIZE];
        uint32_t bn_count[MAX_QUERY_SIZE];
        std::memset(bn_count, 0, sizeof(bn_count));

        bool visited_query[MAX_QUERY_SIZE] = {false};
        visited_query[order[0]] = true;

        for (int depth = 1; depth < max_depth; ++depth)
        {
            VertexId u = order[depth];
            size_t nbr_count;
            const VertexId *nbrs = query_graph_.neighbors(u, nbr_count);

            for (size_t j = 0; j < nbr_count; ++j)
            {
                VertexId nbr = nbrs[j];
                if (visited_query[nbr])
                {
                    bn[depth][bn_count[depth]++] = nbr;
                }
            }
            visited_query[u] = true;
        }

        const auto &all_candidates = te_candidates_.getAllCandidates();
        const VertexId *candidates[MAX_QUERY_SIZE];
        uint32_t candidates_count[MAX_QUERY_SIZE];

        for (int d = 0; d < max_depth; ++d)
        {
            VertexId u = order[d];
            candidates[u] = all_candidates[u].data();
            candidates_count[u] = static_cast<uint32_t>(all_candidates[u].size());
            if (UNLIKELY(candidates_count[u] == 0))
                return 0;
        }

        uint32_t *valid_candidate_idx[MAX_QUERY_SIZE];
        for (int i = 0; i < max_depth; ++i)
        {
            valid_candidate_idx[i] = g_valid_buffer[i];
        }
        uint32_t *temp_buffer = g_temp_buffer;

        std::memset(visited_, false, data_size_);

        uint32_t idx[MAX_QUERY_SIZE];
        uint32_t idx_count[MAX_QUERY_SIZE];
        VertexId embedding[MAX_QUERY_SIZE];
        uint32_t idx_embedding[MAX_QUERY_SIZE];

        const VertexId start_vertex = order[0];
        idx[0] = 0;
        idx_count[0] = candidates_count[start_vertex];

        for (uint32_t i = 0; i < idx_count[0]; ++i)
        {
            valid_candidate_idx[0][i] = i;
        }

        int cur_depth = 0;
        const int last_depth = max_depth - 1;
        size_t embedding_cnt = 0;

        while (true)
        {
            while (idx[cur_depth] < idx_count[cur_depth])
            {
                const uint32_t valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
                const VertexId u = order[cur_depth];
                const VertexId v = candidates[u][valid_idx];

                if (visited_[v])
                {
                    idx[cur_depth]++;
                    continue;
                }

                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_[v] = true;
                idx[cur_depth]++;

                if (cur_depth == last_depth)
                {
                    embedding_cnt++;
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
                    cur_depth++;
                    idx[cur_depth] = 0;

                    const VertexId u_next = order[cur_depth];
                    const VertexId previous_bn = bn[cur_depth][0];
                    const uint32_t previous_index_id = idx_embedding[previous_bn];

                    const auto &edge = te_candidates_.getEdgeTable(previous_bn, u_next);
                    uint32_t valid_count = edge.offset_[previous_index_id + 1] - edge.offset_[previous_index_id];
                    const uint32_t *previous_candidates = edge.edge_ + edge.offset_[previous_index_id];

                    std::memcpy(valid_candidate_idx[cur_depth], previous_candidates, valid_count * sizeof(uint32_t));

                    uint32_t temp_count;
                    for (uint32_t i = 1; i < bn_count[cur_depth]; ++i)
                    {
                        const VertexId current_bn = bn[cur_depth][i];
                        const auto &current_edge = te_candidates_.getEdgeTable(current_bn, u_next);
                        const uint32_t current_index_id = idx_embedding[current_bn];

                        const uint32_t current_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
                        const uint32_t *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

                        temp_count = static_cast<uint32_t>(
                            intersect_adaptive(current_candidates, current_count,
                                               valid_candidate_idx[cur_depth], valid_count,
                                               temp_buffer));

                        std::swap(temp_buffer, valid_candidate_idx[cur_depth]);
                        valid_count = temp_count;
                    }

                    idx_count[cur_depth] = valid_count;
                }
            }

            cur_depth--;
            if (cur_depth < 0)
                break;
            visited_[embedding[order[cur_depth]]] = false;
        }

        count_ = embedding_cnt;
        return count_;
    }

}
