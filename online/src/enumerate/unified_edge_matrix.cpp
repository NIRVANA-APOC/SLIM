#include <slim/enumerate/unified_edge_matrix.h>
#include <algorithm>
#include <cstring>

namespace slim
{

    void UnifiedEdgeMatrix::build(const Graph &data_graph,
                                  const Graph &query_graph,
                                  const CandidateSet &initial_candidates)
    {
        const size_t query_size = query_graph.numVertices();
        const size_t data_size = data_graph.numVertices();

        edge_matrix_.clear();
        edge_matrix_.resize(query_size);
        for (size_t i = 0; i < query_size; ++i)
        {
            edge_matrix_[i].resize(query_size);
        }

        candidates_.resize(query_size);
        std::vector<uint32_t> candidates_count(query_size);
        for (VertexId u = 0; u < query_size; ++u)
        {
            const auto &cands = initial_candidates.get(u);
            candidates_[u].assign(cands.begin(), cands.end());
            candidates_count[u] = static_cast<uint32_t>(cands.size());
        }

        std::vector<uint32_t> flag(data_size, 0);
        std::vector<VertexId> updated_flag;
        updated_flag.reserve(data_size);

        std::vector<VertexId> build_order(query_size);
        for (VertexId i = 0; i < query_size; ++i)
        {
            build_order[i] = i;
        }
        std::sort(build_order.begin(), build_order.end(),
                  [&query_graph](VertexId l, VertexId r)
                  {
                      if (query_graph.degree(l) == query_graph.degree(r))
                      {
                          return l < r;
                      }
                      return query_graph.degree(l) > query_graph.degree(r);
                  });

        std::vector<uint32_t> temp_edges;
        temp_edges.reserve(data_graph.numEdges() * 2);

        for (VertexId u : build_order)
        {
            size_t u_nbrs_count;
            const VertexId *u_nbrs = query_graph.neighbors(u, u_nbrs_count);

            updated_flag.clear();

            for (size_t i = 0; i < u_nbrs_count; ++i)
            {
                VertexId u_nbr = u_nbrs[i];

                if (edge_matrix_[u][u_nbr] != nullptr)
                {
                    continue;
                }

                if (updated_flag.empty())
                {
                    for (uint32_t j = 0; j < candidates_count[u]; ++j)
                    {
                        VertexId v = candidates_[u][j];
                        flag[v] = j + 1;
                        updated_flag.push_back(v);
                    }
                }

                edge_matrix_[u_nbr][u] = std::make_unique<UnifiedEdge>();
                auto &edge_nbr_u = *edge_matrix_[u_nbr][u];
                edge_nbr_u.vertex_count_ = candidates_count[u_nbr];
                edge_nbr_u.offset_ = new uint32_t[candidates_count[u_nbr] + 1];

                edge_matrix_[u][u_nbr] = std::make_unique<UnifiedEdge>();
                auto &edge_u_nbr = *edge_matrix_[u][u_nbr];
                edge_u_nbr.vertex_count_ = candidates_count[u];
                edge_u_nbr.offset_ = new uint32_t[candidates_count[u] + 1];
                std::fill(edge_u_nbr.offset_, edge_u_nbr.offset_ + candidates_count[u] + 1, 0);

                temp_edges.clear();
                uint32_t local_max_degree = 0;

                for (uint32_t j = 0; j < candidates_count[u_nbr]; ++j)
                {
                    VertexId v = candidates_[u_nbr][j];
                    edge_nbr_u.offset_[j] = static_cast<uint32_t>(temp_edges.size());

                    size_t v_nbrs_count;
                    const VertexId *v_nbrs = data_graph.neighbors(v, v_nbrs_count);

                    uint32_t local_degree = 0;
                    for (size_t k = 0; k < v_nbrs_count; ++k)
                    {
                        VertexId v_nbr = v_nbrs[k];
                        if (flag[v_nbr] != 0)
                        {
                            uint32_t position = flag[v_nbr] - 1;
                            temp_edges.push_back(position);
                            edge_u_nbr.offset_[position + 1] += 1;
                            local_degree += 1;
                        }
                    }

                    if (local_degree > local_max_degree)
                    {
                        local_max_degree = local_degree;
                    }
                }

                edge_nbr_u.offset_[candidates_count[u_nbr]] = static_cast<uint32_t>(temp_edges.size());
                edge_nbr_u.max_degree_ = local_max_degree;
                edge_nbr_u.edge_count_ = static_cast<uint32_t>(temp_edges.size());
                edge_nbr_u.edge_ = new uint32_t[edge_nbr_u.edge_count_];
                std::memcpy(edge_nbr_u.edge_, temp_edges.data(), edge_nbr_u.edge_count_ * sizeof(uint32_t));

                edge_u_nbr.edge_count_ = edge_nbr_u.edge_count_;
                edge_u_nbr.edge_ = new uint32_t[edge_u_nbr.edge_count_];

                local_max_degree = 0;
                for (uint32_t j = 1; j <= candidates_count[u]; ++j)
                {
                    if (edge_u_nbr.offset_[j] > local_max_degree)
                    {
                        local_max_degree = edge_u_nbr.offset_[j];
                    }
                    edge_u_nbr.offset_[j] += edge_u_nbr.offset_[j - 1];
                }
                edge_u_nbr.max_degree_ = local_max_degree;

                for (uint32_t j = 0; j < candidates_count[u_nbr]; ++j)
                {
                    for (uint32_t k = edge_nbr_u.offset_[j]; k < edge_nbr_u.offset_[j + 1]; ++k)
                    {
                        uint32_t end = edge_nbr_u.edge_[k];
                        edge_u_nbr.edge_[edge_u_nbr.offset_[end]++] = j;
                    }
                }

                for (uint32_t j = candidates_count[u]; j >= 1; --j)
                {
                    edge_u_nbr.offset_[j] = edge_u_nbr.offset_[j - 1];
                }
                edge_u_nbr.offset_[0] = 0;
            }

            for (VertexId v : updated_flag)
            {
                flag[v] = 0;
            }
        }
    }

}
