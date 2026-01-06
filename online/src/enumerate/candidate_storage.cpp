#include <slim/enumerate/candidate_storage.h>
#include <algorithm>
#include <cstring>

namespace slim
{

    void TreeEdgeCandidates::build(const Graph &data_graph,
                                   const QueryTree &tree,
                                   const CandidateSet &initial_candidates)
    {

        size_t query_size = tree.numVertices();
        size_t data_size = data_graph.numVertices();

        edge_matrix_.clear();
        edge_matrix_.resize(query_size);
        for (size_t i = 0; i < query_size; ++i)
        {
            edge_matrix_[i].resize(query_size);
        }
        parent_.resize(query_size, INVALID_VERTEX);
        candidates_.resize(query_size);

        for (VertexId u = 0; u < query_size; ++u)
        {
            const auto &cands = initial_candidates.get(u);
            candidates_[u].assign(cands.begin(), cands.end());
            parent_[u] = tree.parent(u);
        }

        std::vector<uint32_t> v_to_idx(data_size, UINT32_MAX);
        std::vector<uint32_t> temp_edges;
        temp_edges.reserve(data_size);

        for (VertexId u1 = 0; u1 < query_size; ++u1)
        {
            const auto &u1_cands = candidates_[u1];
            if (u1_cands.empty())
                continue;

            for (VertexId u2 = u1 + 1; u2 < query_size; ++u2)
            {
                const auto &u2_cands = candidates_[u2];
                if (u2_cands.empty())
                    continue;

                for (uint32_t i = 0; i < u2_cands.size(); ++i)
                {
                    v_to_idx[u2_cands[i]] = i;
                }

                EdgeTable &edge12 = edge_matrix_[u1][u2];
                edge12.vertex_count_ = static_cast<uint32_t>(u1_cands.size());
                edge12.offset_ = new uint32_t[edge12.vertex_count_ + 1];

                temp_edges.clear();

                for (uint32_t idx1 = 0; idx1 < u1_cands.size(); ++idx1)
                {
                    VertexId v1 = u1_cands[idx1];
                    edge12.offset_[idx1] = static_cast<uint32_t>(temp_edges.size());

                    size_t nbr_count = 0;
                    const VertexId *nbrs = data_graph.neighbors(v1, nbr_count);

                    size_t start = temp_edges.size();
                    for (size_t i = 0; i < nbr_count; ++i)
                    {
                        VertexId nbr = nbrs[i];
                        if (v_to_idx[nbr] != UINT32_MAX)
                        {
                            temp_edges.push_back(v_to_idx[nbr]);
                        }
                    }
                    std::sort(temp_edges.begin() + start, temp_edges.end());
                }

                edge12.offset_[edge12.vertex_count_] = static_cast<uint32_t>(temp_edges.size());
                edge12.edge_count_ = static_cast<uint32_t>(temp_edges.size());

                if (edge12.edge_count_ > 0)
                {
                    edge12.edge_ = new uint32_t[edge12.edge_count_];
                    std::memcpy(edge12.edge_, temp_edges.data(), edge12.edge_count_ * sizeof(uint32_t));
                }

                for (VertexId v : u2_cands)
                {
                    v_to_idx[v] = UINT32_MAX;
                }
                for (uint32_t i = 0; i < u1_cands.size(); ++i)
                {
                    v_to_idx[u1_cands[i]] = i;
                }

                EdgeTable &edge21 = edge_matrix_[u2][u1];
                edge21.vertex_count_ = static_cast<uint32_t>(u2_cands.size());
                edge21.offset_ = new uint32_t[edge21.vertex_count_ + 1];

                temp_edges.clear();

                for (uint32_t idx2 = 0; idx2 < u2_cands.size(); ++idx2)
                {
                    VertexId v2 = u2_cands[idx2];
                    edge21.offset_[idx2] = static_cast<uint32_t>(temp_edges.size());

                    size_t nbr_count = 0;
                    const VertexId *nbrs = data_graph.neighbors(v2, nbr_count);

                    size_t start = temp_edges.size();
                    for (size_t i = 0; i < nbr_count; ++i)
                    {
                        VertexId nbr = nbrs[i];
                        if (v_to_idx[nbr] != UINT32_MAX)
                        {
                            temp_edges.push_back(v_to_idx[nbr]);
                        }
                    }
                    std::sort(temp_edges.begin() + start, temp_edges.end());
                }

                edge21.offset_[edge21.vertex_count_] = static_cast<uint32_t>(temp_edges.size());
                edge21.edge_count_ = static_cast<uint32_t>(temp_edges.size());

                if (edge21.edge_count_ > 0)
                {
                    edge21.edge_ = new uint32_t[edge21.edge_count_];
                    std::memcpy(edge21.edge_, temp_edges.data(), edge21.edge_count_ * sizeof(uint32_t));
                }

                for (VertexId v : u1_cands)
                {
                    v_to_idx[v] = UINT32_MAX;
                }
            }
        }
    }

    void TreeEdgeCandidates::buildOptimized(const Graph &data_graph,
                                            const Graph &query_graph,
                                            const QueryTree &tree,
                                            const CandidateSet &initial_candidates)
    {

        size_t query_size = query_graph.numVertices();
        size_t data_size = data_graph.numVertices();

        edge_matrix_.clear();
        edge_matrix_.resize(query_size);
        for (size_t i = 0; i < query_size; ++i)
        {
            edge_matrix_[i].resize(query_size);
        }
        parent_.resize(query_size, INVALID_VERTEX);
        candidates_.resize(query_size);

        std::vector<uint32_t> candidates_count(query_size);
        for (VertexId u = 0; u < query_size; ++u)
        {
            const auto &cands = initial_candidates.get(u);
            candidates_[u].assign(cands.begin(), cands.end());
            candidates_count[u] = static_cast<uint32_t>(cands.size());
            parent_[u] = tree.parent(u);
        }

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

        std::vector<uint32_t> flag(data_size, 0);
        std::vector<VertexId> updated_flag;
        updated_flag.reserve(data_size);

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

                if (edge_matrix_[u][u_nbr].offset_ != nullptr)
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

                EdgeTable &edge_nbr_u = edge_matrix_[u_nbr][u];
                edge_nbr_u.vertex_count_ = candidates_count[u_nbr];
                edge_nbr_u.offset_ = new uint32_t[candidates_count[u_nbr] + 1];

                EdgeTable &edge_u_nbr = edge_matrix_[u][u_nbr];
                edge_u_nbr.vertex_count_ = candidates_count[u];
                edge_u_nbr.offset_ = new uint32_t[candidates_count[u] + 1];
                std::fill(edge_u_nbr.offset_, edge_u_nbr.offset_ + candidates_count[u] + 1, 0);

                temp_edges.clear();

                for (uint32_t j = 0; j < candidates_count[u_nbr]; ++j)
                {
                    VertexId v = candidates_[u_nbr][j];
                    edge_nbr_u.offset_[j] = static_cast<uint32_t>(temp_edges.size());

                    size_t v_nbrs_count;
                    const VertexId *v_nbrs = data_graph.neighbors(v, v_nbrs_count);

                    for (size_t k = 0; k < v_nbrs_count; ++k)
                    {
                        VertexId v_nbr = v_nbrs[k];
                        if (flag[v_nbr] != 0)
                        {
                            uint32_t position = flag[v_nbr] - 1;
                            temp_edges.push_back(position);
                            edge_u_nbr.offset_[position + 1] += 1;
                        }
                    }
                }

                edge_nbr_u.offset_[candidates_count[u_nbr]] = static_cast<uint32_t>(temp_edges.size());
                edge_nbr_u.edge_count_ = static_cast<uint32_t>(temp_edges.size());
                edge_nbr_u.edge_ = new uint32_t[edge_nbr_u.edge_count_];
                std::memcpy(edge_nbr_u.edge_, temp_edges.data(), edge_nbr_u.edge_count_ * sizeof(uint32_t));

                edge_u_nbr.edge_count_ = edge_nbr_u.edge_count_;
                edge_u_nbr.edge_ = new uint32_t[edge_u_nbr.edge_count_];

                for (uint32_t j = 1; j <= candidates_count[u]; ++j)
                {
                    edge_u_nbr.offset_[j] += edge_u_nbr.offset_[j - 1];
                }

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
