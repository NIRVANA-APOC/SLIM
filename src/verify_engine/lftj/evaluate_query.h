#ifndef SLIM_VERIFY_ENGINE_LFTJ_EVALUATE_QUERY_H
#define SLIM_VERIFY_ENGINE_LFTJ_EVALUATE_QUERY_H

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <limits>
#include <vector>

#include "computesetintersection.h"
#include "graph.h"

namespace slim_lftj {
class EvaluateQuery {
public:
    static std::size_t LFTJ(const Graph* data_graph, const Graph* query_graph,
                            Edges*** edge_matrix, ui** candidates, ui* candidates_count,
                            ui* order, std::size_t output_limit_num, std::size_t& call_count);

private:
    static void generateBN(const Graph* query_graph, ui* order, ui**& bn, ui*& bn_count);
    static void allocateBuffer(const Graph* data_graph, const Graph* query_graph,
                               ui* candidates_count, ui*& idx, ui*& idx_count, ui*& embedding,
                               ui*& idx_embedding, ui*& temp_buffer, ui**& valid_candidate_idx,
                               bool*& visited_vertices);
    static void generateValidCandidateIndex(ui depth, ui* idx_embedding, ui* idx_count,
                                            ui** valid_candidate_index, Edges*** edge_matrix,
                                            ui** bn, ui* bn_cnt, ui* order, ui*& temp_buffer,
                                            ui* candidates_count);
    static void releaseBuffer(ui query_vertices_num, ui* idx, ui* idx_count, ui* embedding,
                              ui* idx_embedding, ui* temp_buffer, ui** valid_candidate_idx,
                              bool* visited_vertices, ui** bn, ui* bn_count);
};

inline std::size_t EvaluateQuery::LFTJ(const Graph* data_graph, const Graph* query_graph,
                                       Edges*** edge_matrix, ui** candidates, ui* candidates_count,
                                       ui* order, std::size_t output_limit_num, std::size_t& call_count) {
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v]) {
                idx[cur_depth] += 1;
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;

                if (embedding_cnt >= output_limit_num) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count,
                                            valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer, candidates_count);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
            visited_vertices[embedding[u]] = false;
        }
    }

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

inline void EvaluateQuery::generateBN(const Graph* query_graph, ui* order,
                                      ui**& bn, ui*& bn_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

inline void EvaluateQuery::allocateBuffer(const Graph* data_graph, const Graph* query_graph,
                                          ui* candidates_count, ui*& idx, ui*& idx_count, ui*& embedding, ui*& idx_embedding,
                                          ui*& temp_buffer, ui**& valid_candidate_idx, bool*& visited_vertices) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num) {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[query_vertices_num];
    idx_count = new ui[query_vertices_num];
    embedding = new ui[query_vertices_num];
    idx_embedding = new ui[query_vertices_num];
    visited_vertices = new bool[data_vertices_num];
    temp_buffer = new ui[max_candidates_num + kAvx2Slack];
    valid_candidate_idx = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        valid_candidate_idx[i] = new ui[max_candidates_num + kAvx2Slack];
    }

    std::fill(visited_vertices, visited_vertices + data_vertices_num, false);
}

inline void EvaluateQuery::releaseBuffer(ui query_vertices_num, ui* idx,
                                         ui* idx_count, ui* embedding, ui* idx_embedding, ui* temp_buffer,
                                         ui** valid_candidate_idx, bool* visited_vertices, ui** bn, ui* bn_count) {
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    delete[] visited_vertices;
    delete[] bn_count;
    delete[] temp_buffer;
    for (ui i = 0; i < query_vertices_num; ++i) {
        delete[] valid_candidate_idx[i];
        delete[] bn[i];
    }

    delete[] valid_candidate_idx;
    delete[] bn;
}

inline void EvaluateQuery::generateValidCandidateIndex(ui depth, ui* idx_embedding,
                                                       ui* idx_count, ui** valid_candidate_index, Edges*** edge_matrix, ui** bn,
                                                       ui* bn_cnt, ui* order, ui*& temp_buffer, ui* candidates_count) {
    VertexID u = order[depth];

    if (bn_cnt[depth] == 0) {
        const ui count = candidates_count[u];
        for (ui i = 0; i < count; ++i) {
            valid_candidate_index[depth][i] = i;
        }
        idx_count[depth] = count;
        return;
    }

    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;

    Edges& previous_edge = *edge_matrix[previous_bn][u];

    valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i) {
        VertexID current_bn = bn[depth][i];
        Edges& current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                        temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}
}

#endif
