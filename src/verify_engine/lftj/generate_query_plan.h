#ifndef SLIM_VERIFY_ENGINE_LFTJ_GENERATE_QUERY_PLAN_H
#define SLIM_VERIFY_ENGINE_LFTJ_GENERATE_QUERY_PLAN_H

#include <vector>

#include "graph.h"

namespace slim_lftj {
class GenerateQueryPlan {
public:
    static void generateGQLQueryPlan(const Graph* data_graph, const Graph* query_graph,
                                     ui* candidates_count, ui*& order) {
        std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
        std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
        order = new ui[query_graph->getVerticesCount()];

        const VertexID start_vertex = selectGQLStartVertex(query_graph, candidates_count);
        order[0] = start_vertex;
        updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

        for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
            VertexID next_vertex = 0;
            ui min_value = data_graph->getVerticesCount() + 1;
            bool found = false;
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                const VertexID cur_vertex = j;

                if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                    if (candidates_count[cur_vertex] < min_value) {
                        min_value = candidates_count[cur_vertex];
                        next_vertex = cur_vertex;
                        found = true;
                    } else if (candidates_count[cur_vertex] == min_value &&
                               query_graph->getVertexDegree(cur_vertex) >
                                   query_graph->getVertexDegree(next_vertex)) {
                        next_vertex = cur_vertex;
                    }
                }
            }

            if (!found) {
                next_vertex = selectGQLComponentRoot(query_graph, candidates_count,
                                                     visited_vertices);
            }

            updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
            order[i] = next_vertex;
        }
    }

private:
    static VertexID selectGQLStartVertex(const Graph* query_graph, ui* candidates_count) {
        ui start_vertex = 0;

        for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
            const VertexID cur_vertex = i;

            if (candidates_count[cur_vertex] < candidates_count[start_vertex]) {
                start_vertex = cur_vertex;
            } else if (candidates_count[cur_vertex] == candidates_count[start_vertex] &&
                       query_graph->getVertexDegree(cur_vertex) >
                           query_graph->getVertexDegree(start_vertex)) {
                start_vertex = cur_vertex;
            }
        }

        return start_vertex;
    }

    static VertexID selectGQLComponentRoot(const Graph* query_graph, ui* candidates_count,
                                           const std::vector<bool>& visited_vertices) {
        ui start_vertex = 0;
        bool has_start = false;

        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            if (visited_vertices[i]) {
                continue;
            }
            if (!has_start) {
                start_vertex = i;
                has_start = true;
                continue;
            }
            if (candidates_count[i] < candidates_count[start_vertex]) {
                start_vertex = i;
            } else if (candidates_count[i] == candidates_count[start_vertex] &&
                       query_graph->getVertexDegree(i) >
                           query_graph->getVertexDegree(start_vertex)) {
                start_vertex = i;
            }
        }

        return start_vertex;
    }

    static void updateValidVertices(const Graph* query_graph, VertexID query_vertex,
                                    std::vector<bool>& visited, std::vector<bool>& adjacent) {
        visited[query_vertex] = true;
        ui nbr_cnt = 0;
        const ui* nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

        for (ui i = 0; i < nbr_cnt; ++i) {
            adjacent[nbrs[i]] = true;
        }
    }
};
}

#endif
