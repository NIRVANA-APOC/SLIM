

#include <slim/spectral/laplacian.h>
#include <algorithm>
#include <cmath>

namespace slim {

std::vector<std::vector<double>> LaplacianBuilder::buildNeighborhoodLaplacian(
    const Graph& graph,
    VertexId v
) const {

    size_t nbr_count;
    const VertexId* neighbors = graph.neighbors(v, nbr_count);
    const size_t n = nbr_count + 1;

    if (n == 1) {
        double label_value = encoder_.encode(graph.label(v));
        return {{alpha_ * label_value}};
    }

    std::unordered_map<VertexId, size_t> vertex_to_idx;
    vertex_to_idx.reserve(n);
    vertex_to_idx[v] = 0;
    for (size_t i = 0; i < nbr_count; ++i) {
        vertex_to_idx[neighbors[i]] = i + 1;
    }

    std::vector<std::vector<double>> laplacian(n, std::vector<double>(n, 0.0));

    laplacian[0][0] = static_cast<double>(nbr_count) + alpha_ * encoder_.encode(graph.label(v));
    for (size_t i = 0; i < nbr_count; ++i) {
        laplacian[0][i + 1] = -1.0;
        laplacian[i + 1][0] = -1.0;
    }

    for (size_t i = 0; i < nbr_count; ++i) {
        VertexId u = neighbors[i];
        size_t u_idx = i + 1;
        size_t degree_in_neighborhood = 1;

        size_t u_nbr_count;
        const VertexId* u_neighbors = graph.neighbors(u, u_nbr_count);

        for (size_t j = 0; j < u_nbr_count; ++j) {
            VertexId w = u_neighbors[j];
            if (w == v) continue;
            auto it = vertex_to_idx.find(w);
            if (it != vertex_to_idx.end()) {
                size_t w_idx = it->second;
                if (w_idx > u_idx) {
                    laplacian[u_idx][w_idx] = -1.0;
                    laplacian[w_idx][u_idx] = -1.0;
                }
                degree_in_neighborhood++;
            }
        }

        laplacian[u_idx][u_idx] = static_cast<double>(degree_in_neighborhood)
                                 + alpha_ * encoder_.encode(graph.label(u));
    }

    return laplacian;
}

}
