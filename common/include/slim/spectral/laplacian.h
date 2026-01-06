

#ifndef SLIM_LAPLACIAN_H
#define SLIM_LAPLACIAN_H

#include <slim/index.h>
#include <slim/graph.h>
#include <vector>
#include <unordered_map>
#include <stdexcept>

namespace slim {

class LaplacianBuilder {
public:

    LaplacianBuilder(double alpha, const LabelEncoder& encoder)
        : alpha_(alpha), encoder_(encoder) {
        if (alpha <= 0.0) {
            throw std::invalid_argument("Alpha must be positive");
        }
    }

    double alpha() const { return alpha_; }

    const LabelEncoder& encoder() const { return encoder_; }

    std::vector<std::vector<double>> buildNeighborhoodLaplacian(
        const Graph& graph,
        VertexId v
    ) const;

private:
    double alpha_;
    const LabelEncoder& encoder_;
};

}

#endif
