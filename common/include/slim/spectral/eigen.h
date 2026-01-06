

#ifndef SLIM_EIGEN_H
#define SLIM_EIGEN_H

#include <slim/index.h>
#include <slim/graph.h>
#include <slim/spectral/laplacian.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace slim {

std::vector<double> computeEigenvalues(
    const std::vector<std::vector<double>>& matrix
);

SpectralSignature computeNeighborhoodSpectrum(
    const Graph& graph,
    VertexId v,
    const LaplacianBuilder& builder
);

}

#endif
