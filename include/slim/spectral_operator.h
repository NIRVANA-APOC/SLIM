#ifndef SLIM_SPECTRAL_OPERATOR_H
#define SLIM_SPECTRAL_OPERATOR_H

#include <slim/k_hop_neighborhood.h>
#include <slim/label_encoder.h>
#include <slim/labeled_graph.h>
#include <slim/spectral_signature.h>
#include <slim/types.h>

#include <cstddef>
#include <vector>

namespace slim {
class LabelSpectralOperator {
public:
    LabelSpectralOperator(double alpha, const LabelEncoder& encoder)
        : alpha_(alpha), encoder_(encoder) {}

    double alpha() const { return alpha_; }

    double phi_of_label(Label label) const { return encoder_.encode(label); }

    void build_label_spectral_matrix(const LabeledGraph& graph,
                                     const KHopNeighborhood& neighborhood,
                                     const std::vector<VertexId>& hood,
                                     double* const* rows) const;

private:
    double alpha_;
    const LabelEncoder& encoder_;
};

class SpectralEigenSolver {
public:
    void reserve(std::size_t n);

    void build(const LabeledGraph& graph,
               const LabelSpectralOperator& op,
               const KHopNeighborhood& neighborhood,
               const std::vector<VertexId>& hood);

    const double* solve();

    const double* compute(const LabeledGraph& graph,
                          const LabelSpectralOperator& op,
                          const KHopNeighborhood& neighborhood,
                          const std::vector<VertexId>& hood);

    const double* solve_dense(const double* matrix, std::size_t n);

    std::size_t size() const { return n_; }
    const double* eigenvalues() const { return eigenvalues_.data(); }

private:
    static constexpr std::size_t TINY_MAX = 8;

    void reset_matrix();

    std::vector<double> packed_;
    std::vector<double*> rows_;
    std::vector<double> tridiag_;
    std::vector<double> eigenvalues_;
    double tiny_[TINY_MAX * TINY_MAX] = {};
    std::size_t n_ = 0;
    std::size_t capacity_ = 0;
};

std::vector<double> compute_eigenvalues(const std::vector<double>& matrix, std::size_t n);

LocalSpectralSignature compute_local_spectral_signature(
    const LabeledGraph& graph,
    const LabelSpectralOperator& op,
    int k,
    VertexId v,
    KHopNeighborhood& neighborhood,
    const std::vector<VertexId>* known_hood = nullptr,
    SpectralEigenSolver* solver = nullptr);
}

#endif
