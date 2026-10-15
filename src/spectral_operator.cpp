#include <slim/spectral_operator.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace slim {
namespace {
constexpr double CONVERGENCE_TOL = 1e-10;
constexpr int MAX_ITERATIONS = 30;

constexpr std::size_t TINY_MAX = 8;

inline std::size_t packed_offset(std::size_t j, std::size_t n) {
    return j * n - j * (j - 1) / 2;
}

void tridiagonal_ql(double* d, double* e, std::size_t n) {
    const int size = static_cast<int>(n);
    for (int l = 0; l < size; ++l) {
        int iter = 0;
        int m = 0;
        do {
            for (m = l; m < size - 1; ++m) {
                const double dd = std::abs(d[m]) + std::abs(d[m + 1]);
                if (std::abs(e[m]) <= CONVERGENCE_TOL * dd) break;
            }
            if (m != l) {
                if (++iter > MAX_ITERATIONS) break;
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::sqrt(g * g + 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? r : -r));
                double s = 1.0, c = 1.0, p = 0.0;
                int i = 0;
                for (i = m - 1; i >= l; --i) {
                    const double f = s * e[i];
                    const double b = c * e[i];
                    r = std::sqrt(f * f + g * g);
                    e[i + 1] = r;
                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}

template <std::size_t N>
void tridiagonalize(double* const* rows, std::size_t runtime_size, double* __restrict__ v,
                    double* __restrict__ w, double* __restrict__ d, double* __restrict__ e) {
    std::size_t n = runtime_size;
    if constexpr (N != 0) {
        n = N;
    }
    for (std::size_t i = n - 1; i >= 1; --i) {
        double xnorm2 = 0.0;
        for (std::size_t j = 0; j + 1 < i; ++j) {
            const double x = rows[j][i];
            v[j] = x;
            xnorm2 += x * x;
        }
        const double alpha = rows[i - 1][i];
        const double xnorm = std::sqrt(xnorm2);
        double beta = alpha;
        if (xnorm != 0.0) {
            beta = -std::copysign(std::hypot(alpha, xnorm), alpha);
            const double tau = (beta - alpha) / beta;
            const double scale = 1.0 / (alpha - beta);
            for (std::size_t j = 0; j + 1 < i; ++j) {
                v[j] *= scale;
            }
            v[i - 1] = 1.0;
            std::fill(w, w + i, 0.0);
            for (std::size_t j = 0; j < i; ++j) {
                const double* row = rows[j];
                const double vj = v[j];
                double sum = row[j] * vj;
                for (std::size_t k = j + 1; k < i; ++k) {
                    const double a = row[k];
                    sum += a * v[k];
                    w[k] += a * vj;
                }
                w[j] += sum;
            }
            double vw = 0.0;
            for (std::size_t j = 0; j < i; ++j) {
                vw += w[j] * v[j];
            }
            const double shift = -0.5 * tau * tau * vw;
            for (std::size_t j = 0; j < i; ++j) {
                w[j] = tau * w[j] + shift * v[j];
            }
            for (std::size_t j = 0; j < i; ++j) {
                double* row = rows[j];
                const double vj = v[j];
                const double wj = w[j];
                for (std::size_t k = j; k < i; ++k) {
                    row[k] -= vj * w[k] + wj * v[k];
                }
            }
        }
        e[i - 1] = beta;
        d[i] = rows[i][i];
    }
    d[0] = rows[0][0];
}

void eigen2x2(double a, double b, double c, double* eig) {
    const double trace = a + c;
    double disc = trace * trace - 4.0 * (a * c - b * b);
    if (disc < 0.0) disc = 0.0;
    const double root = std::sqrt(disc);
    eig[0] = 0.5 * (trace - root);
    eig[1] = 0.5 * (trace + root);
}

void eigen3x3(const double* m, double* eig) {
    const double a = m[0], b = m[1], c = m[2];
    const double d = m[4], e = m[5], f = m[8];
    constexpr double PI = 3.141592653589793238462643383279502884;

    const double p1 = b * b + c * c + e * e;
    if (p1 < 1e-20) {
        eig[0] = a;
        eig[1] = d;
        eig[2] = f;
        std::sort(eig, eig + 3);
        return;
    }

    const double q = (a + d + f) / 3.0;
    const double p2 = (a - q) * (a - q) + (d - q) * (d - q) + (f - q) * (f - q) + 2.0 * p1;
    const double p = std::sqrt(p2 / 6.0);

    double B[9];
    B[0] = (a - q) / p; B[1] = b / p; B[2] = c / p;
    B[3] = b / p; B[4] = (d - q) / p; B[5] = e / p;
    B[6] = c / p; B[7] = e / p; B[8] = (f - q) / p;

    const double detB = B[0] * (B[4] * B[8] - B[5] * B[7])
                      - B[1] * (B[3] * B[8] - B[5] * B[6])
                      + B[2] * (B[3] * B[7] - B[4] * B[6]);
    double r = detB / 2.0;
    r = std::max(-1.0, std::min(1.0, r));
    const double phi = std::acos(r) / 3.0;

    eig[0] = q + 2.0 * p * std::cos(phi + 2.0 * PI / 3.0);
    eig[2] = q + 2.0 * p * std::cos(phi);
    eig[1] = 3.0 * q - eig[0] - eig[2];
    std::sort(eig, eig + 3);
}

void solve_fixed_rows(std::size_t n, double* const* rows, double* v, double* w,
                      double* d, double* e) {
    switch (n) {
        case 4: tridiagonalize<4>(rows, n, v, w, d, e); break;
        case 5: tridiagonalize<5>(rows, n, v, w, d, e); break;
        case 6: tridiagonalize<6>(rows, n, v, w, d, e); break;
        case 7: tridiagonalize<7>(rows, n, v, w, d, e); break;
        case 8: tridiagonalize<8>(rows, n, v, w, d, e); break;
        default: tridiagonalize<0>(rows, n, v, w, d, e); break;
    }
}
}

void LabelSpectralOperator::build_label_spectral_matrix(
    const LabeledGraph& graph,
    const KHopNeighborhood& neighborhood,
    const std::vector<VertexId>& hood,
    double* const* rows) const {
    const std::size_t n = hood.size();
    for (std::size_t i = 0; i < n; ++i) {
        const VertexId u = hood[i];
        std::size_t count = 0;
        const VertexId* neighbors = graph.neighbors(u, count);
        if (i + 1 < n) {
            std::size_t next_count = 0;
            const VertexId* next_neighbors = graph.neighbors(hood[i + 1], next_count);
            __builtin_prefetch(next_neighbors);
            __builtin_prefetch(next_neighbors + next_count);
        }
        std::size_t degree_in_hood = 0;
        for (std::size_t j = 0; j < count; ++j) {
            const VertexId w = neighbors[j];
            if (!neighborhood.contains(w)) {
                continue;
            }
            ++degree_in_hood;
            const std::size_t wi = static_cast<std::size_t>(neighborhood.local_index(w) - 1);
            if (wi > i) {
                rows[i][wi] = -1.0;
            }
        }
        rows[i][i] = static_cast<double>(degree_in_hood) + alpha_ * encoder_.encode(graph.label(u));
    }
}

void SpectralEigenSolver::reserve(std::size_t n) {
    if (n <= capacity_) {
        return;
    }
    packed_.resize(n * (n + 1) / 2);
    rows_.resize(n);
    tridiag_.resize(3 * n);
    eigenvalues_.resize(n);
    capacity_ = n;
}

void SpectralEigenSolver::reset_matrix() {
    std::fill(packed_.begin(), packed_.begin() + n_ * (n_ + 1) / 2, 0.0);
}

void SpectralEigenSolver::build(const LabeledGraph& graph,
                                const LabelSpectralOperator& op,
                                const KHopNeighborhood& neighborhood,
                                const std::vector<VertexId>& hood) {
    n_ = hood.size();
    reserve(n_);
    if (n_ <= TINY_MAX) {
        std::fill(tiny_, tiny_ + TINY_MAX * TINY_MAX, 0.0);
        for (std::size_t j = 0; j < n_; ++j) {
            rows_[j] = tiny_ + j * TINY_MAX;
        }
    } else {
        reset_matrix();
        for (std::size_t j = 0; j < n_; ++j) {
            rows_[j] = packed_.data() + packed_offset(j, n_) - j;
        }
    }
    op.build_label_spectral_matrix(graph, neighborhood, hood, rows_.data());
}

const double* SpectralEigenSolver::solve() {
    const std::size_t n = n_;
    if (n == 0) {
        return eigenvalues_.data();
    }
    double* d = eigenvalues_.data();
    if (n == 1) {
        d[0] = rows_[0][0];
        return d;
    }
    if (n == 2) {
        eigen2x2(rows_[0][0], rows_[0][1], rows_[1][1], d);
        return d;
    }
    if (n == 3) {
        const double m[9] = {rows_[0][0], rows_[0][1], rows_[0][2],
                             rows_[0][1], rows_[1][1], rows_[1][2],
                             rows_[0][2], rows_[1][2], rows_[2][2]};
        eigen3x3(m, d);
        return d;
    }

    double* v = tridiag_.data();
    double* w = v + n;
    double* e = w + n;
    solve_fixed_rows(n, rows_.data(), v, w, d, e);
    tridiagonal_ql(d, e, n);
    std::sort(d, d + n);
    return d;
}

const double* SpectralEigenSolver::compute(const LabeledGraph& graph,
                                           const LabelSpectralOperator& op,
                                           const KHopNeighborhood& neighborhood,
                                           const std::vector<VertexId>& hood) {
    build(graph, op, neighborhood, hood);
    return solve();
}

const double* SpectralEigenSolver::solve_dense(const double* matrix, std::size_t n) {
    n_ = n;
    reserve(n);
    if (n <= TINY_MAX) {
        std::fill(tiny_, tiny_ + TINY_MAX * TINY_MAX, 0.0);
        for (std::size_t j = 0; j < n; ++j) {
            rows_[j] = tiny_ + j * TINY_MAX;
            for (std::size_t k = j; k < n; ++k) {
                rows_[j][k] = matrix[j * n + k];
            }
        }
    } else {
        reset_matrix();
        for (std::size_t j = 0; j < n; ++j) {
            rows_[j] = packed_.data() + packed_offset(j, n) - j;
            const double* source = matrix + j * n;
            std::memcpy(rows_[j] + j, source + j, (n - j) * sizeof(double));
        }
    }
    return solve();
}

std::vector<double> compute_eigenvalues(const std::vector<double>& matrix, std::size_t n) {
    if (n == 0) {
        return {};
    }
    SpectralEigenSolver solver;
    const double* eigenvalues = solver.solve_dense(matrix.data(), n);
    return std::vector<double>(eigenvalues, eigenvalues + n);
}

LocalSpectralSignature compute_local_spectral_signature(
    const LabeledGraph& graph,
    const LabelSpectralOperator& op,
    int k,
    VertexId v,
    KHopNeighborhood& neighborhood,
    const std::vector<VertexId>* known_hood,
    SpectralEigenSolver* solver) {
    const std::vector<VertexId>& hood =
        (known_hood != nullptr) ? *known_hood : neighborhood.extract(v, k);
    const std::size_t n = hood.size();

    LocalSpectralSignature signature;
    if (n == 1) {
        signature.eigenvalues.push_back(op.alpha() * op.phi_of_label(graph.label(v)));
        return signature;
    }

    static thread_local SpectralEigenSolver fallback;
    SpectralEigenSolver& work = (solver != nullptr) ? *solver : fallback;
    const double* eigenvalues = work.compute(graph, op, neighborhood, hood);
    signature.eigenvalues.assign(eigenvalues, eigenvalues + n);
    return signature;
}
}
