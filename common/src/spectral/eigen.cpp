

#include <slim/spectral/eigen.h>
#include <slim/spectral/laplacian.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstring>

namespace slim {

namespace {

    constexpr double CONVERGENCE_TOL = 1e-10;
    constexpr int MAX_ITERATIONS = 30;

    inline double sign(double a, double b) {
        return b >= 0.0 ? std::abs(a) : -std::abs(a);
    }

    void tred2(double* a, int n, double* d, double* e) {
        for (int i = n - 1; i > 0; --i) {
            int l = i - 1;
            double h = 0.0, scale = 0.0;

            if (l > 0) {
                for (int k = 0; k <= l; ++k) {
                    scale += std::abs(a[i * n + k]);
                }
                if (scale == 0.0) {
                    e[i] = a[i * n + l];
                } else {
                    for (int k = 0; k <= l; ++k) {
                        a[i * n + k] /= scale;
                        h += a[i * n + k] * a[i * n + k];
                    }
                    double f = a[i * n + l];
                    double g = (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
                    e[i] = scale * g;
                    h -= f * g;
                    a[i * n + l] = f - g;
                    f = 0.0;
                    for (int j = 0; j <= l; ++j) {
                        a[j * n + i] = a[i * n + j] / h;
                        g = 0.0;
                        for (int k = 0; k <= j; ++k) {
                            g += a[j * n + k] * a[i * n + k];
                        }
                        for (int k = j + 1; k <= l; ++k) {
                            g += a[k * n + j] * a[i * n + k];
                        }
                        e[j] = g / h;
                        f += e[j] * a[i * n + j];
                    }
                    double hh = f / (h + h);
                    for (int j = 0; j <= l; ++j) {
                        f = a[i * n + j];
                        e[j] = g = e[j] - hh * f;
                        for (int k = 0; k <= j; ++k) {
                            a[j * n + k] -= (f * e[k] + g * a[i * n + k]);
                        }
                    }
                }
            } else {
                e[i] = a[i * n + l];
            }
            d[i] = h;
        }
        d[0] = 0.0;
        e[0] = 0.0;
        for (int i = 0; i < n; ++i) {
            d[i] = a[i * n + i];
        }
    }

    void tqli(double* d, double* e, int n) {
        for (int i = 1; i < n; ++i) {
            e[i - 1] = e[i];
        }
        e[n - 1] = 0.0;

        for (int l = 0; l < n; ++l) {
            int iter = 0;
            int m;
            do {
                for (m = l; m < n - 1; ++m) {
                    double dd = std::abs(d[m]) + std::abs(d[m + 1]);
                    if (std::abs(e[m]) <= CONVERGENCE_TOL * dd) break;
                }
                if (m != l) {
                    if (++iter > MAX_ITERATIONS) break;
                    double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                    double r = std::sqrt(g * g + 1.0);
                    g = d[m] - d[l] + e[l] / (g + sign(r, g));
                    double s = 1.0, c = 1.0, p = 0.0;
                    int i;
                    for (i = m - 1; i >= l; --i) {
                        double f = s * e[i];
                        double b = c * e[i];
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
}

inline void eigen2x2(double a, double b, double d, double& e1, double& e2) {
    double trace = a + d;
    double det = a * d - b * b;
    double disc = trace * trace - 4.0 * det;
    if (disc < 0.0) disc = 0.0;
    double sqrtDisc = std::sqrt(disc);
    e1 = 0.5 * (trace - sqrtDisc);
    e2 = 0.5 * (trace + sqrtDisc);
}

inline void eigen3x3(const double* m, double* eig) {

    double a = m[0], b = m[1], c = m[2];
    double d = m[4], e = m[5], f = m[8];

    double p1 = b*b + c*c + e*e;
    if (p1 < 1e-20) {

        eig[0] = a; eig[1] = d; eig[2] = f;
        if (eig[0] > eig[1]) std::swap(eig[0], eig[1]);
        if (eig[1] > eig[2]) std::swap(eig[1], eig[2]);
        if (eig[0] > eig[1]) std::swap(eig[0], eig[1]);
        return;
    }

    double q = (a + d + f) / 3.0;
    double p2 = (a-q)*(a-q) + (d-q)*(d-q) + (f-q)*(f-q) + 2.0*p1;
    double p = std::sqrt(p2 / 6.0);

    double B[9];
    B[0] = (a - q) / p; B[1] = b / p; B[2] = c / p;
    B[3] = b / p; B[4] = (d - q) / p; B[5] = e / p;
    B[6] = c / p; B[7] = e / p; B[8] = (f - q) / p;

    double detB = B[0]*(B[4]*B[8] - B[5]*B[7])
                - B[1]*(B[3]*B[8] - B[5]*B[6])
                + B[2]*(B[3]*B[7] - B[4]*B[6]);
    double r = detB / 2.0;

    if (r <= -1.0) r = -1.0;
    else if (r >= 1.0) r = 1.0;

    double phi = std::acos(r) / 3.0;

    eig[0] = q + 2.0 * p * std::cos(phi + 2.0 * M_PI / 3.0);
    eig[2] = q + 2.0 * p * std::cos(phi);
    eig[1] = 3.0 * q - eig[0] - eig[2];

    if (eig[0] > eig[1]) std::swap(eig[0], eig[1]);
    if (eig[1] > eig[2]) std::swap(eig[1], eig[2]);
    if (eig[0] > eig[1]) std::swap(eig[0], eig[1]);
}

inline void eigen4x4(double* a, double* eig) {

    double d[4], e[4];

    for (int i = 3; i > 0; --i) {
        int l = i - 1;
        double h = 0.0, scale = 0.0;

        if (l > 0) {
            for (int k = 0; k <= l; ++k) {
                scale += std::abs(a[i * 4 + k]);
            }
            if (scale == 0.0) {
                e[i] = a[i * 4 + l];
            } else {
                for (int k = 0; k <= l; ++k) {
                    a[i * 4 + k] /= scale;
                    h += a[i * 4 + k] * a[i * 4 + k];
                }
                double f = a[i * 4 + l];
                double g = (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
                e[i] = scale * g;
                h -= f * g;
                a[i * 4 + l] = f - g;
                f = 0.0;
                for (int j = 0; j <= l; ++j) {
                    a[j * 4 + i] = a[i * 4 + j] / h;
                    g = 0.0;
                    for (int k = 0; k <= j; ++k) {
                        g += a[j * 4 + k] * a[i * 4 + k];
                    }
                    for (int k = j + 1; k <= l; ++k) {
                        g += a[k * 4 + j] * a[i * 4 + k];
                    }
                    e[j] = g / h;
                    f += e[j] * a[i * 4 + j];
                }
                double hh = f / (h + h);
                for (int j = 0; j <= l; ++j) {
                    f = a[i * 4 + j];
                    e[j] = g = e[j] - hh * f;
                    for (int k = 0; k <= j; ++k) {
                        a[j * 4 + k] -= (f * e[k] + g * a[i * 4 + k]);
                    }
                }
            }
        } else {
            e[i] = a[i * 4 + l];
        }
    }
    e[0] = 0.0;
    for (int i = 0; i < 4; ++i) {
        d[i] = a[i * 4 + i];
    }

    for (int i = 1; i < 4; ++i) {
        e[i - 1] = e[i];
    }
    e[3] = 0.0;

    for (int l = 0; l < 4; ++l) {
        int iter = 0;
        int m;
        do {
            for (m = l; m < 3; ++m) {
                double dd = std::abs(d[m]) + std::abs(d[m + 1]);
                if (std::abs(e[m]) <= 1e-10 * dd) break;
            }
            if (m != l) {
                if (++iter > 30) break;
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::sqrt(g * g + 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? r : -r));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; --i) {
                    double f = s * e[i];
                    double b = c * e[i];
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

    for (int i = 0; i < 4; ++i) eig[i] = d[i];
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (eig[i] > eig[j]) std::swap(eig[i], eig[j]);
        }
    }
}

std::vector<double> computeEigenvalues(
    const std::vector<std::vector<double>>& matrix
) {
    const int n = static_cast<int>(matrix.size());
    if (n == 0) return {};
    if (n == 1) return {matrix[0][0]};

    if (n == 2) {
        double e1, e2;
        eigen2x2(matrix[0][0], matrix[0][1], matrix[1][1], e1, e2);
        return {e1, e2};
    }

    if (n == 3) {
        double m[9] = {
            matrix[0][0], matrix[0][1], matrix[0][2],
            matrix[1][0], matrix[1][1], matrix[1][2],
            matrix[2][0], matrix[2][1], matrix[2][2]
        };
        double eig[3];
        eigen3x3(m, eig);
        return {eig[0], eig[1], eig[2]};
    }

    if (n == 4) {
        double m[16] = {
            matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3],
            matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3],
            matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3],
            matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]
        };
        double eig[4];
        eigen4x4(m, eig);
        return {eig[0], eig[1], eig[2], eig[3]};
    }

    std::vector<double> a(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i * n + j] = matrix[i][j];
        }
    }

    std::vector<double> d(n), e(n);

    tred2(a.data(), n, d.data(), e.data());

    tqli(d.data(), e.data(), n);

    std::sort(d.begin(), d.end());

    return d;
}

SpectralSignature computeNeighborhoodSpectrum(
    const Graph& graph,
    VertexId v,
    const LaplacianBuilder& builder
) {
    auto laplacian = builder.buildNeighborhoodLaplacian(graph, v);
    auto eigenvalues = computeEigenvalues(laplacian);

    SpectralSignature sig;
    sig.eigenvalues = std::move(eigenvalues);
    return sig;
}

}
