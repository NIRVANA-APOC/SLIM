#ifndef SLIM_SPECTRAL_SIGNATURE_H
#define SLIM_SPECTRAL_SIGNATURE_H

#include <slim/types.h>

#include <cmath>
#include <cstddef>
#include <vector>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace slim {
constexpr double DOMINANCE_EPSILON = 1e-6;

struct LocalSpectralSignature {
    std::vector<double> eigenvalues;

    std::size_t dimension() const { return eigenvalues.size(); }

    double lambda_max() const {
        return eigenvalues.empty() ? 0.0 : eigenvalues.back();
    }

    bool dominates(const LocalSpectralSignature& query) const {
        const std::size_t m = query.eigenvalues.size();
        const std::size_t n = eigenvalues.size();
        if (m > n) {
            return false;
        }
        const double* mu = query.eigenvalues.data();
        const double* lambda = eigenvalues.data() + (n - m);

        if (m < 8) {
            for (std::size_t k = 0; k < m; ++k) {
                if (mu[k] > lambda[k] + DOMINANCE_EPSILON) {
                    return false;
                }
            }
            return true;
        }

        std::size_t k = 0;
#if defined(__AVX2__)
        const __m256d eps = _mm256_set1_pd(DOMINANCE_EPSILON);
        for (; k + 4 <= m; k += 4) {
            const __m256d mu_vec = _mm256_loadu_pd(mu + k);
            const __m256d lambda_vec = _mm256_add_pd(_mm256_loadu_pd(lambda + k), eps);
            if (_mm256_movemask_pd(_mm256_cmp_pd(mu_vec, lambda_vec, _CMP_GT_OQ)) != 0) {
                return false;
            }
        }
#endif
        for (; k < m; ++k) {
            if (mu[k] > lambda[k] + DOMINANCE_EPSILON) {
                return false;
            }
        }
        return true;
    }
};

inline bool spectral_dominance(const LocalSpectralSignature& dominated,
                               const LocalSpectralSignature& dominator) {
    return dominator.dominates(dominated);
}

inline double stored_spectral_value(double value) {
    const float rounded = static_cast<float>(value);
    return static_cast<double>(rounded) < value ? static_cast<double>(std::nextafter(rounded, INFINITY))
                                                : static_cast<double>(rounded);
}
}

#endif
