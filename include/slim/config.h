#ifndef SLIM_CONFIG_H
#define SLIM_CONFIG_H

#include <cstdint>

namespace slim {
struct SlimConfig {
    double alpha = 1.0;
    int k = 1;
    double tau_scale = 1.0;
    std::uint64_t neighborhood_size_limit = 2048;
    int threads = 1;

    void validate() const;
};
}

#endif
