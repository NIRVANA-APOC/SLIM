#ifndef SLIM_CONFIG_H
#define SLIM_CONFIG_H

#include <cstdint>
#include <string>

namespace slim {
enum class RadiusMode { Single, Cumulative };

RadiusMode radius_mode_from_string(const std::string& text);
const char* radius_mode_name(RadiusMode mode);

struct SlimConfig {
    double alpha = 1.0;

    int k = 1;

    double tau_scale = 1.0;

    std::uint64_t neighborhood_size_limit = 2048;

    int threads = 1;

    int amsd_ranks = 6;

    int amsd_bins = 32;

    void validate() const;
};
}

#endif
