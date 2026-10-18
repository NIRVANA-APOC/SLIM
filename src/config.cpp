#include <slim/config.h>

#include <stdexcept>

namespace slim {
void SlimConfig::validate() const {
    if (!(alpha > 0.0)) {
        throw std::invalid_argument("alpha must be positive");
    }
    if (k < 1) {
        throw std::invalid_argument("neighborhood radius k must be >= 1");
    }
    if (!(tau_scale > 0.0)) {
        throw std::invalid_argument("tau_scale must be positive");
    }
    if (threads < 1) {
        throw std::invalid_argument("threads must be >= 1");
    }
}
}
