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
    if (amsd_ranks < 0) {
        throw std::invalid_argument("amsd_ranks must be >= 0 (0 disables the rank families)");
    }
    if (!(amsd_budget_ratio >= 1.0)) {
        throw std::invalid_argument(
            "amsd_budget_ratio must be >= 1 (the paper requires B >= |V_G|)");
    }
    if (amsd_bins < 1) {
        throw std::invalid_argument("amsd_bins must be >= 1");
    }
}

RankSelectMode rank_select_mode_from_string(const std::string& text) {
    if (text == "adaptive") {
        return RankSelectMode::Adaptive;
    }
    if (text == "uniform") {
        return RankSelectMode::Uniform;
    }
    if (text == "rank1") {
        return RankSelectMode::Rank1;
    }
    if (text == "all") {
        return RankSelectMode::All;
    }
    throw std::invalid_argument("invalid rank selection mode '" + text +
                                "' (expected adaptive, uniform, rank1 or all)");
}

const char* rank_select_mode_name(RankSelectMode mode) {
    switch (mode) {
        case RankSelectMode::Adaptive:
            return "adaptive";
        case RankSelectMode::Uniform:
            return "uniform";
        case RankSelectMode::Rank1:
            return "rank1";
        case RankSelectMode::All:
            break;
    }
    return "all";
}

RadiusMode radius_mode_from_string(const std::string& text) {
    if (text == "single") {
        return RadiusMode::Single;
    }
    if (text == "cumulative") {
        return RadiusMode::Cumulative;
    }
    throw std::invalid_argument("invalid radius mode '" + text +
                                "' (expected single or cumulative)");
}

const char* radius_mode_name(RadiusMode mode) {
    return mode == RadiusMode::Cumulative ? "cumulative" : "single";
}
}
