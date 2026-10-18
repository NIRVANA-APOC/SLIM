#ifndef SLIM_CONFIG_H
#define SLIM_CONFIG_H

#include <cstdint>
#include <string>

namespace slim {
enum class RadiusMode { Single, Cumulative };

RadiusMode radius_mode_from_string(const std::string& text);
const char* radius_mode_name(RadiusMode mode);

enum class RankSelectMode { Adaptive, Uniform, Rank1, All };

RankSelectMode rank_select_mode_from_string(const std::string& text);
const char* rank_select_mode_name(RankSelectMode mode);

struct SlimConfig {
    double alpha = 1.0;

    int k = 1;

    double tau_scale = 1.0;

    std::uint64_t neighborhood_size_limit = 2048;

    int threads = 1;

    double amsd_budget_ratio = 4.0;

    RankSelectMode amsd_rank_select = RankSelectMode::Adaptive;

    int amsd_ranks = 32;

    int amsd_bins = 32;

    std::string rank_signals_dump_path;
    std::string rank_report_path;

    void validate() const;
};
}

#endif
