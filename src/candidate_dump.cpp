#include <slim/candidate_dump.h>

#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace slim {
void write_candidate_dump(const std::string& path, const CandidateSets& candidates) {
    std::ofstream out(path, std::ios::binary);
    if (!out.is_open()) {
        throw std::runtime_error("cannot create candidate dump: " + path);
    }
    out.write("SLIMCAND1", 9);
    const std::uint32_t query_vertices = static_cast<std::uint32_t>(candidates.size());
    out.write(reinterpret_cast<const char*>(&query_vertices), sizeof(query_vertices));
    for (const auto& set : candidates) {
        const std::uint32_t count = static_cast<std::uint32_t>(set.size());
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
        out.write(reinterpret_cast<const char*>(set.data()),
                  static_cast<std::streamsize>(count) * sizeof(VertexId));
    }
    if (!out) {
        throw std::runtime_error("failed writing candidate dump: " + path);
    }
}
}
