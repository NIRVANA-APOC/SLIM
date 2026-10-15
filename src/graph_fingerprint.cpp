#include <slim/graph_fingerprint.h>

#include <algorithm>
#include <cstdio>
#include <vector>

namespace slim {
namespace {
constexpr std::uint64_t FNV_OFFSET = 14695981039346656037ull;
constexpr std::uint64_t FNV_PRIME = 1099511628211ull;

void mix_byte(std::uint64_t& hash, std::uint8_t value) {
    hash ^= value;
    hash *= FNV_PRIME;
}

void mix_u64(std::uint64_t& hash, std::uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8) {
        mix_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xFFu));
    }
}
}

std::uint64_t data_graph_fingerprint(const LabeledGraph& graph) {
    std::uint64_t hash = FNV_OFFSET;
    mix_byte(hash, 1);
    mix_u64(hash, static_cast<std::uint64_t>(graph.num_vertices()));
    mix_byte(hash, 2);
    mix_u64(hash, static_cast<std::uint64_t>(graph.num_edges()));
    for (VertexId v = 0; v < graph.num_vertices(); ++v) {
        mix_byte(hash, 3);
        mix_u64(hash, static_cast<std::uint64_t>(graph.label(v)));
        std::size_t degree = 0;
        const VertexId* neighbors = graph.neighbors(v, degree);
        mix_u64(hash, static_cast<std::uint64_t>(degree));
        for (std::size_t i = 0; i < degree; ++i) {
            mix_u64(hash, static_cast<std::uint64_t>(neighbors[i]));
        }
    }
    return hash;
}

std::uint64_t label_alphabet_fingerprint(const LabeledGraph& graph) {
    std::vector<Label> alphabet;
    alphabet.reserve(graph.num_vertices());
    for (VertexId v = 0; v < graph.num_vertices(); ++v) {
        alphabet.push_back(graph.label(v));
    }
    std::sort(alphabet.begin(), alphabet.end());
    alphabet.erase(std::unique(alphabet.begin(), alphabet.end()), alphabet.end());

    std::uint64_t hash = FNV_OFFSET;
    mix_byte(hash, 4);
    mix_u64(hash, static_cast<std::uint64_t>(alphabet.size()));
    for (Label label : alphabet) {
        mix_u64(hash, static_cast<std::uint64_t>(label));
    }
    return hash;
}

std::string fingerprint_hex(std::uint64_t fingerprint) {
    char buffer[17];
    std::snprintf(buffer, sizeof(buffer), "%016llx",
                  static_cast<unsigned long long>(fingerprint));
    return std::string(buffer);
}
}
