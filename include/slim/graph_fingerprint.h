#ifndef SLIM_GRAPH_FINGERPRINT_H
#define SLIM_GRAPH_FINGERPRINT_H

#include <slim/labeled_graph.h>

#include <cstdint>
#include <string>

namespace slim {
std::uint64_t data_graph_fingerprint(const LabeledGraph& graph);

std::uint64_t label_alphabet_fingerprint(const LabeledGraph& graph);

std::string fingerprint_hex(std::uint64_t fingerprint);
}

#endif
