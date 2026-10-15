#ifndef SLIM_GRAPH_IO_H
#define SLIM_GRAPH_IO_H

#include <slim/labeled_graph.h>

#include <string>

namespace slim {
LabeledGraph read_graph(const std::string& path);
}

#endif
