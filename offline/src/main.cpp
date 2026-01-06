#include <slim/offline/builder.h>
#include <slim/graph.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <stdexcept>

using namespace slim;
using namespace slim::offline;

struct Args {
    std::string input_path;
    std::string output_path;
    double alpha = 1.0;
    int num_threads = 0;
};

void printUsage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [OPTIONS]\n\n"
              << "Options:\n"
              << "  -i, --input <file>     Input graph file (required)\n"
              << "  -o, --output <file>    Output index file (required)\n"
              << "  -a, --alpha <value>    Alpha parameter (default: 1.0)\n"
              << "  -t, --threads <num>    Number of threads (default: auto)\n"
              << "  -h, --help             Show this help message\n";
}

Args parseArgs(int argc, char* argv[]) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (arg == "-i" || arg == "--input") {
            if (i + 1 < argc) {
                args.input_path = argv[++i];
            } else {
                throw std::runtime_error("Missing value for " + arg);
            }
        } else if (arg == "-o" || arg == "--output") {
            if (i + 1 < argc) {
                args.output_path = argv[++i];
            } else {
                throw std::runtime_error("Missing value for " + arg);
            }
        } else if (arg == "-a" || arg == "--alpha") {
            if (i + 1 < argc) {
                args.alpha = std::stod(argv[++i]);
            } else {
                throw std::runtime_error("Missing value for " + arg);
            }
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                args.num_threads = std::stoi(argv[++i]);
            } else {
                throw std::runtime_error("Missing value for " + arg);
            }
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (args.input_path.empty() || args.output_path.empty()) {
        throw std::runtime_error("Missing required arguments: --input and --output");
    }

    return args;
}

Graph loadGraphFromEdgeList(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open graph file: " + path);
    }

    std::vector<std::pair<VertexId, VertexId>> edges;
    std::unordered_map<VertexId, Label> labels_map;
    VertexId max_vertex = 0;

    std::string line;
    while (std::getline(file, line)) {

        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        VertexId u, v;
        Label label_u, label_v;

        iss >> u >> v;

        if (iss >> label_u >> label_v) {

            labels_map[u] = label_u;
            labels_map[v] = label_v;
        } else {

            labels_map[u] = u;
            labels_map[v] = v;
        }

        edges.push_back({u, v});
        max_vertex = std::max(max_vertex, std::max(u, v));
    }

    size_t num_vertices = max_vertex + 1;
    std::vector<Label> labels(num_vertices);
    for (size_t i = 0; i < num_vertices; ++i) {
        auto it = labels_map.find(i);
        labels[i] = (it != labels_map.end()) ? it->second : static_cast<Label>(i);
    }

    return Graph(num_vertices, edges, labels);
}

int main(int argc, char* argv[]) {
    try {

        Args args = parseArgs(argc, argv);

        std::cout << "SLIM Offline Index Builder\n";

        std::cout << "Loading graph from " << args.input_path << "...\n";
        auto load_start = std::chrono::high_resolution_clock::now();

        Graph graph = loadGraphFromEdgeList(args.input_path);

        auto load_end = std::chrono::high_resolution_clock::now();
        auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            load_end - load_start).count();

        std::cout << "Graph loaded in " << load_duration << " ms\n";

        auto build_start = std::chrono::high_resolution_clock::now();

        IndexBuilder builder(graph, args.alpha);
        if (args.num_threads > 0) {
            builder.setNumThreads(args.num_threads);
        }
        std::unique_ptr<Index> index = builder.build();

        auto build_end = std::chrono::high_resolution_clock::now();
        auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            build_end - build_start).count();

        std::cout << "\nIndex built in " << build_duration << " ms\n";

        std::cout << "\nSaving index to " << args.output_path << "...\n";
        auto save_start = std::chrono::high_resolution_clock::now();

        index->save(args.output_path);

        auto save_end = std::chrono::high_resolution_clock::now();
        auto save_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            save_end - save_start).count();

        std::cout << "Index saved in " << save_duration << " ms\n";

        auto total_duration = load_duration + build_duration + save_duration;
        std::cout << "  Total time: " << total_duration << " ms\n";
        std::cout << "  Index size: " << index->numVertices() << " vertices\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
