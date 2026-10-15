#include <slim/config.h>
#include <slim/graph_io.h>
#include <slim/spectral_dominance_index.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
struct Options {
    std::string input_path;
    std::string output_path;
    slim::SlimConfig config;
    bool verbose = false;
};

void print_usage(const char* program) {
    std::cout
        << "Usage: " << program << " [OPTIONS]\n\n"
        << "SLIM offline stage: builds spectral signatures and the dominance index.\n\n"
        << "  -i, --input <file>                 input graph file (required)\n"
        << "  -o, --output <file>                output index file (required)\n"
        << "  -a, --alpha <value>                label weight alpha > 0 (default: 1.0)\n"
        << "  -k, --k <value>                    neighborhood radius (default: 1)\n"
        << "      --tau-scale <value>            tau = value * lambda_max / sqrt(|V|) (default: 1.0)\n"
        << "      --neighborhood-size-limit <N>  skip eigen-decomposition for neighborhoods\n"
        << "                                      larger than N; 0 disables the guard (default: 2048)\n"
        << "      --threads <N>                  worker threads for the signature pass (default: 1)\n"
        << "  -v, --verbose                      verbose output\n"
        << "  -h, --help                         show this help\n";
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("missing value for ") + name);
            }
            return argv[++i];
        };
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (arg == "-i" || arg == "--input") {
            options.input_path = next(arg.c_str());
        } else if (arg == "-o" || arg == "--output") {
            options.output_path = next(arg.c_str());
        } else if (arg == "-a" || arg == "--alpha") {
            options.config.alpha = std::stod(next(arg.c_str()));
        } else if (arg == "-k" || arg == "--k") {
            options.config.k = std::stoi(next(arg.c_str()));
        } else if (arg == "--tau-scale") {
            options.config.tau_scale = std::stod(next(arg.c_str()));
        } else if (arg == "--neighborhood-size-limit") {
            options.config.neighborhood_size_limit =
                static_cast<std::uint64_t>(std::stoull(next(arg.c_str())));
        } else if (arg == "--threads") {
            options.config.threads = std::stoi(next(arg.c_str()));
        } else if (arg == "-v" || arg == "--verbose") {
            options.verbose = true;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    if (options.input_path.empty() || options.output_path.empty()) {
        throw std::runtime_error("both --input and --output are required");
    }
    options.config.validate();
    return options;
}
}

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);

        const auto start = std::chrono::high_resolution_clock::now();

        const slim::LabeledGraph graph = slim::read_graph(options.input_path);
        if (options.verbose) {
            std::cout << "graph: " << graph.num_vertices() << " vertices, "
                      << graph.num_edges() << " edges\n";
        }

        slim::OfflineIndexBuilder builder(graph, options.config);
        auto index = builder.build();

        if (options.verbose) {
            std::cout << "signatures: " << graph.num_vertices()
                      << " vertices\n";
            std::cout << "tau: " << index->tau() << "\n";
            if (options.config.neighborhood_size_limit > 0) {
                std::cout << "neighborhood-size-limit: "
                          << options.config.neighborhood_size_limit << "\n";
                std::cout << "skipped vertices: " << builder.skipped_vertices() << "\n";
                std::cout << "label-histogram skipped vertices: "
                          << builder.label_skipped_vertices() << "\n";
            }
        }

        index->save(options.output_path);

        const auto end = std::chrono::high_resolution_clock::now();
        const auto milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "index built in " << milliseconds << " ms\n";
        std::cout << "index size: " << graph.num_vertices() << " vertices\n";

        std::ifstream file(options.output_path, std::ios::binary | std::ios::ate);
        if (file.is_open()) {
            std::cout << "index file bytes: " << file.tellg() << "\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
}
