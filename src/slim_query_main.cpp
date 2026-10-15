#include <slim/graph_io.h>
#include <slim/slim_matcher.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Options {
    std::string index_path;
    std::string data_path;
    std::string query_path;
    std::string query_list_path;
    std::string output_csv_path;
    std::size_t limit = 0;
    bool metrics = false;
    bool filter_only = false;
    bool trace_dominance = false;
    bool verbose = false;
};

void print_usage(const char* program) {
    std::cout
        << "Usage: " << program << " [OPTIONS]\n\n"
        << "SLIM online stage: spectral filtering + exact verification.\n\n"
        << "  -i, --index <file>          SLIM index file (required)\n"
        << "  -d, --data <file>           data graph file (required)\n"
        << "  -q, --query <file>          single query graph file\n"
        << "      --query-list <file>     one query path per line\n"
        << "  -l, --limit <N>             stop after N matches; 0 = enumerate all\n"
        << "      --output-csv <file>     write machine-readable metrics CSV\n"
        << "      --metrics               print metrics CSV to stdout\n"
        << "      --filter-only           stop after spectral filtering\n"
        << "      --trace-dominance       print coarse/full dominance counts to stderr\n"
        << "  -v, --verbose               verbose output\n"
        << "  -h, --help                  show this help\n";
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
        } else if (arg == "-i" || arg == "--index") {
            options.index_path = next(arg.c_str());
        } else if (arg == "-d" || arg == "--data") {
            options.data_path = next(arg.c_str());
        } else if (arg == "-q" || arg == "--query") {
            options.query_path = next(arg.c_str());
        } else if (arg == "--query-list") {
            options.query_list_path = next(arg.c_str());
        } else if (arg == "-l" || arg == "--limit") {
            options.limit = std::stoull(next(arg.c_str()));
        } else if (arg == "--output-csv") {
            options.output_csv_path = next(arg.c_str());
            options.metrics = true;
        } else if (arg == "--metrics") {
            options.metrics = true;
        } else if (arg == "--filter-only") {
            options.filter_only = true;
        } else if (arg == "--trace-dominance") {
            options.trace_dominance = true;
        } else if (arg == "-v" || arg == "--verbose") {
            options.verbose = true;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }

    if (options.index_path.empty()) {
        throw std::runtime_error("--index is required");
    }
    if (options.data_path.empty()) {
        throw std::runtime_error("--data is required");
    }
    if (options.query_path.empty() == options.query_list_path.empty()) {
        throw std::runtime_error("use exactly one of --query or --query-list");
    }
    return options;
}

std::vector<std::string> read_query_list(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("cannot open query list: " + path);
    }
    std::vector<std::string> queries;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '#') {
            queries.push_back(line);
        }
    }
    return queries;
}

void write_csv_field(std::ostream& out, const std::string& value) {
    if (value.find_first_of(",\"\n\r") == std::string::npos) {
        out << value;
        return;
    }
    out << '"';
    for (char c : value) {
        if (c == '"') out << "\"\"";
        else out << c;
    }
    out << '"';
}

void write_metrics_header(std::ostream& out) {
    out << "query_path,end_to_end_us,match_core_us,index_load_us,data_load_us,query_load_us,"
        << "query_spectral_us,filter_us,verify_us,"
        << "verify_propagation_us,verify_adjacency_us,verify_enumeration_us,"
        << "candidates_after_filter,candidates_after_propagation,match_count,early_stopped,"
        << "coarse_dominance_checks,full_dominance_checks,"
        << "data_vertices,query_vertices,query_edges\n";
}

void write_metrics_row(std::ostream& out, const std::string& query_path,
                       std::uint64_t end_to_end_us,
                       std::uint64_t index_load_us,
                       std::uint64_t data_load_us,
                       std::uint64_t query_load_us,
                       const slim::MatchStats& stats,
                       const slim::LabeledGraph& data_graph,
                       const slim::LabeledGraph& query_graph) {
    write_csv_field(out, query_path);
    out << ',' << end_to_end_us
        << ',' << stats.total_us
        << ',' << index_load_us
        << ',' << data_load_us
        << ',' << query_load_us
        << ',' << stats.query_spectral_us
        << ',' << stats.filter_us
        << ',' << stats.verify_us
        << ',' << stats.verify_propagation_us
        << ',' << stats.verify_adjacency_us
        << ',' << stats.verify_enumeration_us
        << ',' << stats.candidates_after_filter
        << ',' << stats.candidates_after_propagation
        << ',' << stats.match_count
        << ',' << (stats.early_stopped ? 1 : 0)
        << ',' << stats.coarse_dominance_checks
        << ',' << stats.full_dominance_checks
        << ',' << data_graph.num_vertices()
        << ',' << query_graph.num_vertices()
        << ',' << query_graph.num_edges()
        << '\n';
}
}

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);

        const auto index_start = std::chrono::high_resolution_clock::now();
        slim::SlimMatcher matcher(options.index_path);
        const auto index_end = std::chrono::high_resolution_clock::now();
        const auto data_start = std::chrono::high_resolution_clock::now();
        matcher.set_data_graph(slim::read_graph(options.data_path));
        const auto data_end = std::chrono::high_resolution_clock::now();
        const auto index_load_us = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(index_end - index_start).count());
        const auto data_load_us = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::microseconds>(data_end - data_start).count());

        if (options.verbose) {
            const slim::SpectralDominanceIndex& index = matcher.index();
            std::cerr << "index k=" << index.k()
                      << " alpha=" << index.alpha()
                      << " tau=" << index.tau()
                      << " data-vertices=" << index.num_vertices() << "\n";
        }

        std::ofstream csv_file;
        std::ostream* metrics_out = nullptr;
        if (options.metrics) {
            if (!options.output_csv_path.empty()) {
                csv_file.open(options.output_csv_path);
                if (!csv_file.is_open()) {
                    throw std::runtime_error("cannot create metrics CSV: " + options.output_csv_path);
                }
                metrics_out = &csv_file;
            } else {
                metrics_out = &std::cout;
            }
            write_metrics_header(*metrics_out);
        }

        std::vector<std::string> queries;
        if (!options.query_path.empty()) {
            queries.push_back(options.query_path);
        } else {
            queries = read_query_list(options.query_list_path);
        }

        for (const std::string& query_path : queries) {
            const auto query_start = std::chrono::high_resolution_clock::now();
            const slim::LabeledGraph query_graph = slim::read_graph(query_path);
            const auto query_loaded = std::chrono::high_resolution_clock::now();
            const auto query_load_us = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::microseconds>(query_loaded - query_start).count());
            slim::MatchStats stats;
            const slim::MatchResult result = matcher.match(query_graph, options.limit,
                                                           options.filter_only, &stats);
            const auto query_end = std::chrono::high_resolution_clock::now();
            const auto end_to_end_us = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::microseconds>(
                    query_end - query_start).count());

            if (options.trace_dominance) {
                std::cerr << "dominance query=" << query_path
                          << " coarse=" << stats.coarse_dominance_checks
                          << " full=" << stats.full_dominance_checks << "\n";
            }

            if (metrics_out != nullptr) {
                write_metrics_row(*metrics_out, query_path, end_to_end_us,
                                  index_load_us, data_load_us, query_load_us, stats,
                                  matcher.data_graph(), query_graph);
                continue;
            }

            std::cout << "Query time: " << stats.total_us << " us\n"
                      << "Matches found: " << result.match_count << "\n";
            if (result.early_stopped) {
                std::cout << "Early stopped: 1\n";
            }
        }

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
}
