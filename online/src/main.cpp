

#include <slim/matcher.h>
#include <slim/graph.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <memory>
#include <unordered_map>
#include <cstdlib>
#include <chrono>
#include <numeric>

namespace
{

    struct Args
    {
        std::string index_path;
        std::string query_path;
        std::string data_path;
        size_t limit = 0;
    };

    Args parseArgs(int argc, char *argv[])
    {
        Args args;

        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];

            if (arg == "--index" || arg == "-i")
            {
                if (i + 1 < argc)
                {
                    args.index_path = argv[++i];
                }
            }
            else if (arg == "--query" || arg == "-q")
            {
                if (i + 1 < argc)
                {
                    args.query_path = argv[++i];
                }
            }
            else if (arg == "--data" || arg == "-d")
            {
                if (i + 1 < argc)
                {
                    args.data_path = argv[++i];
                }
            }
            else if (arg == "--limit" || arg == "-l")
            {
                if (i + 1 < argc)
                {
                    args.limit = std::stoull(argv[++i]);
                }
            }
            else if (arg == "--help" || arg == "-h")
            {

                std::exit(0);
            }
        }

        if (args.index_path.empty() || args.query_path.empty())
        {
            throw std::runtime_error("Missing required arguments: --index and --query");
        }

        return args;
    }

    bool isBinaryGraph(const std::string &path)
    {
        std::ifstream file(path, std::ios::binary);
        if (!file.is_open())
            return false;

        uint32_t magic;
        file.read(reinterpret_cast<char *>(&magic), sizeof(magic));
        return magic == 0x47524150;
    }

    slim::Graph loadGraphBinary(const std::string &path)
    {
        std::ifstream file(path, std::ios::binary);
        if (!file.is_open())
        {
            throw std::runtime_error("Failed to open binary graph file: " + path);
        }

        uint32_t magic, num_vertices, num_edges;
        file.read(reinterpret_cast<char *>(&magic), sizeof(magic));
        file.read(reinterpret_cast<char *>(&num_vertices), sizeof(num_vertices));
        file.read(reinterpret_cast<char *>(&num_edges), sizeof(num_edges));

        if (magic != 0x47524150)
        {
            throw std::runtime_error("Invalid binary graph magic number");
        }

        std::vector<slim::Label> labels(num_vertices);
        file.read(reinterpret_cast<char *>(labels.data()), num_vertices * sizeof(slim::Label));

        std::vector<std::pair<slim::VertexId, slim::VertexId>> edges(num_edges);
        file.read(reinterpret_cast<char *>(edges.data()), num_edges * sizeof(std::pair<slim::VertexId, slim::VertexId>));

        return slim::Graph(num_vertices, edges, labels);
    }

    slim::Graph loadGraphText(const std::string &path)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            throw std::runtime_error("Failed to open graph file: " + path);
        }

        std::vector<std::pair<slim::VertexId, slim::VertexId>> edges;
        edges.reserve(10000);

        std::vector<slim::Label> labels;
        slim::VertexId max_vertex = 0;

        char buffer[256];
        std::vector<std::tuple<slim::VertexId, slim::VertexId, slim::Label, slim::Label>> edge_data;
        edge_data.reserve(10000);

        while (file.getline(buffer, sizeof(buffer)))
        {
            if (buffer[0] == '\0' || buffer[0] == '#')
            {
                continue;
            }

            slim::VertexId u, v;
            slim::Label label_u = 0, label_v = 0;

            int parsed = std::sscanf(buffer, "%u %u %u %u", &u, &v, &label_u, &label_v);
            if (parsed >= 2)
            {
                max_vertex = std::max(max_vertex, std::max(u, v));
                edges.push_back({u, v});
                if (parsed >= 4)
                {
                    edge_data.push_back({u, v, label_u, label_v});
                }
                else
                {
                    edge_data.push_back({u, v, u, v});
                }
            }
        }

        size_t num_vertices = max_vertex + 1;
        labels.resize(num_vertices);
        std::vector<bool> label_set(num_vertices, false);

        for (const auto &[u, v, lu, lv] : edge_data)
        {
            if (!label_set[u])
            {
                labels[u] = lu;
                label_set[u] = true;
            }
            if (!label_set[v])
            {
                labels[v] = lv;
                label_set[v] = true;
            }
        }

        for (size_t i = 0; i < num_vertices; ++i)
        {
            if (!label_set[i])
            {
                labels[i] = static_cast<slim::Label>(i);
            }
        }

        return slim::Graph(num_vertices, edges, labels);
    }

    slim::Graph loadGraph(const std::string &path)
    {
        if (isBinaryGraph(path))
        {
            return loadGraphBinary(path);
        }
        else
        {
            return loadGraphText(path);
        }
    }

    void saveGraphBinary(const slim::Graph &graph, const std::string &path)
    {
        std::ofstream file(path, std::ios::binary);
        if (!file.is_open())
        {
            throw std::runtime_error("Failed to create binary graph file: " + path);
        }

        uint32_t magic = 0x47524150;
        uint32_t num_vertices = static_cast<uint32_t>(graph.numVertices());
        uint32_t num_edges = static_cast<uint32_t>(graph.numEdges());

        file.write(reinterpret_cast<const char *>(&magic), sizeof(magic));
        file.write(reinterpret_cast<const char *>(&num_vertices), sizeof(num_vertices));
        file.write(reinterpret_cast<const char *>(&num_edges), sizeof(num_edges));

        std::vector<slim::Label> labels(num_vertices);
        for (uint32_t v = 0; v < num_vertices; ++v)
        {
            labels[v] = graph.label(v);
        }
        file.write(reinterpret_cast<const char *>(labels.data()), num_vertices * sizeof(slim::Label));

        std::vector<std::pair<slim::VertexId, slim::VertexId>> edges;
        edges.reserve(num_edges);
        for (slim::VertexId u = 0; u < num_vertices; ++u)
        {
            size_t nbr_count;
            const slim::VertexId *neighbors = graph.neighbors(u, nbr_count);
            for (size_t i = 0; i < nbr_count; ++i)
            {
                if (u < neighbors[i])
                {
                    edges.push_back({u, neighbors[i]});
                }
            }
        }
        file.write(reinterpret_cast<const char *>(edges.data()), edges.size() * sizeof(std::pair<slim::VertexId, slim::VertexId>));
    }

}

int main(int argc, char *argv[])
{
    try
    {
        auto start_total = std::chrono::high_resolution_clock::now();

        Args args = parseArgs(argc, argv);

        auto start_index = std::chrono::high_resolution_clock::now();
        slim::Matcher matcher(args.index_path);
        auto end_index = std::chrono::high_resolution_clock::now();

        auto start_data = std::chrono::high_resolution_clock::now();
        if (!args.data_path.empty())
        {
            auto data_graph = std::make_unique<slim::Graph>(loadGraph(args.data_path));
            matcher.setDataGraph(std::move(data_graph));
        }
        auto end_data = std::chrono::high_resolution_clock::now();

        auto start_query = std::chrono::high_resolution_clock::now();
        slim::Graph query_graph = loadGraph(args.query_path);
        auto end_query = std::chrono::high_resolution_clock::now();

        auto start_match = std::chrono::high_resolution_clock::now();
        auto [match_count, early_stopped] = matcher.match(query_graph, args.limit);
        auto end_match = std::chrono::high_resolution_clock::now();

        auto query_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                 end_match - start_match)
                                 .count();

        std::cout << "=================================================" << std::endl;

        std::cout << "Query time: " << query_time_us << " us";
        if (query_time_us >= 1000)
        {
            std::cout << " (" << (query_time_us / 1000.0) << " ms)";
        }
        std::cout << std::endl;

        std::cout << "Matches found: " << match_count << std::endl;

        std::cout << "=================================================" << std::endl;

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
