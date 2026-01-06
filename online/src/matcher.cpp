

#include <slim/matcher.h>
#include <slim/filter.h>
#include <iostream>
#include <iomanip>
#include <slim/enumerate/enumerator.h>
#include <slim/enumerate/candidate_storage.h>
#include <chrono>
#include <fstream>

namespace slim
{

    Matcher::Matcher(const std::string &index_path)
    {
        index_ = Index::load(index_path);
        data_graph_ = nullptr;
    }

    std::pair<size_t, bool> Matcher::match(const Graph &query_graph,
                                           size_t limit)
    {
        stats_ = PerformanceStats();

        if (!data_graph_)
        {
            throw std::runtime_error("Data graph not loaded");
        }

        auto start_total = std::chrono::high_resolution_clock::now();

        auto start_filter = std::chrono::high_resolution_clock::now();
        SpectralFilter filter(*index_, *data_graph_);
        std::vector<SpectralSignature> query_signatures(query_graph.numVertices());
        for (VertexId u = 0; u < query_graph.numVertices(); ++u)
        {
            query_signatures[u] = filter.computeQuerySignature(u, query_graph);
        }

        CandidateSet candidates = filter.filterAll(query_graph, query_signatures);

        filter.refineWithNLF(query_graph, candidates);

        auto end_filter = std::chrono::high_resolution_clock::now();
        stats_.filter_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                    end_filter - start_filter)
                                    .count();

        size_t num_query_vertices = query_graph.numVertices();
        size_t num_data_vertices = data_graph_->numVertices();
        size_t label_candidates = 0;
        std::vector<size_t> candidates_count(num_query_vertices);
        for (VertexId u = 0; u < num_query_vertices; ++u)
        {
            candidates_count[u] = candidates.get(u).size();
            stats_.candidate_count += candidates_count[u];
            Label query_label = query_graph.label(u);
            label_candidates += index_->getByLabel(query_label).size();
        }

        QueryTree tree = QueryTree::buildWithCandidates(query_graph, candidates_count);

        size_t total_possible_pairs = num_query_vertices * num_data_vertices;
        double pruning_rate = 0.0;
        if (total_possible_pairs > 0)
        {
            pruning_rate = 1.0 - static_cast<double>(stats_.candidate_count) / total_possible_pairs;
        }

        double spectral_pruning_rate = 0.0;
        if (label_candidates > 0)
        {
            spectral_pruning_rate = 1.0 - static_cast<double>(stats_.candidate_count) / label_candidates;
        }

        auto start_storage = std::chrono::high_resolution_clock::now();

        TreeEdgeCandidates te_candidates;
        te_candidates.buildOptimized(*data_graph_, query_graph, tree, candidates);

        auto end_storage = std::chrono::high_resolution_clock::now();
        stats_.storage_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                     end_storage - start_storage)
                                     .count();

        auto start_enumerate = std::chrono::high_resolution_clock::now();

        Enumerator enumerator(query_graph, *data_graph_, te_candidates);
        size_t match_count = enumerator.enumerate(tree, candidates, limit);
        bool early_stopped = enumerator.wasEarlyStopped();

        auto end_enumerate = std::chrono::high_resolution_clock::now();
        stats_.enumerate_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                       end_enumerate - start_enumerate)
                                       .count();

        stats_.match_count = match_count;
        stats_.early_stopped = early_stopped;

        auto end_total = std::chrono::high_resolution_clock::now();
        stats_.total_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
                                   end_total - start_total)
                                   .count();

        return {match_count, stats_.early_stopped};
    }

}
