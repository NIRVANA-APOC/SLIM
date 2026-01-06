#ifndef SLIM_MATCHER_H
#define SLIM_MATCHER_H

#include <slim/index.h>
#include <slim/graph.h>
#include <slim/filter.h>
#include <slim/enumerate/enumerator.h>
#include <slim/enumerate/candidate_storage.h>
#include <string>
#include <memory>
#include <chrono>

namespace slim
{

    class Matcher
    {
    public:
        explicit Matcher(const std::string &index_path);

        std::pair<size_t, bool> match(const Graph &query_graph,
                                      size_t limit = 0);

        struct PerformanceStats
        {
            size_t filter_time_us = 0;
            size_t storage_time_us = 0;
            size_t enumerate_time_us = 0;
            size_t total_time_us = 0;
            size_t candidate_count = 0;
            size_t match_count = 0;
            bool early_stopped = false;
        };

        PerformanceStats getStats() const { return stats_; }

        void setDataGraph(std::unique_ptr<Graph> data_graph)
        {
            data_graph_ = std::move(data_graph);
        }

        const Graph &dataGraph() const
        {
            if (!data_graph_)
            {
                throw std::runtime_error("Data graph not loaded");
            }
            return *data_graph_;
        }

        const Index &getIndex() const { return *index_; }

    private:
        std::unique_ptr<Index> index_;
        std::unique_ptr<Graph> data_graph_;
        PerformanceStats stats_;
    };

}

#endif
