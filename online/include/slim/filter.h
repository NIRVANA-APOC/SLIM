#ifndef SLIM_FILTER_H
#define SLIM_FILTER_H

#include <slim/index.h>
#include <slim/graph.h>
#include <slim/spectral/laplacian.h>
#include <slim/spectral/eigen.h>
#include <unordered_map>
#include <vector>
#include <cstdint>

namespace slim
{

    class CandidateSet
    {
    public:
        void add(VertexId query_vertex, VertexId data_vertex)
        {
            if (query_vertex >= candidates_.size())
            {
                candidates_.resize(query_vertex + 1);
            }
            candidates_[query_vertex].push_back(data_vertex);
        }

        const std::vector<VertexId> &get(VertexId query_vertex) const
        {
            if (query_vertex >= candidates_.size())
            {
                static const std::vector<VertexId> empty;
                return empty;
            }
            return candidates_[query_vertex];
        }

        void clear() { candidates_.clear(); }

        void reserve(size_t num_vertices)
        {
            candidates_.resize(num_vertices);
        }

        std::vector<VertexId> getQueryVertices() const
        {
            std::vector<VertexId> result;
            for (size_t i = 0; i < candidates_.size(); ++i)
            {
                if (!candidates_[i].empty())
                {
                    result.push_back(static_cast<VertexId>(i));
                }
            }
            return result;
        }

    private:
        std::vector<std::vector<VertexId>> candidates_;
    };

    class SpectralFilter
    {
    public:
        SpectralFilter(const Index &index, const Graph &data_graph)
            : index_(index), data_graph_(data_graph) {}

        CandidateSet filter(VertexId query_vertex, const Graph &query_graph) const;

        CandidateSet filterAll(const Graph &query_graph,
                               const std::vector<SpectralSignature> &query_signatures) const;

        void refineWithNLF(const Graph &query_graph, CandidateSet &candidates) const;

        SpectralSignature computeQuerySignature(
            VertexId query_vertex,
            const Graph &query_graph) const;

    private:
        CandidateSet filterWithHierarchicalIndex(
            VertexId query_vertex,
            Label query_label,
            const SpectralSignature &query_sig) const;

        CandidateSet filterWithLabelIndex(
            VertexId query_vertex,
            Label query_label,
            const SpectralSignature &query_sig) const;

        const Index &index_;
        const Graph &data_graph_;
    };

}

#endif
