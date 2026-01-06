#ifndef SLIM_ENUMERATOR_H
#define SLIM_ENUMERATOR_H

#include <slim/index.h>
#include <slim/graph.h>
#include <slim/enumerate/candidate_storage.h>
#include <cstdint>
#include <cstring>

namespace slim
{

    class Enumerator
    {
    public:
        static constexpr size_t MAX_QUERY_SIZE = 64;

        Enumerator(const Graph &query_graph,
                   const Graph &data_graph,
                   const TreeEdgeCandidates &te_candidates)
            : query_graph_(query_graph),
              data_graph_(data_graph),
              te_candidates_(te_candidates),
              visited_(nullptr),
              data_size_(data_graph.numVertices()),
              count_(0),
              early_stopped_(false)
        {
            visited_ = new bool[data_size_]();
        }

        ~Enumerator()
        {
            delete[] visited_;
        }

        Enumerator(const Enumerator &) = delete;
        Enumerator &operator=(const Enumerator &) = delete;

        size_t enumerate(const QueryTree &tree,
                         const CandidateSet &initial_candidates,
                         size_t limit = 0);

        bool wasEarlyStopped() const { return early_stopped_; }
        size_t count() const { return count_; }

    private:
        void generateValidCandidateIndex(
            int depth,
            uint32_t *idx_embedding,
            uint32_t *idx_count,
            uint32_t **valid_candidate_idx,
            VertexId bn[][MAX_QUERY_SIZE],
            uint32_t *bn_cnt,
            VertexId *order,
            uint32_t *temp_buffer,
            const VertexId **candidates);

        const Graph &query_graph_;
        const Graph &data_graph_;
        const TreeEdgeCandidates &te_candidates_;

        bool *visited_;
        size_t data_size_;
        size_t count_;
        bool early_stopped_;
    };

}

#endif
