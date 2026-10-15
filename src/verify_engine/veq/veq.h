#ifndef SRC_VEQ_H
#define SRC_VEQ_H

#include <slim/labeled_graph.h>

#include <unordered_map>

#include "include/run.h"

namespace VEQ {
void init_veq_data(const slim::LabeledGraph& graph,
                   const std::unordered_map<int, int>& label_ids);

bool run_veq_query(Graph* query, const int* const* candidate_sets,
                   const int* candidate_sizes, std::uint64_t limit,
                   bool self_refinement, long long& embedding_cnt,
                   int64_t& time_total);

void dispose_veq_query(Graph* query);

void dispose_veq_data();
}

#endif
