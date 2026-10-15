#ifndef VEQ_ENTRY_H
#define VEQ_ENTRY_H

#include <slim/labeled_graph.h>
#include <slim/types.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "veq.h"

namespace VEQ {
inline std::unordered_map<int, int> BuildLabelIndex(
    const slim::LabeledGraph& data_graph) {
  std::vector<int> distinct;
  distinct.reserve(64);
  for (std::size_t v = 0; v < data_graph.num_vertices(); ++v) {
    distinct.push_back(static_cast<int>(data_graph.label(static_cast<slim::VertexId>(v))));
  }
  std::sort(distinct.begin(), distinct.end());
  distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
  std::unordered_map<int, int> label_ids;
  label_ids.reserve(distinct.size() * 2);
  for (std::size_t i = 0; i < distinct.size(); ++i) {
    label_ids[distinct[i]] = static_cast<int>(i);
  }
  return label_ids;
}

inline void FillVeqGraph(Graph& out, const slim::LabeledGraph& graph,
                         const std::unordered_map<int, int>& label_ids,
                         int num_labels) {
  const int n = static_cast<int>(graph.num_vertices());
  out.initVertex(n);
  out.initEdge(static_cast<number_type>(graph.num_edges()));
  out.nLabel = num_labels;
  out.maxDegree = 0;
  out.labelFrequency.clear();

  number_type offset = 0;
  for (int v = 0; v < n; ++v) {
    std::size_t degree = 0;
    const slim::VertexId* neighbors =
        graph.neighbors(static_cast<slim::VertexId>(v), degree);
    out.nbrOffset[v] = offset;
    for (std::size_t i = 0; i < degree; ++i) {
      out.nbr[offset + static_cast<number_type>(i)] = static_cast<int>(neighbors[i]);
    }
    offset += static_cast<number_type>(degree);
    out.degree[v] = static_cast<int>(degree);
    out.core[v] = static_cast<int>(degree);
    if (static_cast<int>(degree) > out.maxDegree) out.maxDegree = static_cast<int>(degree);

    const int label =
        label_ids.at(static_cast<int>(graph.label(static_cast<slim::VertexId>(v))));
    out.label[v] = label;
    out.labelFrequency[label] += 1;
  }
  out.nbrOffset[n] = offset;
}

inline void FreeVeqGraphStorage(Graph& graph) {
  delete[] graph.label;
  delete[] graph.nbr;
  delete[] graph.nbrOffset;
  delete[] graph.degree;
  delete[] graph.core;
  delete[] graph.vertex;
  delete[] graph.maxNbrDegree;
  delete[] graph.vertexOffset;
  delete[] graph.NLF;
  graph.label = nullptr;
  graph.nbr = nullptr;
  graph.nbrOffset = nullptr;
  graph.degree = nullptr;
  graph.core = nullptr;
  graph.vertex = nullptr;
  graph.maxNbrDegree = nullptr;
  graph.vertexOffset = nullptr;
  graph.NLF = nullptr;
}

inline void FreeVeqDataTables() {
  delete[] pairV;
  delete[] vCandIndex;
  delete[] flagSumNECCand;
  delete[] vCandInfo;
  if (sumNECCand != nullptr) {
    for (int i = 0; i < nUniqueLabel + 1; ++i) delete[] sumNECCand[i];
    delete[] sumNECCand;
  }
  delete[] nSumNECCand;
  delete[] mapTo;
  delete[] arrayToClean;
  delete[] isVisited;
  delete[] candToPos;
  delete[] answer;
  delete[] cntLabel;
  delete[] positiveLabel;
  pairV = nullptr;
  vCandIndex = nullptr;
  flagSumNECCand = nullptr;
  vCandInfo = nullptr;
  sumNECCand = nullptr;
  nSumNECCand = nullptr;
  mapTo = nullptr;
  arrayToClean = nullptr;
  isVisited = nullptr;
  candToPos = nullptr;
  answer = nullptr;
  cntLabel = nullptr;
  positiveLabel = nullptr;

  for (Graph* graph : dataGraph) {
    FreeVeqGraphStorage(*graph);
    delete graph;
  }
  dataGraph.clear();
}

inline void init_veq_data(const slim::LabeledGraph& graph,
                          const std::unordered_map<int, int>& label_ids) {
  labelMap.clear();
  for (const auto& entry : label_ids) {
    labelMap[std::to_string(entry.first)] = entry.second;
  }

  auto* data = new Graph();
  FillVeqGraph(*data, graph, label_ids, static_cast<int>(label_ids.size()));
  dataGraph.push_back(data);

  ProcessDataGraphs();
  AllocateForDataGraph();

  mpq_init(nMatch);
  mpq_init(nCurrMatch);
  mpq_init(nRemainingMatch);
  mpq_init(nMaxMatch);
  mpq_init(mpqOne);
  mpq_set_d(mpqOne, 1);
  mpq_init(mpqTmp);
  mpq_init(mpqZero);
  mpq_set_d(mpqZero, 0);
}

inline bool run_veq_query(Graph* query, const int* const* candidate_sets,
                          const int* candidate_sizes, std::uint64_t limit,
                          bool self_refinement, long long& embedding_cnt,
                          int64_t& time_total) {
  queryGraph.clear();
  queryGraph.push_back(query);

  injectedCandidates = candidate_sets;
  injectedCandidateSizes = candidate_sizes;
  veqSelfRefinement = self_refinement;

  timeLimit = 1000000000;

  if (limit == 0) {
    mpq_set_d(nMaxMatch, std::numeric_limits<double>::max());
  } else {
    const double requested = (limit > static_cast<std::uint64_t>(9223372036854775807ULL))
                                 ? static_cast<double>(9223372036854775807ULL)
                                 : static_cast<double>(limit);
    mpq_set_d(nMaxMatch, requested);
  }

  mpz_t counter;
  mpz_init(counter);
  const bool exceeded = ProcessQuery(counter, time_total);
  embedding_cnt = mpz_get_si(counter);
  return exceeded;
}

inline void dispose_veq_query(Graph* query) {
  if (query == nullptr) return;
  Deallocate(*query);
  FreeVeqGraphStorage(*query);
  delete query;
  queryGraph.clear();
}

inline void dispose_veq_data() {
  FreeVeqDataTables();
  labelMap.clear();
}
}

#endif
