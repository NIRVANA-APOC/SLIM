#ifndef VEQ_RUN_H
#define VEQ_RUN_H

#include "VEQ.h"

namespace VEQ {
inline bool BuildCS() {
  if (!FilterByCount(*currQ, *currG)) return false;
  root = SelectRootVertex(*currQ, *currG);
  ExtractResidualStructures(*currQ);
  BuildDAG(*currQ);
  if (!InjectCandidateSets(*currQ)) return false;
#ifdef SIGMOD_DAF
  if (!BottomUp(*currQ, *currG)) return false;
  if (!TopDown(*currQ, *currG)) return false;
  ConstructAdjacencyList(*currQ, *currG, false);
#else
#ifdef SUBGRAPH_MATCHING
  if (veqSelfRefinement) {
    if (!filteringWithDAG(*currQ, *currG, bottomUp, childNgb, true))
      return false;
    if (!filteringWithDAG(*currQ, *currG, topDown, parentNgb, true))
      return false;
  }
  if (!ConstructAdjacencyList(*currQ, *currG, veqSelfRefinement)) return false;
#else
  if (veqSelfRefinement &&
      !filteringWithDAG(*currQ, *currG, bottomUp, childNgb, false))
    return false;
  if (veqSelfRefinement &&
      !filteringWithDAG(*currQ, *currG, topDown, parentNgb, true))
    return false;
  if (!ConstructAdjacencyList(*currQ, *currG, veqSelfRefinement)) return false;
#endif
#endif
  return true;
}

inline bool ProcessQuery(mpz_t embedding_cnt, int64_t& time_total) {
  bool isCandidate;
  chrono::high_resolution_clock::time_point startClock, endClock;
  for (int x = 0; x < queryGraph.size(); ++x) {
    exceedTimeLimit = false;
    startPoint = GetClock();
    endTime = GetClockNan();
    endTime += timeLimit * 1000 * 1000;
    filteringTime = 0;
    verificationTime = 0;
    mpq_set_d(nMatch, 0);
    nCandidate = 0;
    recursiveCallCount = 0;
    currQ = queryGraph[x];
    AllocateForQueryGraph(*currQ);
    GetQueryDataStructure(*currQ);
    isTree = true;
    for (int i = 0; i < nQueryVertex; i++) {
      if (currQ->core[i] >= 2) {
        isTree = false;
        break;
      }
    }
#ifdef MULTIPLE_QUERIES
    cout << "*Processing q" << x << "(|V|: " << currQ->nVertex << ")" << endl;
#endif
    for (int y = 0; y < dataGraph.size(); ++y) {
      SetQueryGraphResource(*currQ);
      currG = dataGraph[y];
#ifdef MAPPING_FUNCTION_LOG
      cout << "Processing g" << y << "(|V|: " << currG->nVertex << ")" << endl;
#endif

      startClock = GetClock();
      isCandidate = BuildCS();
      endClock = GetClock();
      filteringTime = std::chrono::duration_cast<std::chrono::nanoseconds>(endClock - startClock).count();
      if (!isCandidate) continue;
      ++nCandidate;
#ifdef SUBGRAPH_MATCHING
      ComputeAuxiliaryDataStructureSize();
#endif
#ifdef FILTERING_ONLY
#else
      startClock = GetClock();
#ifdef SIGMOD_DAF
      if (usePathSizeOrder) InitializeWeight(*currQ);
#endif
      Backtrack(*currQ, *currG, y);
      endClock = GetClock();
      verificationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(endClock - startClock).count();
#endif
#ifdef SPLIT_RECURSIVE_CALLS
      cout << "Recursive calls(g" << y << "): " << recursiveCallCountPerGraph
           << endl;
#endif
    }
#ifdef MULTIPLE_QUERIES
    cout << "Exceeded Time Limit: " << (exceedTimeLimit ? 1 : 0) << endl;
#endif
    mpq_get_num(embedding_cnt, nMatch);
    endPoint = GetClock();
    time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(endPoint - startPoint).count();
    return exceedTimeLimit;
  }
}
}
#endif
