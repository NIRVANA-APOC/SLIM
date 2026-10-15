#ifndef VEQ_PROCESS_H
#define VEQ_PROCESS_H

#include <chrono>
#include <cstring>
#include <iostream>
#include <tuple>
#include <vector>

#include "graph.h"
using namespace std;

namespace VEQ {
CandidateSpace candSpace[MAX_NUM_VERTEX];
#ifdef PRUNING_BY_EQUIVALENCE_SETS
int uI;
int nCandOffset;
int globalCellID = 1;
int* posArray;
int* candArray;
int* candOffset;
#endif

Graph* currG;
int nGraph;
int maxNumCandidate = 0;
int maxDegree = 0;
int maxNumDataVertex = 0;
vector<Graph*> dataGraph;
unordered_map<string, int> labelMap;

Graph* currQ;
int nQueryVertex = 0;
vector<Graph*> queryGraph;
int* mapTo;
bool isMapped[MAX_NUM_VERTEX];
#ifdef PRUNING_BY_EQUIVALENCE_SETS
int __order[MAX_NUM_VERTEX+1];
int* order = __order + 1;
#ifdef SUBGRAPH_MATCHING
int* nActive;
pair<int, int>** active;
#endif
#endif
int extendable[MAX_NUM_VERTEX];
int arrForQuery[MAX_NUM_VERTEX];
int binQuery[MAX_NUM_VERTEX];
int posQuery[MAX_NUM_VERTEX];
int vertQuery[MAX_NUM_VERTEX];
char localFlagQuery[MAX_NUM_VERTEX];
char visitedForQuery[MAX_NUM_VERTEX];
int visitedForQueryInt[MAX_NUM_VERTEX];
vector<pair<int, int>>
    NECElementOffset;
unordered_map<pair<int, int>, int, pairHash> NECChild;
unordered_map<pair<int, int>, int, pairHash> NECSize;

int root;
#ifdef NEIGHBOR_LABEL_FREQUENCY
pair<int, int>* NLFArray[MAX_NUM_VERTEX];
int numNLF[MAX_NUM_VERTEX];
#endif
int* DAG_parent_query[MAX_NUM_VERTEX];
int* DAG_child_query[MAX_NUM_VERTEX];
int* DAG_child_query_parent_index
    [MAX_NUM_VERTEX];
int DAG_parent_query_size[MAX_NUM_VERTEX];
int DAG_child_query_size[MAX_NUM_VERTEX];

int* DAG_ngb_query[MAX_NUM_VERTEX];
int* DAG_ngb_query_ngb_index[MAX_NUM_VERTEX];
int DAG_ngb_query_size[MAX_NUM_VERTEX];

uint64_t* DAG_ancestor[MAX_NUM_VERTEX];

chrono::high_resolution_clock::time_point GetClock() {
  return chrono::high_resolution_clock::now();
}
int64_t GetClockNan() {
  return chrono::high_resolution_clock::now().time_since_epoch().count();
}
double TimeDiffInMilliseconds(
    const chrono::high_resolution_clock::time_point start,
    const chrono::high_resolution_clock::time_point end) {
  chrono::duration<double> durationDiff = end - start;
  return durationDiff.count() * 1000;
}

inline void ProcessDataGraphs() {
  int start, num, maxValue;
  unordered_map<int, int>::const_iterator lit;
  nGraph = dataGraph.size();
  nUniqueLabel = labelMap.size();
#ifdef NEIGHBOR_LABEL_FREQUENCY
#else
  NLFSize = (BITS_PER_LABEL * nUniqueLabel - 1) / UINT64_SIZE + 1;
#endif
  for (int i = 0; i < dataGraph.size(); ++i) {
    Graph& g = *dataGraph[i];
    if (g.maxDegree > maxDegree) maxDegree = g.maxDegree;
    if (g.nVertex > maxNumDataVertex) maxNumDataVertex = g.nVertex;
  }
  positiveLabel = new int[nUniqueLabel];
  cntLabel = new int[nUniqueLabel]();
  int* bin = new int[maxDegree + 1];
  int* pos = new int[maxNumDataVertex];
  int* vert = new int[maxNumDataVertex];
  for (int i = 0; i < dataGraph.size(); ++i) {
    start = 0;
    maxValue = 0;
    memset(bin, 0, sizeof(int) * (maxDegree + 1));
    Graph& g = *dataGraph[i];
    for (lit = g.labelFrequency.begin(); lit != g.labelFrequency.end(); ++lit) {
      if (lit->second > maxValue) {
        g.mostFrequentLabelID = lit->first;
        maxValue = lit->second;
      }
    }
    g.maxNumSameLabelVertex = maxValue;
    if (maxValue > maxNumCandidate) maxNumCandidate = maxValue;
    for (int i = 0; i < g.nVertex; i++) bin[g.core[i]]++;

    for (int d = 0; d <= g.maxDegree; d++) {
      num = bin[d];
      bin[d] = start;
      start += num;
    }
    for (int i = 0; i < g.nVertex; i++) {
      pos[i] = bin[g.core[i]];
      vert[pos[i]] = i;
      bin[g.core[i]]++;
    }
    for (int d = g.maxDegree; d > 0; d--) bin[d] = bin[d - 1];
    bin[0] = 0;

    for (int i = 0; i < g.nVertex; i++) {
      int v = vert[i];
      for (number_type j = g.nbrOffset[v]; j < g.nbrOffset[v + 1]; j++) {
        int u = g.nbr[j];
        if (g.core[u] > g.core[v]) {
          int du = g.core[u];
          int pu = pos[u];
          int pw = bin[du];
          int w = vert[pw];
          if (u != w) {
            pos[u] = pw;
            pos[w] = pu;
            vert[pu] = w;
            vert[pw] = u;
          }
          bin[du]++;
          g.core[u]--;
        }
      }
    }
    g.sortNeighbors();
    g.computeNLF();
  }
  delete[] bin;
  delete[] pos;
  delete[] vert;
}

inline void GetQueryDataStructure(Graph& query) {
  int* pos = posQuery;
  int* vert = vertQuery;
  int* bin = binQuery;
  memset(bin, 0, sizeof(int) * (query.maxDegree + 1));
  for (int i = 0; i < nQueryVertex; i++) bin[query.core[i]]++;

  int num, start = 0;
  for (int d = 0; d <= query.maxDegree; d++) {
    num = bin[d];
    bin[d] = start;
    start += num;
  }
  for (int i = 0; i < nQueryVertex; i++) {
    pos[i] = bin[query.core[i]];
    vert[pos[i]] = i;
    bin[query.core[i]]++;
  }

  for (int d = query.maxDegree; d > 0; d--) bin[d] = bin[d - 1];
  bin[0] = 0;

  for (int i = 0; i < nQueryVertex; i++) {
    int v = vert[i];
    for (number_type j = query.nbrOffset[v]; j < query.nbrOffset[v + 1]; j++) {
      int u = query.nbr[j];
      if (query.core[u] > query.core[v]) {
        int du = query.core[u];
        int pu = pos[u];
        int pw = bin[du];
        int w = vert[pw];
        if (u != w) {
          pos[u] = pw;
          pos[w] = pu;
          vert[pu] = w;
          vert[pw] = u;
        }
        bin[du]++;
        query.core[u]--;
      }
    }
  }
  query.computeNLF();
}
}
#endif
