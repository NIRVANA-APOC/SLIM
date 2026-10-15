#ifndef VEQ_CONFIG_H
#define VEQ_CONFIG_H

namespace VEQ {
#define DYNAMIC_ORDERING
#define FILTERING_BY_NEIGHBOR_SAFETY
#define LEAF_ADAPTIVE_MATCHING
#define N_OPTIMIZATION
#define PRUNING_BY_EQUIVALENCE_SETS
#define SUBGRAPH_MATCHING

typedef long long weight_type;
typedef unsigned short degree_type;
#ifdef HUGE_GRAPH
typedef long long number_type;
#else
typedef int number_type;
#endif
#ifdef SIGMOD_DAF
const int BITS_PER_LABEL =
    1;
#else
#ifdef HUGE_GRAPH
const int BITS_PER_LABEL = 1;
#else
const int BITS_PER_LABEL =
    4;
#endif
#endif
const int Start = 0;
const int LoadNumVertex = 1;
const int LoadVertex = 2;
const int LoadNumEdge = 3;
const int LoadEdge = 4;
const int GFU_FORMAT = 5;
const int IGRAPH_FORMAT = 6;
const int CFL_FORMAT = 7;
const int MAX_NUM_VERTEX = 257;
const int MAX_QUERY_DEGREE = 100;
const int MAX_NUM_LABEL = 48000;
const int UINT64_SIZE = sizeof(uint64_t) * 8;
#ifdef SIGMOD_DAF
const bool usePathSizeOrder = true;
#else
#if defined(QUERY_VERTEX_WITH_MAX) || defined(CANDIDATE_VERTEX_WITH_MAX)
const bool usePathSizeOrder = true;
#else
const bool usePathSizeOrder = false;
#endif
#endif
const bool useFailingSet = true;
bool isTree = false;

int NLFSize = -1;
int nUniqueLabel = 0;
}
#endif
