#ifndef SLIM_CANDIDATE_DUMP_H
#define SLIM_CANDIDATE_DUMP_H

#include <slim/spectral_filter.h>

#include <string>

namespace slim {
void write_candidate_dump(const std::string& path, const CandidateSets& candidates);
}

#endif
