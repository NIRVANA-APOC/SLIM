#ifndef SLIM_BACKTRACKING_VERIFIER_H
#define SLIM_BACKTRACKING_VERIFIER_H

#include <slim/labeled_graph.h>
#include <slim/spectral_filter.h>
#include <slim/types.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace slim {
constexpr std::size_t MAX_QUERY_VERTICES = 256;

enum class VerifierBackend {
    Ive,
    Veq,
    Lftj,
};

class EnumerationBackend {
public:
    struct Result {
        std::size_t count = 0;
        bool early_stopped = false;
        std::uint64_t adjacency_us = 0;
        std::uint64_t enumeration_us = 0;
    };

    virtual ~EnumerationBackend() = default;

    EnumerationBackend(const EnumerationBackend&) = delete;
    EnumerationBackend& operator=(const EnumerationBackend&) = delete;

    virtual Result verify(const LabeledGraph& query_graph,
                          std::vector<std::vector<VertexId>>& candidates,
                          std::size_t limit) = 0;

protected:
    EnumerationBackend() = default;
};

class BacktrackingVerifier {
public:
    struct Result {
        std::size_t count = 0;
        bool early_stopped = false;
        std::size_t candidates_after_propagation = 0;
        std::uint64_t propagation_us = 0;
        std::uint64_t adjacency_us = 0;
        std::uint64_t enumeration_us = 0;
    };

    explicit BacktrackingVerifier(const LabeledGraph& data_graph,
                                  VerifierBackend backend = VerifierBackend::Ive);
    ~BacktrackingVerifier();

    BacktrackingVerifier(const BacktrackingVerifier&) = delete;
    BacktrackingVerifier& operator=(const BacktrackingVerifier&) = delete;

    Result verify(const LabeledGraph& query_graph,
                  CandidateSets& candidates,
                  std::size_t limit);

private:
    void propagate_edge_constraints(const LabeledGraph& query_graph,
                                    std::vector<std::vector<VertexId>>& candidates);

    const LabeledGraph& data_graph_;
    std::vector<std::uint64_t> support_;
    std::vector<char> support_fresh_;

    std::unique_ptr<EnumerationBackend> engine_;
};
}

#endif
