

#include <slim/filter.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <unordered_set>

namespace slim
{

    CandidateSet SpectralFilter::filter(VertexId query_vertex,
                                        const Graph &query_graph) const
    {
        SpectralSignature query_sig = computeQuerySignature(query_vertex, query_graph);
        CandidateSet candidates;
        Label query_label = query_graph.label(query_vertex);

        if (index_.useHierarchical())
        {
            candidates = filterWithHierarchicalIndex(query_vertex, query_label, query_sig);
        }
        else
        {
            candidates = filterWithLabelIndex(query_vertex, query_label, query_sig);
        }

        return candidates;
    }

    struct CompactNLF
    {
        std::vector<std::pair<Label, uint32_t>> entries;

        void add(Label l)
        {
            for (auto &[label, count] : entries)
            {
                if (label == l)
                {
                    count++;
                    return;
                }
            }
            entries.push_back({l, 1});
        }

        void sort()
        {
            std::sort(entries.begin(), entries.end());
        }
    };

    static CompactNLF computeNLF(const Graph &graph, VertexId v)
    {
        CompactNLF nlf;
        size_t nbr_count;
        const VertexId *neighbors = graph.neighbors(v, nbr_count);

        if (nbr_count <= 16)
        {
            for (size_t i = 0; i < nbr_count; ++i)
            {
                nlf.add(graph.label(neighbors[i]));
            }
        }
        else
        {
            std::unordered_map<Label, uint32_t> temp;
            for (size_t i = 0; i < nbr_count; ++i)
            {
                temp[graph.label(neighbors[i])]++;
            }
            nlf.entries.reserve(temp.size());
            for (const auto &[l, c] : temp)
            {
                nlf.entries.push_back({l, c});
            }
        }
        nlf.sort();
        return nlf;
    }

    static bool nlfDominatedBy(const CompactNLF &query_nlf, const CompactNLF &data_nlf)
    {
        size_t qi = 0, di = 0;
        while (qi < query_nlf.entries.size())
        {
            while (di < data_nlf.entries.size() && data_nlf.entries[di].first < query_nlf.entries[qi].first)
            {
                di++;
            }

            if (di >= data_nlf.entries.size() ||
                data_nlf.entries[di].first != query_nlf.entries[qi].first ||
                data_nlf.entries[di].second < query_nlf.entries[qi].second)
            {
                return false;
            }

            qi++;
            di++;
        }
        return true;
    }

    static bool refineCandidate(
        VertexId u, VertexId v,
        const Graph &query_graph,
        const Graph &data_graph,
        const std::vector<std::unordered_set<VertexId>> &candidate_sets)
    {
        size_t u_nbr_count;
        const VertexId *u_neighbors = query_graph.neighbors(u, u_nbr_count);

        size_t v_nbr_count;
        const VertexId *v_neighbors = data_graph.neighbors(v, v_nbr_count);

        for (size_t i = 0; i < u_nbr_count; ++i)
        {
            VertexId u_prime = u_neighbors[i];
            const auto &u_prime_cands = candidate_sets[u_prime];

            bool found = false;
            for (size_t j = 0; j < v_nbr_count; ++j)
            {
                if (u_prime_cands.count(v_neighbors[j]) > 0)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                return false;
            }
        }

        return true;
    }

    static bool refineCandidateByLabelCoverage(
        VertexId u, VertexId v,
        const Graph &query_graph,
        const Graph &data_graph,
        const std::vector<std::unordered_set<VertexId>> &candidate_sets)
    {
        size_t u_nbr_count;
        const VertexId *u_neighbors = query_graph.neighbors(u, u_nbr_count);

        size_t v_nbr_count;
        const VertexId *v_neighbors = data_graph.neighbors(v, v_nbr_count);

        std::unordered_map<Label, uint32_t> required_by_label;
        for (size_t i = 0; i < u_nbr_count; ++i)
        {
            Label l = query_graph.label(u_neighbors[i]);
            required_by_label[l]++;
        }

        std::unordered_map<Label, uint32_t> available_by_label;
        for (size_t j = 0; j < v_nbr_count; ++j)
        {
            VertexId v_prime = v_neighbors[j];
            Label l = data_graph.label(v_prime);

            for (size_t i = 0; i < u_nbr_count; ++i)
            {
                VertexId u_prime = u_neighbors[i];
                if (query_graph.label(u_prime) == l && candidate_sets[u_prime].count(v_prime) > 0)
                {
                    available_by_label[l]++;
                    break;
                }
            }
        }

        for (const auto &[l, required] : required_by_label)
        {
            auto it = available_by_label.find(l);
            uint32_t available = (it != available_by_label.end()) ? it->second : 0;
            if (available < required)
            {
                return false;
            }
        }

        return true;
    }

    static void iterativeRefinement(
        const Graph &query_graph,
        const Graph &data_graph,
        std::vector<std::vector<VertexId>> &candidates,
        std::vector<std::unordered_set<VertexId>> &candidate_sets,
        int max_iterations = 5)
    {
        size_t num_vertices = query_graph.numVertices();

        for (int iter = 0; iter < max_iterations; ++iter)
        {
            bool changed = false;

            for (VertexId u = 0; u < num_vertices; ++u)
            {
                std::vector<VertexId> new_candidates;
                new_candidates.reserve(candidates[u].size());

                for (VertexId v : candidates[u])
                {

                    if (!refineCandidate(u, v, query_graph, data_graph, candidate_sets))
                    {
                        candidate_sets[u].erase(v);
                        changed = true;
                        continue;
                    }

                    if (iter >= 1 && !refineCandidateByLabelCoverage(u, v, query_graph, data_graph, candidate_sets))
                    {
                        candidate_sets[u].erase(v);
                        changed = true;
                        continue;
                    }

                    new_candidates.push_back(v);
                }

                candidates[u] = std::move(new_candidates);
            }

            if (!changed)
            {
                break;
            }
        }
    }

    CandidateSet SpectralFilter::filterAll(const Graph &query_graph,
                                           const std::vector<SpectralSignature> &query_signatures) const
    {
        CandidateSet all_candidates;

        size_t num_vertices = query_graph.numVertices();
        if (query_signatures.size() != num_vertices)
        {
            throw std::runtime_error("Query signatures size mismatch");
        }

        all_candidates.reserve(num_vertices);

        std::vector<std::vector<VertexId>> initial_candidates(num_vertices);

        for (VertexId u = 0; u < num_vertices; ++u)
        {
            const SpectralSignature &query_sig = query_signatures[u];
            Label query_label = query_graph.label(u);

            const auto &label_candidates = index_.getByLabel(query_label);
            const size_t query_dim = query_sig.dimension();
            const double *query_eig = query_sig.eigenvalues.data();
            const size_t query_degree = query_graph.degree(u);

            initial_candidates[u].reserve(label_candidates.size() / 2);

            for (VertexId v : label_candidates)
            {

                if (data_graph_.degree(v) < query_degree)
                {
                    continue;
                }

                const auto &data_sig = index_.signature(v);
                const size_t data_dim = data_sig.dimension();

                if (data_dim < query_dim)
                {
                    continue;
                }

                const size_t delta = data_dim - query_dim;
                const double *data_eig = data_sig.eigenvalues.data();
                bool dominated = true;
                for (size_t k = 0; k < query_dim; ++k)
                {
                    if (query_eig[k] > data_eig[k + delta] + 1e-6)
                    {
                        dominated = false;
                        break;
                    }
                }
                if (!dominated)
                    continue;

                initial_candidates[u].push_back(v);
            }
        }

        for (VertexId u = 0; u < num_vertices; ++u)
        {
            for (VertexId v : initial_candidates[u])
            {
                all_candidates.add(u, v);
            }
        }

        return all_candidates;
    }

    void SpectralFilter::refineWithNLF(const Graph &query_graph, CandidateSet &candidates) const
    {
        size_t num_vertices = query_graph.numVertices();

        std::vector<std::vector<VertexId>> cand_vec(num_vertices);
        std::vector<std::unordered_set<VertexId>> cand_set(num_vertices);

        for (VertexId u = 0; u < num_vertices; ++u)
        {
            const auto &cands = candidates.get(u);
            cand_vec[u].assign(cands.begin(), cands.end());
            cand_set[u].insert(cands.begin(), cands.end());
        }

        for (int iter = 0; iter < 3; ++iter)
        {
            bool changed = false;

            for (VertexId u = 0; u < num_vertices; ++u)
            {
                std::vector<VertexId> new_cands;
                new_cands.reserve(cand_vec[u].size());

                size_t u_nbr_count;
                const VertexId *u_nbrs = query_graph.neighbors(u, u_nbr_count);

                for (VertexId v : cand_vec[u])
                {
                    bool valid = true;

                    for (size_t i = 0; i < u_nbr_count && valid; ++i)
                    {
                        VertexId u_nbr = u_nbrs[i];
                        const auto &u_nbr_cands = cand_set[u_nbr];

                        size_t v_nbr_count;
                        const VertexId *v_nbrs = data_graph_.neighbors(v, v_nbr_count);

                        bool found = false;
                        for (size_t j = 0; j < v_nbr_count && !found; ++j)
                        {
                            if (u_nbr_cands.count(v_nbrs[j]) > 0)
                            {
                                found = true;
                            }
                        }

                        if (!found)
                        {
                            valid = false;
                        }
                    }

                    if (valid)
                    {
                        new_cands.push_back(v);
                    }
                    else
                    {
                        cand_set[u].erase(v);
                        changed = true;
                    }
                }

                cand_vec[u] = std::move(new_cands);
            }

            if (!changed)
                break;
        }

        candidates.clear();
        candidates.reserve(num_vertices);
        for (VertexId u = 0; u < num_vertices; ++u)
        {
            for (VertexId v : cand_vec[u])
            {
                candidates.add(u, v);
            }
        }
    }

    CandidateSet SpectralFilter::filterWithHierarchicalIndex(
        VertexId query_vertex,
        Label query_label,
        const SpectralSignature &query_sig) const
    {
        return filterWithLabelIndex(query_vertex, query_label, query_sig);
    }

    CandidateSet SpectralFilter::filterWithLabelIndex(
        VertexId query_vertex,
        Label query_label,
        const SpectralSignature &query_sig) const
    {
        CandidateSet candidates;
        const auto &label_candidates = index_.getByLabel(query_label);

        const size_t query_dim = query_sig.dimension();
        const double *query_eig = query_sig.eigenvalues.data();

        for (VertexId v : label_candidates)
        {
            const auto &data_sig = index_.signature(v);
            const size_t data_dim = data_sig.dimension();

            if (data_dim < query_dim)
                continue;

            const size_t delta = data_dim - query_dim;
            const double *data_eig = data_sig.eigenvalues.data();
            bool dominated = true;
            for (size_t k = 0; k < query_dim; ++k)
            {
                if (query_eig[k] > data_eig[k + delta] + 1e-6)
                {
                    dominated = false;
                    break;
                }
            }

            if (dominated)
            {
                candidates.add(query_vertex, v);
            }
        }

        return candidates;
    }

    SpectralSignature SpectralFilter::computeQuerySignature(
        VertexId query_vertex,
        const Graph &query_graph) const
    {
        LaplacianBuilder builder(index_.alpha(), index_.labelEncoder());
        return computeNeighborhoodSpectrum(query_graph, query_vertex, builder);
    }

}
