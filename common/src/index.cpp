

#include <slim/index.h>
#include <slim/graph.h>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <unordered_set>

namespace slim {

template<typename T>
void read_binary(std::ifstream& file, T& value) {
    file.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!file) {
        throw std::runtime_error("Failed to read from index file");
    }
}

size_t read_size(std::ifstream& file) {
    uint64_t size;
    read_binary(file, size);
    return static_cast<size_t>(size);
}

bool read_bool(std::ifstream& file) {
    uint8_t value;
    read_binary(file, value);
    return value != 0;
}

void read_vec_u32(std::ifstream& file, std::vector<uint32_t>& vec) {
    size_t size = read_size(file);
    vec.resize(size);
    if (size > 0) {
        file.read(reinterpret_cast<char*>(vec.data()), size * sizeof(uint32_t));
        if (!file) {
            throw std::runtime_error("Failed to read vector from index file");
        }
    }
}

void read_vec_f64(std::ifstream& file, std::vector<double>& vec) {
    size_t size = read_size(file);
    vec.resize(size);
    if (size > 0) {
        file.read(reinterpret_cast<char*>(vec.data()), size * sizeof(double));
        if (!file) {
            throw std::runtime_error("Failed to read vector from index file");
        }
    }
}

struct IndexLoader : public Index {
    IndexLoader() : Index() {}
};

std::unique_ptr<Index> Index::load(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open index file: " + path);
    }

    auto index = std::unique_ptr<Index>(new IndexLoader());

    read_binary(file, index->version_);
    if (index->version_ != INDEX_VERSION) {
        throw std::runtime_error("Unsupported index version: " +
                                std::to_string(index->version_) +
                                " (expected " + std::to_string(INDEX_VERSION) + ")");
    }

    read_binary(file, index->alpha_);

    size_t num_labels = read_size(file);
    index->label_encoder_.num_labels = num_labels;

    for (size_t i = 0; i < num_labels; ++i) {
        Label label;
        double value;
        read_binary(file, label);
        read_binary(file, value);
        index->label_encoder_.encoding[label] = value;
    }

    size_t num_signatures = read_size(file);
    index->signatures_.resize(num_signatures);

    for (size_t i = 0; i < num_signatures; ++i) {
        read_vec_f64(file, index->signatures_[i].eigenvalues);

        std::sort(index->signatures_[i].eigenvalues.begin(),
                  index->signatures_[i].eigenvalues.end());
    }

    std::vector<uint32_t> labels_u32;
    read_vec_u32(file, labels_u32);
    index->labels_.assign(labels_u32.begin(), labels_u32.end());

    size_t label_index_size = read_size(file);

    for (size_t i = 0; i < label_index_size; ++i) {
        Label label;
        read_binary(file, label);
        std::vector<uint32_t> vertices;
        read_vec_u32(file, vertices);
        index->label_index_[label] = std::move(vertices);
    }

    index->use_hierarchical_ = read_bool(file);

    if (index->use_hierarchical_) {
        try {
            index->hierarchical_index_ = HierarchicalIndex::deserialize(file);
            if (!index->hierarchical_index_) {
                index->use_hierarchical_ = false;
            }
        } catch (const std::exception& e) {
            index->hierarchical_index_ = nullptr;
            index->use_hierarchical_ = false;
        }
    }

    return index;
}

std::unique_ptr<HierarchicalIndex> HierarchicalIndex::deserialize(std::ifstream& file) {
    auto index = std::make_unique<HierarchicalIndex>();

    read_binary(file, index->tau_);

    read_binary(file, index->max_layer_);

    read_binary(file, index->max_dimension_);

    std::vector<uint32_t> layers_u32;
    read_vec_u32(file, layers_u32);
    index->layers_ = std::move(layers_u32);

    size_t hierarchical_size = read_size(file);

    for (size_t i = 0; i < hierarchical_size; ++i) {
        HierarchicalKey key;
        read_binary(file, key.label);
        read_binary(file, key.layer);
        read_binary(file, key.dimension);

        std::vector<uint32_t> vertices;
        read_vec_u32(file, vertices);

        index->hierarchical_index_[key] = std::move(vertices);
    }

    size_t label_max_layer_size = read_size(file);

    for (size_t i = 0; i < label_max_layer_size; ++i) {
        Label label;
        uint32_t max_layer;
        read_binary(file, label);
        read_binary(file, max_layer);
        index->label_max_layer_[label] = max_layer;
    }

    size_t label_layer_max_dim_size = read_size(file);

    for (size_t i = 0; i < label_layer_max_dim_size; ++i) {
        Label label;
        uint32_t layer;
        uint32_t max_dim;
        read_binary(file, label);
        read_binary(file, layer);
        read_binary(file, max_dim);
        index->label_layer_max_dim_[{label, layer}] = max_dim;
    }

    return index;
}

QueryTree QueryTree::build(const Graph& query_graph) {
    QueryTree tree;

    if (query_graph.numVertices() == 0) {
        return tree;
    }

    size_t num_vertices = query_graph.numVertices();
    VertexId root = 0;
    size_t max_degree = query_graph.degree(0);

    for (VertexId v = 1; v < num_vertices; ++v) {
        size_t deg = query_graph.degree(v);
        if (deg > max_degree) {
            max_degree = deg;
            root = v;
        }
    }

    tree.root_ = root;
    tree.order_.push_back(root);
    tree.parent_map_[root] = INVALID_VERTEX;

    std::queue<VertexId> queue;
    std::vector<bool> visited(num_vertices, false);
    queue.push(root);
    visited[root] = true;

    while (!queue.empty()) {
        VertexId u = queue.front();
        queue.pop();

        size_t count;
        const VertexId* neighbors = query_graph.neighbors(u, count);

        for (size_t i = 0; i < count; ++i) {
            VertexId v = neighbors[i];
            if (!visited[v]) {
                visited[v] = true;
                tree.order_.push_back(v);
                tree.parent_map_[v] = u;
                queue.push(v);
            }
        }
    }

    return tree;
}

QueryTree QueryTree::buildWithCandidates(const Graph& query_graph,
                                         const std::vector<size_t>& candidates_count) {
    QueryTree tree;

    size_t num_vertices = query_graph.numVertices();
    if (num_vertices == 0) {
        return tree;
    }

    std::vector<bool> visited(num_vertices, false);
    std::vector<bool> adjacent(num_vertices, false);

    VertexId root = 0;
    for (VertexId v = 1; v < num_vertices; ++v) {
        if (candidates_count[v] < candidates_count[root]) {
            root = v;
        } else if (candidates_count[v] == candidates_count[root] &&
                   query_graph.degree(v) > query_graph.degree(root)) {
            root = v;
        }
    }

    tree.root_ = root;
    tree.order_.push_back(root);
    tree.parent_map_[root] = INVALID_VERTEX;
    visited[root] = true;

    size_t nbr_count;
    const VertexId* neighbors = query_graph.neighbors(root, nbr_count);
    for (size_t i = 0; i < nbr_count; ++i) {
        adjacent[neighbors[i]] = true;
    }

    for (size_t i = 1; i < num_vertices; ++i) {
        VertexId next_vertex = INVALID_VERTEX;
        size_t min_cand = SIZE_MAX;

        for (VertexId v = 0; v < num_vertices; ++v) {
            if (!visited[v] && adjacent[v]) {
                if (candidates_count[v] < min_cand) {
                    min_cand = candidates_count[v];
                    next_vertex = v;
                } else if (candidates_count[v] == min_cand &&
                           next_vertex != INVALID_VERTEX &&
                           query_graph.degree(v) > query_graph.degree(next_vertex)) {
                    next_vertex = v;
                }
            }
        }

        if (next_vertex == INVALID_VERTEX) {

            for (VertexId v = 0; v < num_vertices; ++v) {
                if (!visited[v]) {
                    next_vertex = v;
                    break;
                }
            }
        }

        VertexId parent = INVALID_VERTEX;
        const VertexId* next_neighbors = query_graph.neighbors(next_vertex, nbr_count);
        for (size_t j = 0; j < nbr_count; ++j) {
            if (visited[next_neighbors[j]]) {
                parent = next_neighbors[j];
                break;
            }
        }

        tree.order_.push_back(next_vertex);
        tree.parent_map_[next_vertex] = parent;
        visited[next_vertex] = true;

        for (size_t j = 0; j < nbr_count; ++j) {
            adjacent[next_neighbors[j]] = true;
        }
    }

    return tree;
}

template<typename T>
void write_binary(std::ofstream& file, const T& value) {
    file.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!file) {
        throw std::runtime_error("Failed to write to index file");
    }
}

void write_size(std::ofstream& file, size_t size) {
    uint64_t size64 = static_cast<uint64_t>(size);
    write_binary(file, size64);
}

void write_bool(std::ofstream& file, bool value) {
    uint8_t byte = value ? 1 : 0;
    write_binary(file, byte);
}

void write_vec_u32(std::ofstream& file, const std::vector<uint32_t>& vec) {
    write_size(file, vec.size());
    if (!vec.empty()) {
        file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(uint32_t));
        if (!file) {
            throw std::runtime_error("Failed to write vector to index file");
        }
    }
}

void write_vec_f64(std::ofstream& file, const std::vector<double>& vec) {
    write_size(file, vec.size());
    if (!vec.empty()) {
        file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));
        if (!file) {
            throw std::runtime_error("Failed to write vector to index file");
        }
    }
}

std::unique_ptr<Index> Index::create(
    const LabelEncoder& label_encoder,
    double alpha,
    std::vector<SpectralSignature> signatures,
    std::vector<Label> labels,
    std::unordered_map<Label, std::vector<VertexId>> label_index,
    std::unique_ptr<HierarchicalIndex> hierarchical_index) {

    auto index = std::unique_ptr<Index>(new IndexLoader());
    index->version_ = INDEX_VERSION;
    index->label_encoder_ = label_encoder;
    index->alpha_ = alpha;
    index->signatures_ = std::move(signatures);
    index->labels_ = std::move(labels);
    index->label_index_ = std::move(label_index);
    index->hierarchical_index_ = std::move(hierarchical_index);
    index->use_hierarchical_ = (index->hierarchical_index_ != nullptr);
    return index;
}

void Index::save(const std::string& path) const {
    std::ofstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open index file for writing: " + path);
    }

    write_binary(file, version_);

    write_binary(file, alpha_);

    write_size(file, label_encoder_.num_labels);

    std::vector<Label> sorted_labels;
    sorted_labels.reserve(label_encoder_.encoding.size());
    for (const auto& [label, _] : label_encoder_.encoding) {
        sorted_labels.push_back(label);
    }
    std::sort(sorted_labels.begin(), sorted_labels.end());

    for (Label label : sorted_labels) {
        write_binary(file, label);
        write_binary(file, label_encoder_.encoding.at(label));
    }

    write_size(file, signatures_.size());
    for (const auto& sig : signatures_) {
        write_vec_f64(file, sig.eigenvalues);
    }

    write_size(file, labels_.size());
    file.write(reinterpret_cast<const char*>(labels_.data()), labels_.size() * sizeof(Label));
    if (!file) {
        throw std::runtime_error("Failed to write labels to index file");
    }

    write_size(file, label_index_.size());

    std::vector<Label> label_index_keys;
    label_index_keys.reserve(label_index_.size());
    for (const auto& [label, _] : label_index_) {
        label_index_keys.push_back(label);
    }
    std::sort(label_index_keys.begin(), label_index_keys.end());

    for (Label label : label_index_keys) {
        write_binary(file, label);
        const auto& vertices = label_index_.at(label);
        write_size(file, vertices.size());
        file.write(reinterpret_cast<const char*>(vertices.data()), vertices.size() * sizeof(VertexId));
        if (!file) {
            throw std::runtime_error("Failed to write label index to index file");
        }
    }

    write_bool(file, use_hierarchical_);

    if (use_hierarchical_ && hierarchical_index_) {
        hierarchical_index_->serialize(file);
    }
}

void HierarchicalIndex::serialize(std::ofstream& file) const {

    write_binary(file, tau_);

    write_binary(file, max_layer_);

    write_binary(file, max_dimension_);

    write_vec_u32(file, layers_);

    write_size(file, hierarchical_index_.size());

    std::vector<HierarchicalKey> sorted_keys;
    sorted_keys.reserve(hierarchical_index_.size());
    for (const auto& [key, _] : hierarchical_index_) {
        sorted_keys.push_back(key);
    }
    std::sort(sorted_keys.begin(), sorted_keys.end(),
        [](const HierarchicalKey& a, const HierarchicalKey& b) {
            if (a.label != b.label) return a.label < b.label;
            if (a.layer != b.layer) return a.layer < b.layer;
            return a.dimension < b.dimension;
        });

    for (const auto& key : sorted_keys) {
        write_binary(file, key.label);
        write_binary(file, key.layer);
        write_binary(file, key.dimension);

        const auto& vertices = hierarchical_index_.at(key);
        write_size(file, vertices.size());
        file.write(reinterpret_cast<const char*>(vertices.data()), vertices.size() * sizeof(VertexId));
        if (!file) {
            throw std::runtime_error("Failed to write hierarchical index entry");
        }
    }

    write_size(file, label_max_layer_.size());

    std::vector<Label> label_max_layer_keys;
    label_max_layer_keys.reserve(label_max_layer_.size());
    for (const auto& [label, _] : label_max_layer_) {
        label_max_layer_keys.push_back(label);
    }
    std::sort(label_max_layer_keys.begin(), label_max_layer_keys.end());

    for (Label label : label_max_layer_keys) {
        write_binary(file, label);
        write_binary(file, label_max_layer_.at(label));
    }

    write_size(file, label_layer_max_dim_.size());

    std::vector<std::pair<Label, uint32_t>> label_layer_keys;
    label_layer_keys.reserve(label_layer_max_dim_.size());
    for (const auto& [key, _] : label_layer_max_dim_) {
        label_layer_keys.push_back(key);
    }
    std::sort(label_layer_keys.begin(), label_layer_keys.end());

    for (const auto& [label, layer] : label_layer_keys) {
        write_binary(file, label);
        write_binary(file, layer);
        write_binary(file, label_layer_max_dim_.at({label, layer}));
    }
}

}
