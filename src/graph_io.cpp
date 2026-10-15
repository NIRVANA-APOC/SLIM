#include <slim/graph_io.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace slim {
namespace {
struct LineTokenizer {
    const char* p;
    const char* end;

    static bool is_space(char c) {
        return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v';
    }

    bool done() const { return p >= end; }

    bool next_line() {
        while (p < end) {
            while (p < end && *p != '\n' && is_space(*p)) {
                ++p;
            }
            if (p >= end) {
                return false;
            }
            if (*p == '\n') {
                ++p;
                continue;
            }
            if (*p == '#') {
                while (p < end && *p != '\n') {
                    ++p;
                }
                continue;
            }
            return true;
        }
        return false;
    }

    bool next(const char*& begin, const char*& finish) {
        while (p < end && *p != '\n' && is_space(*p)) {
            ++p;
        }
        if (p >= end || *p == '\n') {
            return false;
        }
        begin = p;
        while (p < end && !is_space(*p)) {
            ++p;
        }
        finish = p;
        return true;
    }

    void skip_line() {
        while (p < end && *p != '\n') {
            ++p;
        }
    }

    static bool parse_u64(const char* begin, const char* finish, std::uint64_t& value) {
        if (begin == finish) {
            return false;
        }
        if (*begin == 'L' || *begin == 'l') {
            ++begin;
        }
        if (begin == finish) {
            return false;
        }
        std::uint64_t acc = 0;
        for (const char* q = begin; q != finish; ++q) {
            const unsigned char c = static_cast<unsigned char>(*q);
            if (c < '0' || c > '9') {
                return false;
            }
            acc = acc * 10 + static_cast<std::uint64_t>(c - '0');
        }
        value = acc;
        return true;
    }

    static bool parse_vertex(const char* begin, const char* finish, VertexId& vertex) {
        std::uint64_t value = 0;
        if (!parse_u64(begin, finish, value) ||
            value > static_cast<std::uint64_t>(static_cast<VertexId>(-1))) {
            return false;
        }
        vertex = static_cast<VertexId>(value);
        return true;
    }

    static bool parse_label(const char* begin, const char* finish, Label& label) {
        std::uint64_t value = 0;
        if (!parse_u64(begin, finish, value) ||
            value > static_cast<std::uint64_t>(static_cast<Label>(-1))) {
            return false;
        }
        label = static_cast<Label>(value);
        return true;
    }
};

void assign_label(std::vector<Label>& labels, std::vector<bool>& set,
                  VertexId v, Label label) {
    if (v >= labels.size()) {
        labels.resize(static_cast<std::size_t>(v) + 1, 0);
        set.resize(static_cast<std::size_t>(v) + 1, false);
    }
    labels[v] = label;
    set[v] = true;
}

std::string read_whole_file(const std::string& path) {
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr) {
        throw std::runtime_error("cannot open graph file: " + path);
    }
    std::string buffer;
    char chunk[1 << 20];
    std::size_t got = 0;
    while ((got = std::fread(chunk, 1, sizeof(chunk), file)) > 0) {
        buffer.append(chunk, got);
    }
    const bool failed = std::ferror(file) != 0;
    std::fclose(file);
    if (failed) {
        throw std::runtime_error("failed to read graph file: " + path);
    }
    return buffer;
}
}

LabeledGraph read_graph(const std::string& path) {
    const std::string buffer = read_whole_file(path);

    std::vector<std::pair<VertexId, VertexId>> edges;
    std::vector<Label> vertex_labels;
    std::vector<bool> label_set;
    VertexId max_vertex = 0;
    std::size_t declared_vertices = 0;
    bool saw_vertex_record = false;

    LineTokenizer tokenizer{buffer.data(), buffer.data() + buffer.size()};
    const char* tok[4] = {nullptr, nullptr, nullptr, nullptr};
    const char* tok_end[4] = {nullptr, nullptr, nullptr, nullptr};

    while (tokenizer.next_line()) {
        std::size_t count = 0;
        while (count < 4 && tokenizer.next(tok[count], tok_end[count])) {
            ++count;
        }

        if (count == 0) {
            continue;
        }
        const std::size_t first_len = static_cast<std::size_t>(tok_end[0] - tok[0]);
        if (first_len != 1) {
        } else if (*tok[0] == 't') {
            if (count >= 3) {
                std::uint64_t first = 0;
                std::uint64_t second = 0;
                const bool first_numeric = LineTokenizer::parse_u64(tok[1], tok_end[1], first);
                const bool second_numeric = LineTokenizer::parse_u64(tok[2], tok_end[2], second);
                if (second_numeric) {
                    declared_vertices = static_cast<std::size_t>(first_numeric ? first : second);
                } else if (first_numeric) {
                    declared_vertices = static_cast<std::size_t>(first);
                }
            } else if (count == 2) {
                std::uint64_t only = 0;
                if (LineTokenizer::parse_u64(tok[1], tok_end[1], only)) {
                    declared_vertices = static_cast<std::size_t>(only);
                }
            }
            tokenizer.skip_line();
            continue;
        } else if (*tok[0] == 'v') {
            if (count < 3) {
                throw std::runtime_error("malformed vertex line in " + path);
            }
            VertexId v = 0;
            Label label = 0;
            if (!LineTokenizer::parse_vertex(tok[1], tok_end[1], v) ||
                !LineTokenizer::parse_label(tok[2], tok_end[2], label)) {
                throw std::runtime_error("malformed vertex line in " + path);
            }
            assign_label(vertex_labels, label_set, v, label);
            if (v > max_vertex) max_vertex = v;
            saw_vertex_record = true;
            tokenizer.skip_line();
            continue;
        } else if (*tok[0] == 'e') {
            if (count < 3) {
                throw std::runtime_error("malformed edge line in " + path);
            }
            VertexId u = 0;
            VertexId v = 0;
            if (!LineTokenizer::parse_vertex(tok[1], tok_end[1], u) ||
                !LineTokenizer::parse_vertex(tok[2], tok_end[2], v)) {
                throw std::runtime_error("malformed edge line in " + path);
            }
            edges.push_back({u, v});
            if (u > max_vertex) max_vertex = u;
            if (v > max_vertex) max_vertex = v;
            tokenizer.skip_line();
            continue;
        }

        if (count < 2) {
            continue;
        }
        VertexId u = 0;
        VertexId v = 0;
        if (!LineTokenizer::parse_vertex(tok[0], tok_end[0], u) ||
            !LineTokenizer::parse_vertex(tok[1], tok_end[1], v)) {
            continue;
        }
        Label lu = u;
        Label lv = v;
        if (count >= 4) {
            if (!LineTokenizer::parse_label(tok[2], tok_end[2], lu) ||
                !LineTokenizer::parse_label(tok[3], tok_end[3], lv)) {
                throw std::runtime_error("malformed edge label in " + path);
            }
        }
        edges.push_back({u, v});
        assign_label(vertex_labels, label_set, u, lu);
        assign_label(vertex_labels, label_set, v, lv);
        if (u > max_vertex) max_vertex = u;
        if (v > max_vertex) max_vertex = v;
    }

    const std::size_t num_vertices = std::max<std::size_t>(
        static_cast<std::size_t>(max_vertex) + 1, declared_vertices);
    vertex_labels.resize(num_vertices, 0);
    label_set.resize(num_vertices, false);

    for (VertexId v = 0; v < num_vertices; ++v) {
        if (!label_set[v]) {
            vertex_labels[v] = saw_vertex_record ? 0 : v;
        }
    }

    return LabeledGraph(num_vertices, edges, vertex_labels);
}
}
