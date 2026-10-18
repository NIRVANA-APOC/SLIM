#include <slim/spectral_dominance_index.h>

#include <slim/rank_selection.h>
#include <slim/graph_fingerprint.h>
#include <slim/spectral_operator.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace slim {
namespace {
struct BucketRecord {
    Label label;
    std::uint32_t layer;
    std::uint32_t dimension;
    VertexId vertex;
};

constexpr std::size_t IO_BUFFER = 1u << 20;

constexpr char INDEX_IDENTITY_MAGIC[8] = {'S', 'L', 'I', 'M', 'I', 'D', 'X', '1'};

constexpr char INDEX_AMSD_MAGIC[8] = {'S', 'L', 'I', 'M', 'A', 'M', 'S', 'D'};

void write_rank_signals_dump(const std::string& path, double tau, std::uint32_t max_rank,
                             const std::vector<LocalSpectralSignature>& signatures) {
    std::ofstream dump(path);
    if (!dump) {
        throw std::runtime_error("cannot open rank-signals dump '" + path + "'");
    }
    double max_lambda = 0.0;
    for (const LocalSpectralSignature& signature : signatures) {
        if (signature.dimension() > 0) {
            max_lambda = std::max(max_lambda, stored_spectral_value(signature.lambda_max()));
        }
    }
    std::uint32_t cols = 1;
    if (tau > 0.0 && max_lambda > 0.0) {
        constexpr double kMaxCols = 1u << 21;
        cols = static_cast<std::uint32_t>(std::min(std::floor(max_lambda / tau) + 1.0, kMaxCols));
        cols = std::max<std::uint32_t>(cols, 1);
    }
    dump << "# tau " << std::setprecision(17) << tau << " num_vertices " << signatures.size()
         << " max_rank " << max_rank << " max_lambda " << std::setprecision(17) << max_lambda
         << " cols " << cols << "\n";
    for (std::size_t v = 0; v < signatures.size(); ++v) {
        const std::uint32_t dim = static_cast<std::uint32_t>(signatures[v].dimension());
        dump << v << ' ' << dim;
        for (std::uint32_t rank = 1; rank <= max_rank; ++rank) {
            const double value =
                rank <= dim ? stored_spectral_value(signatures[v].eigenvalues[dim - rank]) : 0.0;
            dump << ' ' << std::setprecision(17) << value;
        }
        dump << '\n';
    }
}

void write_rank_report(const std::string& path, const RankSelectionResult& selection) {
    std::ofstream report(path);
    if (!report) {
        throw std::runtime_error("cannot open rank report '" + path + "'");
    }
    report << "mode " << rank_select_mode_name(selection.mode) << '\n';
    report << "tau " << std::setprecision(17) << selection.tau << '\n';
    report << "num_vertices " << selection.num_vertices << '\n';
    report << "budget " << std::setprecision(17) << selection.budget_entries << '\n';
    report << "max_rank " << selection.max_rank << '\n';
    report << "space " << selection.space_entries << '\n';
    const auto at64 = [&selection](const std::vector<std::uint64_t>& values, std::uint32_t rank) {
        return rank - 1 < values.size() ? values[rank - 1] : 0;
    };
    const auto atd = [&selection](const std::vector<double>& values, std::uint32_t rank) {
        return rank - 1 < values.size() ? values[rank - 1] : 0.0;
    };
    for (std::uint32_t rank = 1; rank <= selection.max_rank; ++rank) {
        report << "rank " << rank << ' ' << at64(selection.c_j, rank) << ' '
               << std::setprecision(17) << atd(selection.h_j, rank) << ' '
               << atd(selection.d_j, rank) << '\n';
    }
    report << "selected";
    for (std::uint32_t rank : selection.selected_ranks) {
        report << ' ' << rank;
    }
    report << '\n';
    report << "order";
    for (std::uint32_t rank : selection.selection_order) {
        report << ' ' << rank;
    }
    report << '\n';
    for (const RankNmiRecord& record : selection.nmi_calls) {
        report << "nmi " << record.rank_a << ' ' << record.rank_b << ' '
               << std::setprecision(17) << record.value << '\n';
    }
}

RankSelectionOptions make_rank_selection_options(const SlimConfig& config, double tau,
                                                 std::size_t num_vertices,
                                                 std::uint32_t max_signature_dim) {
    RankSelectionOptions options;
    options.mode = config.amsd_rank_select;
    options.tau = tau;
    options.num_vertices = num_vertices;
    options.budget_entries = config.amsd_budget_ratio * static_cast<double>(num_vertices);
    options.max_rank = std::min(static_cast<std::uint32_t>(std::max(0, config.amsd_ranks)),
                                max_signature_dim);
    for (std::uint32_t rank = 1; rank <= options.max_rank; ++rank) {
        options.fallback_ranks.push_back(rank);
    }
    return options;
}

constexpr std::uint32_t AMSD_SECTION_VERSION = 2;

class BinaryWriter {
public:
    explicit BinaryWriter(const std::string& path)
        : out_(path, std::ios::binary), buffer_(IO_BUFFER) {
        if (!out_.is_open()) {
            throw std::runtime_error("cannot open index file for writing: " + path);
        }
    }

    ~BinaryWriter() { flush(); }

    void bytes(const char* data, std::size_t count) {
        if (count >= buffer_.size() || used_ + count > buffer_.size()) {
            flush();
            if (count >= buffer_.size()) {
                out_.write(data, static_cast<std::streamsize>(count));
                if (!out_) throw std::runtime_error("failed to write index file");
                return;
            }
        }
        std::memcpy(buffer_.data() + used_, data, count);
        used_ += count;
    }

    template <typename T>
    void pod(const T& value) {
        bytes(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    void varint(std::uint32_t value) {
        char scratch[5];
        std::size_t length = 0;
        do {
            std::uint8_t byte = static_cast<std::uint8_t>(value & 0x7Fu);
            value >>= 7;
            scratch[length++] = static_cast<char>(value != 0 ? (byte | 0x80u) : byte);
        } while (value != 0);
        bytes(scratch, length);
    }

    void flush() {
        if (used_ != 0) {
            out_.write(buffer_.data(), static_cast<std::streamsize>(used_));
            if (!out_) throw std::runtime_error("failed to write index file");
            used_ = 0;
        }
    }

private:
    std::ofstream out_;
    std::vector<char> buffer_;
    std::size_t used_ = 0;
};

class BinaryReader {
public:
    explicit BinaryReader(const std::string& path)
        : in_(path, std::ios::binary), buffer_(IO_BUFFER) {
        if (!in_.is_open()) {
            throw std::runtime_error("cannot open index file: " + path);
        }
    }

    void bytes(char* data, std::size_t count) {
        std::size_t done = 0;
        while (done < count) {
            if (used_ == filled_) {
                refill();
            }
            const std::size_t chunk = std::min(count - done, filled_ - used_);
            std::memcpy(data + done, buffer_.data() + used_, chunk);
            used_ += chunk;
            done += chunk;
        }
    }

    template <typename T>
    T pod() {
        T value{};
        bytes(reinterpret_cast<char*>(&value), sizeof(T));
        return value;
    }

    std::uint32_t varint() {
        std::uint32_t value = 0;
        int shift = 0;
        for (;;) {
            const std::uint8_t byte = pod<std::uint8_t>();
            value |= static_cast<std::uint32_t>(byte & 0x7Fu) << shift;
            if ((byte & 0x80u) == 0) {
                return value;
            }
            shift += 7;
            if (shift > 28) {
                throw std::runtime_error("corrupt index file: bad varint");
            }
        }
    }

private:
    void refill() {
        used_ = 0;
        in_.read(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
        filled_ = static_cast<std::size_t>(in_.gcount());
        if (filled_ == 0) {
            throw std::runtime_error("failed to read index file");
        }
    }

    std::ifstream in_;
    std::vector<char> buffer_;
    std::size_t used_ = 0;
    std::size_t filled_ = 0;
};

void write_dimension(BinaryWriter& out, std::uint32_t dim) {
    if (dim < 0xFFu) {
        out.pod(static_cast<std::uint8_t>(dim));
        return;
    }
    out.pod(static_cast<std::uint8_t>(0xFFu));
    out.pod(dim);
}

std::uint32_t read_dimension(BinaryReader& in) {
    const std::uint8_t code = in.pod<std::uint8_t>();
    return code < 0xFFu ? code : in.pod<std::uint32_t>();
}

float ceil_to_float(double value) {
    return static_cast<float>(stored_spectral_value(value));
}

void write_signature(BinaryWriter& out, const LocalSpectralSignature& signature,
                     std::vector<float>& payload) {
    const std::size_t dim = signature.eigenvalues.size();
    write_dimension(out, static_cast<std::uint32_t>(dim));
    payload.resize(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        payload[i] = ceil_to_float(signature.eigenvalues[i]);
    }
    out.bytes(reinterpret_cast<const char*>(payload.data()), dim * sizeof(float));
}

void write_vertices(BinaryWriter& out, const VertexId* vertices, std::size_t count) {
    std::uint32_t previous = 0;
    for (std::size_t i = 0; i < count; ++i) {
        out.varint(vertices[i] - previous);
        previous = vertices[i];
    }
}

void read_vertices(BinaryReader& in, VertexId* vertices, std::size_t count) {
    std::uint32_t previous = 0;
    for (std::size_t i = 0; i < count; ++i) {
        previous += in.varint();
        vertices[i] = previous;
    }
}
}

const SpectralDominanceIndex::LabelGroup* SpectralDominanceIndex::find_group(Label label) const {
    const auto it = std::lower_bound(
        groups_.begin(), groups_.end(), label,
        [](const LabelGroup& group, Label value) { return group.label < value; });
    if (it == groups_.end() || it->label != label) {
        return nullptr;
    }
    return &*it;
}

std::size_t SpectralDominanceIndex::group_position(Label label) const {
    const auto it = std::lower_bound(
        groups_.begin(), groups_.end(), label,
        [](const LabelGroup& group, Label value) { return group.label < value; });
    if (it == groups_.end() || it->label != label) {
        return groups_.size();
    }
    return static_cast<std::size_t>(it - groups_.begin());
}

const SpectralDominanceIndex::EmptyGroup* SpectralDominanceIndex::find_empty_group(Label label) const {
    const auto it = std::lower_bound(
        empty_groups_.begin(), empty_groups_.end(), label,
        [](const EmptyGroup& group, Label value) { return group.label < value; });
    if (it == empty_groups_.end() || it->label != label) {
        return nullptr;
    }
    return &*it;
}

void SpectralDominanceIndex::save(const std::string& path) const {
    BinaryWriter out(path);

    out.pod(alpha_);
    out.pod(static_cast<std::int32_t>(k_));
    out.pod(tau_);
    out.pod(static_cast<std::uint32_t>(signatures_.size()));

    std::vector<float> payload;
    for (const LocalSpectralSignature& signature : signatures_) {
        write_signature(out, signature, payload);
    }

    out.pod(static_cast<std::uint32_t>(groups_.size()));
    for (const LabelGroup& group : groups_) {
        out.pod(group.label);
        out.pod(static_cast<std::uint32_t>(group.buckets.size()));
        for (const Bucket& bucket : group.buckets) {
            out.varint(bucket.layer);
            out.varint(bucket.dimension);
            out.pod(static_cast<std::uint32_t>(bucket.count));
            write_vertices(out, group.vertices.data() + static_cast<std::size_t>(bucket.begin),
                           static_cast<std::size_t>(bucket.count));
        }
    }

    out.pod(static_cast<std::uint32_t>(empty_groups_.size()));
    for (const EmptyGroup& group : empty_groups_) {
        out.pod(group.label);
        out.pod(static_cast<std::uint32_t>(group.vertices.size()));
        write_vertices(out, group.vertices.data(), group.vertices.size());
    }

    out.pod(num_labels_);
    out.pod(static_cast<std::uint8_t>(has_ball_labels_ ? 1u : 0u));
    if (has_ball_labels_) {
        std::vector<Label> alphabet = alphabet_;
        out.bytes(reinterpret_cast<const char*>(alphabet.data()),
                  alphabet.size() * sizeof(Label));
        out.bytes(reinterpret_cast<const char*>(ball_counts_.data()), ball_counts_.size());
    }

    out.bytes(INDEX_IDENTITY_MAGIC, sizeof(INDEX_IDENTITY_MAGIC));
    out.pod(identity_.format_version);
    out.pod(identity_.data_fingerprint);
    out.pod(identity_.label_fingerprint);
    out.pod(identity_.neighborhood_size_limit);
    out.pod(identity_.data_vertices);
    out.pod(identity_.data_edges);
    out.pod(identity_.tau_scale);

    if (amsd_.present) {
        out.bytes(INDEX_AMSD_MAGIC, sizeof(INDEX_AMSD_MAGIC));
        out.pod(AMSD_SECTION_VERSION);
        out.pod(amsd_.configured_ranks);
        out.pod(amsd_.configured_bins);
        out.pod(amsd_.rank_select_mode);
        out.pod(amsd_.budget_ratio);
        out.pod(amsd_.budget_entries);
        out.pod(amsd_.space_entries);
        out.pod(static_cast<std::uint32_t>(amsd_.selected_ranks.size()));
        for (std::uint32_t rank : amsd_.selected_ranks) {
            out.pod(rank);
        }
        out.pod(static_cast<std::uint32_t>(amsd_.selection_order.size()));
        for (std::uint32_t rank : amsd_.selection_order) {
            out.pod(rank);
        }
        out.pod(static_cast<std::uint32_t>(amsd_.degrees.size()));
        out.bytes(reinterpret_cast<const char*>(amsd_.degrees.data()),
                  amsd_.degrees.size() * sizeof(std::uint32_t));
        out.pod(static_cast<std::uint32_t>(amsd_.groups.size()));
        for (const AmsdGroup& group : amsd_.groups) {
            out.pod(group.label);
            out.pod(group.size);
            out.pod(static_cast<std::uint32_t>(group.ranks.size()));
            for (const AmsdRankIndex& rank : group.ranks) {
                out.pod(rank.rank);
                out.pod(static_cast<std::uint32_t>(rank.rungs.size()));
                for (double rung : rank.rungs) {
                    out.pod(rung);
                }
                for (std::uint32_t cardinality : rank.suffix_sets) {
                    out.pod(cardinality);
                }
                const std::size_t words = (static_cast<std::size_t>(group.size) + 63) / 64;
                for (const BitVector& bits : rank.suffix) {
                    out.bytes(reinterpret_cast<const char*>(bits.word_data()),
                              words * sizeof(std::uint64_t));
                }
            }
        }
    }
    out.flush();
}

std::unique_ptr<SpectralDominanceIndex> SpectralDominanceIndex::load(const std::string& path) {
    BinaryReader in(path);

    auto index = std::unique_ptr<SpectralDominanceIndex>(new SpectralDominanceIndex());
    index->alpha_ = in.pod<double>();
    index->k_ = in.pod<std::int32_t>();
    index->tau_ = in.pod<double>();

    const std::uint32_t num_vertices = in.pod<std::uint32_t>();
    index->signatures_.resize(static_cast<std::size_t>(num_vertices));
    std::vector<float> payload;
    for (std::size_t v = 0; v < num_vertices; ++v) {
        const std::uint32_t dim = read_dimension(in);
        LocalSpectralSignature& signature = index->signatures_[v];
        signature.eigenvalues.resize(dim);
        payload.resize(dim);
        in.bytes(reinterpret_cast<char*>(payload.data()), dim * sizeof(float));
        for (std::size_t i = 0; i < dim; ++i) {
            signature.eigenvalues[i] = static_cast<double>(payload[i]);
        }
    }

    const std::uint32_t num_groups = in.pod<std::uint32_t>();
    index->groups_.resize(static_cast<std::size_t>(num_groups));
    for (std::size_t g = 0; g < num_groups; ++g) {
        LabelGroup& group = index->groups_[g];
        group.label = in.pod<Label>();
        const std::uint32_t num_buckets = in.pod<std::uint32_t>();
        group.buckets.resize(static_cast<std::size_t>(num_buckets));
        std::vector<VertexId> vertices;
        for (std::size_t b = 0; b < num_buckets; ++b) {
            Bucket& bucket = group.buckets[b];
            bucket.layer = in.varint();
            bucket.dimension = in.varint();
            bucket.count = in.pod<std::uint32_t>();
            bucket.begin = vertices.size();
            const std::size_t old_size = vertices.size();
            vertices.resize(old_size + static_cast<std::size_t>(bucket.count));
            read_vertices(in, vertices.data() + old_size,
                          static_cast<std::size_t>(bucket.count));
        }
        group.vertices = std::move(vertices);
    }

    const std::uint32_t num_empty_groups = in.pod<std::uint32_t>();
    index->empty_groups_.resize(static_cast<std::size_t>(num_empty_groups));
    for (std::size_t g = 0; g < num_empty_groups; ++g) {
        EmptyGroup& group = index->empty_groups_[g];
        group.label = in.pod<Label>();
        const std::uint32_t count = in.pod<std::uint32_t>();
        group.vertices.resize(static_cast<std::size_t>(count));
        read_vertices(in, group.vertices.data(), group.vertices.size());
    }

    try {
        const std::uint32_t num_labels = in.pod<std::uint32_t>();
        const std::uint8_t present = in.pod<std::uint8_t>();
        if (present != 0) {
            std::vector<Label> alphabet(num_labels);
            in.bytes(reinterpret_cast<char*>(alphabet.data()),
                     alphabet.size() * sizeof(Label));
            std::vector<std::uint8_t> counts(
                static_cast<std::size_t>(num_vertices) * num_labels);
            in.bytes(reinterpret_cast<char*>(counts.data()), counts.size());
            index->num_labels_ = num_labels;
            index->alphabet_ = std::move(alphabet);
            index->ball_counts_ = std::move(counts);
            index->has_ball_labels_ = true;
        }
    } catch (const std::exception&) {
        index->num_labels_ = 0;
        index->alphabet_.clear();
        index->ball_counts_.clear();
        index->has_ball_labels_ = false;
    }

    try {
        char magic[sizeof(INDEX_IDENTITY_MAGIC)];
        in.bytes(magic, sizeof(magic));
        if (std::memcmp(magic, INDEX_IDENTITY_MAGIC, sizeof(magic)) == 0) {
            index->identity_.format_version = in.pod<std::uint32_t>();
            index->identity_.data_fingerprint = in.pod<std::uint64_t>();
            index->identity_.label_fingerprint = in.pod<std::uint64_t>();
            index->identity_.neighborhood_size_limit = in.pod<std::uint64_t>();
            index->identity_.data_vertices = in.pod<std::uint64_t>();
            index->identity_.data_edges = in.pod<std::uint64_t>();
            index->identity_.tau_scale = in.pod<double>();
            index->has_identity_ = true;
        }
    } catch (const std::exception&) {
        index->has_identity_ = false;
    }

    try {
        char magic[sizeof(INDEX_AMSD_MAGIC)];
        in.bytes(magic, sizeof(magic));
        if (std::memcmp(magic, INDEX_AMSD_MAGIC, sizeof(magic)) == 0) {
            const std::uint32_t version = in.pod<std::uint32_t>();
            if (version != 1 && version != AMSD_SECTION_VERSION) {
                throw std::runtime_error("unknown AMSD section version");
            }
            index->amsd_.present = true;
            index->amsd_.configured_ranks = in.pod<std::uint32_t>();
            index->amsd_.configured_bins = in.pod<std::uint32_t>();
            if (version >= 2) {
                index->amsd_.rank_select_mode = in.pod<std::uint32_t>();
                index->amsd_.budget_ratio = in.pod<double>();
                index->amsd_.budget_entries = in.pod<double>();
                index->amsd_.space_entries = in.pod<std::uint64_t>();
                const std::uint32_t num_selected = in.pod<std::uint32_t>();
                index->amsd_.selected_ranks.resize(num_selected);
                for (std::uint32_t& rank : index->amsd_.selected_ranks) {
                    rank = in.pod<std::uint32_t>();
                }
                const std::uint32_t num_order = in.pod<std::uint32_t>();
                index->amsd_.selection_order.resize(num_order);
                for (std::uint32_t& rank : index->amsd_.selection_order) {
                    rank = in.pod<std::uint32_t>();
                }
            }
            const std::uint32_t num_vertices_with_degrees = in.pod<std::uint32_t>();
            index->amsd_.degrees.resize(static_cast<std::size_t>(num_vertices_with_degrees));
            in.bytes(reinterpret_cast<char*>(index->amsd_.degrees.data()),
                     index->amsd_.degrees.size() * sizeof(std::uint32_t));
            const std::uint32_t num_groups = in.pod<std::uint32_t>();
            index->amsd_.groups.resize(static_cast<std::size_t>(num_groups));
            for (std::size_t g = 0; g < num_groups; ++g) {
                AmsdGroup& group = index->amsd_.groups[g];
                group.label = in.pod<Label>();
                group.size = in.pod<std::uint32_t>();
                const std::uint32_t num_ranks = in.pod<std::uint32_t>();
                group.ranks.resize(static_cast<std::size_t>(num_ranks));
                const std::size_t words = (static_cast<std::size_t>(group.size) + 63) / 64;
                for (std::size_t r = 0; r < num_ranks; ++r) {
                    AmsdRankIndex& rank = group.ranks[r];
                    rank.rank = in.pod<std::uint32_t>();
                    const std::uint32_t num_rungs = in.pod<std::uint32_t>();
                    rank.rungs.resize(static_cast<std::size_t>(num_rungs));
                    for (std::size_t b = 0; b < num_rungs; ++b) {
                        rank.rungs[b] = in.pod<double>();
                    }
                    rank.suffix_sets.resize(static_cast<std::size_t>(num_rungs));
                    for (std::size_t b = 0; b < num_rungs; ++b) {
                        rank.suffix_sets[b] = in.pod<std::uint32_t>();
                    }
                    rank.suffix.resize(static_cast<std::size_t>(num_rungs));
                    for (std::size_t b = 0; b < num_rungs; ++b) {
                        rank.suffix[b].resize_bits(static_cast<std::size_t>(group.size));
                        in.bytes(reinterpret_cast<char*>(rank.suffix[b].word_data()),
                                 words * sizeof(std::uint64_t));
                    }
                    index->amsd_.bitmap_bytes +=
                        static_cast<std::uint64_t>(num_rungs) *
                        static_cast<std::uint64_t>(words) * sizeof(std::uint64_t);
                }
            }
            if (version == 1) {
                for (const AmsdGroup& group : index->amsd_.groups) {
                    for (const AmsdRankIndex& rank : group.ranks) {
                        index->amsd_.selected_ranks.push_back(rank.rank);
                    }
                }
                std::sort(index->amsd_.selected_ranks.begin(),
                          index->amsd_.selected_ranks.end());
                index->amsd_.selected_ranks.erase(
                    std::unique(index->amsd_.selected_ranks.begin(),
                                index->amsd_.selected_ranks.end()),
                    index->amsd_.selected_ranks.end());
                index->amsd_.selection_order = index->amsd_.selected_ranks;
                index->amsd_.rank_select_mode = 3;
                std::uint64_t space = 0;
                for (std::uint32_t rank : index->amsd_.selected_ranks) {
                    std::uint64_t count = 0;
                    for (const LocalSpectralSignature& signature : index->signatures_) {
                        if (signature.dimension() >= rank) {
                            ++count;
                        }
                    }
                    space += count;
                }
                index->amsd_.space_entries = space;
            }
        }
    } catch (const std::exception&) {
        index->amsd_ = AmsdIndex();
    }

    return index;
}

void SpectralDominanceIndex::rebuild_amsd(const SlimConfig& config) {
    config.validate();
    if (!has_amsd() || amsd_.degrees.size() != signatures_.size()) {
        throw std::runtime_error(
            "--rebuild-amsd needs a base index that carries the AMSD section (its degree array)");
    }
    const std::size_t num_vertices = signatures_.size();
    const std::vector<std::uint32_t> degrees = amsd_.degrees;
    std::uint32_t max_signature_dim = 0;
    for (const LocalSpectralSignature& signature : signatures_) {
        max_signature_dim =
            std::max(max_signature_dim, static_cast<std::uint32_t>(signature.dimension()));
    }
    const RankSelectionOptions selection_options =
        make_rank_selection_options(config, tau_, num_vertices, max_signature_dim);
    const RankSelectionResult selection = select_rank_set(signatures_, selection_options);

    std::vector<AmsdBuildGroup> groups;
    groups.reserve(groups_.size());
    for (const LabelGroup& group : groups_) {
        groups.push_back({group.label, group.vertices.data(), group.vertices.size()});
    }
    amsd_ = build_amsd_index(groups, signatures_, degrees, selection.selected_ranks,
                             static_cast<std::uint32_t>(config.amsd_bins));
    amsd_.configured_ranks = selection.max_rank;
    amsd_.rank_select_mode = static_cast<std::uint32_t>(config.amsd_rank_select);
    amsd_.budget_ratio = config.amsd_budget_ratio;
    amsd_.budget_entries = selection.budget_entries;
    amsd_.space_entries = selection.space_entries;
    amsd_.selection_order = selection.selection_order;
    if (!config.rank_report_path.empty()) {
        write_rank_report(config.rank_report_path, selection);
    }
    if (!config.rank_signals_dump_path.empty()) {
        write_rank_signals_dump(config.rank_signals_dump_path, tau_, selection.max_rank,
                                signatures_);
    }
}

OfflineIndexBuilder::OfflineIndexBuilder(const LabeledGraph& graph, SlimConfig config)
    : graph_(graph), config_(std::move(config)) {
    config_.validate();
    if (graph_.num_vertices() == 0) {
        throw std::runtime_error("cannot build an index for an empty graph");
    }
}

std::unique_ptr<SpectralDominanceIndex> OfflineIndexBuilder::build() {
    const std::size_t num_vertices = graph_.num_vertices();

    const auto preprocess_start = std::chrono::steady_clock::now();
    const LabelEncoder encoder = LabelEncoder::from_graph(graph_);
    const LabelSpectralOperator op(config_.alpha, encoder);

    std::vector<LocalSpectralSignature> signatures(num_vertices);
    std::vector<Label> labels(num_vertices);
    std::vector<std::pair<Label, VertexId>> empty_records;

    double lambda_max_global = 0.0;
    skipped_vertices_ = 0;

    const auto signature_pass = [&](KHopNeighborhood& neighborhood, SpectralEigenSolver& solver,
                                    VertexId begin, VertexId end, double& lambda_local,
                                    std::vector<std::pair<Label, VertexId>>& empty_local,
                                    std::size_t& skipped_local) {
        for (VertexId v = begin; v < end; ++v) {
            labels[v] = graph_.label(v);

            const std::vector<VertexId>* known_hood = nullptr;
            if (config_.neighborhood_size_limit > 0) {
                known_hood = &neighborhood.extract(v, config_.k);
                if (known_hood->size() > config_.neighborhood_size_limit) {
                    ++skipped_local;
                    empty_local.push_back({labels[v], v});
                    continue;
                }
            }

            const std::vector<VertexId>& hood =
                known_hood != nullptr ? *known_hood : neighborhood.extract(v, config_.k);
            const double* eigenvalues = solver.compute(graph_, op, neighborhood, hood);
            const std::size_t dim = solver.size();
            if (dim > 0 && eigenvalues[dim - 1] > lambda_local) {
                lambda_local = eigenvalues[dim - 1];
            }
            signatures[v].eigenvalues.assign(eigenvalues, eigenvalues + dim);
        }
    };

    int num_threads = config_.threads;
    if (num_vertices < static_cast<std::size_t>(num_threads)) {
        num_threads = static_cast<int>(num_vertices);
    }
    const unsigned int hardware = std::thread::hardware_concurrency();
    if (hardware > 0 && static_cast<unsigned int>(num_threads) > hardware) {
        num_threads = static_cast<int>(hardware);
    }
    if (num_threads <= 1) {
        KHopNeighborhood neighborhood(graph_);
        SpectralEigenSolver solver;
        signature_pass(neighborhood, solver, 0, static_cast<VertexId>(num_vertices),
                       lambda_max_global, empty_records, skipped_vertices_);
    } else {
        std::vector<double> lambda_local(static_cast<std::size_t>(num_threads), 0.0);
        std::vector<std::vector<std::pair<Label, VertexId>>> empty_local(
            static_cast<std::size_t>(num_threads));
        std::vector<std::size_t> skipped_local(static_cast<std::size_t>(num_threads), 0);

        std::atomic<VertexId> next_vertex{0};
        constexpr VertexId CHUNK = 32;
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(num_threads - 1));
        for (int t = 0; t < num_threads; ++t) {
            const auto body = [&, t]() {
                KHopNeighborhood neighborhood(graph_);
                SpectralEigenSolver solver;
                for (;;) {
                    const VertexId begin = next_vertex.fetch_add(CHUNK);
                    if (begin >= num_vertices) {
                        break;
                    }
                    const VertexId end =
                        std::min<VertexId>(static_cast<VertexId>(num_vertices), begin + CHUNK);
                    signature_pass(neighborhood, solver, begin, end,
                                   lambda_local[static_cast<std::size_t>(t)],
                                   empty_local[static_cast<std::size_t>(t)],
                                   skipped_local[static_cast<std::size_t>(t)]);
                }
            };
            if (t == num_threads - 1) {
                body();
            } else {
                workers.emplace_back(body);
            }
        }
        for (std::thread& worker : workers) {
            worker.join();
        }
        for (int t = 0; t < num_threads; ++t) {
            lambda_max_global = std::max(lambda_max_global, lambda_local[static_cast<std::size_t>(t)]);
            skipped_vertices_ += skipped_local[static_cast<std::size_t>(t)];
            empty_records.insert(empty_records.end(),
                                 std::make_move_iterator(empty_local[static_cast<std::size_t>(t)].begin()),
                                 std::make_move_iterator(empty_local[static_cast<std::size_t>(t)].end()));
        }
    }

    std::vector<Label> alphabet;
    std::vector<std::uint8_t> ball_counts;
    {
        alphabet.reserve(num_vertices / 4 + 8);
        for (VertexId v = 0; v < num_vertices; ++v) {
            alphabet.push_back(graph_.label(v));
        }
        std::sort(alphabet.begin(), alphabet.end());
        alphabet.erase(std::unique(alphabet.begin(), alphabet.end()), alphabet.end());
        const std::size_t num_labels = alphabet.size();
        ball_counts.assign(num_vertices * num_labels, 0);

        const std::size_t label_limit = config_.neighborhood_size_limit == 0
                                            ? 0
                                            : 8 * config_.neighborhood_size_limit;

        const std::size_t slice =
            (num_vertices + static_cast<std::size_t>(num_threads) - 1) /
            static_cast<std::size_t>(num_threads);
        std::vector<std::size_t> thread_skipped(static_cast<std::size_t>(num_threads), 0);

        const auto label_body = [&](int t) {
            const std::size_t t_index = static_cast<std::size_t>(t);
            const VertexId begin = static_cast<VertexId>(
                std::min<std::size_t>(num_vertices, slice * t_index));
            const VertexId end = static_cast<VertexId>(
                std::min<std::size_t>(num_vertices, slice * (t_index + 1)));
            KHopNeighborhood neighborhood(graph_);
            std::vector<std::uint32_t> counts(num_labels, 0);
            std::vector<std::uint32_t> stamps(num_labels, 0);
            std::uint32_t generation = 0;
            std::vector<std::uint32_t> touched;
            for (VertexId v = begin; v < end; ++v) {
                std::uint8_t* row = ball_counts.data() + static_cast<std::size_t>(v) * num_labels;
                const std::vector<VertexId>& hood = neighborhood.extract(v, config_.k);
                if (label_limit > 0 && hood.size() > label_limit) {
                    std::fill(row, row + num_labels, 255);
                    ++thread_skipped[t_index];
                    continue;
                }
                if (generation == 0) {
                    std::fill(stamps.begin(), stamps.end(), 0);
                }
                ++generation;
                touched.clear();
                for (VertexId w : hood) {
                    const Label label = graph_.label(w);
                    const auto it = std::lower_bound(alphabet.begin(), alphabet.end(), label);
                    const std::size_t id = static_cast<std::size_t>(it - alphabet.begin());
                    if (stamps[id] != generation) {
                        stamps[id] = generation;
                        counts[id] = 0;
                        touched.push_back(static_cast<std::uint32_t>(id));
                    }
                    if (counts[id] < 255) {
                        ++counts[id];
                    }
                }
                for (std::uint32_t id : touched) {
                    row[id] = static_cast<std::uint8_t>(counts[id]);
                }
            }
        };

        if (num_threads <= 1) {
            label_body(0);
        } else {
            std::vector<std::thread> workers;
            workers.reserve(static_cast<std::size_t>(num_threads - 1));
            for (int t = 0; t < num_threads; ++t) {
                if (t == num_threads - 1) {
                    label_body(t);
                } else {
                    workers.emplace_back(label_body, t);
                }
            }
            for (std::thread& worker : workers) {
                worker.join();
            }
        }
        for (std::size_t t = 0; t < thread_skipped.size(); ++t) {
            label_skipped_vertices_ += thread_skipped[t];
        }
    }

    const auto preprocess_end = std::chrono::steady_clock::now();
    timings_.preprocess_ms = static_cast<std::uint64_t>(std::chrono::duration_cast<
        std::chrono::milliseconds>(preprocess_end - preprocess_start).count());

    const double tau = config_.tau_scale * lambda_max_global /
                       std::sqrt(static_cast<double>(num_vertices));

    if (!config_.rank_signals_dump_path.empty()) {
        std::uint32_t max_rank = 0;
        for (const LocalSpectralSignature& signature : signatures) {
            max_rank = std::max(max_rank, static_cast<std::uint32_t>(signature.dimension()));
        }
        max_rank = std::min(max_rank, static_cast<std::uint32_t>(std::max(0, config_.amsd_ranks)));
        write_rank_signals_dump(config_.rank_signals_dump_path, tau, max_rank, signatures);
    }

    auto kappa_of = [tau](const LocalSpectralSignature& signature)
        -> SpectralDominanceIndex::DominanceKey {
        return {static_cast<std::uint32_t>(std::floor(signature.lambda_max() / tau)),
                static_cast<std::uint32_t>(signature.dimension())};
    };

    std::vector<BucketRecord> records;
    records.reserve(num_vertices - skipped_vertices_);
    if (tau > 0.0) {
        for (VertexId v = 0; v < num_vertices; ++v) {
            if (signatures[v].eigenvalues.empty()) {
                continue;
            }
            const SpectralDominanceIndex::DominanceKey key = kappa_of(signatures[v]);
            records.push_back({labels[v], key.layer, key.dimension, v});
        }
    }

    std::sort(records.begin(), records.end(),
              [](const BucketRecord& a, const BucketRecord& b) {
                  if (a.label != b.label) return a.label < b.label;
                  if (a.layer != b.layer) return a.layer < b.layer;
                  if (a.dimension != b.dimension) return a.dimension < b.dimension;
                  return a.vertex < b.vertex;
              });

    std::sort(empty_records.begin(), empty_records.end());
    std::vector<SpectralDominanceIndex::EmptyGroup> empty_groups;
    for (const auto& [label, vertex] : empty_records) {
        if (empty_groups.empty() || empty_groups.back().label != label) {
            empty_groups.push_back({label, {}});
        }
        empty_groups.back().vertices.push_back(vertex);
    }

    std::vector<SpectralDominanceIndex::LabelGroup> groups;
    for (const BucketRecord& record : records) {
        if (groups.empty() || groups.back().label != record.label) {
            groups.push_back({record.label, {}, {}});
        }
        SpectralDominanceIndex::LabelGroup& group = groups.back();
        if (group.buckets.empty() ||
            group.buckets.back().layer != record.layer ||
            group.buckets.back().dimension != record.dimension) {
            group.buckets.push_back({record.layer, record.dimension,
                                     static_cast<std::uint64_t>(group.vertices.size()), 0});
        }
        SpectralDominanceIndex::Bucket& bucket = group.buckets.back();
        group.vertices.push_back(record.vertex);
        ++bucket.count;
    }

    auto index = std::unique_ptr<SpectralDominanceIndex>(new SpectralDominanceIndex());
    index->alpha_ = config_.alpha;
    index->k_ = config_.k;
    index->tau_ = tau;
    index->signatures_ = std::move(signatures);
    index->groups_ = std::move(groups);
    index->empty_groups_ = std::move(empty_groups);
    index->num_labels_ = static_cast<std::uint32_t>(alphabet.size());
    index->alphabet_ = std::move(alphabet);
    index->ball_counts_ = std::move(ball_counts);
    index->has_ball_labels_ = true;

    const auto assemble_end = std::chrono::steady_clock::now();
    timings_.assemble_ms = static_cast<std::uint64_t>(std::chrono::duration_cast<
        std::chrono::milliseconds>(assemble_end - preprocess_end).count());

    const auto amsd_start = std::chrono::steady_clock::now();
    std::uint32_t max_signature_dim = 0;
    for (const LocalSpectralSignature& signature : index->signatures_) {
        max_signature_dim =
            std::max(max_signature_dim, static_cast<std::uint32_t>(signature.dimension()));
    }
    const RankSelectionOptions selection_options =
        make_rank_selection_options(config_, tau, num_vertices, max_signature_dim);
    const RankSelectionResult selection =
        select_rank_set(index->signatures_, selection_options);

    std::vector<AmsdBuildGroup> amsd_groups;
    amsd_groups.reserve(index->groups_.size());
    for (const SpectralDominanceIndex::LabelGroup& group : index->groups_) {
        amsd_groups.push_back({group.label, group.vertices.data(), group.vertices.size()});
    }
    std::vector<std::uint32_t> degrees(num_vertices, 0);
    for (VertexId v = 0; v < num_vertices; ++v) {
        degrees[v] = static_cast<std::uint32_t>(graph_.degree(v));
    }
    index->amsd_ = build_amsd_index(amsd_groups, index->signatures_, degrees,
                                    selection.selected_ranks,
                                    static_cast<std::uint32_t>(config_.amsd_bins));
    index->amsd_.configured_ranks = selection.max_rank;
    index->amsd_.rank_select_mode = static_cast<std::uint32_t>(config_.amsd_rank_select);
    index->amsd_.budget_ratio = config_.amsd_budget_ratio;
    index->amsd_.budget_entries = selection.budget_entries;
    index->amsd_.space_entries = selection.space_entries;
    index->amsd_.selection_order = selection.selection_order;

    if (!config_.rank_report_path.empty()) {
        write_rank_report(config_.rank_report_path, selection);
    }
    const auto amsd_end = std::chrono::steady_clock::now();
    timings_.amsd_ms = static_cast<std::uint64_t>(std::chrono::duration_cast<
        std::chrono::milliseconds>(amsd_end - amsd_start).count());

    index->identity_.format_version = INDEX_FORMAT_VERSION;
    index->identity_.data_fingerprint = data_graph_fingerprint(graph_);
    index->identity_.label_fingerprint = label_alphabet_fingerprint(graph_);
    index->identity_.neighborhood_size_limit = config_.neighborhood_size_limit;
    index->identity_.data_vertices = static_cast<std::uint64_t>(num_vertices);
    index->identity_.data_edges = static_cast<std::uint64_t>(graph_.num_edges());
    index->identity_.tau_scale = config_.tau_scale;
    index->has_identity_ = true;

    return index;
}
}
