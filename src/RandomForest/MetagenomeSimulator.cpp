// SPDX-License-Identifier: GPL-2.0-only
#include "MetagenomeSimulator.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <condition_variable>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../Utilities/Benchmark.h"
#include "../Utilities/Zstd.h"
#include "../IO/Bgzf.h"
#include "ThreadedGzStream.h"

namespace fs = std::filesystem;

namespace protal::sim {

static std::string lowercase(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static std::string trim_ws(std::string text) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    text.erase(text.begin(), std::find_if(text.begin(), text.end(), not_space));
    text.erase(std::find_if(text.rbegin(), text.rend(), not_space).base(), text.end());
    return text;
}

static std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        auto pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.emplace_back(trim_ws(line.substr(start)));
            break;
        }
        fields.emplace_back(trim_ws(line.substr(start, pos - start)));
        start = pos + 1;
    }
    return fields;
}

static std::vector<std::string> split_semicolon(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        auto pos = line.find(';', start);
        if (pos == std::string::npos) {
            fields.emplace_back(line.substr(start));
            break;
        }
        fields.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

static bool has_rank_prefix(const std::string& token) {
    if (token.size() < 3 || token[1] != '_' || token[2] != '_') {
        return false;
    }
    switch (token[0]) {
        case 'd':
        case 'k':
        case 'p':
        case 'c':
        case 'o':
        case 'f':
        case 'g':
        case 's':
            return true;
        default:
            return false;
    }
}

static std::string strip_rank_prefix(const std::string& token) {
    if (has_rank_prefix(token) && token.size() > 3) {
        return token.substr(3);
    }
    return token;
}

static bool taxonomy_has_unclassified(const std::string& taxonomy) {
    for (const auto& token : split_semicolon(taxonomy)) {
        if (token.empty()) {
            continue;
        }
        if (strip_rank_prefix(token) == "Unclassified") {
            return true;
        }
    }
    return false;
}

static std::string make_unclassified_lineage(const std::string& genome_name) {
    return "k__" + genome_name + ";p__" + genome_name + ";c__" + genome_name + ";o__" + genome_name +
           ";f__" + genome_name + ";g__" + genome_name + ";s__" + genome_name;
}

std::vector<GenomeRecord> read_genome_table(const fs::path& tsv_path) {
    std::cout << "Load genome table" << std::endl;
    std::ifstream in(tsv_path);
    if (!in) {
        throw std::runtime_error("Unable to open genome table: " + tsv_path.string());
    }
    std::vector<GenomeRecord> genomes;
    int name_idx = 0;
    int tax_idx = 1;
    int path_idx = 2;
    int len_idx = -1;
    bool require_length = false;
    bool header_checked = false;
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') {
            continue;
        }

        auto fields = split_tab(line);
        if (!header_checked) {
            int header_name = -1;
            int header_tax = -1;
            int header_path = -1;
            int header_len = -1;
            for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
                auto lf = lowercase(fields[i]);
                if (header_name == -1 &&
                    (lf.find("name") != std::string::npos || lf.find("genome") != std::string::npos ||
                     lf.find("accession") != std::string::npos)) {
                    header_name = i;
                }
                if (header_tax == -1 &&
                    (lf.find("tax") != std::string::npos || lf.find("taxonomy") != std::string::npos)) {
                    header_tax = i;
                }
                if (header_path == -1 &&
                    (lf.find("path") != std::string::npos || lf.find("fasta") != std::string::npos ||
                     lf.find("file") != std::string::npos)) {
                    header_path = i;
                }
                if (header_len == -1 && lf.find("length") != std::string::npos) {
                    header_len = i;
                }
            }
            const bool looks_like_header = header_name != -1 && header_tax != -1 && header_path != -1;
            if (looks_like_header) {
                name_idx = header_name;
                tax_idx = header_tax;
                path_idx = header_path;
                len_idx = header_len;
                require_length = header_len != -1;
                header_checked = true;
                continue;  // header row
            }
            len_idx = static_cast<int>(fields.size()) > 3 ? 3 : -1;  // fallback: fourth column if present
            header_checked = true;
        }

        if (fields.size() <= static_cast<std::size_t>(std::max({name_idx, tax_idx, path_idx}))) {
            throw std::runtime_error("Malformed line " + std::to_string(line_no) + " in " + tsv_path.string());
        }

        std::string name = fields[name_idx];
        std::string taxonomy = fields[tax_idx];
        if (taxonomy_has_unclassified(taxonomy)) {
            taxonomy = make_unclassified_lineage(name);
        }
        std::string fasta_path = fields[path_idx];
        std::optional<std::uint64_t> provided_length;
        if (len_idx >= 0 && static_cast<std::size_t>(len_idx) < fields.size()) {
            if (!fields[len_idx].empty()) {
                try {
                    provided_length = std::stoull(fields[len_idx]);
                } catch (const std::exception&) {
                    throw std::runtime_error("Invalid genome_length on line " + std::to_string(line_no) + " in " +
                                             tsv_path.string());
                }
            } else if (require_length) {
                throw std::runtime_error("Missing genome_length on line " + std::to_string(line_no) + " in " +
                                         tsv_path.string());
            }
        } else if (require_length) {
            throw std::runtime_error("Missing genome_length column on line " + std::to_string(line_no) + " in " +
                                     tsv_path.string());
        }
        genomes.push_back(GenomeRecord{
            std::move(name),
            std::move(taxonomy),
            fs::path(std::move(fasta_path)),
            provided_length});
    }
    if (genomes.empty()) {
        throw std::runtime_error("Genome table is empty: " + tsv_path.string());
    }
    return genomes;
}

std::vector<double> parse_strain_probabilities(const std::string& text) {
    std::vector<double> probs;
    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = text.substr(start, end - start);
        token.erase(token.begin(), std::find_if(token.begin(), token.end(), [](unsigned char c) { return !std::isspace(c); }));
        token.erase(std::find_if(token.rbegin(), token.rend(), [](unsigned char c) { return !std::isspace(c); }).base(), token.end());
        if (!token.empty()) {
            probs.push_back(std::stod(token));
        }
        start = end + 1;
    }
    return probs;
}

std::vector<std::string> parse_species_list(const std::string& text) {
    std::vector<std::string> species;
    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = text.substr(start, end - start);
        token.erase(token.begin(),
                    std::find_if(token.begin(), token.end(), [](unsigned char c) { return !std::isspace(c); }));
        token.erase(std::find_if(token.rbegin(), token.rend(), [](unsigned char c) { return !std::isspace(c); }).base(),
                    token.end());
        if (!token.empty()) {
            species.push_back(token);
        }
        start = end + 1;
    }
    return species;
}

std::unordered_map<std::string, std::size_t> parse_genus_selection(const std::string& text) {
    std::unordered_map<std::string, std::size_t> genus_counts;
    auto trim = [](std::string s) {
        s.erase(s.begin(),
                std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isspace(c); }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c) { return !std::isspace(c); }).base(), s.end());
        return s;
    };

    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = trim(text.substr(start, end - start));
        if (!token.empty()) {
            auto sep = token.find(':');
            if (sep == std::string::npos) {
                throw std::runtime_error("Invalid --genus entry (expected genus:count): " + token);
            }
            std::string genus = trim(token.substr(0, sep));
            if (genus.rfind("g__", 0) == 0 && genus.size() > 3) {
                genus = genus.substr(3);  // allow g__ prefix from GTDB-style names
            }
            std::string count_str = trim(token.substr(sep + 1));
            if (genus.empty() || count_str.empty()) {
                throw std::runtime_error("Invalid --genus entry (missing genus or count): " + token);
            }
            std::size_t count = 0;
            try {
                count = static_cast<std::size_t>(std::stoull(count_str));
            } catch (const std::exception&) {
                throw std::runtime_error("Invalid --genus count for entry: " + token);
            }
            if (count > 0) {
                genus_counts[genus] += count;
            }
        }
        start = end + 1;
    }
    return genus_counts;
}

std::unordered_map<std::string, std::size_t> parse_taxon_selection(const std::string& text) {
    std::unordered_map<std::string, std::size_t> taxon_counts;
    auto trim = [](std::string s) {
        s.erase(s.begin(),
                std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isspace(c); }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c) { return !std::isspace(c); }).base(), s.end());
        return s;
    };

    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = trim(text.substr(start, end - start));
        if (!token.empty()) {
            auto sep = token.find(':');
            if (sep == std::string::npos) {
                throw std::runtime_error("Invalid --taxon entry (expected taxon:count): " + token);
            }
            std::string taxon = trim(token.substr(0, sep));
            std::string count_str = trim(token.substr(sep + 1));
            if (taxon.empty() || count_str.empty()) {
                throw std::runtime_error("Invalid --taxon entry (missing taxon or count): " + token);
            }
            std::size_t count = 0;
            try {
                count = static_cast<std::size_t>(std::stoull(count_str));
            } catch (const std::exception&) {
                throw std::runtime_error("Invalid --taxon count for entry: " + token);
            }
            if (count > 0) {
                taxon_counts[taxon] += count;
            }
        }
        start = end + 1;
    }
    return taxon_counts;
}

void parse_congener_groups(const std::string& text, ProfileDesignOptions& options) {
    if (text.empty()) return;
    auto const fail = [&text]() {
        throw std::runtime_error("Invalid --congener_groups (expected SHARE:MIN-MAX, a share of the sample's species "
                                 "0-1 and 2 <= MIN <= MAX species per genus, e.g. 0.25:2-5): " + text);
    };
    auto const colon = text.find(':');
    auto const dash = text.find('-', colon == std::string::npos ? 0 : colon);
    if (colon == std::string::npos || dash == std::string::npos) fail();
    double share = -1;
    std::size_t lo = 0, hi = 0;
    try {
        std::size_t used = 0;
        share = std::stod(text.substr(0, colon), &used);
        if (used != colon) fail();
        std::string const lo_text = text.substr(colon + 1, dash - colon - 1), hi_text = text.substr(dash + 1);
        if (lo_text.empty() || hi_text.empty() ||
            lo_text.find_first_not_of("0123456789") != std::string::npos ||
            hi_text.find_first_not_of("0123456789") != std::string::npos) fail();
        lo = static_cast<std::size_t>(std::stoull(lo_text));
        hi = static_cast<std::size_t>(std::stoull(hi_text));
    } catch (const std::invalid_argument&) {
        fail();
    } catch (const std::out_of_range&) {
        fail();
    }
    if (!(share >= 0 && share <= 1) || lo < 2 || hi < lo) fail();
    options.congener_share = share;
    options.congener_min = lo;
    options.congener_max = hi;
}

// fasta_path and art_seed are appended after the historical columns, never inserted
// between them: downstream consumers read manifests by column name, so growing the
// header on the right keeps every existing reader working.
static constexpr const char* kManifestHeader =
    "sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance"
    "\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\n";

static void write_manifest_row(std::ofstream& out, const SampleOutput& sample, const GenomeAssignment& assignment) {
    out << sample.sample_name << '\t' << assignment.genome.name << '\t' << assignment.species << '\t'
        << assignment.genome.taxonomy << '\t' << assignment.genome_length << '\t' << assignment.read_pairs << '\t'
        << assignment.vertical_coverage << '\t' << assignment.relative_abundance << '\t'
        << sample.read1_path.string() << '\t' << sample.read2_path.string() << '\t'
        << assignment.genome.fasta_path.string() << '\t';
    if (assignment.art_seed) {
        out << *assignment.art_seed;
    }
    out << '\n';
}

void write_sample_manifest(const SampleOutput& sample, const fs::path& manifest_path) {
    if (!manifest_path.parent_path().empty()) {
        fs::create_directories(manifest_path.parent_path());
    }
    std::ofstream out(manifest_path);
    if (!out) {
        throw std::runtime_error("Unable to write manifest: " + manifest_path.string());
    }
    out << kManifestHeader;
    for (const auto& assignment : sample.assignments) {
        write_manifest_row(out, sample, assignment);
    }
}

void write_combined_manifest(const std::vector<SampleOutput>& samples, const fs::path& manifest_path) {
    if (!manifest_path.parent_path().empty()) {
        fs::create_directories(manifest_path.parent_path());
    }
    std::ofstream out(manifest_path);
    if (!out) {
        throw std::runtime_error("Unable to write manifest: " + manifest_path.string());
    }
    out << kManifestHeader;
    for (const auto& sample : samples) {
        for (const auto& assignment : sample.assignments) {
            write_manifest_row(out, sample, assignment);
        }
    }
}

MetagenomeSimulator::MetagenomeSimulator(
    std::vector<GenomeRecord> genomes, ArtIlluminaOptions art_options, std::uint64_t seed)
    : genomes_(std::move(genomes)),
      art_(std::move(art_options)),
      designer_(genomes_),
      rng_(seed) {}

static std::uint64_t read_genome_length(const fs::path& fasta_path) {
    auto accumulate_length = [](auto&& getter, auto&& handle) -> std::uint64_t {
        std::string line;
        std::uint64_t total = 0;
        while (getter(handle, line)) {
            if (!line.empty() && line[0] == '>') {
                continue;
            }
            for (char c : line) {
                if (std::isalpha(static_cast<unsigned char>(c))) {
                    ++total;
                }
            }
        }
        return total;
    };

    if (fasta_path.extension() == ".gz") {
        // BGZF with libdeflate, other gzip with zlib-ng (ThreadedGzStream.h).
        protal::ThreadedGzIstream input(fasta_path.string().c_str());
        if (!input.rdbuf()->is_open()) {
            throw std::runtime_error("Unable to open compressed fasta: " + fasta_path.string());
        }
        auto getter = [](protal::ThreadedGzIstream& file, std::string& out) -> bool {
            if (!std::getline(file, out)) return false;
            if (!out.empty() && out.back() == '\r') out.pop_back();
            return true;
        };
        std::uint64_t len = accumulate_length(getter, input);
        if (input.rdbuf()->read_failed()) {
            throw std::runtime_error("The compressed fasta " + fasta_path.string() + " is truncated or corrupt (" +
                                     input.rdbuf()->read_error_message() + ")");
        }
        return len;
    }

    std::ifstream in(fasta_path);
    if (!in) {
        throw std::runtime_error("Unable to open fasta: " + fasta_path.string());
    }
    auto getter = [](std::ifstream& file, std::string& out) -> bool {
        return static_cast<bool>(std::getline(file, out));
    };
    return accumulate_length(getter, in);
}

static std::unordered_map<std::string, std::uint64_t> build_length_cache(const std::vector<GenomeRecord>& genomes) {
    std::unordered_map<std::string, std::uint64_t> lengths;
    lengths.reserve(genomes.size());
    size_t read_genomes = 0;
    for (const auto& genome : genomes) {
        if ((read_genomes % 100) == 0) {
            std::cout << "genomes processed: " << read_genomes << std::endl;
        }
        const auto len = genome.genome_length ? *genome.genome_length : read_genome_length(genome.fasta_path);
        lengths.emplace(genome.name, len);
        ++read_genomes;
    }
    return lengths;
}

//static void append_fastq(const fs::path& src, std::ofstream& dst) {
//    std::ifstream in(src, std::ios::binary);
//    if (!in) {
//        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
//    }
//    dst << in.rdbuf();
//}

//static void append_fastq(const fs::path& src, std::ofstream& dst) {
//    std::ifstream in(src, std::ios::binary);
//    if (!in) {
//        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
//    }
//
//    dst << in.rdbuf();
//
//    if (!dst) {
//        throw std::runtime_error("Write failed while appending " + src.string());
//    }
//
//    dst.clear();  // ← THIS IS THE CRITICAL FIX
//}

// A sample's read file as its reads arrive: BGZF (libdeflate) or one zstd frame at level 3 with a checksum, without a
// long window (FASTQ gains nothing from one). Error() says why writing failed; Close() finishes the file.
class ReadsWriter {
public:
    ReadsWriter(std::string const& path, ReadsCompression compression) {
        if (compression == ReadsCompression::Zstd) {
            m_zstd = std::make_unique<protal::zstd::OStream>(path, protal::zstd::Params{3, 0, 1, 0});
        } else {
            m_bgzf = std::make_unique<protal::bgzf::Writer>(path);
        }
    }

    void Write(char const* data, std::size_t size) {
        if (m_zstd) m_zstd->write(data, static_cast<std::streamsize>(size));
        else m_bgzf->Write(data, size);
    }

    std::string Error() const {
        if (!m_zstd) return m_bgzf->Error();
        if (!m_zstd->Buffer().Error().empty()) return m_zstd->Buffer().Error();
        return m_zstd->fail() ? "writing the zstd stream failed" : "";
    }

    bool Close() { return m_zstd ? m_zstd->Close() : m_bgzf->Close(); }

private:
    std::unique_ptr<protal::bgzf::Writer> m_bgzf;
    std::unique_ptr<protal::zstd::OStream> m_zstd;
};

// The file name suffix of a sample's reads: .fq.gz (BGZF) or .fq.zst.
static std::string ReadsSuffix(ReadsCompression compression) {
    return compression == ReadsCompression::Zstd ? ".fq.zst" : ".fq.gz";
}

static void append_fastq(const fs::path& src, ReadsWriter& dst)
{
    std::ifstream in(src, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
    }

    constexpr std::size_t bufsize = 1 << 20; // 1 MB
    std::vector<char> buffer(bufsize);

    while (in) {
        in.read(buffer.data(), buffer.size());
        std::streamsize n = in.gcount();
        if (n > 0) {
            dst.Write(buffer.data(), static_cast<std::size_t>(n));
        }
    }
    if (in.bad()) {
        throw std::runtime_error("Read failed while appending " + src.string());
    }
    if (!dst.Error().empty()) {
        throw std::runtime_error("Write failed while appending " + src.string() + ": " + dst.Error());
    }
}


static void normalize_relative_abundance(SampleOutput& sample) {
    double coverage_sum = 0.0;
    for (const auto& assignment : sample.assignments) {
        coverage_sum += assignment.relative_abundance;
    }
    for (auto& assignment : sample.assignments) {
        assignment.relative_abundance =
            coverage_sum > 0.0 ? assignment.relative_abundance / coverage_sum : 0.0;
    }
}

SampleOutput MetagenomeSimulator::simulate_single(
        const ProfileDesignOptions& profile_options,
        const std::string& sample_name,
        const std::unordered_map<std::string, std::uint64_t>& genome_lengths,
        std::uint64_t paired_read_length)
{
    auto assignments = designer_.design_profile(profile_options, rng_);
    for (auto& assignment : assignments) {
        auto it_len = genome_lengths.find(assignment.genome.name);
        if (it_len == genome_lengths.end()) {
            throw std::runtime_error("Missing genome length for " + assignment.genome.name);
        }
        assignment.genome_length = it_len->second;
    }
    SampleOutput sample{sample_name, {}, {}, std::move(assignments)};
    prepare_sample(sample, paired_read_length);
    return sample;
}

// Runs ART over assignments whose genome, read_pairs and genome_length are already
// fixed — by the profile designer on a fresh run, or by a manifest on a replay.
// Assignments whose genome, read_pairs and genome_length are already fixed — by the profile
// designer on a fresh run, or by a manifest on a replay.
void MetagenomeSimulator::prepare_sample(SampleOutput& sample, std::uint64_t paired_read_length)
{
    for (auto& assignment : sample.assignments) {
            const auto genome_len = assignment.genome_length;
            if (genome_len == 0) {
                throw std::runtime_error("Missing genome length for " + assignment.genome.name);
            }
            // Drawn even when reads are skipped, so that a --test design and the
            // corresponding real run consume the RNG identically. ART reads -rs as a signed
            // 32-bit number, so larger seeds would alias; the draw keeps 31 bits.
            const auto drawn_seed = static_cast<unsigned int>(rng_() & 0x7fffffffu);
            if (auto forced_seed = art_.seed_override()) {
                // An -rs in extra_art_args beats both the drawn seed and a replayed one.
                static bool warned = false;
                if (*forced_seed > 0x7fffffffu && !warned) {
                    std::cerr << "Warning: -rs " << *forced_seed << " is above 2^31-1, which ART reads "
                                 "as a different (signed 32-bit) seed" << std::endl;
                    warned = true;
                }
                assignment.art_seed = *forced_seed;
            } else if (!assignment.art_seed) {
                assignment.art_seed = drawn_seed;
            }
            const double bases = static_cast<double>(assignment.read_pairs * paired_read_length);
            const double coverage = bases / static_cast<double>(genome_len);
            assignment.vertical_coverage = coverage;
            assignment.relative_abundance = coverage;  // normalized later in simulate_samples
    }
}

void MetagenomeSimulator::write_reads(
        SampleOutput& sample,
        const fs::path& output_dir,
        bool skip_reads,
        bool keep_tmp,
        std::size_t threads) const
{
    const std::string& sample_name = sample.sample_name;
    fs::path reads_dir = output_dir / "reads";
    fs::create_directories(reads_dir);
    const fs::path sample_prefix = reads_dir / sample_name;
    const fs::path r1_gz = sample_prefix.string() + "_R1" + ReadsSuffix(reads_compression_);
    const fs::path r2_gz = sample_prefix.string() + "_R2" + ReadsSuffix(reads_compression_);
    sample.read1_path = r1_gz;
    sample.read2_path = r2_gz;
    if (skip_reads) {
        // Create empty placeholder files for manifests to point to.
        std::ofstream placeholder1(r1_gz, std::ios::binary);
        std::ofstream placeholder2(r2_gz, std::ios::binary);
        (void)placeholder1;
        (void)placeholder2;
        return;
    }

    const fs::path temp_dir = output_dir / (sample_name + "_tmp");
    fs::create_directories(temp_dir);
    // Truncated: a leftover file of an interrupted run must not be extended. Compressed as the reads arrive, so that
    // a deep sample's uncompressed reads are never on disk; replays give byte-identical files.
    ReadsWriter r1_out(r1_gz.string(), reads_compression_);
    ReadsWriter r2_out(r2_gz.string(), reads_compression_);
    if (!r1_out.Error().empty() || !r2_out.Error().empty()) {
        throw std::runtime_error("Unable to create output FASTQ files for " + sample_name + ": " + r1_out.Error() + r2_out.Error());
    }
    // ART runs for the genomes in the sample's order, each in a folder of its own (its place in the sample), on up
    // to `threads` threads (this one included) and at most 2 x threads genomes ahead of the next to append; this
    // thread appends them in that order, so the files are the same for any number of threads.
    auto const& assignments = sample.assignments;
    std::size_t const n = assignments.size();
    std::size_t const helpers = std::min(std::max<std::size_t>(1, threads), std::max<std::size_t>(1, n)) - 1;
    std::size_t const window = 2 * (helpers + 1);
    std::vector<std::pair<fs::path, fs::path>> made(n);
    std::vector<char> ready(n, 0);
    std::vector<std::exception_ptr> failed(n);
    std::size_t next = 0, appended = 0;
    bool stop = false;
    std::mutex mutex;
    std::condition_variable changed;
    auto genome_dir = [&](std::size_t i) { return temp_dir / std::to_string(i); };
    auto simulate = [&](std::size_t i) {  // without the lock
        std::pair<fs::path, fs::path> files;
        std::exception_ptr error;
        try {
            auto const& assignment = assignments[i];
            fs::create_directories(genome_dir(i));
            files = art_.simulate_read_pairs(
                assignment.genome, assignment.read_pairs, assignment.genome_length,
                genome_dir(i) / assignment.genome.name, static_cast<unsigned int>(*assignment.art_seed), genome_dir(i));
        } catch (...) {
            error = std::current_exception();
        }
        std::lock_guard<std::mutex> lock(mutex);
        made[i] = std::move(files);
        failed[i] = error;
        ready[i] = 1;
        changed.notify_all();
    };
    auto may_start = [&] { return !stop && next < n && next < appended + window; };
    auto helper = [&] {
        for (;;) {
            std::size_t i;
            {
                std::unique_lock<std::mutex> lock(mutex);
                changed.wait(lock, [&] { return stop || next >= n || may_start(); });
                if (!may_start()) return;
                i = next++;
            }
            simulate(i);
        }
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < helpers; t++) pool.emplace_back(helper);
    auto finish_pool = [&] {
        {
            std::lock_guard<std::mutex> lock(mutex);
            stop = true;
            changed.notify_all();
        }
        for (auto& thread : pool) thread.join();
        pool.clear();
    };
    try {
        for (std::size_t j = 0; j < n; j++) {
            for (;;) {  // until genome j is simulated; this thread simulates the next genome itself meanwhile
                std::unique_lock<std::mutex> lock(mutex);
                if (ready[j]) break;
                if (may_start()) {
                    std::size_t const i = next++;
                    lock.unlock();
                    simulate(i);
                    continue;
                }
                changed.wait(lock, [&] { return ready[j] || may_start(); });
            }
            if (failed[j]) std::rethrow_exception(failed[j]);
            append_fastq(made[j].first, r1_out);
            append_fastq(made[j].second, r2_out);
            if (!keep_tmp) {
                // This genome's reads, and its decompressed copy (ArtIlluminaWrapper::ensure_fasta), are not needed again.
                std::error_code ec;
                fs::remove_all(genome_dir(j), ec);
            }
            std::lock_guard<std::mutex> lock(mutex);
            appended = j + 1;
            changed.notify_all();
        }
    } catch (...) {
        finish_pool();
        throw;
    }
    finish_pool();
    if (!r1_out.Close() || !r2_out.Close()) {
        throw std::runtime_error("Unable to finish writing FASTQ files for " + sample_name + ": " + r1_out.Error() + r2_out.Error());
    }
    if (!keep_tmp) {
        fs::remove_all(temp_dir);
    }
}

void MetagenomeSimulator::write_all_reads(
        std::vector<SampleOutput>& samples,
        const fs::path& output_dir,
        bool skip_reads,
        bool keep_tmp) const
{
    // A sample per thread; threads beyond the samples go to their genomes (write_reads), the first samples taking
    // the remainder.
    std::size_t const threads = static_cast<std::size_t>(std::max(1, art_.options().threads));
    std::size_t const workers = std::min<std::size_t>(samples.size(), threads);
    std::atomic<std::size_t> next{0};
    std::vector<std::exception_ptr> errors(samples.size());
    std::mutex print;
    auto work = [&](std::size_t worker) {
        std::size_t const genome_threads = workers ? threads / workers + (worker < threads % workers ? 1 : 0) : 1;
        for (std::size_t i = next++; i < samples.size(); i = next++) {
            protal::Benchmark timer("write the reads of " + samples[i].sample_name);
            timer.Start();
            try {
                write_reads(samples[i], output_dir, skip_reads, keep_tmp, genome_threads);
            } catch (...) {
                errors[i] = std::current_exception();
                next = samples.size();  // no new samples: the run fails
            }
            timer.Stop();
            std::lock_guard<std::mutex> lock(print);
            timer.PrintResults();
        }
    };
    std::vector<std::thread> pool;
    for (std::size_t w = 1; w < workers; w++) pool.emplace_back(work, w);
    work(0);
    for (auto& thread : pool) thread.join();
    for (auto const& error : errors) {
        if (error) std::rethrow_exception(error);
    }
}

std::vector<SampleOutput> MetagenomeSimulator::simulate_samples(
    const ProfileDesignOptions& profile_options,
    std::size_t sample_count,
    const std::string& sample_prefix,
    const fs::path& output_dir,
    bool skip_reads,
    bool keep_tmp)
{
    if (sample_count == 0) {
        return {};
    }
    const auto genome_lengths = build_length_cache(genomes_);
    const auto paired_read_length = static_cast<std::uint64_t>(art_.options().read_length) * 2ULL;

    // Pre-assign strains across samples for all strain_sharing specs.
    auto strain_assignments = designer_.assign_strains_across_samples(
        profile_options.strain_sharing, sample_count,
        profile_options.total_read_pairs, paired_read_length,
        genome_lengths, rng_);

    // Every sample's design and seeds first, in the order a sample-by-sample run draws them from rng_, then the
    // samples' reads on threads (write_all_reads): the same samples for any thread count.
    std::vector<SampleOutput> outputs;
    outputs.reserve(sample_count);
    for (std::size_t i = 0; i < sample_count; ++i) {
        std::ostringstream name;
        name << sample_prefix << "_" << (i + 1);

        ProfileDesignOptions per_sample_opts = profile_options;
        per_sample_opts.pln_sigma = SigmaForSample(profile_options, i);  // the samples' sigmas in turn
        per_sample_opts.total_read_pairs = ReadPairsForSample(profile_options, i);  // and depths
        if (profile_options.species_per_sample_min > 0 &&
            profile_options.species_per_sample_min < profile_options.species_per_sample) {
            per_sample_opts.species_per_sample = std::uniform_int_distribution<std::size_t>(
                profile_options.species_per_sample_min,
                profile_options.species_per_sample)(rng_);
        }
        per_sample_opts.forced_strains      = strain_assignments[i].forced_strains;
        per_sample_opts.species_min_abundance = strain_assignments[i].species_min_abundance;

        auto sample = simulate_single(per_sample_opts, name.str(), genome_lengths, paired_read_length);
        normalize_relative_abundance(sample);
        outputs.push_back(std::move(sample));
    }
    write_all_reads(outputs, output_dir, skip_reads, keep_tmp);
    return outputs;
}

std::vector<SampleOutput> MetagenomeSimulator::replay_samples(
    std::vector<SampleOutput> design,
    const fs::path& output_dir,
    bool skip_reads,
    bool keep_tmp)
{
    const auto paired_read_length = static_cast<std::uint64_t>(art_.options().read_length) * 2ULL;
    for (auto& sample : design) {  // rng_ in the order of the samples, then their reads on threads
        prepare_sample(sample, paired_read_length);
        normalize_relative_abundance(sample);
    }
    write_all_reads(design, output_dir, skip_reads, keep_tmp);
    return design;
}

std::vector<SampleOutput> read_manifest(const fs::path& manifest_path) {
    std::ifstream in(manifest_path);
    if (!in) {
        throw std::runtime_error("Unable to open manifest: " + manifest_path.string());
    }

    std::unordered_map<std::string, int> col;
    std::vector<SampleOutput> samples;
    std::unordered_map<std::string, std::size_t> sample_index;
    std::string line;
    std::size_t line_no = 0;
    bool header_seen = false;

    auto column = [&](const std::vector<std::string>& fields, const char* name, bool required) -> std::string {
        auto it = col.find(name);
        if (it == col.end() || static_cast<std::size_t>(it->second) >= fields.size() || fields[it->second].empty()) {
            if (required) {
                throw std::runtime_error("Manifest line " + std::to_string(line_no) + ": missing column \"" +
                                         std::string(name) + "\" in " + manifest_path.string());
            }
            return {};
        }
        return fields[it->second];
    };

    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        auto fields = split_tab(line);
        if (!header_seen) {
            for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
                col.emplace(lowercase(fields[i]), i);
            }
            for (const char* required : {"sample", "genome", "species", "taxonomy", "genome_length", "read_pairs"}) {
                if (!col.count(required)) {
                    throw std::runtime_error("Manifest " + manifest_path.string() +
                                             " is missing the required column \"" + std::string(required) + "\"");
                }
            }
            header_seen = true;
            continue;
        }

        const std::string sample_name = column(fields, "sample", true);
        auto [it, inserted] = sample_index.emplace(sample_name, samples.size());
        if (inserted) {
            samples.push_back(SampleOutput{sample_name, {}, {}, {}});
        }

        GenomeAssignment assignment;
        assignment.genome.name = column(fields, "genome", true);
        assignment.genome.taxonomy = column(fields, "taxonomy", true);
        assignment.genome.fasta_path = fs::path(column(fields, "fasta_path", false));
        assignment.species = column(fields, "species", true);
        try {
            assignment.genome_length = std::stoull(column(fields, "genome_length", true));
            assignment.read_pairs = std::stoull(column(fields, "read_pairs", true));
        } catch (const std::exception&) {
            throw std::runtime_error("Manifest line " + std::to_string(line_no) +
                                     ": genome_length and read_pairs must be integers");
        }
        assignment.genome.genome_length = assignment.genome_length;
        // Recorded so callers can recover the read length the run used; recomputed on replay.
        const std::string vcov_field = column(fields, "vertical_coverage", false);
        if (!vcov_field.empty()) {
            assignment.vertical_coverage = std::stod(vcov_field);
        }
        const std::string seed_field = column(fields, "art_seed", false);
        if (!seed_field.empty()) {
            assignment.art_seed = std::stoull(seed_field);
        }
        samples[it->second].assignments.push_back(std::move(assignment));
    }

    if (samples.empty()) {
        throw std::runtime_error("Manifest has no rows: " + manifest_path.string());
    }
    return samples;
}

void resolve_manifest_fasta_paths(std::vector<SampleOutput>& samples, const std::vector<GenomeRecord>& genomes) {
    std::unordered_map<std::string, fs::path> by_name;
    by_name.reserve(genomes.size());
    for (const auto& genome : genomes) {
        by_name.emplace(genome.name, genome.fasta_path);
    }
    std::vector<std::string> unresolved;
    for (auto& sample : samples) {
        for (auto& assignment : sample.assignments) {
            if (!assignment.genome.fasta_path.empty()) {
                continue;
            }
            auto it = by_name.find(assignment.genome.name);
            if (it == by_name.end()) {
                unresolved.push_back(assignment.genome.name);
                continue;
            }
            assignment.genome.fasta_path = it->second;
        }
    }
    if (!unresolved.empty()) {
        std::sort(unresolved.begin(), unresolved.end());
        unresolved.erase(std::unique(unresolved.begin(), unresolved.end()), unresolved.end());
        std::string message = "Manifest has no fasta_path for " + std::to_string(unresolved.size()) +
                              " genome(s) and they are absent from the genome table: ";
        for (std::size_t i = 0; i < unresolved.size() && i < 5; ++i) {
            message += (i ? ", " : "") + unresolved[i];
        }
        if (unresolved.size() > 5) {
            message += ", ...";
        }
        throw std::runtime_error(message);
    }
}

void write_abundance_matrix(const std::vector<SampleOutput>& samples, const fs::path& matrix_path) {
    if (samples.empty()) {
        return;
    }
    if (!matrix_path.parent_path().empty()) {
        fs::create_directories(matrix_path.parent_path());
    }

    // Collect all species
    std::unordered_map<std::string, std::vector<double>> species_to_samples;
    const std::size_t sample_count = samples.size();
    for (const auto& sample : samples) {
        for (const auto& assignment : sample.assignments) {
            species_to_samples[assignment.species].resize(sample_count, 0.0);
        }
    }
    // Sum relative abundances per species per sample
    for (std::size_t idx = 0; idx < samples.size(); ++idx) {
        for (const auto& assignment : samples[idx].assignments) {
            species_to_samples[assignment.species][idx] += assignment.relative_abundance;
        }
    }

    std::ofstream out(matrix_path);
    if (!out) {
        throw std::runtime_error("Unable to write abundance matrix: " + matrix_path.string());
    }
    out << "species";
    for (const auto& sample : samples) {
        out << '\t' << sample.sample_name;
    }
    out << '\n';
    for (const auto& [species, values] : species_to_samples) {
        out << species;
        for (double val : values) {
            out << '\t' << val;
        }
        out << '\n';
    }
}

}  // namespace protal::sim
