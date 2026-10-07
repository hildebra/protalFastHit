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
    std::vector<GenomeRecord> genomes, IlluminaOptions illumina_options, std::uint64_t seed)
    : genomes_(std::move(genomes)),
      illumina_(std::move(illumina_options)),
      seed_(seed),
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
        // BGZF or other gzip, inflated with ISA-L (ThreadedGzStream.h).
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

static std::string ReadsSuffix(ReadsCompression compression) {
    return compression == ReadsCompression::Zstd ? ".fq.zst" : ".fq.gz";
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

// Prepares assignments whose genome, read_pairs and genome_length are already fixed — by the profile
// designer on a fresh run, or by a manifest on a replay.
void MetagenomeSimulator::prepare_sample(SampleOutput& sample, std::uint64_t paired_read_length)
{
    for (auto& assignment : sample.assignments) {
            const auto genome_len = assignment.genome_length;
            if (genome_len == 0) {
                throw std::runtime_error("Missing genome length for " + assignment.genome.name);
            }
            // Drawn even when reads are skipped, so that a --test design and the corresponding real run consume the
            // RNG identically (31 bits, as when ART read it as a signed 32-bit number); a replayed seed is kept.
            const auto drawn_seed = static_cast<unsigned int>(rng_() & 0x7fffffffu);
            if (!assignment.art_seed) {
                assignment.art_seed = drawn_seed;
            }
            const double bases = static_cast<double>(assignment.read_pairs * paired_read_length);
            const double coverage = bases / static_cast<double>(genome_len);
            assignment.vertical_coverage = coverage;
            assignment.relative_abundance = coverage;  // normalized later in simulate_samples
    }
}

static bool is_named_pipe(const fs::path& path) {
    std::error_code ec;
    return fs::is_fifo(path, ec);
}

void MetagenomeSimulator::write_all_reads(
        std::vector<SampleOutput>& samples,
        const fs::path& output_dir,
        bool skip_reads) const
{
    fs::path reads_dir = output_dir / "reads";
    fs::create_directories(reads_dir);
    std::vector<PairedSample> paired;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        auto& sample = samples[i];
        const fs::path prefix = reads_dir / sample.sample_name;
        sample.read1_path = prefix.string() + "_R1" + ReadsSuffix(reads_compression_);
        sample.read2_path = prefix.string() + "_R2" + ReadsSuffix(reads_compression_);
        if (skip_reads) {
            for (const auto& path : {sample.read1_path, sample.read2_path}) {
                if (!is_named_pipe(path)) std::ofstream placeholder(path, std::ios::binary);  // for manifests to name
            }
            continue;
        }
        PairedSample reads;
        reads.name = sample.sample_name;
        reads.r1 = sample.read1_path;
        if (!illumina_.first_reads_only) reads.r2 = sample.read2_path;
        for (const auto& assignment : sample.assignments) {
            reads.genomes.push_back({assignment.genome.name, assignment.genome.fasta_path, assignment.read_pairs,
                                     assignment.art_seed.value_or(0)});
        }
        if (!illumina_.host_pairs.empty()) reads.host_pairs = illumina_.host_pairs[i % illumina_.host_pairs.size()];
        reads.host_seed = MixSeed(seed_, 0x686f7374ULL + i);   // "host"
        reads.run_seed = MixSeed(seed_, 0x72756e00ULL + i);    // "run"
        paired.push_back(std::move(reads));
    }
    if (skip_reads) return;
    PairedOptions options;
    options.setup.profile = IlluminaProfile::Named(illumina_.sequencer);
    options.setup.read_length = illumina_.read_length;
    options.setup.fragment_mean = illumina_.fragment_mean;
    options.setup.fragment_sd = illumina_.fragment_stdev;
    options.setup.mean_quality = illumina_.mean_quality;
    options.host = illumina_.host_folder;
    options.threads = std::max(1, illumina_.threads);
    protal::Benchmark timer("write the reads of " + std::to_string(samples.size()) + " samples");
    timer.Start();
    SimulatePairs(paired, options);
    timer.Stop();
    timer.PrintResults();
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
    const auto paired_read_length = static_cast<std::uint64_t>(illumina_.read_length) * 2ULL;

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
    write_all_reads(outputs, output_dir, skip_reads);
    return outputs;
}

std::vector<SampleOutput> MetagenomeSimulator::replay_samples(
    std::vector<SampleOutput> design,
    const fs::path& output_dir,
    bool skip_reads,
    bool keep_tmp)
{
    const auto paired_read_length = static_cast<std::uint64_t>(illumina_.read_length) * 2ULL;
    for (auto& sample : design) {  // rng_ in the order of the samples, then their reads on threads
        prepare_sample(sample, paired_read_length);
        normalize_relative_abundance(sample);
    }
    write_all_reads(design, output_dir, skip_reads);
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
