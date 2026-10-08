// SPDX-License-Identifier: GPL-2.0-only
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <memory>
#include <thread>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <cxxopts.hpp>
#include "RandomForest/LongReadSimulator.h"
#include "RandomForest/MetagenomeSimulator.h"
#include "RandomForest/ReadPipeline.h"
#include "IO/Bgzf.h"
#include "Utilities/BuildInfo.h"
#include "Utilities/Zstd.h"

namespace fs = std::filesystem;
using protal::sim::AbundanceDistribution;
using protal::sim::IlluminaOptions;
using protal::sim::MetagenomeSimulator;
using protal::sim::ProfileDesignOptions;

static std::optional<fs::path> find_executable_path(const std::string& arg0) {
    std::error_code ec;
    fs::path direct = arg0;
    if (!direct.empty() && fs::exists(direct, ec)) {
        fs::path canonical = fs::canonical(direct, ec);
        if (ec) {
            canonical = fs::absolute(direct);
        }
        return canonical;
    }
    const char* path_env = std::getenv("PATH");
    if (!path_env) {
        return std::nullopt;
    }
    std::string path_list(path_env);
    std::size_t start = 0;
    while (start <= path_list.size()) {
        auto end = path_list.find(':', start);
        std::string dir = path_list.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (!dir.empty()) {
            fs::path candidate = fs::path(dir) / arg0;
            if (fs::exists(candidate, ec)) {
                fs::path canonical = fs::canonical(candidate, ec);
                if (ec) {
                    canonical = fs::absolute(candidate);
                }
                return canonical;
            }
        }
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    return std::nullopt;
}

static std::string derive_prefix_from_r1(const fs::path& r1) {
    std::string base = r1.filename().string();
    auto strip_suffix = [&](const std::string& suf) {
        if (base.size() >= suf.size() && base.compare(base.size() - suf.size(), suf.size(), suf) == 0) {
            base.erase(base.size() - suf.size());
        }
    };
    strip_suffix(".gz");
    strip_suffix(".zst");
    strip_suffix(".fastq");
    strip_suffix(".fq");
    strip_suffix(".fasta");
    strip_suffix(".fa");
    strip_suffix(".FASTQ");
    strip_suffix(".FQ");
    if (base.size() > 3 && base.compare(base.size() - 3, 3, "_R1") == 0) {
        base.erase(base.size() - 3);
    } else if (base.size() > 2 && base.compare(base.size() - 2, 2, "_1") == 0) {
        base.erase(base.size() - 2);
    }
    return base;
}

static void write_protal_metafile(
    const std::vector<protal::sim::SampleOutput>& samples,
    const fs::path& output_dir,
    const fs::path& reads_dir,
    const fs::path& metafile_path,
    const std::vector<fs::path>& profile_truth_paths,
    const std::optional<fs::path>& output_dir_override) {
    if (samples.empty()) {
        return;
    }
    if (!metafile_path.parent_path().empty()) {
        fs::create_directories(metafile_path.parent_path());
    }
    auto to_abs = [](const fs::path& p) -> fs::path {
        std::error_code ec;
        auto c = fs::canonical(p, ec);
        return ec ? fs::absolute(p) : c;
    };
    fs::path out_abs = to_abs(output_dir_override ? *output_dir_override : output_dir);
    fs::path reads_abs = to_abs(reads_dir);

    std::ofstream out(metafile_path);
    if (!out) {
        throw std::runtime_error("Unable to write protal metafile: " + metafile_path.string());
    }
    out << "#OUTPUT_DIR\t" << out_abs.string() << "\n";
    out << "#INPUT_DIR\t" << reads_abs.string() << "\n";
    out << "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n";
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        auto first = sample.read1_path.filename().string();
        auto second = sample.read2_path.empty() ? std::string("-") : sample.read2_path.filename().string();  // first reads only
        auto prefix = derive_prefix_from_r1(sample.read1_path);
        fs::path truth_abs = to_abs(profile_truth_paths[i]);
        out << sample.sample_name << '\t' << first << '\t' << second << '\t' << prefix << ".sam.zst"
            << '\t' << prefix << '\t' << prefix << ".profile" << '\t' << truth_abs.string() << '\n';
    }
}

struct CliOptions {
    fs::path genome_table;
    fs::path output_dir;
    std::size_t samples{1};
    std::string sample_prefix{"sample"};
    std::uint64_t total_read_pairs{100'000};
    std::vector<std::uint64_t> total_read_pairs_per_sample;  // --total_read_pairs with several values: the samples' in turn
    std::size_t species_per_sample{10};
    std::size_t species_per_sample_min{0};
    std::vector<std::size_t> species_per_sample_list;  // --species_per_sample with several values: the samples' in turn
    AbundanceDistribution distribution{AbundanceDistribution::Lognormal};
    double alpha{2.0};
    int nb_r{5};
    double nb_p{0.5};
    double pln_mu{0.0};
    double pln_sigma{1.3};
    double abundance_floor{0.001};
    std::vector<double> pln_sigmas;  // --pln_sigma with several values: the samples' in turn
    std::string strain_probabilities;
    std::string include_species;
    std::string genus_counts;
    std::string taxon_counts;
    std::string congener_groups;
    IlluminaOptions illumina;
    std::uint64_t illumina_report{0};  // --illumina_report PAIRS: a profile's statistics instead of samples
    std::optional<std::uint64_t> seed;
    bool plot_png{false};
    bool test_mode{false};
    bool keep_tmp{false};
    protal::sim::ReadsCompression reads_compression{protal::sim::ReadsCompression::Bgzf};
    bool pick_random_demand_if_fail{false};
    int threads{1};
    std::optional<fs::path> protal_metafile_output_dir;
    fs::path strain_sharing_file;
    std::optional<fs::path> from_manifest;
    // --long_samples: long or Ultima reads of given communities instead (LongReadSimulator.h)
    fs::path long_samples, long_genomes, long_model, long_stats, long_templates, long_out;
    std::string long_setup;
    fs::path genome_store;    // --genome_store
    std::string tiles;        // --tiles LENGTH:STRIDE (tile_genomes, tile_fasta)
    fs::path tile_out, tile_exclude;
    fs::path tile_fasta;      // --tile_fasta: the records of a FASTA file instead of the genome table's genomes
    std::size_t tile_per_header{0};  // --tile_per_header: at most this many records of each header (0: all)
    bool plain_pipes{false};  // --plain_pipes
};

static std::vector<protal::sim::StrainSharingSpec> parse_strain_sharing_file(const fs::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open strain sharing file: " + path.string());

    std::vector<protal::sim::StrainSharingSpec> specs;
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;

        // tab-split
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            auto pos = line.find('\t', start);
            std::string tok = (pos == std::string::npos) ? line.substr(start) : line.substr(start, pos - start);
            // trim whitespace
            auto nb = [](unsigned char c) { return !std::isspace(c); };
            tok.erase(tok.begin(), std::find_if(tok.begin(), tok.end(), nb));
            tok.erase(std::find_if(tok.rbegin(), tok.rend(), nb).base(), tok.end());
            fields.push_back(tok);
            if (pos == std::string::npos) break;
            start = pos + 1;
        }

        if (fields.size() < 4) {
            throw std::runtime_error("strain_sharing_file line " + std::to_string(line_no) +
                ": expected at least 4 tab-separated columns (SPECIES, SAMPLE_FRACTION, N_STRAINS, MIN_OCCURRENCE)");
        }

        // Catch common mistake: decimal point in an integer column usually means a tab was
        // replaced by a space, shifting all subsequent columns by one.
        auto require_integer_field = [&](const std::string& val, const char* col_name) -> std::size_t {
            if (val.find('.') != std::string::npos) {
                throw std::runtime_error(
                    "strain_sharing_file line " + std::to_string(line_no) + ": " + col_name +
                    " must be an integer, got \"" + val + "\" (decimal point found — "
                    "check for spaces instead of tabs or a missing column)");
            }
            return static_cast<std::size_t>(std::stoull(val));
        };

        protal::sim::StrainSharingSpec spec;
        spec.species         = fields[0];
        spec.sample_fraction = std::stod(fields[1]);
        spec.n_strains       = require_integer_field(fields[2], "N_STRAINS");
        spec.min_occurrence  = require_integer_field(fields[3], "MIN_OCCURRENCE");
        spec.min_vcov        = (fields.size() >= 5 && !fields[4].empty()) ? std::stod(fields[4]) : 0.0;
        if (fields.size() >= 6 && !fields[5].empty()) {
            spec.conspecific_strain_probabilities = protal::sim::parse_strain_probabilities(fields[5]);
        }

        if (spec.sample_fraction < 0.0 || spec.sample_fraction > 1.0) {
            throw std::runtime_error("strain_sharing_file line " + std::to_string(line_no) +
                ": SAMPLE_FRACTION must be in [0,1]");
        }
        specs.push_back(std::move(spec));
    }
    return specs;
}

static cxxopts::Options build_cxxopts() {
    cxxopts::Options options("simulate_metagenomes", "Simulate shotgun metagenomes from a genome table");

    options.add_options("I/O")
        ("genome_table", "TSV with genome name, GTDB taxonomy, FASTA path (.gz ok)", cxxopts::value<std::string>())
        ("o,output_dir",  "Output directory for FASTQs and manifest", cxxopts::value<std::string>())
        ("from_manifest", "Replay a previous run from its manifest.tsv (combined or per-sample) instead of "
                          "designing a new community. All --distribution/--species_per_sample/--seed style "
                          "sampling options are ignored. Manifests with art_seed, run_seed and host_seed columns "
                          "(since 2026-10-07) reproduce the reads exactly, given the original read settings "
                          "(--read_length, --sequencer, ..., --host_folder, --host_pairs; checked against its "
                          "run_params.tsv); older ones reproduce the composition with fresh reads. "
                          "--genome_table is only needed if the manifest has no fasta_path column.",
                          cxxopts::value<std::string>());

    options.add_options("Sampling")
        ("n,samples",           "Number of metagenome samples", cxxopts::value<std::size_t>()->default_value("1"))
        ("sample_prefix",       "Prefix for sample names", cxxopts::value<std::string>()->default_value("sample"))
        ("total_read_pairs",    "Read pairs per sample; several comma-separated values are given to the samples in turn (sample 1 the first, sample 2 the second, ...), so that one run's samples differ in depth", cxxopts::value<std::string>()->default_value("100000"))
        ("species_per_sample",  "Number of species per sample; an inclusive range e.g. 20-80 (each sample's drawn from it); or several comma-separated numbers, given to the samples in turn (as --total_read_pairs)", cxxopts::value<std::string>()->default_value("10"))
        ("distribution",        "Abundance model: lognormal (a continuous long tail, never below --abundance_floor of its median), "
                                "power_law, negative_binomial or poisson_lognormal (Poisson counts + 1 of a lognormal mean: "
                                "at --pln_mu 0 about 40-45% of the species at the lowest weight)",
                                cxxopts::value<std::string>()->default_value("lognormal"))
        ("alpha",               "Power law alpha", cxxopts::value<double>()->default_value("2.0"))
        ("nb_r",                "Negative binomial r", cxxopts::value<int>()->default_value("5"))
        ("nb_p",                "Negative binomial p", cxxopts::value<double>()->default_value("0.5"))
        ("pln_mu",              "Lognormal (and Poisson-lognormal) mean (log-scale)", cxxopts::value<double>()->default_value("0.0"))
        ("pln_sigma",           "Lognormal (and Poisson-lognormal) sigma (log-scale); several comma-separated values are given to the samples in turn (sample 1 the first, sample 2 the second, ...), so that one run mixes abundance distributions", cxxopts::value<std::string>()->default_value("1.3"))
        ("abundance_floor",     "lognormal: the lowest weight of a species, as a share of the distribution's median (e^pln_mu), "
                                "so that the tail's species keep some abundance (every species gets a read pair or more in "
                                "any case); 0: none", cxxopts::value<double>()->default_value("0.001"))
        ("strains_per_species", "Probabilities for adding 2nd, 3rd, ... strains per species, e.g. \"0.4,0.2,0.1\"", cxxopts::value<std::string>()->default_value(""))
        ("include_species",     "Comma-separated species to force-include in each sample", cxxopts::value<std::string>()->default_value(""))
        ("genus",               "Comma-separated genus:count pairs, e.g. \"g__A:10,g__B:2\"", cxxopts::value<std::string>()->default_value(""))
        ("taxon",               "Comma-separated taxon:count pairs, e.g. \"d__Archaea:10\"", cxxopts::value<std::string>()->default_value(""))
        ("pick_random_demand_if_fail", "If --genus/--taxon demand more species than available, cap and fill randomly instead of failing")
        ("congener_groups",     "SHARE:MIN-MAX, e.g. 0.25:2-5: about SHARE of each sample's species come in groups of MIN "
                                "to MAX species of one genus, genera drawn per sample among those with MIN species or "
                                "more (after --genus and --taxon). Default: none, every species drawn from all species, "
                                "so that congeners hardly ever share a sample", cxxopts::value<std::string>()->default_value(""))
        ("strain_sharing_file", "TSV file for cross-sample strain sharing. Columns (tab-separated): "
                                "SPECIES  SAMPLE_FRACTION  N_STRAINS  MIN_OCCURRENCE  MIN_VCOV  CONSPECIFIC_STRAINS. "
                                "Lines starting with # are ignored. MIN_VCOV and CONSPECIFIC_STRAINS are optional. "
                                "CONSPECIFIC_STRAINS uses the same probability format as --strains_per_species.",
                                cxxopts::value<std::string>()->default_value(""))
        ("seed",                "RNG seed (default: random)", cxxopts::value<std::uint64_t>());

    options.add_options("Illumina")
        ("read_length",     "Read length", cxxopts::value<int>()->default_value("150"))
        ("fragment_mean",   "Fragment mean", cxxopts::value<int>()->default_value("350"))
        ("fragment_stdev",  "Fragment stdev", cxxopts::value<int>()->default_value("50"))
        ("sequencer",       "The instrument the reads model (IlluminaSimulator.h): HS20 (HiSeq 2000), HS25 (HiSeq 2500), "
                            "HSXt (HiSeq X Ten / 4000), NovaSeq (NovaSeq 6000), MSv3 (MiSeq v3)",
                            cxxopts::value<std::string>()->default_value("HS25"))
        ("mean_quality",    "Shift each read's qualities so that they average this (Phred), e.g. 35; default: the "
                            "instrument's", cxxopts::value<double>())
        ("host_folder",     "A host genome (scenarios.prepare_host's folder) for --host_pairs", cxxopts::value<std::string>())
        ("host_pairs",      "Read pairs of the host per sample, after its community's; several comma-separated values are "
                            "given to the samples in turn", cxxopts::value<std::string>()->default_value(""))
        ("first_reads_only", "Write only the _R1 files (the same reads as a run with both)")
        ("illumina_report", "Instead of samples: statistics of the --sequencer's reads at --read_length, --fragment_mean "
                            "and --mean_quality (per cycle mean quality, Q30 share, error rates), from this many pairs "
                            "of a random genome", cxxopts::value<std::uint64_t>())
        ("art_path",        "Unused: ART is no longer run (kept so that older commands fail clearly)", cxxopts::value<std::string>())
        ("extra_art_args",  "Unused: ART is no longer run", cxxopts::value<std::string>());

    options.add_options("Long reads")
        ("long_samples",    "Long (PacBio HiFi, Nanopore) or Ultima reads of given communities instead of a design: a TSV of "
                            "sample, out (its reads: .fq.zst or .fq.gz), bases, seed. Templates are drawn from the genomes "
                            "by weight (a gamma length of the setup, a uniform place where it fits on the genome's contigs, "
                            "cut at a contig's end only if longer than every contig; either strand) until each sample's "
                            "bases, and each made into one read named g<i>x_<n> (i: the genome's "
                            "place in the sample's list), written compressed as they are made; the same reads for any "
                            "number of threads", cxxopts::value<std::string>())
        ("long_genomes",    "With --long_samples: a TSV of sample, genome, fasta, weight (relative abundance x length), "
                            "host (1: fasta is a host folder of scenarios.prepare_host, read by memory map)",
                            cxxopts::value<std::string>())
        ("long_setup",      "With --long_samples: hifi:LENGTH_MEAN:LENGTH_SD:Q_SD (HiFi reads, quality by length), "
                            "ultima:LENGTH_MEAN:LENGTH_SD:Q_MEAN:Q_SD (flow-based reads of a mean quality) or "
                            "qshmm:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN[:SUB/INS/DEL] (pbsim3's qshmm model, Nanopore)",
                            cxxopts::value<std::string>())
        ("long_model",      "With a qshmm setup: pbsim3's model file (e.g. QSHMM-ONT-HQ.model)", cxxopts::value<std::string>())
        ("long_templates",  "Instead of --long_samples: one read of each sequence of this FASTA by --long_setup's model, "
                            "named after it, into --long_out (as hifi_reads.py and pbsim3 --strategy templ; --seed)",
                            cxxopts::value<std::string>())
        ("long_out",        "With --long_templates: the reads (.fq.zst or .fq.gz)", cxxopts::value<std::string>())
        ("long_stats",      "With --long_samples: write sample, reads, template_bases, read_bases, errors, rounds here "
                            "(default: standard output)", cxxopts::value<std::string>());

    options.add_options("Tiles")
        ("tiles",           "Instead of samples: error-free reads of LENGTH bases every STRIDE bases of every contig of every "
                            "genome of --genome_table (LENGTH:STRIDE, e.g. 150:500; forward strand, quality I), named "
                            "<genome>:<contig number>:<start>, into --tile_out, in table order: a scan of the genomes against a "
                            "database (scripts/foreign_rates.py)", cxxopts::value<std::string>())
        ("tile_out",        "With --tiles: the reads (.fq.zst or .fq.gz; a named pipe with --plain_pipes gets plain FASTQ)",
                            cxxopts::value<std::string>())
        ("tile_exclude",    "With --tiles: species to leave out, one per line in the first column (s__Genus species, as "
                            "heldout_species.txt)", cxxopts::value<std::string>())
        ("tile_fasta",      "With --tiles: the records of this FASTA file (plain or zstd; a database's full_reference.fna.zst: "
                            "every genome's marker genes, >taxid_geneid) instead of the genome table's genomes, the reads named "
                            "<header>:<record number>:<start>; <tile_out>.sources.tsv lists every header's records and tiles "
                            "(scripts/foreign_rates.py)", cxxopts::value<std::string>())
        ("tile_per_header", "With --tile_fasta: at most this many records of each header (0: all), so that a species with "
                            "thousands of genomes gives no more reads of a gene than one with this many",
                            cxxopts::value<std::size_t>()->default_value("0"));

    options.add_options("General")
        ("t,threads",       "Threads: the reads are made in work items of a few hundred kB on all of them, several samples at a time when one leaves threads idle; the same files for any number", cxxopts::value<int>()->default_value("1"))
        ("pigz_path",       "Unused: the reads are compressed in process (kept so that older commands still run)", cxxopts::value<std::string>()->default_value(""))
        ("reads_compression", "How the read files are written: bgzf (_R1.fq.gz, ISA-L level 1) or zstd (_R1.fq.zst, level 3: smaller files; protal reads both)", cxxopts::value<std::string>()->default_value("bgzf"))
        ("plain_pipes",     "Outputs that are named pipes get plain FASTQ, whatever their names say: nothing is compressed only "
                            "for the reader to inflate it again (protal takes plain FASTQ as it is); regular files stay "
                            "compressed")
        ("genome_store",    "A folder of genomes decoded once (2 bits a base, memory-mapped): a genome's FASTA is read and "
                            "parsed the first time any run needs it, its file written here, and every later sample, round "
                            "and run maps that file instead; a FASTA newer than its file is read again. The same reads as "
                            "without it", cxxopts::value<std::string>())
        ("protal_metafile", "Write a Protal meta file (output_dir/protal.meta) but set OUTPUT_DIR to <path>", cxxopts::value<std::string>())
        ("test",            "Generate profiles/manifests but skip read simulation (fast dry run)")
        ("keep_tmp",        "Unused: no temporary files are written (kept so that older commands still run)")
        ("plot_png",        "Generate barplot PNG of species abundances")
        ("v,version",       "Print the version (and the commit it was built from).")
        ("h,help",          "Print help.");

    return options;
}

static CliOptions parse_cli(int argc, char** argv) {
    auto cxx = build_cxxopts();

    if (argc <= 1) {
        std::cout << cxx.help({"I/O", "Sampling", "Illumina", "Long reads", "Tiles", "General"}) << std::endl;
        std::exit(0);
    }

    auto result = cxx.parse(argc, argv);

    if (result.count("help")) {
        std::cout << cxx.help({"I/O", "Sampling", "Illumina", "Long reads", "Tiles", "General"}) << std::endl;
        std::exit(0);
    }
    if (result.count("version")) {
        std::cout << "simulate_metagenomes " << protal::VersionText() << std::endl;
        std::exit(0);
    }

    if (result.count("illumina_report")) {  // a profile's statistics: no design, no genome table
        CliOptions opts;
        opts.illumina_report = result["illumina_report"].as<std::uint64_t>();
        opts.illumina.read_length = result["read_length"].as<int>();
        opts.illumina.fragment_mean = result["fragment_mean"].as<int>();
        opts.illumina.fragment_stdev = result["fragment_stdev"].as<int>();
        opts.illumina.sequencer = result["sequencer"].as<std::string>();
        if (result.count("mean_quality")) opts.illumina.mean_quality = result["mean_quality"].as<double>();
        if (result.count("seed")) opts.seed = result["seed"].as<std::uint64_t>();
        return opts;
    }
    if (result.count("tiles")) {  // tiles of the genome table's genomes, or of a FASTA file's records: no design
        CliOptions opts;
        opts.tiles = result["tiles"].as<std::string>();
        if (!result.count("tile_out") || (!result.count("genome_table") && !result.count("tile_fasta"))) {
            throw std::runtime_error("--tiles needs --tile_out and --genome_table or --tile_fasta");
        }
        opts.tile_out = result["tile_out"].as<std::string>();
        if (result.count("genome_table")) opts.genome_table = result["genome_table"].as<std::string>();
        if (result.count("tile_fasta")) opts.tile_fasta = result["tile_fasta"].as<std::string>();
        opts.tile_per_header = result["tile_per_header"].as<std::size_t>();
        if (result.count("tile_exclude")) opts.tile_exclude = result["tile_exclude"].as<std::string>();
        opts.threads = result["threads"].as<int>();
        if (result.count("genome_store")) opts.genome_store = result["genome_store"].as<std::string>();
        opts.plain_pipes = result.count("plain_pipes") > 0;
        return opts;
    }
    if (result.count("long_samples") || result.count("long_templates")) {  // long reads: no design, no genome table
        CliOptions opts;
        if (result.count("long_templates")) {
            opts.long_templates = result["long_templates"].as<std::string>();
            if (!result.count("long_out") || !result.count("long_setup")) {
                throw std::runtime_error("--long_templates needs --long_out and --long_setup");
            }
            opts.long_out = result["long_out"].as<std::string>();
            opts.seed = result.count("seed") ? result["seed"].as<std::uint64_t>() : 1;
        } else {
            opts.long_samples = result["long_samples"].as<std::string>();
            if (!result.count("long_genomes") || !result.count("long_setup")) {
                throw std::runtime_error("--long_samples needs --long_genomes and --long_setup");
            }
            opts.long_genomes = result["long_genomes"].as<std::string>();
        }
        opts.long_setup = result["long_setup"].as<std::string>();
        if (result.count("long_model")) opts.long_model = result["long_model"].as<std::string>();
        if (result.count("long_stats")) opts.long_stats = result["long_stats"].as<std::string>();
        opts.threads = result["threads"].as<int>();
        if (result.count("genome_store")) opts.genome_store = result["genome_store"].as<std::string>();
        opts.plain_pipes = result.count("plain_pipes") > 0;
        return opts;
    }

    const bool replay = result.count("from_manifest") > 0;
    if ((!result.count("genome_table") && !replay) || !result.count("output_dir")) {
        std::cerr << "Error: --output_dir is required, as is --genome_table unless --from_manifest is given.\n\n";
        std::cout << cxx.help({"I/O", "Sampling", "Illumina", "Long reads", "Tiles", "General"}) << std::endl;
        std::exit(1);
    }

    CliOptions opts;
    if (replay) {
        opts.from_manifest = result["from_manifest"].as<std::string>();
    }
    if (result.count("genome_table")) {
        opts.genome_table  = result["genome_table"].as<std::string>();
    }
    opts.output_dir        = result["output_dir"].as<std::string>();
    opts.samples           = result["samples"].as<std::size_t>();
    opts.sample_prefix     = result["sample_prefix"].as<std::string>();
    {
        const std::string pairs_arg = result["total_read_pairs"].as<std::string>();
        std::stringstream ss(pairs_arg);
        std::string item;
        while (std::getline(ss, item, ',')) {
            if (item.empty()) continue;
            std::size_t used = 0;
            const unsigned long long pairs = std::stoull(item, &used);
            if (used != item.size() || pairs == 0) {
                throw std::runtime_error("--total_read_pairs takes positive whole numbers, comma-separated: " + pairs_arg);
            }
            opts.total_read_pairs_per_sample.push_back(pairs);
        }
        if (opts.total_read_pairs_per_sample.empty()) throw std::runtime_error("--total_read_pairs needs a value");
        opts.total_read_pairs = opts.total_read_pairs_per_sample.front();
        if (opts.total_read_pairs_per_sample.size() == 1) opts.total_read_pairs_per_sample.clear();
    }
    {
        const std::string sps_arg = result["species_per_sample"].as<std::string>();
        const auto dash = sps_arg.find('-');
        if (sps_arg.find(',') != std::string::npos) {
            std::stringstream ss(sps_arg);
            std::string item;
            while (std::getline(ss, item, ',')) {
                std::size_t used = 0;
                const unsigned long long count = item.empty() ? 0 : std::stoull(item, &used);
                if (item.empty() || used != item.size() || count == 0) {
                    throw std::runtime_error("--species_per_sample takes N, MIN-MAX or positive whole numbers, "
                                             "comma-separated: " + sps_arg);
                }
                opts.species_per_sample_list.push_back(count);
            }
            opts.species_per_sample = *std::max_element(opts.species_per_sample_list.begin(),
                                                        opts.species_per_sample_list.end());
        } else if (dash != std::string::npos) {
            opts.species_per_sample_min = std::stoull(sps_arg.substr(0, dash));
            opts.species_per_sample     = std::stoull(sps_arg.substr(dash + 1));
            if (opts.species_per_sample_min > opts.species_per_sample) {
                throw std::runtime_error("--species_per_sample range min (" +
                    std::to_string(opts.species_per_sample_min) + ") exceeds max (" +
                    std::to_string(opts.species_per_sample) + ")");
            }
        } else {
            opts.species_per_sample = std::stoull(sps_arg);
        }
    }
    opts.alpha             = result["alpha"].as<double>();
    opts.nb_r              = result["nb_r"].as<int>();
    opts.nb_p              = result["nb_p"].as<double>();
    opts.pln_mu            = result["pln_mu"].as<double>();
    {
        const std::string sigma_arg = result["pln_sigma"].as<std::string>();
        std::stringstream ss(sigma_arg);
        std::string item;
        while (std::getline(ss, item, ',')) {
            if (item.empty()) continue;
            opts.pln_sigmas.push_back(std::stod(item));
        }
        if (opts.pln_sigmas.empty()) throw std::runtime_error("--pln_sigma needs a value");
        for (double s : opts.pln_sigmas) {
            if (!(s > 0)) throw std::runtime_error("--pln_sigma values must be positive: " + sigma_arg);
        }
        opts.pln_sigma = opts.pln_sigmas.front();
        if (opts.pln_sigmas.size() == 1) opts.pln_sigmas.clear();
    }
    opts.strain_probabilities = result["strains_per_species"].as<std::string>();
    opts.include_species   = result["include_species"].as<std::string>();
    opts.genus_counts      = result["genus"].as<std::string>();
    opts.taxon_counts      = result["taxon"].as<std::string>();
    opts.congener_groups   = result["congener_groups"].as<std::string>();
    {
        ProfileDesignOptions check;  // fails before any genome is read
        protal::sim::parse_congener_groups(opts.congener_groups, check);
    }
    opts.threads          = result["threads"].as<int>();
    opts.test_mode         = result.count("test") > 0;
    opts.keep_tmp          = result.count("keep_tmp") > 0;
    {
        std::string const compression = result["reads_compression"].as<std::string>();
        if (compression == "zstd") opts.reads_compression = protal::sim::ReadsCompression::Zstd;
        else if (compression != "bgzf") throw std::runtime_error("--reads_compression: bgzf or zstd, not " + compression);
    }
    opts.plot_png          = result.count("plot_png") > 0;
    opts.pick_random_demand_if_fail = result.count("pick_random_demand_if_fail") > 0;

    if (result.count("seed")) {
        opts.seed = result["seed"].as<std::uint64_t>();
    }
    if (result.count("protal_metafile")) {
        opts.protal_metafile_output_dir = result["protal_metafile"].as<std::string>();
    }
    {
        const std::string ssf = result["strain_sharing_file"].as<std::string>();
        if (!ssf.empty()) opts.strain_sharing_file = ssf;
    }

    const std::string dist_str = result["distribution"].as<std::string>();
    if (dist_str == "power_law") {
        opts.distribution = AbundanceDistribution::PowerLaw;
    } else if (dist_str == "negative_binomial") {
        opts.distribution = AbundanceDistribution::NegativeBinomial;
    } else if (dist_str == "poisson_lognormal") {
        opts.distribution = AbundanceDistribution::PoissonLognormal;
    } else if (dist_str == "lognormal") {
        opts.distribution = AbundanceDistribution::Lognormal;
    } else {
        throw std::runtime_error("Unknown distribution: " + dist_str);
    }
    opts.abundance_floor = result["abundance_floor"].as<double>();
    if (!(opts.abundance_floor >= 0 && opts.abundance_floor < 1)) {
        throw std::runtime_error("--abundance_floor is a share of the median, from 0 to below 1");
    }

    if (result.count("art_path") || result.count("extra_art_args")) {
        throw std::runtime_error("--art_path and --extra_art_args: ART is no longer run; the reads are made in process "
                                 "(--sequencer, --mean_quality)");
    }
    opts.illumina.read_length   = result["read_length"].as<int>();
    opts.illumina.fragment_mean = result["fragment_mean"].as<int>();
    opts.illumina.fragment_stdev = result["fragment_stdev"].as<int>();
    opts.illumina.sequencer     = result["sequencer"].as<std::string>();
    protal::sim::IlluminaProfile::Named(opts.illumina.sequencer);  // fails before any genome is read
    if (result.count("mean_quality")) opts.illumina.mean_quality = result["mean_quality"].as<double>();
    if (result.count("host_folder")) opts.illumina.host_folder = result["host_folder"].as<std::string>();
    {
        const std::string pairs_arg = result["host_pairs"].as<std::string>();
        std::stringstream ss(pairs_arg);
        std::string item;
        while (std::getline(ss, item, ',')) {
            if (item.empty()) continue;
            std::size_t used = 0;
            const unsigned long long pairs = std::stoull(item, &used);
            if (used != item.size()) throw std::runtime_error("--host_pairs takes whole numbers, comma-separated: " + pairs_arg);
            opts.illumina.host_pairs.push_back(pairs);
        }
        bool const any = std::any_of(opts.illumina.host_pairs.begin(), opts.illumina.host_pairs.end(),
                                     [](auto p) { return p > 0; });
        if (any && opts.illumina.host_folder.empty()) throw std::runtime_error("--host_pairs needs --host_folder");
    }
    opts.illumina.first_reads_only = result.count("first_reads_only") > 0;
    if (result.count("genome_store")) opts.illumina.genome_store = result["genome_store"].as<std::string>();
    opts.illumina.plain_pipes = result.count("plain_pipes") > 0;

    return opts;
}

// --host_pairs as given: comma-separated.
static std::string HostPairsText(const std::vector<std::uint64_t>& pairs) {
    std::string text;
    for (std::size_t i = 0; i < pairs.size(); ++i) text += (i ? "," : "") + std::to_string(pairs[i]);
    return text;
}

// A manifest pins the community; this pins how it was produced. Written on every run
// so that neither the seed nor the read settings live only in shell history.
static void write_run_params(
    const CliOptions& cli, std::uint64_t seed, int argc, char** argv, const fs::path& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Unable to write run parameters: " + path.string());
    }
    out << "key\tvalue\n";
    std::string command;
    for (int i = 0; i < argc; ++i) {
        command += (i ? " " : "");
        command += argv[i];
    }
    out << "command\t" << command << '\n';
    out << "mode\t" << (cli.from_manifest ? "replay" : "design") << '\n';
    if (cli.from_manifest) {
        out << "from_manifest\t" << cli.from_manifest->string() << '\n';
    }
    out << "seed\t" << seed << '\n';
    out << "genome_table\t" << cli.genome_table.string() << '\n';
    out << "read_length\t" << cli.illumina.read_length << '\n';
    out << "fragment_mean\t" << cli.illumina.fragment_mean << '\n';
    out << "fragment_stdev\t" << cli.illumina.fragment_stdev << '\n';
    out << "sequencer\t" << cli.illumina.sequencer << '\n';
    out << "mean_quality\t" << (cli.illumina.mean_quality ? std::to_string(*cli.illumina.mean_quality) : "") << '\n';
    out << "host_folder\t" << cli.illumina.host_folder.string() << '\n';
    out << "host_pairs\t" << HostPairsText(cli.illumina.host_pairs) << '\n';
    out << "reads\t" << "simulate_metagenomes IlluminaSimulator (no ART)" << '\n';
}

static std::vector<protal::sim::SampleOutput> design_and_simulate(
    CliOptions& cli, std::vector<protal::sim::GenomeRecord> genomes) {
    ProfileDesignOptions profile{};
    profile.species_per_sample     = cli.species_per_sample;
    profile.species_per_sample_min = cli.species_per_sample_min;
    profile.species_per_sample_list = cli.species_per_sample_list;
    profile.distribution = cli.distribution;
    profile.powerlaw_alpha = cli.alpha;
    profile.negative_binomial_r = cli.nb_r;
    profile.negative_binomial_p = cli.nb_p;
    profile.pln_mu = cli.pln_mu;
    profile.pln_sigma = cli.pln_sigma;
    profile.pln_sigmas = cli.pln_sigmas;
    profile.abundance_floor = cli.abundance_floor;
    profile.strain_probabilities = protal::sim::parse_strain_probabilities(cli.strain_probabilities);
    profile.include_species = protal::sim::parse_species_list(cli.include_species);
    profile.genus_species_counts = protal::sim::parse_genus_selection(cli.genus_counts);
    profile.taxon_species_counts = protal::sim::parse_taxon_selection(cli.taxon_counts);
    protal::sim::parse_congener_groups(cli.congener_groups, profile);
    profile.total_read_pairs = cli.total_read_pairs;
    profile.total_read_pairs_per_sample = cli.total_read_pairs_per_sample;
    profile.pick_random_demand_if_fail = cli.pick_random_demand_if_fail;
    if (!cli.strain_sharing_file.empty()) {
        profile.strain_sharing = parse_strain_sharing_file(cli.strain_sharing_file);
        std::cerr << "Loaded " << profile.strain_sharing.size()
                  << " strain sharing spec(s) from " << cli.strain_sharing_file << '\n';
    }

    cli.illumina.threads = std::max(1, cli.threads);
    MetagenomeSimulator simulator(std::move(genomes), cli.illumina, *cli.seed);
    simulator.set_reads_compression(cli.reads_compression);

    return simulator.simulate_samples(
        profile, cli.samples, cli.sample_prefix, cli.output_dir, cli.test_mode, cli.keep_tmp);
}

// The manifest records the community and the reads' seeds but not the read settings, which come from
// this command line. Compare them with the run_params.tsv of the original run (next to manifest.tsv,
// or one level up for manifests/<sample>.tsv) and warn about every difference: with other settings
// the replayed reads differ.
static void check_replay_read_settings(const CliOptions& cli) {
    const fs::path dir = cli.from_manifest->parent_path();
    fs::path params = dir / "run_params.tsv";
    if (!fs::exists(params)) params = dir.parent_path() / "run_params.tsv";
    if (!fs::exists(params)) {
        std::cerr << "Note: no run_params.tsv next to the manifest, so the read settings cannot be checked "
                     "against the original run.\n";
        return;
    }

    std::unordered_map<std::string, std::string> recorded;
    std::ifstream in(params);
    std::string line;
    while (std::getline(in, line)) {
        auto tab = line.find('\t');
        if (tab != std::string::npos) recorded[line.substr(0, tab)] = line.substr(tab + 1);
    }

    const std::vector<std::pair<std::string, std::string>> current = {
        {"read_length", std::to_string(cli.illumina.read_length)},
        {"fragment_mean", std::to_string(cli.illumina.fragment_mean)},
        {"fragment_stdev", std::to_string(cli.illumina.fragment_stdev)},
        {"sequencer", cli.illumina.sequencer},
        {"mean_quality", cli.illumina.mean_quality ? std::to_string(*cli.illumina.mean_quality) : ""},
        {"host_folder", cli.illumina.host_folder.string()},
        {"host_pairs", HostPairsText(cli.illumina.host_pairs)},
    };
    for (const auto& [key, value] : current) {
        auto it = recorded.find(key);
        if (it != recorded.end() && it->second != value) {
            std::cerr << "Warning: --" << key << " is '" << value << "' but was '" << it->second
                      << "' in the original run (" << params.string() << "); replayed reads will differ.\n";
        }
    }
}

static std::vector<protal::sim::SampleOutput> replay_from_manifest(
    CliOptions& cli, const std::vector<protal::sim::GenomeRecord>& genomes) {
    check_replay_read_settings(cli);
    auto design = protal::sim::read_manifest(*cli.from_manifest);
    protal::sim::resolve_manifest_fasta_paths(design, genomes);

    std::size_t rows = 0;
    std::size_t seeded = 0;
    std::size_t samples_seeded = 0;
    for (const auto& sample : design) {
        samples_seeded += sample.run_seed && sample.host_seed;
        for (const auto& assignment : sample.assignments) {
            ++rows;
            seeded += assignment.art_seed.has_value();
        }
    }
    std::cerr << "Replaying " << design.size() << " sample(s), " << rows << " genome assignment(s) from "
              << cli.from_manifest->string() << '\n';
    if (seeded == rows && samples_seeded == design.size()) {
        std::cerr << "Every row carries an art_seed and every sample its run_seed and host_seed: the reads are "
                     "reproduced exactly (with the original read settings).\n";
    } else {
        if (seeded < rows) {
            std::cerr << "Warning: " << (rows - seeded) << " of " << rows
                      << " rows have no art_seed. Composition and depth are reproduced exactly, but those "
                         "reads are fresh realizations.\n";
        }
        if (samples_seeded < design.size()) {
            std::cerr << "Warning: " << (design.size() - samples_seeded) << " of " << design.size()
                      << " samples have no run_seed or host_seed (a manifest of before 2026-10-07): their qualities, "
                         "and so their errors, and their host reads are this run's (--seed), not the original's.\n";
        }
    }

    // The manifest does not store the read length, but it is recoverable from any row:
    // vertical_coverage = read_pairs * 2 * read_length / genome_length. Replaying at a
    // different --read_length silently changes every depth, so check rather than trust.
    for (const auto& sample : design) {
        for (const auto& assignment : sample.assignments) {
            if (assignment.read_pairs == 0 || assignment.vertical_coverage <= 0.0) {
                continue;
            }
            const double implied = assignment.vertical_coverage * static_cast<double>(assignment.genome_length) /
                                   (2.0 * static_cast<double>(assignment.read_pairs));
            if (std::abs(implied - static_cast<double>(cli.illumina.read_length)) > 1.0) {
                std::cerr << "Warning: manifest implies read_length ~" << static_cast<int>(implied + 0.5)
                          << " but --read_length is " << cli.illumina.read_length
                          << "; depths will not match the original run.\n";
            }
            break;
        }
        break;
    }

    // The profile designer is unused on a replay, but the simulator still wants a
    // genome set; the manifest's own genomes are it.
    std::vector<protal::sim::GenomeRecord> replay_genomes;
    std::unordered_set<std::string> seen_genomes;
    for (const auto& sample : design) {
        for (const auto& assignment : sample.assignments) {
            if (seen_genomes.insert(assignment.genome.name).second) {
                replay_genomes.push_back(assignment.genome);
            }
        }
    }

    // The seed makes only what the manifest does not record: art_seed of rows without one, run_seed and host_seed of
    // samples without them.
    cli.illumina.threads = std::max(1, cli.threads);
    MetagenomeSimulator simulator(std::move(replay_genomes), cli.illumina, *cli.seed);
    simulator.set_reads_compression(cli.reads_compression);

    return simulator.replay_samples(std::move(design), cli.output_dir, cli.test_mode, cli.keep_tmp);
}

// --illumina_report PAIRS: per read, the mean written quality, the share at Q30 or more, the errors per base and the
// written qualities' shares; then the mean written quality by cycle against the profile's (published) curve.
static int illumina_report(const CliOptions& cli) {
    protal::sim::IlluminaSetup setup;
    setup.profile = protal::sim::IlluminaProfile::Named(cli.illumina.sequencer);
    setup.read_length = cli.illumina.read_length;
    setup.fragment_mean = cli.illumina.fragment_mean;
    setup.fragment_sd = cli.illumina.fragment_stdev;
    setup.mean_quality = cli.illumina.mean_quality;
    protal::sim::IlluminaModel const model(setup);
    auto const report = protal::sim::ReportProfile(setup, cli.illumina_report, cli.seed.value_or(1));
    std::cout << "# " << setup.profile.name << " (" << setup.profile.description << "), " << setup.read_length
              << " bp, fragments " << setup.fragment_mean << " (SD " << setup.fragment_sd << "), "
              << cli.illumina_report << " pairs\n";
    std::cout << "read\tmean_Q\tQ30\tsubstitutions\tinsertions\tdeletions\tNs\twritten qualities (share)\n";
    for (int r = 0; r < 2; ++r) {
        double mean = 0;
        for (double q : report.mean_by_cycle[r]) mean += q;
        mean /= static_cast<double>(report.mean_by_cycle[r].size());
        std::cout << "R" << r + 1 << '\t' << mean << '\t' << report.q30[r] << '\t' << report.substitutions[r] << '\t'
                  << report.insertions[r] << '\t' << report.deletions[r] << '\t' << report.ns[r] << '\t';
        for (int q = 0; q < 94; ++q) {
            if (report.quality_share[r][q] >= 0.0005) std::cout << 'Q' << q << ':' << report.quality_share[r][q] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << "cycle\tR1\tR2\tR1_published\tR2_published\n";
    for (int c = 0; c < setup.read_length; ++c) {
        std::cout << c + 1 << '\t' << report.mean_by_cycle[0][c] << '\t' << report.mean_by_cycle[1][c] << '\t'
                  << model.TargetMean(0)[c] << '\t' << model.TargetMean(1)[c] << '\n';
    }
    return 0;
}

// --tiles LENGTH:STRIDE: error-free reads of LENGTH bases every STRIDE bases of every contig of every genome of the table
// (forward strand, quality 'I'), named <genome>:<contig number>:<start> (0-based), into --tile_out, the genomes in table
// order; those of the species --tile_exclude names (first column: s__Genus species, as heldout_species.txt) left out. A scan
// of the genomes at hand against a database (scripts/foreign_rates.py): which gene copies other species' reads reach.
// The genomes are read on --threads threads, a batch at a time, and written in order: the same file for any number.
static std::pair<std::uint64_t, std::uint64_t> parse_tiles(std::string const& text) {
    auto const colon = text.find(':');
    std::uint64_t length = 0, stride = 0;
    try {
        if (colon == std::string::npos) throw std::invalid_argument("no colon");
        length = std::stoull(text.substr(0, colon));
        stride = std::stoull(text.substr(colon + 1));
    } catch (std::exception const&) {
        throw std::runtime_error("--tiles takes LENGTH:STRIDE (e.g. 150:500), not '" + text + "'");
    }
    if (length == 0 || stride == 0) throw std::runtime_error("--tiles: LENGTH and STRIDE must be positive");
    return { length, stride };
}

// The reads of a tiling, packed for --tile_out (.fq.zst, .fq.gz or plain) and written in order; open until Close().
class TileWriter {
public:
    TileWriter(fs::path const& out, bool plain_pipes, std::size_t threads)
        : m_out(out), m_pipe(fs::is_fifo(out)),
          m_packing(m_pipe && plain_pipes ? protal::sim::pipeline::Packing::Plain : protal::sim::pipeline::PackingOf(out)),
          m_target(m_pipe ? out : fs::path(out.string() + ".partial")),
          m_file(std::fopen(m_target.c_str(), "wb"), &std::fclose), m_threads(std::max<std::size_t>(1, threads)) {
        if (!m_file) throw std::runtime_error("cannot write " + m_target.string());
    }

    protal::sim::pipeline::Packing Packing() const { return m_packing; }

    // Packs these chunks of FASTQ text on the threads and writes them in order.
    void Write(std::vector<std::string>& chunks) {
        std::vector<std::string> packed(chunks.size());
        std::atomic<std::size_t> next{0};
        auto work = [&] {
            for (std::size_t c; (c = next++) < chunks.size();) packed[c] = protal::sim::pipeline::Pack(chunks[c], m_packing);
        };
        std::vector<std::thread> pool;
        for (std::size_t t = 1; t < std::min(m_threads, chunks.size()); t++) pool.emplace_back(work);
        work();
        for (auto& t : pool) t.join();
        for (auto const& bytes : packed) {
            if (std::fwrite(bytes.data(), 1, bytes.size(), m_file.get()) != bytes.size()) {
                throw std::runtime_error("writing " + m_target.string() + " failed");
            }
        }
        chunks.clear();
    }

    void Close() {
        if (m_packing == protal::sim::pipeline::Packing::Bgzf &&
            std::fwrite(protal::bgzf::kEof, 1, sizeof(protal::bgzf::kEof), m_file.get()) != sizeof(protal::bgzf::kEof)) {
            throw std::runtime_error("writing " + m_target.string() + " failed");
        }
        if (std::fclose(m_file.release()) != 0) throw std::runtime_error("writing " + m_target.string() + " failed");
        if (!m_pipe) fs::rename(m_target, m_out);
    }

private:
    fs::path m_out;
    bool m_pipe;
    protal::sim::pipeline::Packing m_packing;
    fs::path m_target;
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> m_file;
    std::size_t m_threads;
};

// --tiles with --tile_fasta: error-free reads of LENGTH bases every STRIDE bases of every record of a FASTA file (plain or
// zstd, read sequentially; a database's full_reference.fna.zst: every genome's marker genes under >taxid_geneid), named
// <header>:<record>:<start> (the header's first word, the record's number in the file from 0, the tile's 0-based start),
// into --tile_out. --tile_per_header N keeps the first N records of each header within one gene's run of records (the
// full reference is gene by gene, the gene the header's part after '_'): a species with thousands of genomes then gives
// no more reads of a gene than one with N. <tile_out>.sources.tsv lists, per header in file order (a header again when
// its gene's records are not contiguous), the records seen and kept and the tiles made, so that scripts/foreign_rates.py
// gives every gene copy a row, reached by a read or not. The records are read on one thread, the reads packed on
// --threads threads a batch at a time and written in order: the same file for any number of threads.
static int tile_fasta(const CliOptions& cli) {
    auto const [length, stride] = parse_tiles(cli.tiles);
    protal::zstd::InputFile input(cli.tile_fasta.string());
    if (!input.IsOpen()) throw std::runtime_error("cannot read " + cli.tile_fasta.string());
    std::istream& is = input.Stream();
    std::size_t const threads = static_cast<std::size_t>(std::max(1, cli.threads));
    TileWriter writer(cli.tile_out, cli.plain_pipes, threads);
    fs::path const sources_path(cli.tile_out.string() + ".sources.tsv");
    std::ofstream sources(sources_path.string() + ".partial");
    if (!sources) throw std::runtime_error("cannot write " + sources_path.string());
    sources << "header\trecords\tkept\ttiles\n";
    struct Source { std::uint64_t seen = 0, kept = 0, tiles = 0; };
    std::unordered_map<std::string, Source> current;  // the headers of the gene being read
    std::vector<std::string> order;                   // in order of first appearance
    std::string gene;                                 // the current records' gene (the header after '_'; "" without one)
    std::uint64_t headers = 0, records = 0, kept = 0, reads = 0;
    auto flush_sources = [&] {
        for (auto const& key : order) {
            auto const& s = current.at(key);
            sources << key << '\t' << s.seen << '\t' << s.kept << '\t' << s.tiles << '\n';
        }
        headers += order.size();
        current.clear();
        order.clear();
    };
    std::string const quality(length, 'I');
    std::size_t const chunk_bytes = std::size_t{4} << 20;
    std::vector<std::string> chunks;
    std::string chunk;
    auto take = [&](std::string const& key, std::string const& seq) {
        auto const us = key.find('_');
        std::string group = us == std::string::npos ? std::string() : key.substr(us + 1);
        if (group != gene) {
            flush_sources();
            gene = std::move(group);
        }
        auto& s = current[key];
        if (s.seen++ == 0) order.push_back(key);
        if (cli.tile_per_header > 0 && s.kept >= cli.tile_per_header) return;
        s.kept++;
        kept++;
        for (std::uint64_t pos = 0; pos + length <= seq.size(); pos += stride) {
            chunk += '@';
            chunk += key;
            chunk += ':';
            chunk += std::to_string(records);
            chunk += ':';
            chunk += std::to_string(pos);
            chunk += '\n';
            chunk.append(seq, pos, length);
            chunk += "\n+\n";
            chunk += quality;
            chunk += '\n';
            s.tiles++;
            reads++;
        }
        if (chunk.size() >= chunk_bytes) {
            chunks.push_back(std::move(chunk));
            chunk.clear();
            if (chunks.size() >= threads * 4) writer.Write(chunks);
        }
    };
    std::string line, key, seq;
    bool in_record = false;
    while (std::getline(is, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            if (in_record) {
                take(key, seq);
                records++;
            }
            auto const end = line.find_first_of(" \t", 1);
            key = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
            seq.clear();
            in_record = true;
        } else if (in_record) {
            seq += line;
        }
    }
    if (is.bad() || (!is.eof() && is.fail())) throw std::runtime_error("cannot read " + cli.tile_fasta.string() + " to its end");
    if (in_record) {
        take(key, seq);
        records++;
    }
    flush_sources();
    if (!chunk.empty()) chunks.push_back(std::move(chunk));
    if (!chunks.empty()) writer.Write(chunks);
    writer.Close();
    sources.close();
    if (!sources) throw std::runtime_error("writing " + sources_path.string() + " failed");
    fs::rename(sources_path.string() + ".partial", sources_path);
    std::cout << "Tiles: " << reads << " reads of " << length << " bases every " << stride << " bases of " << kept << " of "
              << records << " records (" << headers << " headers"
              << (cli.tile_per_header > 0 ? ", at most " + std::to_string(cli.tile_per_header) + " records each" : "")
              << ") of " << cli.tile_fasta.string() << " into " << cli.tile_out.string() << "; their headers' records and tiles: "
              << sources_path.string() << std::endl;
    return 0;
}

static int tile_genomes(const CliOptions& cli) {
    auto const [length, stride] = parse_tiles(cli.tiles);
    std::unordered_set<std::string> excluded;
    if (!cli.tile_exclude.empty()) {
        std::ifstream in(cli.tile_exclude);
        if (!in) throw std::runtime_error("cannot read " + cli.tile_exclude.string());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::string name = line.substr(0, line.find('\t'));
            while (!name.empty() && (name.back() == '\r' || name.back() == ' ')) name.pop_back();
            if (!name.empty()) excluded.insert(name);
        }
    }
    auto species_of = [](std::string const& taxonomy) {
        auto const at = taxonomy.rfind(';');
        std::string s = at == std::string::npos ? taxonomy : taxonomy.substr(at + 1);
        while (!s.empty() && s.front() == ' ') s.erase(s.begin());
        return s;
    };
    auto const genomes = protal::sim::read_genome_table(cli.genome_table);
    std::vector<std::size_t> chosen;
    for (std::size_t i = 0; i < genomes.size(); i++) {
        if (!excluded.contains(species_of(genomes[i].taxonomy))) chosen.push_back(i);
    }
    bool const pipe = fs::is_fifo(cli.tile_out);
    auto const packing = pipe && cli.plain_pipes ? protal::sim::pipeline::Packing::Plain
                                                 : protal::sim::pipeline::PackingOf(cli.tile_out);
    fs::path const target = pipe ? cli.tile_out : fs::path(cli.tile_out.string() + ".partial");
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> file(std::fopen(target.c_str(), "wb"), &std::fclose);
    if (!file) throw std::runtime_error("cannot write " + target.string());
    std::string const quality(length, 'I');
    std::size_t const threads = static_cast<std::size_t>(std::max(1, cli.threads));
    std::size_t const batch = threads * 4;
    std::uint64_t reads = 0, skipped = 0;
    for (std::size_t from = 0; from < chosen.size(); from += batch) {
        std::size_t const to = std::min(chosen.size(), from + batch);
        std::vector<std::string> packed(to - from);
        std::vector<std::uint64_t> made(to - from, 0);
        std::vector<std::string> errors(to - from);
        std::atomic<std::size_t> next{from};
        auto work = [&] {
            std::string fastq, seq;
            for (std::size_t c; (c = next++) < to;) {
                auto const& genome = genomes[chosen[c]];
                std::shared_ptr<protal::sim::Contigs const> contigs;
                try {
                    contigs = protal::sim::LoadContigs(genome.name, genome.fasta_path, static_cast<std::uint32_t>(length), cli.genome_store);
                } catch (std::exception const& ex) {
                    errors[c - from] = ex.what();
                    continue;
                }
                fastq.clear();
                for (std::size_t k = 0; k < contigs->names.size(); k++) {
                    for (std::uint64_t pos = 0; pos + length <= contigs->lengths[k]; pos += stride) {
                        contigs->Extract(k, pos, length, seq);
                        fastq += '@' + genome.name + ':' + std::to_string(k) + ':' + std::to_string(pos) + '\n' + seq + "\n+\n" +
                                 quality + '\n';
                        made[c - from]++;
                        if (fastq.size() >= (4u << 20)) {
                            packed[c - from] += protal::sim::pipeline::Pack(fastq, packing);
                            fastq.clear();
                        }
                    }
                }
                if (!fastq.empty()) packed[c - from] += protal::sim::pipeline::Pack(fastq, packing);
            }
        };
        std::vector<std::thread> pool;
        for (std::size_t t = 1; t < threads; t++) pool.emplace_back(work);
        work();
        for (auto& t : pool) t.join();
        for (std::size_t c = from; c < to; c++) {
            if (!errors[c - from].empty()) {
                std::cerr << "--tiles: " << genomes[chosen[c]].name << " left out: " << errors[c - from] << '\n';
                skipped++;
                continue;
            }
            auto const& bytes = packed[c - from];
            if (std::fwrite(bytes.data(), 1, bytes.size(), file.get()) != bytes.size()) {
                throw std::runtime_error("writing " + target.string() + " failed");
            }
            reads += made[c - from];
        }
    }
    if (packing == protal::sim::pipeline::Packing::Bgzf &&
        std::fwrite(protal::bgzf::kEof, 1, sizeof(protal::bgzf::kEof), file.get()) != sizeof(protal::bgzf::kEof)) {
        throw std::runtime_error("writing " + target.string() + " failed");
    }
    if (std::fclose(file.release()) != 0) throw std::runtime_error("writing " + target.string() + " failed");
    if (!pipe) fs::rename(target, cli.tile_out);
    std::cout << "Tiles: " << reads << " reads of " << length << " bases every " << stride << " bases of "
              << chosen.size() - skipped << " genomes (" << genomes.size() - chosen.size() << " of excluded species left out, "
              << skipped << " unreadable) into " << cli.tile_out.string() << std::endl;
    return 0;
}

// --long_samples: the samples' long or Ultima reads (LongReadSimulator), and a line of counts per sample.
static int simulate_long_reads(const CliOptions& cli) {
    protal::sim::LongReadOptions options;
    options.setup = protal::sim::LongReadSetup::Parse(cli.long_setup);
    options.model = cli.long_model;
    options.threads = std::max(1, cli.threads);
    options.genome_store = cli.genome_store;
    options.plain_pipes = cli.plain_pipes;
    if (options.setup.method == protal::sim::LongReadSetup::Method::Qshmm && options.model.empty()) {
        throw std::runtime_error("a qshmm setup needs pbsim3's model file: --long_model");
    }
    std::vector<protal::sim::LongSampleResult> results;
    if (!cli.long_templates.empty()) {
        results.push_back(protal::sim::MutateTemplates(cli.long_templates, cli.long_out, options, cli.seed.value_or(1)));
    } else {
        results = protal::sim::SimulateLongReads(protal::sim::ReadLongSamples(cli.long_samples, cli.long_genomes), options);
    }
    std::ofstream file;
    if (!cli.long_stats.empty()) {
        file.open(cli.long_stats);
        if (!file) throw std::runtime_error("cannot write " + cli.long_stats.string());
    }
    std::ostream& out = cli.long_stats.empty() ? std::cout : file;
    out << "sample\treads\ttemplate_bases\tread_bases\terrors\trounds\n";
    for (const auto& r : results) {
        out << r.name << '\t' << r.reads << '\t' << r.template_bases << '\t' << r.read_bases << '\t' << r.errors << '\t'
            << r.rounds << '\n';
    }
    out.flush();
    if (!out) throw std::runtime_error("writing the counts failed");
    return 0;
}

int main(int argc, char** argv) {
    CliOptions cli;
    try {
        cli = parse_cli(argc, argv);
    } catch (const std::exception& ex) {
        std::cerr << "Error parsing arguments: " << ex.what() << '\n';
        return 1;
    }

    if (cli.illumina_report > 0) {
        try {
            return illumina_report(cli);
        } catch (const std::exception& ex) {
            std::cerr << "Report failed: " << ex.what() << '\n';
            return 1;
        }
    }
    if (!cli.tiles.empty()) {
        try {
            return cli.tile_fasta.empty() ? tile_genomes(cli) : tile_fasta(cli);
        } catch (const std::exception& ex) {
            std::cerr << "Tiling failed: " << ex.what() << '\n';
            return 1;
        }
    }
    if (!cli.long_samples.empty() || !cli.long_templates.empty()) {
        try {
            return simulate_long_reads(cli);
        } catch (const std::exception& ex) {
            std::cerr << "Simulation failed: " << ex.what() << '\n';
            return 1;
        }
    }

    try {
        std::vector<protal::sim::GenomeRecord> genomes;
        if (!cli.genome_table.empty()) {
            genomes = protal::sim::read_genome_table(cli.genome_table);
        }
        // Always resolve the seed here, so an unseeded run still records the seed it
        // actually used and can be repeated.
        if (!cli.seed) {
            cli.seed = std::random_device{}();
            std::cerr << "Using random seed: " << *cli.seed << '\n';
        } else {
            std::cerr << "Using fixed seed: " << *cli.seed << '\n';
        }

        auto samples = cli.from_manifest ? replay_from_manifest(cli, genomes)
                                         : design_and_simulate(cli, std::move(genomes));

        fs::create_directories(cli.output_dir);
        write_run_params(cli, *cli.seed, argc, argv, cli.output_dir / "run_params.tsv");

        auto combined_manifest_path = cli.output_dir / "manifest.tsv";
        protal::sim::write_combined_manifest(samples, combined_manifest_path);

        fs::path manifests_dir = cli.output_dir / "manifests";
        fs::create_directories(manifests_dir);
        for (const auto& sample : samples) {
            protal::sim::write_sample_manifest(sample, manifests_dir / (sample.sample_name + ".tsv"));
        }

        protal::sim::write_abundance_matrix(samples, cli.output_dir / "abundance_matrix.tsv");

        std::cout << "Wrote " << samples.size() << " samples to " << cli.output_dir << '\n';

        if (cli.plot_png) {
            fs::path plots_dir = cli.output_dir / "plots";
            fs::create_directories(plots_dir);
            fs::path script_path = fs::path("scripts/plot_abundances.R");
            if (!fs::exists(script_path)) {
                // Try locating next to the executable.
                if (auto exe_path = find_executable_path(argv[0])) {
                    script_path = exe_path->parent_path() / "plot_abundances.R";
                }
            }
            if (!fs::exists(script_path)) {
                throw std::runtime_error("plot_abundances.R not found; please run from repo root");
            }
            for (const auto& sample : samples) {
                fs::path manifest_path = manifests_dir / (sample.sample_name + ".tsv");
                fs::path plot_out = plots_dir / (sample.sample_name + ".png");
                std::stringstream cmd;
                cmd << "Rscript \"" << script_path.string() << "\" \"" << manifest_path.string() << "\" \""
                    << plot_out.string() << "\"";
                int rc = std::system(cmd.str().c_str());
                if (rc != 0) {
                    throw std::runtime_error("Plot generation failed for " + sample.sample_name + " with exit code " +
                                             std::to_string(rc));
                }
                std::cout << "Generated plot: " << plot_out << '\n';
            }
        }

        if (cli.protal_metafile_output_dir) {
            fs::path reads_dir = cli.output_dir / "reads";
            fs::path goldstd_dir = cli.output_dir / "protal_goldstd";
            fs::create_directories(goldstd_dir);

            std::vector<fs::path> truth_paths;
            truth_paths.reserve(samples.size());
            for (const auto& sample : samples) {
                fs::path truth_path = goldstd_dir / (sample.sample_name + ".profile_truth");
                truth_paths.push_back(truth_path);
                std::unordered_set<std::string> seen;
                std::ofstream truth_out(truth_path);
                if (!truth_out) {
                    throw std::runtime_error("Unable to write profile truth: " + truth_path.string());
                }
                for (const auto& assignment : sample.assignments) {
                    if (seen.insert(assignment.genome.taxonomy).second) {
                        truth_out << assignment.genome.taxonomy << '\n';
                    }
                }
            }

            fs::path metafile_path = cli.output_dir / "protal.meta";
            write_protal_metafile(
                samples, cli.output_dir, reads_dir, metafile_path, truth_paths, cli.protal_metafile_output_dir);
            std::cout << "Wrote Protal metafile: " << metafile_path << '\n';
        }
    } catch (const std::exception& ex) {
        std::cerr << "Simulation failed: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
