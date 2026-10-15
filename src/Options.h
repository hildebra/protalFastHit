//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <cerrno>
#include <cstring>
#include <unistd.h>
#include <glob.h>
#include <filesystem>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <utility>
#include "LineSplitter.h"
#include <fstream>
#include <regex>
#include <cctype>
#include <algorithm>
#include <numeric>
#include "Utilities.h"
#include "Zstd.h"
#include "Database.h"
#include "SamHandler.h"
#include "GeneIncongruence.h"
#include "BuildSettings.h"
#include "SpeciesPriors.h"
#include "NcbiTaxids.h"
#include "SpeciesNeighbours.h"
#include "CongenerGapsTable.h"
#include "ForeignRatesTable.h"
#include "StrainAlleles.h"
#include "ColumnWeights.h"
#include "AncestrySites.h"
#include "UniqueCores.h"
#include "SamFile.h"
#include "ReadType.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/ReadTypeDetection.h"
#include "SequenceUtils/GeneConservation.h"
#include "SequenceUtils/GeneNeighbours.h"
#include "SequenceUtils/GeneTableFile.h"
#include "gzstream/gzstream.h"


namespace protal {
    // Reads longer than this in a sample of short single-end reads are taken for long reads.
    static const size_t MAX_SHORT_READ_LENGTH = 1000;

    static const size_t DEFAULT_THREADS = 1;
    // --build writes a zstd-compressed database unless --no_compress.
    static const int DEFAULT_COMPRESS_LEVEL = 19;
    static const int DEFAULT_COMPRESS_WINDOW_LOG = 27;
    static const int DEFAULT_COMPRESS_FRAME_MB = 64;  // independent frames: loading uses -t threads
    static const size_t DEFAULT_ALIGN_TOP = 3;
    // Anchors a divergent read may try beyond --align_top on congeners of its best alignment that its seeds could not
    // tell apart (SimpleAlignmentHandler::operator(), --adaptive_candidates).
    static const size_t DEFAULT_ADAPTIVE_CANDIDATES = 7;
    static const double DEFAULT_MAX_SCORE_ANI = 0.9;
    static const size_t DEFAULT_MSA_MIN_HCOV = 1000;
    // Reads a position needs to be written in a strain MSA (a mixture needs --snp_min_cov per allele).
    static const size_t DEFAULT_MSA_MIN_DEPTH = 1;
    static const size_t DEFAULT_MIN_SUCCESSFUL_LOOKUPS = 4;
    static const size_t DEFAULT_X_DROP = 1000;
    static const size_t DEFAULT_MAX_KEY_UBIQUITY = 256;
    static const size_t DEFAULT_MAX_SEED_SIZE = 128;
    static const size_t DEFAULT_MAX_OUT = 1;
    static const bool DEFAULT_NO_STRAIN = false;
    static const std::string PROTAL_DB_ENV_VARIABLE = "PROTAL_DB_PATH";

    static const size_t DEFAULT_MIN_SNP_COV = 2;
    static const size_t DEFAULT_MIN_SNP_PHRED_SUM = 90;
    // Below it, an allele is sequencing noise or a read of a relative rather than a strain: at 0,
    // clean 50x samples carried ~180 IUPAC codes each (strain audit, 2026-09-29).
    static const double DEFAULT_MIN_SNP_AF = 0.15;
    static const size_t DEFAULT_MIN_SNP_MEAN_QUAL = 15;
    static const bool   DEFAULT_SNP_REQUIRE_STRAND = true; // disabled via --snp_no_strand
    static const size_t DEFAULT_SNP_MAX_ALLELES = 3;
    // protal emits a RAW MSA (all observed genes, all detected samples). Gene/sample
    // coverage filtering, multi-allelicity outlier removal and site cleanup all live
    // in the qcmsa post-filter (default on), which owns their thresholds and can be
    // re-run on the raw MSA without re-running protal. Reach them via --qcmsa_args.

    // The options' definitions (CxxOptions), the parsing and the checks before a run: Options.cpp.

    struct OptionsData {
        // flags
        bool build = false;
        bool no_profile = false;
        bool preload_genomes = false;
        bool benchmark_alignment = false;
        bool show_help = false;
        bool show_help_dev = false;
        bool show_map_help = false;
        bool show_version = false;
        bool no_strains = false;
        bool no_mate_guidance = false;
        bool no_gene_neighbours = false;
        bool keep_foreign_genes = false;
        bool keep_suspect_copies = false;
        double suspect_copy_distance = gene_incongruence::kDefaultSuspectDistance;
        size_t strain_alleles = 16;
        double allele_genome_share = 1.0;
        size_t column_weights = 10;
        std::string column_weights_alignment = "auto";
        std::string grid_choice = "scored";
        bool no_ancestry_sites = false;
        bool no_unique_cores = false;
        bool no_site_weights = false;
        bool weighted_site_shift = false;
        double index_alleles = 1;
        bool no_phasing = false;
        std::string scratch_dir;   // --scratch_dir (or the older --strain_spill): where the run's scratch folder goes, or none
        bool sams_named = false;   // a map's SAM column named the SAM files (not protal)
        bool fastalign = false;
        bool profile_only = false;
        bool force = false;
        bool verbose = false;
        bool mapq_debug_out = false;
        bool whole_read_alignment = false;
        bool no_alignment_screen = false;
        bool no_indel_bound = false;
        bool no_single_seed_anchors = false;
        bool compact_keys = true;
        size_t long_read_budget = 0;
        bool sequential_load = false;
        SparseFirstGates sparse_first;                // --sparse_first per read type; nothing: auto (the database's)
        std::optional<size_t> sparse_first_score;     // --sparse_first_score; nothing: the database's, else 14
        std::optional<size_t> sparse_first_hits;      // --sparse_first_hits; nothing: the database's, else 1
        SparseFirstPicks sparse_first_pick;           // --sparse_first_pick per read type; nothing: auto (the database's)
        bool stage1_probe = false;
        bool sparse_seeding = false;                  // --sparse_seeding
        std::string sparse_index = "auto";            // --sparse_index: auto, on or off
        bool taxon_statistics = false;
        bool full_sam_header = false;
        bool write_unmapped_reads = false;
        bool serial_index_passes = false;
        bool profile_ahead = false;
        bool no_unknown_share = false;
        bool ncbi_taxids = false;
        size_t index_batch_kb = 1024;
        std::array<bool, kReadTypeCount> flex_indels{};  // --flex_indels per read type

        // build
        std::vector<uint8_t> build_gene_mask;
        bool compress = true;
        bool bundle = true;
        bool compress_db = false;
        bool decompress_db = false;
        bool unpack_db = false;
        int compress_level = DEFAULT_COMPRESS_LEVEL;
        int compress_window_log = DEFAULT_COMPRESS_WINDOW_LOG;
        int compress_frame_mb = DEFAULT_COMPRESS_FRAME_MB;
        std::string unpack_dir;

        // the command line, for messages that show how to rerun protal
        std::vector<std::string> command_line;

        // paths
        std::string sequence_file;
        std::string full_sequence_file;
        std::string database_path;
        std::string output_dir;
        std::string strain_output_dir;
        std::string misc_output_dir;
        std::string map_file;
        std::string benchmark_alignment_output;

        // sample lists
        std::vector<std::string> first_list;
        std::vector<std::string> second_list;
        std::vector<std::string> prefix_list;
        std::vector<std::string> samplename_list;
        std::vector<std::string> sam_list;
        std::vector<std::string> profile_list;
        std::vector<std::string> profile_truth_list;
        std::vector<std::string> read_type_list;  // per sample: a ReadType token, or empty: pe or se by the second file
        std::vector<std::string> unmapped_reads_list;  // per sample, a map's UNMAPPED_READS: write, count or '-'
        std::vector<std::string> sparse_first_list;    // per sample, a map's SPARSE_FIRST: C or '-'
        std::vector<std::string> sparse_first_pick_list;  // per sample, a map's SPARSE_FIRST_PICK: hash, stride or '-'
        std::vector<std::string> flex_indels_list;     // per sample, a map's FLEX_INDELS: 1, 0 or '-'
        std::vector<size_t> range;
        // What reading the arguments found (patterns without a SAM, how samples were named), reported with the
        // other checks' errors and notes (PrepareAndCheckValidity).
        std::vector<std::string> input_errors;
        std::vector<std::string> input_notes;

        // profiling
        std::string profile_truth;
        std::string model;
        std::string model_se;
        std::string model_pb;
        std::string model_ont;
        std::string read_type;  // --read_type as given (empty if not): the samples' (read_type_list) and --add_model's
        std::string add_model;
        std::string model_gate;  // --model_gate as given (with --add_model)
        std::string model_kind = "auto";  // --model_kind as given
        std::string add_tables;  // --add_tables
        std::string write_species_neighbours;  // --write_species_neighbours
        double knob = 0.5;
        bool knob_given = false;  // --knob on the command line: the models' depth knobs are not used
        std::optional<double> fdr;  // --fdr; none: a model's calibrated calls are not used
        size_t singleton_congener = 0;  // --singleton_congener; 0: no rule
        double depth_identity_margin = 0.08;
        double msa_identity_margin = 0.04;
        std::string gene_conservation = "none";  // --gene_conservation: none, db (the database's), or a file

        // alignment
        size_t threads = DEFAULT_THREADS;
        size_t align_top = DEFAULT_ALIGN_TOP;
        size_t adaptive_candidates = DEFAULT_ADAPTIVE_CANDIDATES;
        bool no_allele_scores = false;
        double max_score_ani = DEFAULT_MAX_SCORE_ANI;
        bool max_score_ani_given = false;  // -a given: it applies to all read types
        size_t x_drop = DEFAULT_X_DROP;
        std::optional<WfaTraceback> wfa_traceback;  // --wfa_traceback; nothing: auto
        std::optional<FragmentCount> fragment_count;  // --fragment_count; nothing: auto
        size_t max_key_ubiquity = DEFAULT_MAX_KEY_UBIQUITY;
        size_t max_seed_size = DEFAULT_MAX_SEED_SIZE;
        size_t min_successful_lookups = DEFAULT_MIN_SUCCESSFUL_LOOKUPS;
        size_t max_out = DEFAULT_MAX_OUT;

        // strains / MSA
        size_t msa_min_hcov = DEFAULT_MSA_MIN_HCOV;
        size_t msa_min_depth = DEFAULT_MSA_MIN_DEPTH;
        std::optional<double> msa_knob;  // --msa_knob; none: --knob
        std::vector<std::string> msa_species;
        size_t snp_min_cov = DEFAULT_MIN_SNP_COV;
        size_t snp_min_phred_sum = DEFAULT_MIN_SNP_PHRED_SUM;
        double snp_min_af = DEFAULT_MIN_SNP_AF;
        bool snp_min_af_given = false;  // --snp_min_af given: it applies to all read types
        size_t snp_min_mean_qual = DEFAULT_MIN_SNP_MEAN_QUAL;
        bool   snp_require_strand = DEFAULT_SNP_REQUIRE_STRAND;
        size_t snp_max_alleles = DEFAULT_SNP_MAX_ALLELES;
        bool   run_qcmsa = true;
        std::string qcmsa_script;
        std::string qcmsa_args;
    };

    class Options {
    private:
        bool m_build = false;
        bool m_no_profile = false;
        bool m_preload_genomes = false;
        bool m_benchmark_alignment = false;

        bool m_show_help = false;
        bool m_show_help_dev = false;
        bool m_show_map_help = false;
        bool m_show_version = false;

        bool m_no_strains = false;
        bool m_no_mate_guidance = false;
        bool m_no_gene_neighbours = false;
        bool m_keep_foreign_genes = false;
        bool m_keep_suspect_copies = false;  // --keep_suspect_copies
        double m_suspect_copy_distance = gene_incongruence::kDefaultSuspectDistance;  // --suspect_copy_distance (--build)
        size_t m_strain_alleles = 16;  // --strain_alleles (--build)
        double m_allele_genome_share = 1.0;  // --allele_genome_share (--build)
        double m_index_alleles = 1;  // --index_alleles (--build): the least divergence of an allele whose k-mers are indexed
        size_t m_column_weights = 10;  // --column_weights (--build): genera per family that vote; 0: no table
        std::string m_column_weights_alignment = "auto";  // --column_weights_alignment (--build): auto, msa, aa or nt
        std::string m_grid_choice = "scored";  // --grid_choice (--build): scored, even or phase (GridLevels.h Choice)
        bool m_no_ancestry_sites = false;  // --no_ancestry_sites (--build): no ancestry_sites.bin
        bool m_no_unique_cores = false;  // --no_unique_cores: no unique_cores.bin (--build), or the table left unread (a run)
        bool m_no_site_weights = false;  // --no_site_weights: the table left unread
        bool m_weighted_site_shift = false;  // --weighted_site_shift (off by default since 2026-10-12): the shift discounts
                                             // differences at variable columns
        bool m_no_phasing = false;
        std::string m_scratch_dir;     // --scratch_dir
        std::string m_scratch_folder;  // the run's folder in it (UseScratchFolder), once made
        bool m_sams_named = false;     // a map's SAM column named the SAM files
        std::vector<bool> m_scratch_sams;  // by sample: its SAM is in the run's scratch folder (MoveSamsToScratch)
        bool m_fastalign = false;
        bool m_profile_only = false;
        bool m_force = false;
        bool m_verbose = false;

        bool m_mapq_debug_out = false;
        bool m_whole_read_alignment = false;
        bool m_no_alignment_screen = false;
        bool m_no_indel_bound = false;
        bool m_no_single_seed_anchors = false;
        bool m_compact_keys = true;  // the key map of the used cores only (CompactKeys.h); --direct_keys: false
        size_t m_long_read_budget = 0;
        bool m_sequential_load = false;
        bool m_stage1_probe = false;  // --stage1_probe
        bool m_sparse_seeding = false;  // --sparse_seeding
        std::string m_sparse_index = "auto";  // --sparse_index
        // --sparse_first as given (per read type; nothing: auto), and as the run gates each read type (C, 0 off): given,
        // else the database's (SetSparseFirst; --build: as given, 0 where not). Its score and hits likewise.
        SparseFirstGates m_sparse_first_given;
        std::array<size_t, kReadTypeCount> m_sparse_first{};
        std::optional<size_t> m_sparse_first_score_given;
        std::optional<size_t> m_sparse_first_hits_given;
        size_t m_sparse_first_score = 14;
        size_t m_sparse_first_hits = 1;
        // --sparse_first_pick as given (per read type; nothing: auto), and as the run picks for each read type: given,
        // else the database's (SetSparseFirstPick; --build: as given, stride where not).
        SparseFirstPicks m_sparse_first_pick_given;
        std::array<SparseFirstPick, kReadTypeCount> m_sparse_first_pick{};
        // The gates each read type has models of, and the database's own gates (SetModelGates); with a kind of models asked
        // for a read type (--model_kind), the database's member of that kind at each of its gates.
        std::array<std::vector<size_t>, kReadTypeCount> m_model_gates;
        std::array<size_t, kReadTypeCount> m_database_gates{};
        std::array<std::map<size_t, std::string>, kReadTypeCount> m_kind_models;
        // --model_kind of a run: each read type's kind of models, empty for auto (the database's default files).
        std::array<std::string, kReadTypeCount> m_model_kinds;
        bool m_taxon_statistics = false;
        bool m_full_sam_header = false;
        bool m_write_unmapped_reads = false;
        bool m_serial_index_passes = false;
        bool m_profile_ahead = false;
        bool m_no_unknown_share = false;
        bool m_ncbi_taxids = false;
        size_t m_index_batch_kb = 1024;
        std::array<bool, kReadTypeCount> m_flex_indels{};  // --flex_indels per read type

        size_t m_current_index = 0;

        std::vector<uint8_t> m_build_gene_mask;
        bool m_compress = true;
        bool m_bundle_db = true;
        bool m_compress_db = false;
        bool m_decompress_db = false;
        bool m_unpack_db = false;
        int m_compress_level = DEFAULT_COMPRESS_LEVEL;
        int m_compress_window_log = DEFAULT_COMPRESS_WINDOW_LOG;
        int m_compress_frame_mb = DEFAULT_COMPRESS_FRAME_MB;
        std::string m_unpack_dir;
        std::vector<std::string> m_command_line;

        std::string m_sequence_file;
        std::string m_full_sequence_file;
        std::string m_database_path;
        // Where m_database_path points (Database.h); set by PrepareAndCheckValidity. m_bundle is the
        // opened single-file database, if it is one.
        db::Location m_location;
        std::shared_ptr<db::Bundle const> m_bundle;
        std::string m_output_dir;
        std::string m_strain_output_dir;
        std::string m_misc_output_dir;

        std::vector<std::string> m_first_list;
        std::vector<std::string> m_second_list;
        std::vector<std::string> m_prefix_list;
        std::vector<std::string> m_sam_list;
        std::vector<std::string> m_profile_list;
        std::vector<std::string> m_sampleid_list;
        std::vector<std::string> m_profile_truth_list;
        // Per sample: the read type given (a token, or empty: by the second file), and its kind (see
        // ResolveReadTypes).
        std::vector<std::string> m_read_type_list;
        std::vector<ReadType> m_read_types;
        // Per sample, from a map's UNMAPPED_READS column: whether its unaligned reads get an unmapped record each
        // (empty without the column: --write_unmapped_reads for every sample).
        std::vector<char> m_unmapped_reads_list;
        // Per sample, from a map's SPARSE_FIRST column: its gate (C, 0 off), or nothing: as its read type's (empty
        // without the column).
        std::vector<std::optional<size_t>> m_sparse_first_list;
        // Per sample, from a map's SPARSE_FIRST_PICK column: its gate's pick, or nothing: as its read type's (empty
        // without the column).
        std::vector<std::optional<SparseFirstPick>> m_sparse_first_pick_list;
        // Per sample, from a map's FLEX_INDELS column: 1 or 0, or nothing: as its read type's (empty without the column).
        std::vector<std::optional<bool>> m_flex_indels_list;
        std::vector<std::string> m_input_errors;  // OptionsData::input_errors
        std::vector<std::string> m_input_notes;

        std::vector<size_t> m_range;

        std::string m_profile_truth;
        std::string m_model;
        std::string m_model_se;
        std::string m_model_pb;
        std::string m_model_ont;
        std::string m_read_type;  // --read_type as given; empty if not
        std::string m_add_model;  // --add_model
        std::string m_model_gate;  // --model_gate
        std::string m_model_kind;  // --model_kind as given (ModelKindsOf parses it)
        std::string m_add_tables;  // --add_tables
        std::string m_write_species_neighbours;  // --write_species_neighbours
        double m_knob = 0.5;
        bool m_knob_given = false;  // --knob
        std::optional<double> m_fdr;  // --fdr
        size_t m_singleton_congener = 0;  // --singleton_congener; 0: no rule
        double m_depth_identity_margin = 0.08;
        double m_msa_identity_margin = 0.04;
        std::string m_gene_conservation = "none";  // --gene_conservation

        size_t m_threads = DEFAULT_THREADS;

        std::string m_benchmark_alignment_output;

        size_t m_align_top = DEFAULT_ALIGN_TOP;
        size_t m_adaptive_candidates = DEFAULT_ADAPTIVE_CANDIDATES;
        bool m_no_allele_scores = false;  // --no_allele_scores
        double m_max_score_ani = DEFAULT_MAX_SCORE_ANI;
        bool m_max_score_ani_given = false;
        size_t m_x_drop = DEFAULT_X_DROP;
        std::optional<WfaTraceback> m_wfa_traceback_given;            // --wfa_traceback; nothing: auto
        WfaTraceback m_wfa_traceback = WfaTraceback::Singletrack;      // what aligns: given, else --build's singletrack, or the database's (SetWfaTraceback)
        std::optional<FragmentCount> m_fragment_count_given;          // --fragment_count; nothing: auto
        FragmentCount m_fragment_count = FragmentCount::Mapq;          // given, else --build's mapq, or the database's (SetFragmentCount)
        size_t m_max_key_ubiquity = DEFAULT_MAX_KEY_UBIQUITY;
        size_t m_max_seed_size = DEFAULT_MAX_SEED_SIZE;
        size_t m_min_successful_lookups = DEFAULT_MIN_SUCCESSFUL_LOOKUPS;
        size_t m_max_out = DEFAULT_MAX_OUT;
        size_t m_msa_min_hcov = DEFAULT_MSA_MIN_HCOV;
        size_t m_msa_min_depth = DEFAULT_MSA_MIN_DEPTH;
        std::optional<double> m_msa_knob;
        std::vector<std::string> m_msa_species;

        size_t m_snp_min_cov = DEFAULT_MIN_SNP_COV;
        size_t m_snp_min_phred_sum = DEFAULT_MIN_SNP_PHRED_SUM;
        double m_snp_min_af = DEFAULT_MIN_SNP_AF;
        bool m_snp_min_af_given = false;
        size_t m_snp_min_mean_qual = DEFAULT_MIN_SNP_MEAN_QUAL;
        bool   m_snp_require_strand = DEFAULT_SNP_REQUIRE_STRAND;
        size_t m_snp_max_alleles = DEFAULT_SNP_MAX_ALLELES;
        bool   m_run_qcmsa = true;
        std::string m_qcmsa_script;
        std::string m_qcmsa_args;

    public:
        static inline const std::string PROTAL_INDEX_FILE = "index.prx";
        static inline const std::string PROTAL_SEQUENCE_FILE = "reference.fna";
        static inline const std::string PROTAL_SEQUENCE_MAP_FILE = "reference.map";
        static inline const std::string PROTAL_HITTABLE_GENES_FILE = "species_gene_mask.tsv";
        static inline const std::string PROTAL_UNIQUE_KMER_FILE = "unique_kmers.tsv";
        static inline const std::string PROTAL_GENE_CONSERVATION_FILE = gene_conservation::kFileName;
        static inline const std::string PROTAL_SUSPECT_COPIES_FILE = gene_incongruence::kFileName;
        static inline const std::string PROTAL_SPECIES_PRIORS_FILE = species_priors::kFileName;
        static inline const std::string PROTAL_NCBI_TAXIDS_FILE = ncbi_taxids::kFileName;
        static inline const std::string PROTAL_SPECIES_NEIGHBOURS_FILE = species_neighbours::kFileName;
        static inline const std::string PROTAL_CONGENER_GAPS_FILE = congener_gaps::kFileName;
        static inline const std::string PROTAL_STRAIN_ALLELES_FILE = strain_alleles::kFileName;
        static inline const std::string PROTAL_STRAIN_ALLELES_BIN_FILE = strain_alleles::kBinaryFileName;
        static inline const std::string PROTAL_COLUMN_WEIGHTS_FILE = column_weights::kFileName;
        static inline const std::string PROTAL_COLUMN_WEIGHTS_BIN_FILE = column_weights::kBinaryFileName;
        static inline const std::string PROTAL_ANCESTRY_SITES_FILE = ancestry::kFileName;
        static inline const std::string PROTAL_UNIQUE_CORES_FILE = unique_cores::kFileName;
        static inline const std::string PROTAL_FOREIGN_RATES_FILE = foreign_rates::kFileName;
        static inline const std::string PROTAL_GENE_NEIGHBOURS_FILE = gene_neighbours::kFileName;
        static inline const std::string PROTAL_GENE_POSITIONS_FILE = gene_neighbours::kPositionsFileName;
        static inline const std::string PROTAL_TAXONOMY_FILE = "internal_taxonomy.dmp";
        static inline const std::string PROTAL_BUILD_SETTINGS_FILE = build_settings::kFileName;
        // The notes a database carries for its users, which no run reads (NoteFiles): build_gtdb_database.py's account of
        // the build, and the converter's names of the genes and genomes (gtdb_to_protal_db.py).
        static inline const std::string PROTAL_BUILD_METADATA_FILE = "build_metadata.tsv";
        static inline const std::string PROTAL_GENE_IDS_FILE = "gene2geneid.tsv";
        static inline const std::string PROTAL_GENOME_IDS_FILE = "genome2tiid.tsv";
        // The presence models of paired-end and of single-end reads.

        static inline const std::string MAP_SAMPLEID = "#SAMPLEID";
        static inline const std::string MAP_FIRST_READ = "FIRST";
        static inline const std::string MAP_SECOND_READ = "SECOND";
        static inline const std::string MAP_NO_SECOND_READ = "-";  // SECOND of a single-end sample
        // A sample's read files, read one after the other as one sample (its lanes or runs): in a FIRST or SECOND
        // cell, and with --one_sample in -1/-2, separated by commas; so kept in the lists of first and second files.
        static constexpr char kReadFileSeparator = ',';
        static inline const std::string MAP_READ_TYPE = "READ_TYPE";
        // A sample's reads that seeded on taxa but aligned nowhere: an unmapped record each (write) or counted per taxon
        // in the SAM header (count); '-': as --write_unmapped_reads says for every sample.
        static inline const std::string MAP_UNMAPPED_READS = "UNMAPPED_READS";
        static inline const std::string MAP_UNMAPPED_WRITE = "write";
        static inline const std::string MAP_UNMAPPED_COUNT = "count";
        static inline const std::string MAP_SPARSE_FIRST = "SPARSE_FIRST";
        static inline const std::string MAP_SPARSE_FIRST_PICK = "SPARSE_FIRST_PICK";
        static inline const std::string MAP_FLEX_INDELS = "FLEX_INDELS";
        static inline const std::string MAP_SAM = "SAM";
        static inline const std::string MAP_PROFILE = "PROFILE";
        static inline const std::string MAP_PROFILE_TRUTH = "PROFILE_TRUTH";
        static inline const std::string MAP_HEADER_TO_SPECIES = "HEADER_TO_SPECIES";
        static inline const std::string MAP_PREFIX = "PREFIX";

        static inline const std::string MAP_VAR_INPUT_DIR = "#INPUT_DIR";
        static inline const std::string MAP_VAR_OUTPUT_DIR = "#OUTPUT_DIR";
        static inline const std::string MAP_VAR_SAM_OUTPUT_DIR = "#SAM_OUTPUT_DIR";
        static inline const std::string MAP_VAR_PROFILE_OUTPUT_DIR = "#PROFILE_OUTPUT_DIR";
        static inline const std::string MAP_VAR_STRAIN_OUTPUT_DIR = "#STRAIN_OUTPUT_DIR";
        static inline const std::string MAP_VAR_MISC_OUTPUT_DIR = "#MISC_OUTPUT_DIR";

        static inline const std::string MAP_VAR_DEFAULT_SAM_OUTPUT_DIR = "alignments";
        static inline const std::string MAP_VAR_DEFAULT_PROFILE_OUTPUT_DIR = "profiles";
        static inline const std::string MAP_VAR_DEFAULT_STRAIN_OUTPUT_DIR = "strains";
        static inline const std::string MAP_VAR_DEFAULT_MISC_OUTPUT_DIR = "misc";


        const size_t MAP_SAMPLE_ID_COL = 0;

        Options(bool show_help, bool show_map_help, bool show_version, bool show_help_dev) :
                m_show_help(show_help), m_show_map_help(show_map_help),
                m_show_version(show_version), m_show_help_dev(show_help_dev) {}

        explicit Options(OptionsData d) :
                m_build(d.build),
                m_no_profile(d.no_profile),
                m_profile_only(d.profile_only),
                m_no_strains(d.no_strains),
                m_no_mate_guidance(d.no_mate_guidance),
                m_no_gene_neighbours(d.no_gene_neighbours),
                m_keep_foreign_genes(d.keep_foreign_genes),
                m_keep_suspect_copies(d.keep_suspect_copies),
                m_suspect_copy_distance(d.suspect_copy_distance),
                m_strain_alleles(d.strain_alleles),
                m_allele_genome_share(d.allele_genome_share),
                m_index_alleles(d.index_alleles),
                m_column_weights(d.column_weights),
                m_column_weights_alignment(std::move(d.column_weights_alignment)),
                m_grid_choice(std::move(d.grid_choice)),
                m_no_ancestry_sites(d.no_ancestry_sites),
                m_no_unique_cores(d.no_unique_cores),
                m_no_site_weights(d.no_site_weights),
                m_weighted_site_shift(d.weighted_site_shift),
                m_no_phasing(d.no_phasing),
                m_preload_genomes(d.preload_genomes),
                m_show_help(d.show_help),
                m_show_help_dev(d.show_help_dev),
                m_show_map_help(d.show_map_help),
                m_show_version(d.show_version),
                m_benchmark_alignment(d.benchmark_alignment),
                m_benchmark_alignment_output(std::move(d.benchmark_alignment_output)),
                m_mapq_debug_out(d.mapq_debug_out),
                m_whole_read_alignment(d.whole_read_alignment),
                m_no_alignment_screen(d.no_alignment_screen),
                m_no_indel_bound(d.no_indel_bound),
                m_no_single_seed_anchors(d.no_single_seed_anchors),
                m_compact_keys(d.compact_keys),
                m_long_read_budget(d.long_read_budget),
                m_sequential_load(d.sequential_load),
                m_stage1_probe(d.stage1_probe),
                m_sparse_seeding(d.sparse_seeding),
                m_sparse_index(d.sparse_index),
                m_sparse_first_given(d.sparse_first),
                m_sparse_first(GatesOrOff(d.sparse_first)),
                m_sparse_first_score_given(d.sparse_first_score),
                m_sparse_first_hits_given(d.sparse_first_hits),
                m_sparse_first_score(d.sparse_first_score.value_or(14)),
                m_sparse_first_hits(d.sparse_first_hits.value_or(1)),
                m_sparse_first_pick_given(d.sparse_first_pick),
                m_sparse_first_pick(PicksOrStride(d.sparse_first_pick)),
                m_taxon_statistics(d.taxon_statistics),
                m_full_sam_header(d.full_sam_header),
                m_write_unmapped_reads(d.write_unmapped_reads),
                m_serial_index_passes(d.serial_index_passes),
                m_profile_ahead(d.profile_ahead),
                m_no_unknown_share(d.no_unknown_share),
                m_ncbi_taxids(d.ncbi_taxids),
                m_index_batch_kb(d.index_batch_kb),
                m_flex_indels(d.flex_indels),
                m_fastalign(d.fastalign),
                m_force(d.force),
                m_verbose(d.verbose),
                m_build_gene_mask(std::move(d.build_gene_mask)),
                m_compress(d.compress),
                m_bundle_db(d.bundle),
                m_compress_db(d.compress_db),
                m_decompress_db(d.decompress_db),
                m_unpack_db(d.unpack_db),
                m_compress_level(d.compress_level),
                m_compress_window_log(d.compress_window_log),
                m_compress_frame_mb(d.compress_frame_mb),
                m_unpack_dir(std::move(d.unpack_dir)),
                m_command_line(std::move(d.command_line)),
                m_sequence_file(std::move(d.sequence_file)),
                m_full_sequence_file(std::move(d.full_sequence_file)),
                m_database_path(std::move(d.database_path)),
                m_output_dir(std::move(d.output_dir)),
                m_strain_output_dir(std::move(d.strain_output_dir)),
                m_misc_output_dir(std::move(d.misc_output_dir)),
                m_first_list(std::move(d.first_list)),
                m_second_list(std::move(d.second_list)),
                m_prefix_list(std::move(d.prefix_list)),
                m_sam_list(std::move(d.sam_list)),
                m_profile_list(std::move(d.profile_list)),
                m_profile_truth_list(std::move(d.profile_truth_list)),
                m_range(std::move(d.range)),
                m_profile_truth(std::move(d.profile_truth)),
                m_read_type_list(std::move(d.read_type_list)),
                m_unmapped_reads_list(UnmappedReadsFlags(d.unmapped_reads_list, d.write_unmapped_reads)),
                m_sparse_first_list(SparseFirstList(d.sparse_first_list)),
                m_sparse_first_pick_list(SparseFirstPickList(d.sparse_first_pick_list)),
                m_flex_indels_list(FlexIndelsList(d.flex_indels_list)),
                m_model(std::move(d.model)),
                m_model_se(std::move(d.model_se)),
                m_model_pb(std::move(d.model_pb)),
                m_model_ont(std::move(d.model_ont)),
                m_read_type(std::move(d.read_type)),
                m_add_model(std::move(d.add_model)),
                m_model_gate(std::move(d.model_gate)),
                m_model_kind(std::move(d.model_kind)),
                m_add_tables(std::move(d.add_tables)),
                m_write_species_neighbours(std::move(d.write_species_neighbours)),
                m_knob(d.knob),
                m_knob_given(d.knob_given),
                m_fdr(d.fdr),
                m_singleton_congener(d.singleton_congener),
                m_depth_identity_margin(d.depth_identity_margin),
                m_msa_identity_margin(d.msa_identity_margin),
                m_gene_conservation(std::move(d.gene_conservation)),
                m_threads(d.threads),
                m_align_top(d.align_top),
                m_adaptive_candidates(d.adaptive_candidates),
                m_no_allele_scores(d.no_allele_scores),
                m_max_score_ani(d.max_score_ani),
                m_max_score_ani_given(d.max_score_ani_given),
                m_x_drop(d.x_drop),
                m_wfa_traceback_given(d.wfa_traceback),
                m_wfa_traceback(d.wfa_traceback.value_or(WfaTraceback::Singletrack)),
                m_fragment_count_given(d.fragment_count),
                m_fragment_count(d.fragment_count.value_or(FragmentCount::Mapq)),
                m_max_key_ubiquity(d.max_key_ubiquity),
                m_max_seed_size(d.max_seed_size),
                m_min_successful_lookups(d.min_successful_lookups),
                m_max_out(d.max_out),
                m_msa_min_hcov(d.msa_min_hcov),
                m_msa_min_depth(d.msa_min_depth),
                m_msa_knob(d.msa_knob),
                m_msa_species(std::move(d.msa_species)),
                m_snp_min_cov(d.snp_min_cov),
                m_snp_min_phred_sum(d.snp_min_phred_sum),
                m_snp_min_af(d.snp_min_af),
                m_snp_min_af_given(d.snp_min_af_given),
                m_snp_min_mean_qual(d.snp_min_mean_qual),
                m_snp_require_strand(d.snp_require_strand),
                m_snp_max_alleles(d.snp_max_alleles),
                m_run_qcmsa(d.run_qcmsa),
                m_qcmsa_script(std::move(d.qcmsa_script)),
                m_qcmsa_args(std::move(d.qcmsa_args)) {
            m_input_errors = std::move(d.input_errors);
            m_input_notes = std::move(d.input_notes);
            m_scratch_dir = std::move(d.scratch_dir);
            m_sams_named = d.sams_named;
            if (d.samplename_list.empty()) {
                // A sample is named after its prefix's file name, not its path.
                for (auto const& prefix : m_prefix_list) {
                    m_sampleid_list.emplace_back(std::filesystem::path(prefix).filename().string());
                }
            } else {
                m_sampleid_list = std::move(d.samplename_list);
            }
            if (m_profile_list.empty()) {
                for (auto& prefix : m_prefix_list) {
                    m_profile_list.emplace_back(prefix + ".profile");
                    std::cout << m_profile_list.size() << " " << prefix + ".profile" << std::endl;
                }
            }
        };

        std::string ToString() const;

        bool BuildMode() const {
            return m_build;
        }

        bool NoProfile() const {
            return m_no_profile;
        }

        bool ProfileOnly() const {
            return m_profile_only;
        }

        std::string& ProfileTruthFile() {
            return m_profile_truth;
        }

        double GetKnob() const {
            return m_knob;
        }

        // Whether --knob was given: then it applies to every sample; otherwise a model's knob for the sample's depth
        // (its depth knobs, machine_learning_cmdline.py --depth-knobs) where it has one, else --knob's default.
        bool KnobGiven() const {
            return m_knob_given;
        }

        // Whether --fdr was given, and its value: the expected share of false calls for a model with calibrated calls
        // (ProfileWrapper, context::FalseCallKnob), which are used only when it is given and above 0.
        bool FdrGiven() const {
            return m_fdr.has_value();
        }

        double GetFdr() const {
            return m_fdr.value_or(0);
        }

        // --singleton_congener: a taxon of one fragment beside a congener of at least this many is not reported; 0: no rule.
        size_t GetSingletonCongener() const {
            return m_singleton_congener;
        }

        // Whether --msa_knob was given: otherwise each sample's MSA knob is its profile's knob.
        bool MSAKnobGiven() const {
            return m_msa_knob.has_value();
        }

        // The model probability a sample's taxon needs to enter its strain MSA: --msa_knob, else --knob (ProfileWrapper:
        // without --msa_knob, each sample's own knob, which may be its model's for its depth).
        double GetMSAKnob() const {
            return m_msa_knob.value_or(m_knob);
        }

        double GetDepthIdentityMargin() const {
            return m_depth_identity_margin;
        }

        double GetMSAIdentityMargin() const {
            return m_msa_identity_margin;
        }

        // The table files --add_tables stores (comma-separated as given); empty without it.
        std::vector<std::string> AddTables() const {
            std::vector<std::string> files;
            std::string item;
            std::istringstream is(m_add_tables);
            while (std::getline(is, item, ',')) {
                if (!item.empty()) files.push_back(item);
            }
            return files;
        }

        // The file --write_species_neighbours writes the species' nearest congeners to (Build.h WriteSpeciesNeighboursOnly);
        // empty without it.
        std::string const& WriteSpeciesNeighboursFile() const {
            return m_write_species_neighbours;
        }

        // The tables --add_tables may store, by file name: those runs read, and the notes (NoteFiles).
        static std::vector<std::string> AddableTables() {
            std::vector<std::string> tables = { PROTAL_CONGENER_GAPS_FILE, PROTAL_FOREIGN_RATES_FILE, PROTAL_SPECIES_PRIORS_FILE,
                                                PROTAL_NCBI_TAXIDS_FILE, PROTAL_STRAIN_ALLELES_FILE, PROTAL_COLUMN_WEIGHTS_FILE };
            for (auto const& note : NoteFiles()) tables.push_back(note);
            return tables;
        }

        // The files a database carries that no run reads (since 2026-10-13): what it was built from and with, and the names
        // of its genes and genomes, so that database.protal holds everything of a database, not only what runs read.
        // --build and --compress_db pack them where the folder has them and leave the files there (build inputs, read by
        // scripts); --add_tables stores them; --unpack_db writes them back. build_metadata.tsv goes last, after the models
        // (Build.h ModelsLast), so that a build's last step replaces it in place with them.
        static std::vector<std::string> NoteFiles() {
            return { PROTAL_GENE_IDS_FILE, PROTAL_GENOME_IDS_FILE, PROTAL_BUILD_METADATA_FILE };
        }

        static bool IsNoteFile(std::string const& name) {
            auto const notes = NoteFiles();
            return std::find(notes.begin(), notes.end(), name) != notes.end();
        }

        // --no_allele_scores: candidates scored against the reference alone, not their species' strain alleles.
        bool NoAlleleScores() const {
            return m_no_allele_scores;
        }

        // --adaptive_candidates: anchors a divergent read may try beyond --align_top (0: none).
        size_t GetAdaptiveCandidates() const {
            return m_adaptive_candidates;
        }

        // --add_model as given: one PMML file, or several, comma-separated (AddModels).
        std::string const& GetAddModel() const {
            return m_add_model;
        }

        // A model --add_model stores: its file, the reads it is for, the gate its training runs were gated at
        // (--model_gate; nothing: the database's own), and its kind (--model_kind; nothing: the read type's default file).
        struct AddedModel {
            std::string file;
            ReadType type;
            std::optional<size_t> gate;
            std::optional<std::string> kind = std::nullopt;
        };

        // The models --add_model stores: its files, comma-separated, with --read_type's tokens, --model_gate's gates and
        // --model_kind's kinds at the same places (a single model: pe if --read_type is not given; the database's gate where
        // --model_gate gives none or '-'; the default file where --model_kind gives none, '-' or auto). Checked in
        // PrepareAndCheckValidity; a read type that is no token is paired-end here.
        std::vector<AddedModel> AddModels() const {
            auto const files = SplitCommas(m_add_model);
            auto const types = SplitCommas(m_read_type);
            auto const gates = SplitCommas(m_model_gate);
            auto const kinds = m_model_kind == "auto" ? std::vector<std::string>() : SplitCommas(m_model_kind);
            std::vector<AddedModel> models;
            for (size_t i = 0; i < files.size(); i++) {
                models.push_back({files[i], i < types.size() ? ReadTypeFromToken(types[i]).value_or(ReadType::Paired)
                                                             : ReadType::Paired,
                                  i < gates.size() && gates[i] != "-" ? ParseSparseFirstLevel(gates[i]) : std::nullopt,
                                  i < kinds.size() && kinds[i] != "-" ? std::optional(kinds[i]) : std::nullopt});
            }
            return models;
        }

        // Comma-separated items, empty ones kept ("a,,b": a, "", b).
        static std::vector<std::string> SplitCommas(std::string const& text) {
            std::vector<std::string> items;
            std::string item;
            std::istringstream is(text);
            while (std::getline(is, item, ',')) items.push_back(item);
            return items;
        }

        // --model_kind of a run: each read type's kind of models (ReadType.h kModelKinds), empty for auto; nullopt with a
        // message in error if `text` is none of auto, a kind, or read types and theirs (pe:forest,pb:gbm).
        static std::optional<std::array<std::string, kReadTypeCount>> ModelKindsOf(std::string const& text, std::string& error) {
            std::array<std::string, kReadTypeCount> kinds{};
            if (text == "auto") return kinds;
            if (IsModelKind(text)) {
                kinds.fill(text);
                return kinds;
            }
            auto const items = SplitCommas(text);
            for (auto const& item : items.empty() ? std::vector<std::string>{ "" } : items) {
                auto const colon = item.find(':');
                auto const type = colon == std::string::npos ? std::nullopt : ReadTypeFromToken(item.substr(0, colon));
                std::string const kind = colon == std::string::npos ? "" : item.substr(colon + 1);
                if (!type || !(IsModelKind(kind) || kind == "auto")) {
                    error = "--model_kind is auto, " + ModelKindNames() + ", or read types and their kinds (pe:forest,pb:gbm; "
                            "read types " + ReadTypeTokens() + "), not '" + text + "'";
                    return std::nullopt;
                }
                kinds[static_cast<size_t>(*type)] = kind == "auto" ? "" : kind;
            }
            return kinds;
        }

        // The kind of models a run profiles reads of `type` with (--model_kind), empty for the database's default files.
        std::string const& ModelKind(ReadType type) const {
            return m_model_kinds[static_cast<size_t>(type)];
        }

        bool PreloadGenomes() const {
            return m_preload_genomes;
        }

        bool Help() const {
            return m_show_help;
        }

        bool HelpDev() const {
            return m_show_help_dev;
        }

        bool ShowMapHelp() const {
            return m_show_map_help;
        }

        bool ShowVersion() const {
            return m_show_version;
        }

        bool NoStrains() const {
            return m_no_strains;
        }

        // Paired-end reads: whether a mate sure of its alignment guides the other (MateGuidance.h).
        bool MateGuidance() const {
            return !m_no_mate_guidance;
        }

        // Whether the database's gene neighbours are used (GeneNeighbours.h), if it has any.
        bool GeneNeighbours() const {
            return !m_no_gene_neighbours;
        }

        // Whether foreign genes (profiler::Taxon::ForeignGene) are left out of their taxon's depth and strain MSAs.
        bool DropForeignGenes() const {
            return !m_keep_foreign_genes;
        }

        // Whether records on the database's suspect gene copies (GeneIncongruence.h) are left out of the evidence.
        bool DropSuspectCopies() const {
            return !m_keep_suspect_copies;
        }

        // --build: the k-mer distance within which another genus's copy makes a gene copy suspect; 0: no scan.
        double GetSuspectCopyDistance() const {
            return m_suspect_copy_distance;
        }

        // --strain_alleles (--build): alleles kept per species and gene; 0: none.
        size_t GetStrainAlleles() const {
            return m_strain_alleles;
        }

        // --allele_genome_share (--build): the share of the genomes that may give alleles (strain_alleles::AlleleGenome).
        double GetAlleleGenomeShare() const {
            return m_allele_genome_share;
        }

        // --index_alleles (--build): the k-mers of the strain alleles at least this far from the representative are
        // indexed (StrainAllelesBuild.h AlleleSeeds); 1 or more (the default): none.
        double IndexAlleles() const {
            return m_index_alleles;
        }

        // --grid_choice (--build, a developer option): how the grid levels choose each window's core: scored, even or phase
        // (grid_levels::Choice).
        std::string const& GridChoice() const {
            return m_grid_choice;
        }

        // --column_weights (--build): genera per family that vote on the column weights; 0: no table.
        size_t GetColumnWeights() const {
            return m_column_weights;
        }

        // --column_weights_alignment (--build): auto (msa where the database folder has GTDB's alignments, else aa), msa, aa
        // or nt.
        std::string const& ColumnWeightsAlignment() const {
            return m_column_weights_alignment;
        }

        // --no_ancestry_sites (--build): no table of every copy's ancestry sites; a run computes them as it needs them.
        bool NoAncestrySites() const {
            return m_no_ancestry_sites;
        }

        // --no_unique_cores: --build writes no unique_cores.bin (UniqueCores.h); a run counts ZS, ZU and ZT from its anchors'
        // seeds, as before 2026-10-14, though its database has the table.
        bool NoUniqueCores() const {
            return m_no_unique_cores;
        }

        // --no_site_weights: the database's column weights left unread.
        bool NoSiteWeights() const {
            return m_no_site_weights;
        }

        // --weighted_site_shift (off by default; the default 2026-10-11 to 2026-10-12): the site shift discounts a
        // candidate's differences at variable columns (StrainAlleles.h).
        bool WeightedSiteShift() const {
            return m_weighted_site_shift;
        }

        // Whether a long-read sample's strain MSA row is split into its strains' (Haplotypes.h).
        bool Phasing() const {
            return !m_no_phasing;
        }

        // --scratch_dir: the folder protal makes the run's scratch folder in (ScratchSpace.h), or empty (none).
        std::string const& ScratchDir() const {
            return m_scratch_dir;
        }

        // The run's scratch folder (UseScratchFolder), or empty without --scratch_dir: the strain stage's spill files,
        // the SAMs of a run without strain MSAs (IsScratchSam), the records of a SAM before its header, qcMSA's files.
        std::string const& ScratchFolder() const {
            return m_scratch_folder;
        }

        // Whether sample `index`'s SAM is in the run's scratch folder (UseScratchFolder): removed once it is profiled.
        bool IsScratchSam(size_t index) const {
            return index < m_scratch_sams.size() && m_scratch_sams[index];
        }

        // The run's scratch folder, made in --scratch_dir (RunProtal: scratch::Space). With --no_strains a run that
        // profiles reads its SAMs only to profile them: the SAM files that protal names (not a map's SAM column) of the
        // samples it aligns (in the range, without an earlier run's SAM in their place, or with --force) are written
        // into the folder, under their file names (with the sample's number in front where two share one). Returns
        // how many.
        size_t UseScratchFolder(std::string const& folder) {
            m_scratch_folder = folder;
            m_scratch_sams.assign(m_sam_list.size(), false);
            if (folder.empty() || !m_no_strains || m_no_profile || m_profile_only || m_build || m_sams_named) return 0;
            std::set<std::string> names;
            size_t moved = 0;
            for (auto const index : m_range) {
                if (index >= m_sam_list.size()) continue;
                auto const& sam = m_sam_list[index];
                if (!m_force && (std::filesystem::exists(sam) || std::filesystem::exists(UncompressedSamName(sam)))) continue;
                std::string name = std::filesystem::path(sam).filename().string();
                while (!names.insert(name).second) name = std::to_string(index + 1) + "_" + name;
                m_sam_list[index] = (std::filesystem::path(folder) / name).string();
                m_scratch_sams[index] = true;
                moved++;
            }
            return moved;
        }

        bool BenchmarkAlignment() const {
            return m_benchmark_alignment;
        }

        bool FastAlign() const {
            return m_fastalign;
        }

        std::string GetSequenceFilePath() const {
            return m_sequence_file;
        }

        std::string GetFullSequenceFilePath() const {
            return m_full_sequence_file;
        }

        // The raw file names; the build writes the index as GetIndexFile() + ".zst" unless
        // --no_compress. Readers use the Resolved* variants, which pick the file that exists.
        std::string GetIndexFile() const {
            return m_database_path + "/" + PROTAL_INDEX_FILE;
        }

        std::string ResolvedIndexFile() const {
            return zstd::Resolve(GetIndexFile());
        }

        std::string ResolvedSequenceFile() const {
            return zstd::Resolve(GetSequenceFile());
        }

        bool Compress() const {
            return m_compress;
        }

        bool CompressDbMode() const {
            return m_compress_db;
        }

        bool DecompressDbMode() const {
            return m_decompress_db;
        }

        zstd::Params CompressionParams() const {
            return { m_compress_level, m_compress_window_log, static_cast<int>(std::max<size_t>(m_threads, 1)),
                     static_cast<uint64_t>(std::max(m_compress_frame_mb, 0)) << 20 };
        }

        std::string GetInternalTaxonomyFile() const {
            return m_database_path + "/" + PROTAL_TAXONOMY_FILE;
        }

        std::string GetSequenceFile() const {
            return m_database_path + "/" + PROTAL_SEQUENCE_FILE;
        }

        std::string GetSequenceMapFile() const {
            return m_database_path + "/" + PROTAL_SEQUENCE_MAP_FILE;
        }

        std::string GetHittableGenesMap() const {
            return m_database_path + "/" + PROTAL_HITTABLE_GENES_FILE;
        }

        std::string GetUniqueKmersFile() const {
            return m_database_path + "/" + PROTAL_UNIQUE_KMER_FILE;
        }

        std::string GetGeneConservationFile() const {
            return m_database_path + "/" + PROTAL_GENE_CONSERVATION_FILE;
        }

        std::string GetGeneNeighboursFile() const {
            return m_database_path + "/" + PROTAL_GENE_NEIGHBOURS_FILE;
        }

        // The database's suspect gene copies (--build, GeneIncongruence.h), packed by --build.
        std::string GetSuspectCopiesFile() const {
            return m_database_path + "/" + PROTAL_SUSPECT_COPIES_FILE;
        }

        // The species' priors (gtdb_to_protal_db.py, SpeciesPriors.h), packed by --build.
        std::string GetSpeciesPriorsFile() const {
            return m_database_path + "/" + PROTAL_SPECIES_PRIORS_FILE;
        }

        // Each species' NCBI species taxid (gtdb_to_protal_db.py, NcbiTaxids.h), packed by --build.
        std::string GetNcbiTaxidsFile() const {
            return m_database_path + "/" + PROTAL_NCBI_TAXIDS_FILE;
        }

        // The settings --build used that runs should match (BuildSettings.h: the WFA2 traceback), packed by --build.
        std::string GetBuildSettingsFile() const {
            return m_database_path + "/" + PROTAL_BUILD_SETTINGS_FILE;
        }

        // Each species' nearest congeners (SpeciesNeighbours.h), written and packed by --build.
        std::string GetSpeciesNeighboursFile() const {
            return m_database_path + "/" + PROTAL_SPECIES_NEIGHBOURS_FILE;
        }

        // Each gene copy's gap to its congeners' copies (CongenerGaps.h), written and packed by --build.
        std::string GetCongenerGapsFile() const {
            return m_database_path + "/" + PROTAL_CONGENER_GAPS_FILE;
        }

        // Each gene copy's strain alleles (StrainAlleles.h), written and packed by --build (--strain_alleles).
        std::string GetStrainAllelesFile() const {
            return m_database_path + "/" + PROTAL_STRAIN_ALLELES_FILE;
        }

        // Each family's column weights (ColumnWeights.h), written and packed by --build (--column_weights).
        std::string GetColumnWeightsFile() const {
            return m_database_path + "/" + PROTAL_COLUMN_WEIGHTS_FILE;
        }

        // Every gene copy's ancestry sites (AncestrySites.h), written and packed by --build (unless --no_ancestry_sites).
        std::string GetAncestrySitesFile() const {
            return m_database_path + "/" + PROTAL_ANCESTRY_SITES_FILE;
        }

        // Which of every gene copy's index cores are unique (UniqueCores.h), written and packed by --build (unless --no_unique_cores).
        std::string GetUniqueCoresFile() const {
            return m_database_path + "/" + PROTAL_UNIQUE_CORES_FILE;
        }

        // Each gene copy's reads of a tiled scan and their foreign share (ForeignRatesTable.h): packed by --build if the
        // folder has it, or by --add_tables.
        std::string GetForeignRatesFile() const {
            return m_database_path + "/" + PROTAL_FOREIGN_RATES_FILE;
        }

        // Where gene_neighbours.py placed each gene in each genome: packed by --build, never read by a run.
        std::string GetGenePositionsFile() const {
            return m_database_path + "/" + PROTAL_GENE_POSITIONS_FILE;
        }

        bool HittableGenesMapExists() const {
            return Utils::exists(GetHittableGenesMap());
        }

        bool UniqueKmersFileExists() const {
            return UniqueKmersDbFile().Exists();
        }

        // ---- The database as protal reads it (Database.h): a member of the single-file database, or
        // the file in the database directory (the Get*File names above, .zst siblings resolved).

        bool IsBundle() const {
            return m_bundle != nullptr;
        }

        std::shared_ptr<db::Bundle const> const& GetBundle() const {
            return m_bundle;
        }

        db::Location const& GetLocation() const {
            return m_location;
        }

        db::DbFile DbFileNamed(std::string const& name, std::string const& path_in_dir) const {
            return m_bundle ? db::DbFile::InBundle(*m_bundle, name) : db::DbFile::OnDisk(path_in_dir);
        }

        db::DbFile IndexDbFile() const {
            return DbFileNamed(PROTAL_INDEX_FILE, ResolvedIndexFile());
        }

        db::DbFile SequenceDbFile() const {
            return DbFileNamed(PROTAL_SEQUENCE_FILE, ResolvedSequenceFile());
        }

        db::DbFile SequenceMapDbFile() const {
            return DbFileNamed(PROTAL_SEQUENCE_MAP_FILE, GetSequenceMapFile());
        }

        db::DbFile TaxonomyDbFile() const {
            return DbFileNamed(PROTAL_TAXONOMY_FILE, GetInternalTaxonomyFile());
        }

        db::DbFile UniqueKmersDbFile() const {
            return DbFileNamed(PROTAL_UNIQUE_KMER_FILE, GetUniqueKmersFile());
        }

        // The single-file database's binary gene table (GeneTableFile.h), if it has one; a folder's text tables are read as they are.
        std::optional<db::DbFile> GeneTableDbFile() const {
            if (!m_bundle) return std::nullopt;
            auto file = db::DbFile::InBundle(*m_bundle, gene_table_file::kFileName);
            return file.Exists() ? std::optional<db::DbFile>(std::move(file)) : std::nullopt;
        }

        // The database's gene_conservation.tsv (--build); Exists() is false if it has none.
        db::DbFile DatabaseGeneConservationDbFile() const {
            return DbFileNamed(PROTAL_GENE_CONSERVATION_FILE, GetGeneConservationFile());
        }

        // The database's suspect_copies.tsv (--build); Exists() is false if it has none.
        db::DbFile SuspectCopiesDbFile() const {
            return DbFileNamed(PROTAL_SUSPECT_COPIES_FILE, GetSuspectCopiesFile());
        }

        // The database's species_priors.tsv (the converter's); Exists() is false if it has none.
        db::DbFile SpeciesPriorsDbFile() const {
            return DbFileNamed(PROTAL_SPECIES_PRIORS_FILE, GetSpeciesPriorsFile());
        }

        // The database's ncbi_taxids.tsv (the converter's since 2026-10-11); Exists() is false if it has none.
        db::DbFile NcbiTaxidsDbFile() const {
            return DbFileNamed(PROTAL_NCBI_TAXIDS_FILE, GetNcbiTaxidsFile());
        }

        // The database's build_settings.tsv (--build's since 2026-10-11); Exists() is false if it has none.
        db::DbFile BuildSettingsDbFile() const {
            return DbFileNamed(PROTAL_BUILD_SETTINGS_FILE, GetBuildSettingsFile());
        }

        // The database's species_neighbours.tsv (--build's since 2026-10-06); Exists() is false if it has none.
        db::DbFile SpeciesNeighboursDbFile() const {
            return DbFileNamed(PROTAL_SPECIES_NEIGHBOURS_FILE, GetSpeciesNeighboursFile());
        }

        // The database's congener_gaps.tsv (--build's since 2026-10-07); Exists() is false if it has none.
        db::DbFile CongenerGapsDbFile() const {
            return DbFileNamed(PROTAL_CONGENER_GAPS_FILE, GetCongenerGapsFile());
        }

        // The database's strain_alleles.tsv (--build's since 2026-10-08); Exists() is false if it has none.
        db::DbFile StrainAllelesDbFile() const {
            return DbFileNamed(PROTAL_STRAIN_ALLELES_FILE, GetStrainAllelesFile());
        }

        // The database's column_weights.tsv (--build's since 2026-10-10); Exists() is false if it has none.
        db::DbFile ColumnWeightsDbFile() const {
            return DbFileNamed(PROTAL_COLUMN_WEIGHTS_FILE, GetColumnWeightsFile());
        }

        // The single-file database's binary tables (BinaryTable.h: column_weights.bin, strain_alleles.bin), if it has them; a
        // folder's text tables are read as they are.
        std::optional<db::DbFile> ColumnWeightsBinDbFile() const {
            return BundleMember(PROTAL_COLUMN_WEIGHTS_BIN_FILE);
        }
        std::optional<db::DbFile> StrainAllelesBinDbFile() const {
            return BundleMember(PROTAL_STRAIN_ALLELES_BIN_FILE);
        }
        std::optional<db::DbFile> BundleMember(std::string const& name) const {
            if (!m_bundle) return std::nullopt;
            auto file = db::DbFile::InBundle(*m_bundle, name);
            return file.Exists() ? std::optional<db::DbFile>(std::move(file)) : std::nullopt;
        }

        // The database's ancestry_sites.bin (--build's since 2026-10-11); Exists() is false if it has none.
        db::DbFile AncestrySitesDbFile() const {
            return DbFileNamed(PROTAL_ANCESTRY_SITES_FILE, GetAncestrySitesFile());
        }

        // The database's unique_cores.bin (--build's since 2026-10-14); Exists() is false if it has none.
        db::DbFile UniqueCoresDbFile() const {
            return DbFileNamed(PROTAL_UNIQUE_CORES_FILE, GetUniqueCoresFile());
        }

        // The database's foreign_rates.tsv (scripts/foreign_rates.py, --add_tables); Exists() is false if it has none.
        db::DbFile ForeignRatesDbFile() const {
            return DbFileNamed(PROTAL_FOREIGN_RATES_FILE, GetForeignRatesFile());
        }

        // The database's gene_neighbours.tsv (scripts/mini_db/gene_neighbours.py, packed by --build); Exists() is
        // false if it has none.
        db::DbFile GeneNeighboursDbFile() const {
            return DbFileNamed(PROTAL_GENE_NEIGHBOURS_FILE, GetGeneNeighboursFile());
        }

        // The genes' conservation factors: --gene_conservation's file, else the database's (for none, the default,
        // and db). They give the model's conservation features (conserved_fast_depth_ratio, conserved_hit_share)
        // either way, and scale the depth identity margin with db or a file (ScaleDepthMarginByConservation).
        db::DbFile GeneConservationDbFile() const {
            if (m_gene_conservation.empty() || m_gene_conservation == "none" || m_gene_conservation == "db") {
                return DatabaseGeneConservationDbFile();
            }
            return db::DbFile::OnDisk(m_gene_conservation);
        }

        // Whether the depth identity margin is scaled per gene (--gene_conservation db or a file).
        bool ScaleDepthMarginByConservation() const {
            return !m_gene_conservation.empty() && m_gene_conservation != "none";
        }

        // The model given for reads of `type`: by the type's option (--model_se, --model_pb,
        // --model_ont), else by --model; empty if neither is given.
        std::string const& ModelNamed(ReadType type) const {
            std::string const& own = type == ReadType::Single ? m_model_se : type == ReadType::PacBio ? m_model_pb :
                                     type == ReadType::ONT ? m_model_ont : m_model;
            return own.empty() ? m_model : own;
        }

        // The PMML model of samples with reads of `type`: the model named by the type's option
        // (--model_se, --model_pb, --model_ont) or else by --model, as an existing file or else as <name> in the
        // database (<name>.xml without an extension); by default the first ModelCandidates file
        // (ReadType.h) the database has. Exists() is false if it has none.
        db::DbFile ModelDbFile(ReadType type = ReadType::Paired) const {
            std::string const& name = ModelNamed(type);
            if (!name.empty()) {
                if (std::filesystem::exists(name)) return db::DbFile::OnDisk(name);
                std::string const file = std::filesystem::path(name).extension().empty() ? name + ".xml" : name;
                return DbFileNamed(file, (std::filesystem::path(m_database_path) / file).string());
            }
            auto const candidates = ModelCandidates(type);
            for (auto const& file : candidates) {
                auto model = DbFileNamed(file, (std::filesystem::path(m_database_path) / file).string());
                if (model.Exists()) return model;
            }
            return DbFileNamed(candidates.front(), (std::filesystem::path(m_database_path) / candidates.front()).string());
        }

        // The database's model of reads of `type` trained at another gate than its own (GateModelFile); Exists() is false
        // if it has none.
        db::DbFile GateModelDbFile(ReadType type, size_t gate) const {
            std::string const file = GateModelFile(type, gate);
            return DbFileNamed(file, (std::filesystem::path(m_database_path) / file).string());
        }

        // A file of the database by its name: a member of database.protal, or the file in the database's folder.
        db::DbFile DatabaseFileNamed(std::string const& name) const {
            return DbFileNamed(name, (std::filesystem::path(m_database_path) / name).string());
        }

        // The PMML model sample `index` is profiled with: the one given for its reads (ModelNamed), else the database's
        // model of its read type trained at the gate nearest its own (ModelGate): with a kind asked for (--model_kind), the
        // database's model of that kind there (SetModelGates); else the type's model at the database's gate (ModelDbFile),
        // its model at another gate (GateModelDbFile) there.
        db::DbFile SampleModelDbFile(size_t index) const {
            auto const type = GetReadType(index);
            size_t const gate = ModelGate(index);
            if (ModelNamed(type).empty()) {
                auto const& kinds = m_kind_models[static_cast<size_t>(type)];
                if (auto const it = kinds.find(gate); it != kinds.end()) return DatabaseFileNamed(it->second);
            }
            if (!ModelNamed(type).empty() || gate == m_database_gates[static_cast<size_t>(type)]) return ModelDbFile(type);
            return GateModelDbFile(type, gate);
        }

        // The kind of reads of sample `index` (see ResolveReadTypes).
        ReadType GetReadType(size_t index) const {
            return index < m_read_types.size() ? m_read_types[index] : ReadType::Paired;
        }

        // Whether any sample in the range has reads of `type`.
        bool AnySample(ReadType type) const {
            return std::any_of(m_range.begin(), m_range.end(), [&](size_t i) { return GetReadType(i) == type; });
        }

        bool WriteBundle() const {
            return m_bundle_db && m_compress;
        }

        bool UnpackDbMode() const {
            return m_unpack_db;
        }

        // Where --unpack_db writes the files: --unpack_dir, or the folder database.protal is in.
        std::string UnpackDir() const {
            return m_unpack_dir.empty() ? m_location.dir : m_unpack_dir;
        }

        // The name to show for protal in commands: as it was called.
        std::string ProgramName() const {
            if (m_command_line.empty()) return "protal";
            return m_command_line.front();
        }

        // A word as a POSIX shell reads it back.
        static std::string ShellWord(std::string const& word) {
            bool const plain = !word.empty() && std::all_of(word.begin(), word.end(), [](unsigned char c) {
                return std::isalnum(c) || std::strchr("_-./:=,+@%", c) != nullptr;
            });
            if (plain) return word;
            std::string quoted = "'";
            for (char c : word) quoted += c == '\'' ? std::string("'\\''") : std::string(1, c);
            return quoted + "'";
        }

        // This run's command line with --db set to db (added if the database came from $PROTAL_DB_PATH).
        std::string CommandWithDb(std::string const& db) const {
            std::string command = ShellWord(ProgramName());
            bool replaced = false;
            for (size_t i = 1; i < m_command_line.size(); i++) {
                std::string const& arg = m_command_line[i];
                if (arg == "--db" && i + 1 < m_command_line.size()) {
                    command += " --db " + ShellWord(db);
                    replaced = true;
                    i++;
                } else if (arg.rfind("--db=", 0) == 0) {
                    command += " --db=" + ShellWord(db);
                    replaced = true;
                } else {
                    command += " " + ShellWord(arg);
                }
            }
            if (!replaced) command += " --db " + ShellWord(db);
            return command;
        }

        void SetCurrentIndex(size_t i) {
            if (i >= m_prefix_list.size()) {
                std::cerr << "SetCurrentIndex to " << i << " not possible." << std::endl;
                std::cerr << "Max is " << m_prefix_list.size() << std::endl;
                exit(9);
            }
            m_current_index = i;
        }

        size_t GetCurrentIndex() const {
            return m_current_index;
        }

        std::vector<size_t> GetRange() const {
            return m_range;
        }

        bool HasBuildGeneSubset() const {
            return !m_build_gene_mask.empty();
        }

        bool BuildGeneAllowed(size_t gene_id) const {
            if (m_build_gene_mask.empty()) return true;
            if (gene_id >= m_build_gene_mask.size()) return false;
            return m_build_gene_mask[gene_id] != 0;
        }

//        std::string GetOutputPrefix() const {
//            return m_output_prefix;
//        }
//
//        std::string GetFirstFile() const {
//            return m_first;
//        }
//
//        std::string GetSecondFile() const {
//            return m_second;
//        }

        size_t GetFileCount() const {
            return m_prefix_list.size();
        }

        // A sample's read files as kept (kReadFileSeparator between them), one by one; none for an empty entry (the
        // second files of single-end reads).
        static std::vector<std::string> SplitReadFiles(std::string const& files) {
            std::vector<std::string> list;
            if (files.empty()) return list;
            size_t begin = 0;
            while (true) {
                size_t const end = files.find(kReadFileSeparator, begin);
                list.push_back(files.substr(begin, end == std::string::npos ? std::string::npos : end - begin));
                if (end == std::string::npos) break;
                begin = end + 1;
            }
            return list;
        }

        static std::string JoinReadFiles(std::vector<std::string> const& files) {
            std::string joined;
            for (auto const& file : files) joined += (joined.empty() ? "" : std::string(1, kReadFileSeparator)) + file;
            return joined;
        }

        // The name a sample of several read files (its lanes or runs) gets without --prefix: the longest beginning
        // their names share, without a lane number cut short (_L00 of _L001 and _L002) or separators at its end;
        // empty if they share none.
        static std::string LaneSamplePrefix(std::vector<std::string> const& names) {
            if (names.empty()) return "";
            std::string prefix = names.front();
            for (auto const& name : names) {
                size_t i = 0;
                while (i < prefix.size() && i < name.size() && prefix[i] == name[i]) i++;
                prefix.resize(i);
            }
            auto trim_separators = [&prefix] {
                while (!prefix.empty() && (prefix.back() == '_' || prefix.back() == '.' || prefix.back() == '-')) prefix.pop_back();
            };
            trim_separators();
            size_t digits = prefix.size();
            while (digits > 0 && std::isdigit(static_cast<unsigned char>(prefix[digits - 1]))) digits--;
            if (digits >= 2 && prefix[digits - 1] == 'L' && prefix[digits - 2] == '_') {
                prefix.resize(digits - 2);
                trim_separators();
            }
            return prefix;
        }

        // The first-in-pair read files of sample `index` (its single-end or long reads), one or several (its lanes),
        // read one after the other.
        std::vector<std::string> GetFirstFiles(int index) const {
            if (index >= m_first_list.size()) {
                std::cerr << "Cannot access index " << index << " of first files (Length: " << m_first_list.size() << ")" << std::endl;
                exit(33);
            }
            return SplitReadFiles(m_first_list[index]);
        }

        // The second-in-pair read files of sample `index`, paired with its first by position; none for single-end
        // and long reads.
        std::vector<std::string> GetSecondFiles(int index) const {
            if (index >= m_second_list.size()) {
                std::cerr << "Cannot access index " << index << " of second files (Length: " << m_second_list.size() << ")" << std::endl;
                exit(33);
            }
            return SplitReadFiles(m_second_list[index]);
        }

        // The read files of sample `index` for messages: its first files, then its second ones.
        std::string ReadFilesText(int index) const {
            std::string text;
            for (auto const& files : { GetFirstFiles(index), GetSecondFiles(index) }) {
                for (auto const& file : files) text += (text.empty() ? "" : ", ") + file;
            }
            return text;
        }

        // The SAM file of sample `index`, and whether its name asks for compression (".gz" or ".zst",
        // see SamFile.h); with strip_compression the name without that extension.
        std::pair<std::string, bool> SamFile(int index, bool strip_compression=false) const {
            if (index >= m_sam_list.size()) {
                std::cerr << "Cannot access index " << index << " of sam files (Length: " << m_sam_list.size() << ")" << std::endl;
                exit(33);
            }
            auto sam = m_sam_list[index];
            bool compressed = SamCompressionOf(sam) != SamCompression::None;
            if (strip_compression && compressed) {
                sam = UncompressedSamName(sam);
                compressed = false;
            }
            return {sam, compressed};
        }

        // Points sample `index` at the uncompressed name of its SAM file (without ".gz" or ".zst").
        void UseUncompressedSamFile(int index) {
            m_sam_list[index] = UncompressedSamName(m_sam_list[index]);
        }

        std::vector<std::string> SamFiles() {
            return m_sam_list;
        }

        std::string ProfileFile(int index, bool gzipped = false) const {
            if (index >= m_profile_list.size()) {
                std::cerr << "Cannot access index " << index << " of profile files (Length: " << m_profile_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_profile_list[index] + (gzipped ? ".gz" : "");
        }

        std::string ProfileTruthFile(int index) const {
            if (index >= m_profile_truth_list.size()) {
                std::cerr << "Cannot access index " << index << " of profile files (Length: " << m_profile_truth_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_profile_truth_list[index];
        }

        std::string GetPrefix(int index) const {
            if (index > m_prefix_list.size()) {
                std::cerr << "Cannot access index " << index << " of prefix files (Length: " << m_prefix_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_prefix_list[index];
        }

        std::string GetSampleId(int index) const {
            if (index > m_sampleid_list.size()) {
                std::cerr << "Cannot access index " << index << " of sample names (Length: " << m_sampleid_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_sampleid_list[index];
        }

        std::string GetOutputDir() const {
            return m_output_dir;
        }

        std::string GetStrainOutputDir() const {
            return m_strain_output_dir;
        }

        std::string GetMiscOutputDir() const {
            return m_misc_output_dir;
        }

        std::string GetSimilarityMatrixOutput(std::string species_name) const {
            return m_output_dir + '/' + species_name + ".tsv";
        }

        // protal's native (pre-qcmsa) MSA. The qcmsa post-filter writes the final
        // <species>.msa.fna from this; see docs/strains.md.
        std::string GetMSAOutput(std::string species_name) const {
            return m_strain_output_dir + '/' + species_name + ".raw.msa.fna";
        }

        std::string GetSpeciesMetaOutput(std::string species_name) const {
            return m_strain_output_dir + '/' + species_name + ".meta.tsv";
        }

        std::string GetMSAPartitionOutput(std::string species_name) const {
            return m_strain_output_dir + '/' + species_name + ".raw.partition.txt";
        }

        // How a species' long-read samples were phased into strain rows (Haplotypes.h).
        std::string GetHaplotypesOutput(std::string species_name) const {
            return m_strain_output_dir + '/' + species_name + ".haplotypes.tsv";
        }

        std::string GetMSAStatsOutput(std::string species_name) const {
            return m_strain_output_dir + '/' + species_name + ".snp_stats.tsv";
        }

        // The species this run wrote strain outputs for, with their MSA files.
        std::string GetStrainSpeciesListOutput() const {
            return m_strain_output_dir + "/species.tsv";
        }

        std::string GetBenchmarkAlignmentOutputFile() const {
            return m_benchmark_alignment_output;
        }

        const std::string& GetBenchmarkAlignmentOutputFile() {
            return m_benchmark_alignment_output;
        }

        bool Force() const {
            return m_force;
        }

        bool Verbose() const {
            return m_verbose;
        }

        size_t GetThreads() const {
            return m_threads;
        }

        bool GetMAPQDebugOut() const {
            return m_mapq_debug_out;
        }

        // --whole_read_alignment: short reads aligned as a whole into their window, not from their anchors.
        bool WholeReadAlignment() const {
            return m_whole_read_alignment;
        }

        // --no_alignment_screen: every candidate aligned with WFA2, without the k-mer screen (AlignmentScreen.h).
        bool NoAlignmentScreen() const {
            return m_no_alignment_screen;
        }

        // --no_indel_bound: ONT reads screened by their k-mers, not by the indel bound (AlignmentScreen.h IndelBound).
        bool NoIndelBound() const {
            return m_no_indel_bound;
        }

        // --no_single_seed_anchors: with the grid pick, short reads anchor on two seeds of a gene, as other runs, not on one
        // (ChainAnchorFinder::SetSingleSeedAnchors).
        bool NoSingleSeedAnchors() const {
            return m_no_single_seed_anchors;
        }

        // The loaded key map held as CompactKeys (the used cores' starts only; Seedmap::SetCompactLoad, CompactKeymap): the
        // default since 2026-10-14; --direct_keys keeps the direct one.
        bool CompactKeys() const {
            return m_compact_keys;
        }

        // --flex_indels: whether sample `index`'s lookups score their flex cells with an indel per flank allowed
        // (ChainAnchorFinder::SetFlexIndels): the map's FLEX_INDELS where it gives 1 or 0, else its read type's.
        bool FlexIndels(size_t index) const {
            if (index < m_flex_indels_list.size() && m_flex_indels_list[index]) return *m_flex_indels_list[index];
            return FlexIndels(GetReadType(index));
        }
        bool FlexIndels(ReadType type) const {
            return m_flex_indels[static_cast<size_t>(type)];
        }

        // --flex_indels' value: none, all, or read type tokens comma-separated (pe,se,pb,ont); nothing if it is none of them.
        static std::optional<std::array<bool, kReadTypeCount>> ParseFlexIndels(std::string const& text) {
            std::array<bool, kReadTypeCount> on{};
            if (text == "none") return on;
            if (text == "all") {
                on.fill(true);
                return on;
            }
            std::vector<std::string> tokens;
            LineSplitter::Split(text, ",", tokens);
            if (tokens.empty()) return std::nullopt;
            for (auto const& token : tokens) {
                auto const type = ReadTypeFromToken(token);
                if (!type) return std::nullopt;
                on[static_cast<size_t>(*type)] = true;
            }
            return on;
        }
        // "none", "all" or the read types on ("pb,ont"), for the options summary.
        static std::string FlexIndelsText(std::array<bool, kReadTypeCount> const& on) {
            if (std::none_of(on.begin(), on.end(), [](bool b) { return b; })) return "none";
            if (std::all_of(on.begin(), on.end(), [](bool b) { return b; })) return "all";
            std::string text;
            for (auto const& info : kReadTypes) {
                if (on[static_cast<size_t>(info.type)]) text += (text.empty() ? "" : ",") + info.token;
            }
            return text;
        }
        // A map's FLEX_INDELS values: 1 or 0, or nothing for '-' (LoadFromMap admits no other value).
        static std::vector<std::optional<bool>> FlexIndelsList(std::vector<std::string> const& values) {
            std::vector<std::optional<bool>> list;
            list.reserve(values.size());
            for (auto const& value : values) list.push_back(value == "-" ? std::nullopt : std::optional<bool>(value == "1"));
            return list;
        }

        // --long_read_budget: edits past a segment's best hit its other candidates may cost (0: the ANI floor's budget).
        size_t GetLongReadBudget() const {
            return m_long_read_budget;
        }

        // --sequential_load: the database's parts loaded one after another (RunProtal::Run).
        bool SequentialLoad() const {
            return m_sequential_load;
        }

        // --stage1_probe: what a sparse first stage would keep, tallied while aligning (Stage1Probe.h, RunProtal).
        bool Stage1Probe() const {
            return m_stage1_probe;
        }

        // --sparse_seeding: the reads a gate keeps are seeded from its cores alone (ChainAnchorFinder::SetSeedLevel).
        bool SparseSeeding() const {
            return m_sparse_seeding;
        }
        // --sparse_index: auto, on or off (RunProtal::SparseIndexLevel).
        std::string const& SparseIndex() const {
            return m_sparse_index;
        }

        // --sparse_first: the gate of sample `index` (C: 1/C of a read's cores looked up first, 0 off): the map's
        // SPARSE_FIRST where it gives one, else its read type's (GetSparseFirst(ReadType)); its score and hits
        // (ChainAnchorFinder::SetSparseFirst, Classify's short-read loops).
        size_t GetSparseFirst(size_t index) const {
            if (index < m_sparse_first_list.size() && m_sparse_first_list[index]) return *m_sparse_first_list[index];
            return GetSparseFirst(GetReadType(index));
        }
        // The gate of a read type: --sparse_first's, else the database's (SetSparseFirst; --build: 0 where not given).
        size_t GetSparseFirst(ReadType type) const {
            return m_sparse_first[static_cast<size_t>(type)];
        }
        size_t GetSparseFirstScore() const {
            return m_sparse_first_score;
        }
        size_t GetSparseFirstHits() const {
            return m_sparse_first_hits;
        }
        // --sparse_first, --sparse_first_score and --sparse_first_hits as given (nothing: auto); what --build records.
        SparseFirstGates const& SparseFirstGiven() const {
            return m_sparse_first_given;
        }
        std::optional<size_t> SparseFirstScoreGiven() const {
            return m_sparse_first_score_given;
        }
        std::optional<size_t> SparseFirstHitsGiven() const {
            return m_sparse_first_hits_given;
        }
        // The run's gates by read type, score and hits, once the database's are known (RunProtal: ApplyBuildSettings).
        void SetSparseFirst(std::array<size_t, kReadTypeCount> const& gates, size_t score, size_t hits) {
            m_sparse_first = gates;
            m_sparse_first_score = score;
            m_sparse_first_hits = hits;
        }
        // Whether a sample of the run is gated (map rows included).
        bool AnySparseFirst() const {
            return std::any_of(m_range.begin(), m_range.end(), [&](size_t i) { return GetSparseFirst(i) > 0; });
        }

        // --sparse_first per read type with 0 where not given (--build, and a run until the database's are known).
        static std::array<size_t, kReadTypeCount> GatesOrOff(SparseFirstGates const& gates) {
            std::array<size_t, kReadTypeCount> levels{};
            for (size_t t = 0; t < kReadTypeCount; t++) levels[t] = gates[t].value_or(0);
            return levels;
        }
        // A map's SPARSE_FIRST values as gates: C, or nothing for '-' (LoadFromMap admits no other value).
        static std::vector<std::optional<size_t>> SparseFirstList(std::vector<std::string> const& values) {
            std::vector<std::optional<size_t>> list;
            list.reserve(values.size());
            for (auto const& value : values) list.push_back(value == "-" ? std::nullopt : ParseSparseFirstLevel(value));
            return list;
        }

        // --sparse_first_pick: the pick of sample `index`'s gate (which 1/C of a read's cores it looks up first): the
        // map's SPARSE_FIRST_PICK where it gives one, else its read type's (ChainAnchorFinder::SetSparseFirst).
        SparseFirstPick GetSparseFirstPick(size_t index) const {
            if (index < m_sparse_first_pick_list.size() && m_sparse_first_pick_list[index]) return *m_sparse_first_pick_list[index];
            return GetSparseFirstPick(GetReadType(index));
        }
        // The pick of a read type: --sparse_first_pick's, else the database's (SetSparseFirstPick; --build: stride where
        // not given).
        SparseFirstPick GetSparseFirstPick(ReadType type) const {
            return m_sparse_first_pick[static_cast<size_t>(type)];
        }
        // --sparse_first_pick as given (nothing: auto); what --build records.
        SparseFirstPicks const& SparseFirstPickGiven() const {
            return m_sparse_first_pick_given;
        }
        // The run's picks by read type, once the database's are known (RunProtal: ApplySparseFirst).
        void SetSparseFirstPick(std::array<SparseFirstPick, kReadTypeCount> const& picks) {
            m_sparse_first_pick = picks;
        }
        // --sparse_first_pick per read type with stride where not given (--build, and a run until the database's are
        // known).
        static std::array<SparseFirstPick, kReadTypeCount> PicksOrStride(SparseFirstPicks const& picks) {
            std::array<SparseFirstPick, kReadTypeCount> out{};
            for (size_t t = 0; t < kReadTypeCount; t++) out[t] = picks[t].value_or(SparseFirstPick::Stride);
            return out;
        }

        // The gates each read type has models of: the database's own gate (its model_file) and those of its models at
        // other gates (GateModelFile), sorted; set once the database's gates are known (RunProtal: ApplySparseFirst). With a
        // kind asked for (--model_kind), the gates it has a model of that kind at, and `kinds` that model's member at each.
        void SetModelGates(std::array<std::vector<size_t>, kReadTypeCount> gates, std::array<size_t, kReadTypeCount> own,
                           std::array<std::map<size_t, std::string>, kReadTypeCount> kinds = {}) {
            m_model_gates = std::move(gates);
            m_database_gates = own;
            m_kind_models = std::move(kinds);
        }
        std::vector<size_t> const& ModelGates(ReadType type) const {
            return m_model_gates[static_cast<size_t>(type)];
        }
        // The gate whose model profiles sample `index`: of the gates its read type has models of, the one nearest its
        // gate (NearestGate); its read type's database gate before SetModelGates, or with a model given (--model, ...).
        size_t ModelGate(size_t index) const {
            auto const type = GetReadType(index);
            auto const& gates = ModelGates(type);
            if (gates.empty() || !ModelNamed(type).empty()) return m_database_gates[static_cast<size_t>(type)];
            return NearestGate(GetSparseFirst(index), gates);
        }
        // A map's SPARSE_FIRST_PICK values: hash or stride, or nothing for '-' (LoadFromMap admits no other value).
        static std::vector<std::optional<SparseFirstPick>> SparseFirstPickList(std::vector<std::string> const& values) {
            std::vector<std::optional<SparseFirstPick>> list;
            list.reserve(values.size());
            for (auto const& value : values) list.push_back(value == "-" ? std::nullopt : ParseSparseFirstPick(value));
            return list;
        }

        // --taxon_statistics: one misc/<taxon>.statistics.tsv per taxon with reads (off by default).
        bool TaxonStatistics() const {
            return m_taxon_statistics;
        }

        // --full_sam_header: every gene of the database in the SAM header, not only those aligned to.
        bool FullSamHeader() const {
            return m_full_sam_header;
        }

        // --write_unmapped_reads (or --full_sam_header, whose header is written first): an unmapped record for every read that
        // seeded on taxa but aligned nowhere, rather than their counts per taxon in the SAM header.
        bool WriteUnmappedReads() const {
            return m_write_unmapped_reads || m_full_sam_header;
        }

        // The same for sample `index`: a map's UNMAPPED_READS (write or count) where the map has the column, else as
        // above; --full_sam_header, whose header is written before the reads, writes them in any case.
        bool WriteUnmappedReads(size_t index) const {
            if (m_full_sam_header) return true;
            return index < m_unmapped_reads_list.size() ? m_unmapped_reads_list[index] != 0 : WriteUnmappedReads();
        }

        // A map's UNMAPPED_READS values as flags: write 1, count 0, '-' `all` (--write_unmapped_reads). LoadFromMap admits
        // no other value.
        static std::vector<char> UnmappedReadsFlags(std::vector<std::string> const& values, bool all) {
            std::vector<char> flags;
            flags.reserve(values.size());
            for (auto const& value : values) {
                flags.push_back(value == MAP_UNMAPPED_WRITE ? 1 : value == MAP_UNMAPPED_COUNT ? 0 : static_cast<char>(all));
            }
            return flags;
        }

        // --serial_index_passes: the index build's passes and value pointers on one thread (build::Run).
        bool SerialIndexPasses() const {
            return m_serial_index_passes;
        }

        // --profile_ahead: samples profiled while the next ones are aligned (ProfilingAhead in RunProtal.h); off by default.
        bool ProfileAhead() const {
            return m_profile_ahead;
        }

        // --no_unknown_share: the profile's abundances add up to 1 over the called species, without the "?" line
        // (Composition.h).
        bool NoUnknownShare() const {
            return m_no_unknown_share;
        }

        // --ncbi_taxids: the profile's lines end with the species' NCBI taxid and name (the database's ncbi_taxids.tsv).
        bool NcbiTaxids() const {
            return m_ncbi_taxids;
        }

        // --index_batch_kb: bytes of the reference per batch in the parallel index passes.
        size_t IndexBatchBytes() const {
            return m_index_batch_kb * 1024;
        }

        bool HasProfileTruths() const {
            return m_profile_truth_list.size() == m_profile_list.size();
        }

        auto GetMSAMinHCOV() {
            return m_msa_min_hcov;
        }
        auto GetMSAMinDepth() const { return m_msa_min_depth; }
        const std::vector<std::string>& GetMSASpecies() const {
            return m_msa_species;
        }

        auto GetSNPMinCov() const { return m_snp_min_cov; }
        auto GetSNPMinPhredSum() const { return m_snp_min_phred_sum; }
        auto GetSNPMinAF() const { return m_snp_min_af; }
        // --snp_min_af for reads of `type`: the read type's default unless it is given (ReadTypeInfo).
        double GetSNPMinAF(ReadType type) const {
            return !m_snp_min_af_given && Info(type).snp_min_af ? *Info(type).snp_min_af : m_snp_min_af;
        }
        auto GetSNPMinMeanQual() const { return m_snp_min_mean_qual; }
        auto GetSNPRequireStrand() const { return m_snp_require_strand; }
        auto GetSNPMaxAlleles() const { return m_snp_max_alleles; }
        bool GetRunQCMSA() const { return m_run_qcmsa; }
        const std::string& GetQCMSAScript() const { return m_qcmsa_script; }
        const std::string& GetQCMSAArgs() const { return m_qcmsa_args; }

        size_t GetAlignTop() const {
            return m_align_top;
        }

        size_t GetMaxOut() const {
            return m_max_out;
        }

        size_t GetMaxSeedSize() const {
            return m_max_seed_size;
        }

        size_t GetXDrop() const {
            return m_x_drop;
        }

        // --wfa_traceback as given; nothing for auto (the database's).
        std::optional<WfaTraceback> WfaTracebackGiven() const {
            return m_wfa_traceback_given;
        }

        // The traceback the aligners use: the one given, else (--build: singletrack, since 0.8.2) the database's once a run has set it.
        WfaTraceback GetWfaTraceback() const {
            return m_wfa_traceback;
        }

        void SetWfaTraceback(WfaTraceback traceback) {
            m_wfa_traceback = traceback;
        }

        // --fragment_count as given; nothing for auto (the database's).
        std::optional<FragmentCount> FragmentCountGiven() const {
            return m_fragment_count_given;
        }

        // How the model's `fragments` feature counts: given, else --build's mapq, or the database's (SetFragmentCount).
        FragmentCount GetFragmentCount() const {
            return m_fragment_count;
        }

        void SetFragmentCount(FragmentCount count) {
            m_fragment_count = count;
        }

        size_t GetMaxKeyUbiquity() const {
            return m_max_key_ubiquity;
        }

        size_t GetMinSuccessfulLookups() const {
            return m_min_successful_lookups;
        }

        double GetMaxScoreAni() const {
            return m_max_score_ani;
        }

        // -a for reads of `type`: the read type's default unless it is given (ReadTypeInfo).
        double GetMaxScoreAni(ReadType type) const {
            return !m_max_score_ani_given && Info(type).max_score_ani ? *Info(type).max_score_ani : m_max_score_ani;
        }

        void PrintHelp(bool show_dev=false);

        void PrintMapHelp(std::ostream& os = std::cout) {
            os << R"(
The map file helps you organise input files spanning different folders without having to use complicated wildcard terms.
Header lines are indicated with a # and you can define individual output-base folders for strain output, regular output, sam output,
and profile output. You could leave them empty and always specify the full path in the sample rows (not starting with a hashtag),
but this can get complicated quite quickly. A SAM name ending in .sam.zst (the default when the SAM column is left out) is
written zstd-compressed, .sam.gz gzip-compressed, any other name plain. A relative OUTPUT_DIR (or -o) and INPUT_DIR are
relative to the folder protal runs in; a relative SAM_OUTPUT_DIR, PROFILE_OUTPUT_DIR, STRAIN_OUTPUT_DIR or MISC_OUTPUT_DIR
lies in OUTPUT_DIR, and without them those folders are OUTPUT_DIR/alignments, profiles, strains and misc.
				
#OUTPUT_DIR	/path-to-your-results-dir/					
#SAM_OUTPUT_DIR	/path-to-your-results-dir/alignments					
#PROFILE_OUTPUT_DIR	/path-to-your-results-dir/profiles
#STRAIN_OUTPUT_DIR	/path-to-your-results-dir/strains	
#MISC_OUTPUT_DIR	/path-to-your-results-dir/misc
#INPUT_DIR	/path-to-input-dir/
#SAMPLEID	FIRST	SECOND	SAM	PREFIX	PROFILE
SAMPLE1	sample1/reads_1.fq	sample1/reads_2.fq	1.sam	AIR1	1.profile
SAMPLE2	sample2/reads_1.fq	sample2/reads_2.fq	2.sam	AIR2	2.profile
SAMPLE3	sample3/reads_1.fq	sample3/reads_2.fq	3.sam	AIR3	3.profile
SAMPLE4	sample4/reads.fq	-	4.sam	AIR4	4.profile

FIRST and PREFIX are mandatory. SECOND holds the second-in-pair files of paired-end reads; a
sample with single-end reads (in FIRST) has '-' there, as SAMPLE4 above. Without a SECOND
column, all samples are single-end. Single-end samples are profiled with the database's
single-end model (model_se.xml, or --model_se).
A sample sequenced in several lanes or runs lists its files in one cell, comma-separated, read
one after the other as one sample: FIRST s_L001_R1.fq.gz,s_L002_R1.fq.gz and SECOND
s_L001_R2.fq.gz,s_L002_R2.fq.gz (pairing up by position; a lane whose two files hold
different numbers of reads fails the sample). A file name cannot contain a comma.
An optional READ_TYPE column names each sample's reads: pe (paired-end), se (single-end), pb
(PacBio) or ont (Oxford Nanopore long reads; pb and ont with '-' as SECOND), each profiled with
its model (model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml). Without the column,
--read_type applies to all samples; without either, or with '-' as READ_TYPE, a sample is pe
with a SECOND file and se without.
The first column, #SAMPLEID, names the sample in the outputs (MSA rows, logs, statistics) and in
file names, so it must be safe as a file name on Linux and macOS: no '/' or ':', no control
characters, not starting with '-', not '.' or '..', UTF-8, at most 200 bytes; protal stops on
reading the map otherwise. With strain MSAs it holds no whitespace either.
An optional UNMAPPED_READS column says per sample what becomes of the reads that seeded on taxa
but aligned nowhere: write (an unmapped record each, with the taxa in its ZF tag, as
--write_unmapped_reads writes them for all samples) or count (counted per taxon in the SAM
header, the default); '-' leaves it to --write_unmapped_reads.
An optional SPARSE_FIRST column gives a sample its own sparse-first gate (0 off, 1, 2, 4, 8, 16 or
32; see --sparse_first); '-' leaves it to its read type's.
An optional SPARSE_FIRST_PICK column gives a sample its gate's pick (hash or stride; see
--sparse_first_pick); '-' leaves it to its read type's.
An optional FLEX_INDELS column turns --flex_indels on (1) or off (0) for a sample; '-' leaves it to
its read type's.
SAM and PROFILE are optional and default to <PREFIX>.sam.zst (--sam_format) and <PREFIX>.profile,
in OUTPUT_DIR. Every sample needs its own SAM and PROFILE file; protal stops if two samples share
one. A sample whose SAM exists is not aligned again (unless --force), and a SAM given as an
absolute path may lie anywhere, e.g. in another run's folder: a map of existing SAMs (protal_map_utils
merge writes one from the runs' maps) builds strain MSAs over the samples of several runs.)" << std::endl;
        }


        // The ending of the SAM files protal names, by --sam_format (zst, gz, sam); empty if unknown.
        static std::string SamEnding(std::string const& format) {
            if (format == "zst") return ".sam.zst";
            if (format == "gz") return ".sam.gz";
            if (format == "sam") return ".sam";
            return "";
        }
        static inline const std::string kDefaultSamEnding = ".sam.zst";

        // sam_ending: of the SAM names given to samples without a SAM column (<prefix><sam_ending>); sams_named (if
        // given): set to whether the map's SAM column named them.
        static bool LoadFromMap(std::string map_path, std::string& output_dir, std::string& strain_output_dir, std::string& misc_output_dir,
                                std::vector<std::string>& prefix_list,
                                std::vector<std::string>& first_list, std::vector<std::string>& second_list,
                                std::vector<std::string>& sam_list, std::vector<std::string>& profile_list,
                                std::vector<std::string>& samplenames_list, std::vector<std::string>& profile_truth_list,
                                std::vector<std::string>& read_type_list, std::vector<std::string>& unmapped_reads_list,
                                std::string const& sam_ending = kDefaultSamEnding, bool* sams_named = nullptr,
                                std::vector<std::string>* sparse_first_list = nullptr,
                                std::vector<std::string>* flex_indels_list = nullptr,
                                std::vector<std::string>* sparse_first_pick_list = nullptr);

        // Finds the database --db names (db::Locate) and opens it if it is a single file. --build always
        // writes into a folder.
        void ResolveDatabase(std::vector<std::string>& error_log) {
            m_location = db::Locate(m_database_path);
            m_bundle.reset();
            if (m_database_path.empty()) {
                error_log.emplace_back("No database given: give --db, or set $" + PROTAL_DB_ENV_VARIABLE);
                return;
            }
            if (!m_location.missing.empty()) {
                std::string where;
                std::error_code ec;
                auto const cwd = std::filesystem::current_path(ec);
                if (std::filesystem::path(m_database_path).is_relative() && !ec) where = " (relative to the working directory " + cwd.string() + ")";
                error_log.emplace_back("--db " + m_database_path + " " + m_location.missing + where);
                return;
            }
            if (m_build) {
                if (m_location.bundle == m_database_path && !m_database_path.empty()) {
                    error_log.emplace_back("--build writes a database into a folder: give --db the folder with reference.fna, "
                                           "reference.map and internal_taxonomy.dmp, not the file " + m_database_path);
                }
                return;
            }
            if (m_location.bundle.empty()) return;
            std::string error;
            auto bundle = db::Bundle::Open(m_location.bundle, error);
            std::string const journal = m_location.bundle + db::kJournalExtension;
            if (!bundle && std::filesystem::exists(journal)) {
                // --add_model stopped while it replaced the models in place (db::ReplaceTail); the journal has the old bytes.
                if (error.empty()) error = "it does not read as one";
                if (m_add_model.empty()) {
                    error += " (an interrupted protal --add_model left " + journal + ": run --add_model again, which first "
                             "writes the database's old content back from it)";
                } else if (std::string restore_error; db::RestoreFromJournal(m_location.bundle, restore_error)) {
                    std::cout << "Wrote the old content of " << m_location.bundle << " back from " << journal
                              << ", which an interrupted --add_model left" << std::endl;
                    error.clear();
                    bundle = db::Bundle::Open(m_location.bundle, error);
                } else {
                    error += "; writing its old content back from " + journal + " (an interrupted --add_model left it) failed: " +
                             restore_error;
                }
            }
            if (!bundle) {
                error_log.emplace_back(error.empty() ? "--db " + m_database_path + " is neither a folder nor a single-file protal database (" +
                                                       db::kFileName + ")"
                                                     : "Cannot read the database " + m_location.bundle + ": " + error);
                return;
            }
            m_bundle = std::make_shared<db::Bundle const>(std::move(*bundle));
        }

        // Why --preload_genomes_off cannot use a single-file database, with the commands that unpack it
        // and rerun this protal command on the unpacked files.
        std::string PreloadOffNeedsFilesMessage() const {
            std::string const threads = m_threads > 1 ? " -t " + std::to_string(m_threads) : "";
            return "--preload_genomes_off reads genes one by one from an uncompressed reference.fna, so it needs the database "
                   "as separate files, but " + m_bundle->Path() + " is a single-file database. Unpack it into " + m_location.dir +
                   " (protal then uses the separate files there; " + m_bundle->Path() + " is kept, and can be removed once the "
                   "files are unpacked):\n"
                   "    " + ShellWord(ProgramName()) + " --unpack_db --db " + ShellWord(m_bundle->Path()) + threads + "\n"
                   "and then rerun protal on the unpacked database:\n"
                   "    " + CommandWithDb(m_location.dir);
        }

        // The kind of reads of each sample: the one its READ_TYPE (or --read_type) names, else
        // paired-end with a second read file, and without one what the first reads of the file show
        // (ReadTypeDetection.h: single-end for short reads, PacBio or ONT for long ones, by their names
        // or qualities), noted in note_log for samples of this run, or warned of in warning_log where
        // the reads' qualities alone, or their lack, decided (DescribeReadType). With --profile_only, from
        // the SAM: the kind its header names (protal writes it), else single-end if its alignments
        // are unpaired (no 0x1); a SAM without usable alignments, or that cannot be read (reported
        // when it is profiled), counts as paired-end. A kind given for a SAM wins over its own, e.g. to
        // profile it with another read type's model, with a warning if they differ. A rerun that
        // profiles the SAM an earlier run wrote (it exists, no --force) takes the kind its header
        // names in the same way.
        void ResolveReadTypes(std::vector<std::string>& warning_log, std::vector<std::string>* note_log = nullptr) {
            m_read_types.assign(m_prefix_list.size(), ReadType::Paired);
            for (size_t i = 0; i < m_read_types.size(); i++) {
                auto const given = i < m_read_type_list.size() ? ReadTypeFromToken(m_read_type_list[i]) : std::nullopt;
                if (!m_profile_only) {
                    bool const single = i < m_second_list.size() && m_second_list[i].empty();
                    m_read_types[i] = given ? *given : single ? ReadType::Single : ReadType::Paired;
                    bool const in_run = std::find(m_range.begin(), m_range.end(), i) != m_range.end();
                    if (!given && single && in_run && i < m_first_list.size() && !m_first_list[i].empty()) {
                        auto const first_file = SplitReadFiles(m_first_list[i]).front();  // of several (lanes), the first
                        auto const guess = GuessReadType(SampleReads(first_file), MAX_SHORT_READ_LENGTH);
                        if (guess) {
                            m_read_types[i] = guess->type;
                            auto const message = DescribeReadType(*guess, i < m_sampleid_list.size() ? m_sampleid_list[i] : m_prefix_list[i],
                                                                  first_file);
                            if (message && message->warning) warning_log.push_back(message->text);
                            else if (message && note_log) note_log->push_back(message->text);
                        }
                    }
                    if (m_force || i >= m_sam_list.size()) continue;
                    // As RunProtal picks the SAM it skips the alignment for.
                    std::string sam = m_sam_list[i];
                    if (!std::filesystem::exists(sam)) sam = UncompressedSamName(sam);
                    if (!std::filesystem::exists(sam)) continue;
                    std::optional<ReadType> own;
                    try {
                        SamInput input(sam);
                        own = ReadsOfSam(input.Stream()).declared;
                    } catch (SamFormatError const&) {}
                    if (!own || *own == m_read_types[i]) continue;
                    if (given) {
                        warning_log.emplace_back(sam + " holds " + ReadTypeName(*own) + " reads; profiled as " +
                                                 ReadTypeName(*given) + " reads (" + Info(*given).token +
                                                 ", --read_type or READ_TYPE; --force aligns them again)");
                    } else {
                        warning_log.emplace_back(sam + " holds " + ReadTypeName(*own) + " reads (its header says), so they are "
                                                 "profiled as such, not as " + ReadTypeName(m_read_types[i]) +
                                                 " reads; --force aligns them again");
                        m_read_types[i] = *own;
                    }
                    continue;
                }
                std::optional<ReadType> own;  // the SAM's, if it tells
                if (i < m_sam_list.size() && std::filesystem::exists(m_sam_list[i])) {
                    SamInput input(m_sam_list[i]);
                    try {
                        auto const reads = ReadsOfSam(input.Stream());
                        if (reads.declared) {
                            own = *reads.declared;
                        } else if (reads.paired.has_value()) {
                            own = *reads.paired ? ReadType::Paired : ReadType::Single;
                        }
                    } catch (SamFormatError const&) {}
                }
                m_read_types[i] = given ? *given : own.value_or(ReadType::Paired);
                if (given && own && *given != *own) {
                    warning_log.emplace_back(m_sam_list[i] + " holds " + ReadTypeName(*own) + " reads; profiled as " +
                                             ReadTypeName(*given) + " reads (" + Info(*given).token + ", --read_type or READ_TYPE)");
                }
            }
        }

        // Why a read file cannot be read, or an empty string if it can. Checked with access(), not by
        // opening it: opening a FIFO would wait for its writer, and closing it would end its stream.
        static std::string Unreadable(std::string const& path) {
            std::error_code ec;
            if (std::filesystem::is_directory(path, ec)) return "it is a directory";
            if (::access(path.c_str(), R_OK) != 0) return std::strerror(errno);
            return {};
        }

        // The longest of the first `records` reads of a FASTQ/FASTA file (plain, gzip or zstd); 0 if it
        // cannot be read.
        static size_t LongestRead(std::string const& path, size_t records) {
            if (!std::filesystem::is_regular_file(path)) return 0;
            ThreadedGzIstream is(path.c_str());
            SeqReaderSE reader(is);
            FastxRecord record;
            size_t longest = 0;
            for (size_t n = 0; n < records && reader(record); n++) longest = std::max(longest, record.sequence.size());
            return longest;
        }

        bool PrepareAndCheckValidity(bool force_read_check=false);

        // The longest sample ID: it names files (misc/<sample>_runtime.tsv), and Linux and macOS file names hold 255 bytes.
        static constexpr size_t kMaxSampleIdBytes = 200;

        // Why a sample ID cannot name a file safely on Linux and macOS, or an empty string if it can. Sample IDs name
        // files (misc/<sample>_runtime.tsv) and the rows of the strain MSAs; whitespace is checked with the MSAs.
        static std::string SampleIdProblem(std::string const& id) {
            if (id.empty()) return "it is empty";
            if (id == "." || id == "..") return "'" + id + "' names a folder";
            if (id.find('/') != std::string::npos) return "it contains '/', which separates folders";
            if (id.find(':') != std::string::npos) return "it contains ':', which macOS shows as '/' and its older file systems refuse";
            if (std::any_of(id.begin(), id.end(), [](unsigned char c) { return c < 0x20 || c == 0x7f; })) {
                return "it contains a control character";
            }
            if (id.front() == '-') return "it starts with '-', which commands take for an option";
            if (id.size() > kMaxSampleIdBytes) {
                return "it is " + std::to_string(id.size()) + " bytes long; protal adds to it in file names, which hold 255 bytes "
                       "(at most " + std::to_string(kMaxSampleIdBytes) + ")";
            }
            // macOS (APFS) file names must be UTF-8.
            for (size_t i = 0; i < id.size();) {
                auto const c = static_cast<unsigned char>(id[i]);
                size_t const length = c < 0x80 ? 1 : (c >> 5) == 0x6 ? 2 : (c >> 4) == 0xe ? 3 : (c >> 3) == 0x1e ? 4 : 0;
                bool valid = length > 0 && i + length <= id.size();
                for (size_t k = 1; valid && k < length; k++) valid = (static_cast<unsigned char>(id[i + k]) >> 6) == 0x2;
                if (!valid) return "it is not UTF-8 text, which macOS file names must be";
                i += length;
            }
            return {};
        }

        // Whether a file name is a SAM's: .sam, .sam.gz or .sam.zst.
        static bool IsSamName(std::string const& name) {
            return UncompressedSamName(name).ends_with(".sam");
        }

        // The SAM files of --profile_only (its comma-separated items, then the arguments an unquoted pattern became): an
        // item that names a file is taken as it is; one with a wildcard (*, ? or [...]) stands for the SAM files it
        // matches (glob(3); .sam, .sam.gz and .sam.zst files only, so a run's .err, profile and .partial files are
        // passed over), sorted. A pattern that matches no SAM, and a folder, are errors; notes say what each pattern
        // matched. A missing file is left in the list, for the existence check to report with the others.
        static std::vector<std::string> ExpandSamFiles(std::vector<std::string> const& items, std::vector<std::string>& errors,
                                                       std::vector<std::string>& notes) {
            std::vector<std::string> sams;
            for (auto const& item : items) {
                if (item.empty()) continue;
                std::error_code ec;
                if (std::filesystem::is_directory(item, ec)) {
                    errors.emplace_back("--profile_only " + item + " is a folder: give its SAM files, e.g. '" +
                                        (std::filesystem::path(item) / "*.sam.zst").string() + "'");
                    continue;
                }
                if (std::filesystem::exists(item, ec) || item.find_first_of("*?[") == std::string::npos) {
                    sams.push_back(item);
                    continue;
                }
                glob_t matches{};
                int flags = 0;
#ifdef GLOB_TILDE
                flags |= GLOB_TILDE;  // ~/... in a quoted pattern
#endif
                int const status = glob(item.c_str(), flags, nullptr, &matches);
                std::vector<std::string> found;
                size_t others = 0;
                if (status == 0) {
                    for (size_t i = 0; i < matches.gl_pathc; i++) {
                        std::string const path = matches.gl_pathv[i];
                        if (IsSamName(path) && std::filesystem::is_regular_file(path, ec)) {
                            found.push_back(path);
                        } else {
                            others++;
                        }
                    }
                }
                globfree(&matches);
                std::sort(found.begin(), found.end());
                std::string const passed_over = others ? " (and " + std::to_string(others) + " other file(s) or folder(s), passed over)" : "";
                if (status != 0 && status != GLOB_NOMATCH) {
                    errors.emplace_back("--profile_only " + item + ": its folders cannot be read");
                } else if (found.empty()) {
                    errors.emplace_back("--profile_only " + item + " matches no SAM file (.sam, .sam.gz or .sam.zst)" + passed_over);
                } else {
                    notes.emplace_back("--profile_only " + item + " matches " + std::to_string(found.size()) + " SAM file(s)" + passed_over);
                }
                sams.insert(sams.end(), found.begin(), found.end());
            }
            return sams;
        }

        // The names of --profile_only's SAM files as samples: each file's name without .sam, .sam.gz or .sam.zst, unless
        // two files share one (aln.sam.zst in a folder per sample, or two studies' sa.sam.zst). Then each is named by the
        // folders in which the files' paths differ, joined by '_', followed by its file name unless all files share one:
        // /d/S1/aln.sam.zst and /d/S2/aln.sam.zst are S1 and S2; /p/study1/alignments/sa.sam.zst,
        // /p/study2/alignments/sa.sam.zst and /p/study2/alignments/sb.sam.zst are study1_sa, study2_sa and study2_sb.
        // The same file given twice keeps its name (the duplicate check reports it).
        static std::vector<std::string> SamSampleNames(std::vector<std::string> const& sams) {
            std::vector<std::string> stems;
            for (auto const& sam : sams) {
                std::string name = std::filesystem::path(UncompressedSamName(sam)).filename().string();
                if (name.ends_with(".sam")) name.resize(name.size() - 4);
                stems.push_back(name);
            }
            if (std::set<std::string>(stems.begin(), stems.end()).size() == stems.size()) return stems;

            std::vector<std::vector<std::string>> folders;  // of each file's absolute path
            size_t shortest = SIZE_MAX, longest = 0;
            for (auto const& sam : sams) {
                std::error_code ec;
                auto const absolute = std::filesystem::absolute(sam, ec).lexically_normal();
                std::vector<std::string> parts;
                for (auto const& part : absolute.parent_path()) {
                    if (!part.empty() && part != part.root_directory()) parts.push_back(part.string());
                }
                shortest = std::min(shortest, parts.size());
                longest = std::max(longest, parts.size());
                folders.push_back(std::move(parts));
            }
            std::vector<size_t> differ;  // the folder levels at which the paths differ
            for (size_t level = 0; level < longest; level++) {
                bool const differs = level >= shortest ||
                                     std::any_of(folders.begin(), folders.end(), [&](auto const& f) { return f[level] != folders.front()[level]; });
                if (differs) differ.push_back(level);
            }
            bool const one_stem = std::all_of(stems.begin(), stems.end(), [&](auto const& s) { return s == stems.front(); });
            std::vector<std::string> names;
            for (size_t i = 0; i < sams.size(); i++) {
                std::string name;
                for (auto level : differ) {
                    if (level < folders[i].size()) name += (name.empty() ? "" : "_") + folders[i][level];
                }
                if (!one_stem || name.empty()) name += (name.empty() ? "" : "_") + stems[i];
                names.push_back(name);
            }
            return names;
        }

        static std::vector<size_t> ProcessRange(std::string const& s, size_t const& total) {
            std::regex pattern("^\\d+(-\\d+)?(,\\d+(-(\\d+)?)?)*$");
            if (!std::regex_match(s, pattern)) {
                std::cout << "String " << s << " does not match the regex." << std::endl;
                exit(9);
            }

            std::vector<size_t> samples;
            std::vector<std::string> tokens;
            std::vector<std::string> subtokens;

            Utils::split(tokens, s, ",");

            for (auto& token : tokens) {
                if (token.empty()) continue;
                if (token.find('-') != std::string::npos) {
                    Utils::split(subtokens, token, "-");

                    auto begin = stoll(subtokens[0]);
                    auto end = !subtokens[1].empty() ? stoll(subtokens[1]) + 1 : total + 1;

                    if (begin >= total || end - 1 > total) {
                        std::cout << "There are only " << total << " samples." << std::endl;
                        exit(9);
                    }

                    for (auto i = begin; i < end; i++) {
                        samples.emplace_back(i-1);
                    }
                } else {
                    samples.emplace_back(stoll(token)-1);
                }
            }

            std::unordered_set<size_t> as_set(samples.begin(), samples.end());
            if (as_set.size() != samples.size()) {
                std::cout << "Ranges and numbers may not overlap." << std::endl;
                exit(9);
            }

            return samples;
        }

        static Options OptionsFromArguments(int argc, char *argv[]);
    };


}
