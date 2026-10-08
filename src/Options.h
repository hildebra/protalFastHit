//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <cxxopts.hpp>
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
#include "SpeciesPriors.h"
#include "SpeciesNeighbours.h"
#include "CongenerGapsTable.h"
#include "ForeignRatesTable.h"
#include "StrainAlleles.h"
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

    static cxxopts::Options CxxOptions() {
        cxxopts::Options options(
                "protal",
                "Protal help text");

        options.positional_help("Help");


        // I/O related options
        options.add_options("I/O")
                ("db", "Path to the protal database: a single-file database (database.protal, as --build writes it), or a folder holding one or the database's separate files. If not given, it is taken from the environment variable $" + PROTAL_DB_ENV_VARIABLE + ".", cxxopts::value<std::string>())
                ("1,first", "Comma separated list of read files: the first-in-pair files of paired-end reads (the second-in-pair files go to -2/--second), or single-end reads when -2/--second is not given.", cxxopts::value<std::string>()->default_value(""))
                ("2,second", "Comma separated list of second-in-pair read files, one per file given via -1/--first. Leave it out for single-end reads.", cxxopts::value<std::string>()->default_value(""))
                ("read_type", "The reads of -1/--first: pe (paired-end reads, with -2/--second), se (single-end reads), pb (PacBio long reads, e.g. HiFi), ont (Oxford Nanopore long reads). Without it, pe with -2/--second; with one read file, its first reads decide, and protal says what it found: se for short reads, for long reads pb or ont by their names (PacBio movie/ZMW, the UUIDs of MinKNOW and dorado), else by their median read quality (pb from Q25 on, ont below; pb without qualities). A map gives it per sample in a READ_TYPE column. With --profile_only, the SAM's header names its reads, and --read_type replaces that (e.g. to profile a SAM with another read type's model); with --add_model, it names the model's reads (pe if not given).", cxxopts::value<std::string>()->default_value(""))
                ("prefix", "Comma separated list of output prefixes (optional). If not specified, output file prefixes are generated from the input file names: the longest common prefix of the two files of paired-end reads (which must then be in the same folder), the file name without its FASTQ/FASTA and compression extensions for single-end reads.", cxxopts::value<std::string>()->default_value(""))
                ("o,outdir", "Overwrites #OUTPUT_DIR in map and needs to be defined if #OUTPUT_DIR is not defined in the map. If not otherwise specified in the map file, sam files, profiles, msas, and other miscellaneous files will be stored in the subfolders to this directory 'alignments', 'profiles', 'strains', and 'misc'. With -1/-2: the folder of the SAM files and profiles, with 'strains' and 'misc' in it (default: the current folder).", cxxopts::value<std::string>())
                
                ("map", "For larger datasets you can define parameters -1, -2, --prefix and -o in a tsv-file.", cxxopts::value<std::string>()->default_value(""))
                ("map_range", "If you specified a map file with --map you can also pass a range to protal to run protal only on a subset. The first entry is 1, the end is inclusive. e.g.: 1-10. If the end open or larger than the number of entries in the map file, the last entry in the map file is selected as end.", cxxopts::value<std::string>()->default_value(""))
                ("profile_only", "Comma separated list of existing sam files to profile without re-running the alignment, e.g. to build strain MSAs over the samples of several runs. An item with a wildcard, such as 'runs/*/alignments/*.sam.zst' (quoted: protal expands it; unquoted, the shell does, which works as well), stands for the .sam, .sam.gz and .sam.zst files it matches, sorted. Read files given via -1/-2 are ignored. Output prefixes are either given via --prefix (one per sam file) or derived from the sam file names (where two files share a name, as aln.sam.zst in a folder per sample, from the folders in which their paths differ); the outputs then go to -o if it is given, else the profiles next to each sam file and 'strains' and 'misc' into the current folder.", cxxopts::value<std::string>()->default_value(""))
                ("sam_format", "Format of the SAM files that protal names (with -1/-2, or a map without a SAM column): zst (zstd-compressed, <prefix>.sam.zst), gz (gzip, .sam.gz) or sam (plain). A map's SAM names keep their own ending.", cxxopts::value<std::string>()->default_value("zst"));

        // Alignment / algorithm options
        options.add_options("Alignment")
                ("c,align_top", "After seeding, anchor are sorted by quality passed to alignment. <take_top> specifies how many anchors should be aligned starting with the most promising anchor.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_ALIGN_TOP)))
                ("adaptive_candidates", "A short read whose best alignment is divergent (identity below 0.99) and whose seeds fit more taxa equally well than --align_top took (ZN) tries up to this many more of those anchors, of taxa of the best alignment's genus: a strain's own species that its seeds ranked below its congeners'. 0: never. Needs the taxonomy (it is loaded for profiling).", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_ADAPTIVE_CANDIDATES)))
                ("no_allele_scores", "Short reads: score each candidate alignment against the reference alone. By default, with a database that has strain alleles (strain_alleles.tsv, --build --strain_alleles), a candidate's score counts the read's differences from the best of its species' known alleles where those explain them (StrainAlleles.h), so a read of a known strain scores on its species as on the strain's own gene; the records keep their alignment to the reference. A database's models are trained with its alleles' scores: use this to measure them, not on a production run.")
                ("m,max_out", "Maximum alignments that should be outputted", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_OUT)))
                ("no_mate_guidance", "Paired-end reads: do not let a mate that is sure of its alignment guide the other one when they did not align together (to the guiding mate's taxon from the other's own anchor, or on its gene where the fragment can reach, partly if it runs past the gene's end).")
                ("no_gene_neighbours", "Do not use the database's gene neighbours (gene_neighbours.tsv: how often each marker gene end faces which other in a clade's genomes): no looking for a mate past the end of its guiding mate's gene, no pairs of mates on two neighbouring genes, no looking on a long read for the genes next to its genes, adjacent_expected_share and adjacent_unlikely_share 0 and adjacent_support 0.5.")
                ("keep_suspect_copies", "Count reads on the database's suspect gene copies (suspect_copies.tsv: a species' copy of a marker gene that is near-identical to the copy of a species of another genus, family, order, class, phylum or domain while its own congeners' copies are farther, which --build finds; a contaminating contig or a transferred gene) as evidence of the species. By default their records are left out as if the reads had not aligned: every present organism with such a gene puts a perfect read on the copy, a fifth of the false species calls at GTDB r226 (docs/claude/2026-10-03-false-positive-anatomy).")
                ("keep_foreign_genes", "Keep foreign genes in their taxon's depth and strain MSAs. A gene is foreign when the genes next to it on its reads (a pair's mates on two genes, a long read's consecutive genes) are mostly unlikely neighbours in its taxon's clade by the database's gene neighbours (4 or more such links judged, more than half of them unlikely): its reads come from another genome. By default foreign genes are left out of both; .profile.genes.log lists them either way.")
                ("u,max_key_ubiquity", "Max key ubiquity. Best matching Flexkey count for seed must be lower or equal", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_KEY_UBIQUITY)))
                ("s,max_seed_size", "Max seed size after which seeding is stopped.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_SEED_SIZE)))
                ("w,min_successful_lookups", "If the number of seeds is >=max_seed_size and the number of successful core-mer lookups is >= min_successful_lookups, stop looking for further seeds.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MIN_SUCCESSFUL_LOOKUPS)))
                ("a,max_score_ani", "A max score makes an alignment stop if the alignment diverges too much. This parameter estimates the score for a given ani and is a tradeoff between speed/accuracy. Given, it applies to all read types; else ONT reads take 0.85 (their indels count twice).", cxxopts::value<double>()->default_value(std::to_string(DEFAULT_MAX_SCORE_ANI)))
                ("x,x_drop", "X-drop of the alignment of short reads (WFA2), on top of its adaptive pruning: alignment branches whose score falls this far behind are cut. 0 turns it off. The default changes no 150 bp alignment in tests; small values (50) lose alignments and change MAPQs. Long reads (pb, ont) are aligned without it: over their gene-long windows even 1000 lost alignments.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_X_DROP)));

        // Profiling options
        options.add_options("Profiling")
                ("no_profile", "Do NOT perform taxonomic profiling, only output alignments.")
                ("knob", "Prediction threshold, 0 to 1: taxa whose model probability is at least this are reported. Lower finds more of the taxa present, higher reports fewer absent ones. How much a change matters depends on the model and the samples, so choose it on data like yours. Default 0.5, or, for a model with knobs by depth (machine_learning_cmdline.py --depth-knobs; build_gtdb_database.py gives them to its PacBio and ONT models), the model's knob for each sample's depth; a --knob given applies to every sample. A model's calibrated calls (machine_learning_cmdline.py --fdr-calls) are used only with --fdr.", cxxopts::value<double>()->default_value("0.5"))
                ("fdr", "Expected share of false calls, 0 to 1, for a model with calibrated calls (machine_learning_cmdline.py --fdr-calls; build_gtdb_database.py --call-mode fdr): each sample reports its highest-scoring taxa as long as the mean of their 1 - probability stays at or below this, the model's scores made probabilities by its calibration and adjusted to the share of present taxa among the sample's candidates (a deep sample, with many more absent candidates, needs higher scores). Off by default: without --fdr the model's knob curve or --knob applies even to a model with calibrated calls (at GTDB r226 they called 0.001-0.007 F1 below the curve for every read type). --fdr F uses them at F; 0: off. Not with --knob.", cxxopts::value<double>())
                ("singleton_congener", "A taxon of a single fragment beside a species of its genus with at least this many fragments in the sample is not reported, whatever its score, if its read looks like that species': the abundance-weighted assignment of the reads' alternatives leaves it less than half of it, or its identity is below 0.95. 0 (default): no such rule; on GTDB r226 training data the model called none of the taxa the rule would have vetoed, so it is off (docs/claude/2026-10-03-false-positive-anatomy). The training dump's prediction column and the trainer (machine_learning_cmdline.py --singleton-congener) apply the same rule.", cxxopts::value<size_t>()->default_value("0"))
                ("depth_identity_margin", "Reads count towards a species' abundance when their identity is at most this far below that of its best-matching reads (98th percentile). Reads below that, e.g. of a relative the database lacks, still count for detection. The margin is the same on every gene unless --gene_conservation scales it. The default, 0.08, counts the reads of strains up to about 5% from the reference, of which 0.04 dropped up to 40%. 1 lets every read count.", cxxopts::value<double>()->default_value("0.08"))
                ("gene_conservation", "Scale --depth_identity_margin per gene by how fast each gene diverges within species: db for the database's factors (gene_conservation.tsv, which --build estimates from --full_reference and stores in the database), or a file of them (geneid, factor, species; 1 for a gene of typical conservation). A gene's margin is then 0.03 for read errors plus the rest times its factor (0.08: 0.05 at factor 0.4, 0.10 at 1.4). none (default): the same margin on every gene, which did best summed over three simulated worlds, among them one with many congeners missing from the database. The factors (the file's, else the database's) give the model's conservation features (conserved_fast_depth_ratio, conserved_hit_share) whatever this says.", cxxopts::value<std::string>()->default_value("none"))
                ("model", "PMML model file: an existing path is used as is, otherwise <name> in the database (<name>.xml without an extension). Default: the database's model of each sample's read type: model_pe.xml (or, in older databases, model.xml) for paired-end, model_se.xml for single-end, model_PB.xml for PacBio and model_ONT.xml for ONT samples; --model replaces all of them unless --model_se, --model_pb or --model_ont is given.", cxxopts::value<std::string>()->default_value(""))
                ("model_se", "PMML model file for single-end samples, given as --model. Default: --model if given, else the database's model_se.xml.", cxxopts::value<std::string>()->default_value(""))
                ("model_pb", "PMML model file for PacBio samples, given as --model. Default: --model if given, else the database's model_PB.xml.", cxxopts::value<std::string>()->default_value(""))
                ("model_ont", "PMML model file for ONT samples, given as --model. Default: --model if given, else the database's model_ONT.xml.", cxxopts::value<std::string>()->default_value(""))
                ("profile_dir", "Override profile output directory. Takes precedence over the directory specified in the map file.", cxxopts::value<std::string>()->default_value(""))
                ("no_unknown_share", "Report the called species' abundances as shares of their summed depth (adding up to 1), as protal did before 2026-10-08. By default, where the database gives genome sizes (species_priors.tsv) and the SAM says how many reads were scanned, the profile ends with a line '?' of the genomes the called species do not explain (the reads they leave over, at their average genome size), and the abundances are shares of all genomes sequenced. <profile>.composition gives the numbers either way.");

        // Strain / SNP options
        options.add_options("Strains")
                ("no_strains", "Stay on species level: do not write strain MSAs or SNP tables. Variants are still called, as the model uses them, so profiles are the same with or without this flag.")
                ("no_phasing", "Write one strain MSA row per long-read sample (PacBio, ONT), the consensus of its reads. By default a sample whose reads show two or more strains of a species (on 3 genes or more) gets a row per strain, <sample>_hap1, _hap2, ..., the most abundant first, each called from its strain's reads: the reads are clustered by their alleles across the genes each one covers, and the blocks of genes that no read links are joined by the strains' shares of the reads where those tell which strain is which. <species>.haplotypes.tsv says how each block was phased.")
                ("strain_spill", "A folder (best on a local disk) for what each sample keeps for the strain MSAs: written there once its profile is written and read back one species at a time, so that the memory of many samples (e.g. --profile_only over several runs) does not grow with them. A file per sample, removed at the end; protal keeps a sample's in memory if its file cannot be written.", cxxopts::value<std::string>()->default_value(""))
                ("snp_min_cov", "Minimum number of reads supporting a variant to call a SNP.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MIN_SNP_COV)))
                ("snp_min_phred_sum", "Minimum cumulative phred score (sum of base qualities) across all supporting reads. Combined with --snp_min_mean_qual via OR: a variant passes quality if phred_sum >= snp_min_phred_sum OR mean_qual >= snp_min_mean_qual.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MIN_SNP_PHRED_SUM)))
                ("snp_min_mean_qual", "Minimum mean base quality across supporting reads. Combined with --snp_min_phred_sum via OR: a variant passes quality if mean_qual >= snp_min_mean_qual OR phred_sum >= snp_min_phred_sum. Note: at low coverage, --snp_min_cov is the binding constraint regardless.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MIN_SNP_MEAN_QUAL)))
                ("snp_min_af", "Minimum allele frequency for an allele (its reads / the reads with a base at the position), so also the least share of reads a second strain needs to show as an IUPAC code. Interacts with --snp_min_cov: below coverage = snp_min_cov/snp_min_af, the count filter is stricter. Given, it applies to all read types; else ONT reads take 0.2 (their errors put low-frequency alleles at many positions). Profiles do not depend on it.", cxxopts::value<double>()->default_value(std::to_string(DEFAULT_MIN_SNP_AF)))
                ("snp_no_strand", "Disable the strand-bias filter. By default an allele (the reference included) that is seen on one strand only fails where that is unlikely given the strands of all reads at the position (p < 0.05): with reads on both strands, an allele on just one of them is an artefact. Where the reads are from one strand, as often at low depth, it passes.")
                ("msa_min_hcov", "Minimum non-N/non-'-' bases required per sequence to keep it in the MSA.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MSA_MIN_HCOV)))
                ("msa_knob", "Model probability, 0 to 1, a sample's taxon needs for its reads to enter the taxon's strain MSA. Default: the sample's knob (--knob, or the model's knob for its depth), so the MSA holds the samples whose profile reports the taxon.", cxxopts::value<double>())
                ("msa_min_depth", "Reads a position needs to be written in a strain MSA; with fewer it is '-'. Where all its reads show one allele, that many suffice; a second allele (an IUPAC code) needs --snp_min_cov reads of its own, and a site whose reads disagree otherwise is N.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MSA_MIN_DEPTH)))
                ("msa_identity_margin", "Reads enter a species' strain MSA rows when their identity is at most this far below that of its best-matching reads (98th percentile), as for its abundance (--depth_identity_margin), but by default more strictly: a relative's reads within 0.08 of the best add false alleles (10 times the false calls before qcmsa, twice after). 1 lets every read in.", cxxopts::value<double>()->default_value("0.04"))
                ("msa_species", "Restrict MSAs to a single species (s__Genus_species) or a comma-separated list.", cxxopts::value<std::string>()->default_value(""))
                ("snp_max_alleles", "Maximum number of alleles at a position to encode as an IUPAC ambiguity code in the MSA. 1 = only the top allele (standard), 2 = encode two-allele mixtures (e.g. R,Y), 3 = also encode three-allele mixtures (e.g. B,H). Alleles are ranked by observation count; ties go to higher-quality allele.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_SNP_MAX_ALLELES)));

        // MSA post-filter. protal emits a raw MSA and hands it to qcmsa, which does
        // all gene/sample coverage gating, multi-allelicity outlier removal and site
        // cleanup. qcmsa owns these thresholds; protal only forwards.
        options.add_options("qcmsa")
                ("no_qcmsa", "Disable the qcmsa post-filter, leaving only protal's raw MSA. By default protal runs qcmsa (coverage gating, adaptive multi-allelicity outlier removal, site cleanup) on each species' raw MSA; it requires python3 and fails gracefully with a warning if python3 or the script cannot be found.")
                ("qcmsa_args", "Extra arguments forwarded verbatim to qcmsa, e.g. --qcmsa_args \"--preset strict --gene-min-hcov 0.5\". This is how you reach every qcmsa knob (aggressiveness preset, coverage thresholds, site cleanup) -- protal has no equivalents of its own and keeps qcmsa's defaults. Appended last, so they override what protal passes itself (--reapply-hcov). Run 'qcmsa --help' for the full list; protal does not validate them.", cxxopts::value<std::string>()->default_value(""))
                ("qcmsa_script", "Path to (or name of) the qcmsa executable. If empty, protal looks for 'qcmsa' next to its own binary and then on $PATH, honours the PROTAL_QCMSA_SCRIPT environment variable, and finally falls back to the qcmsa.py script of a source checkout.", cxxopts::value<std::string>()->default_value(""));

        
        // Advanced / benchmarking / build
        options.add_options("DevOptions")
                ("mapq_debug_output", "Output mapq debug info to stderr")
                ("whole_read_alignment", "Align each short read as a whole into its gene window, as protal did before it aligned from the anchor's exact matches (slower; the results differ in a few alignments). Long reads are always aligned as a whole.")
                ("no_alignment_screen", "Align every candidate with WFA2, without the k-mer screen that refuses a candidate whose read and gene window share too few k-mers for any alignment within the score budget to exist (AlignmentScreen.h). The screen changes no alignment; this is for measuring it.")
                ("long_read_budget", "Long reads: once a candidate of a read's gene has aligned, the gene's other candidates are aligned with a budget of this many edits (as mismatches) more than the best so far; one that would cost more fails there and counts as a failed candidate (ZF). Saves WFA2 work on far relatives; MAPQ and the listed alternatives of such genes change, so the models should be retrained with it. 0 (default): every candidate gets the ANI floor's budget.", cxxopts::value<size_t>()->default_value("0"))
                ("sequential_load", "Load the database's parts one after another, as protal did before 0.7.6, rather than the index beside the genome preload and the single-threaded tables (taxonomy, models, gene conservation, suspect copies, species priors, gene neighbours) each on a thread of its own. For measuring; the run is the same either way.")
                ("taxon_statistics", "Write misc/<taxon>.statistics.tsv for every taxon with reads: its coverage, reads, ANI and MAPQ in each sample, and whether it is reported. Off by default: a sample on a GTDB-sized database has reads on thousands of taxa, and that many small files took 15 s on a network file system. The profile files (<prefix>.profile, .profile.log, .profile.genes.log) hold the same per sample.")
                ("full_sam_header", "List every gene of the database in the SAM header (@SQ), as protal did before; by default only the genes that alignments name are listed.")
                ("write_unmapped_reads", "Write an unmapped record (FLAG 4, no sequence; its ZF tag names the taxa the read seeded on) for every read that seeded on taxa but aligned nowhere, as protal did up to 0.7.6. By default these reads are counted per taxon in the SAM header instead (@CO protal failed candidates of unaligned reads), which is all the profiler takes from them: on a GTDB-sized database they were 95% of a sample's records. --full_sam_header, whose header is written before the reads, writes them too. A map's UNMAPPED_READS column (write or count) chooses per sample.")
                ("serial_index_passes", "With --build: count and place the reference's k-mers and compute the value pointers on one thread, as protal did before these ran in -t threads (slower; the index is the same).")
                ("profile_ahead", "With several samples: profile each sample whose SAM is complete while the next sample's reads are aligned, on a worker with a quarter of the threads beside the alignment's; the profiling stage after the alignment takes the rest on all threads. The profiles are the same either way. Measured on a 6-core laptop it gained nothing (the worker's CPU time came out of the alignment's); it may pay on a node where the profiling stage leaves cores idle.")
                ("index_batch_kb", "With --build: KB of the reference per batch in the parallel index passes (smaller batches are for testing).", cxxopts::value<size_t>()->default_value("1024"))
                ("build", "Build index from reference file with header format ()")
                ("no_compress", "With --build: write the database as separate, uncompressed files (index.prx, reference.fna, ...). By default --build writes the single-file database database.protal (zstd-compressed; see --no_bundle). protal reads every form.")
                ("no_bundle", "With --build or --compress_db: keep the database as separate compressed files (index.prx.zst, reference.fna.zst, reference.map, internal_taxonomy.dmp, unique_kmers.tsv, gene_conservation.tsv, suspect_copies.tsv, species_priors.tsv, species_neighbours.tsv, congener_gaps.tsv, strain_alleles.tsv, foreign_rates.tsv, gene_neighbours.tsv, gene_positions.tsv, the models model_*.xml) instead of packing them into database.protal.")
                ("compress_level", "With --build: zstd compression level (1-22). Higher levels compress more but more slowly (level 19: ~3 MB/s per thread, -t threads are used); decompression speed barely depends on it.", cxxopts::value<int>()->default_value(std::to_string(DEFAULT_COMPRESS_LEVEL)))
                ("compress_window_log", "With --build: zstd long-distance matching window, as log2 bytes (27 = 128 MB, capped at the frame size); finds repeats between distant related sequences. 0 turns it off.", cxxopts::value<int>()->default_value(std::to_string(DEFAULT_COMPRESS_WINDOW_LOG)))
                ("compress_frame_mb", "With --build or --compress_db: size of the independent zstd frames in MB (1-4095). protal loads a database with -t threads, one frame per thread at a time. 0 writes a single frame, which loads with one thread.", cxxopts::value<int>()->default_value(std::to_string(DEFAULT_COMPRESS_FRAME_MB)))
                ("compress_db", "Compress the database folder --db in place, without rebuilding it: index.prx (raw or compressed in an older way) in protal's column format, reference.fna as seekable zstd (see --compress_level, --compress_frame_mb, -t), all packed into database.protal (--no_bundle: kept as index.prx.zst, reference.fna.zst, ...). Everything is read back and compared before the old files are removed. Needs the index in memory. On a single-file database --db: adds the binary gene table (gene_table.bin) a run loads instead of reference.map and unique_kmers.tsv, if it has no current one.")
                ("decompress_db", "Write the database --db as separate raw files (index.prx, reference.fna, ...) and remove its compressed files (index.prx.zst, reference.fna.zst, or database.protal), e.g. for older protal versions. zstd -d does not give a raw index.prx from protal's column format.")
                ("unpack_db", "Write the files of the single-file database --db into --unpack_dir (default: the folder it is in): index.prx.zst, an uncompressed reference.fna (as --preload_genomes_off needs), reference.map, internal_taxonomy.dmp, unique_kmers.tsv, gene_conservation.tsv, suspect_copies.tsv, species_priors.tsv, species_neighbours.tsv, congener_gaps.tsv, strain_alleles.tsv, foreign_rates.tsv, gene_neighbours.tsv and gene_positions.tsv (if it has them) and the models it has (model_pe.xml or model.xml, model_se.xml, model_PB.xml, model_ONT.xml). database.protal is kept; protal uses the separate files when both are there.")
                ("unpack_dir", "With --unpack_db: the folder to write the files into (default: the one database.protal is in).", cxxopts::value<std::string>()->default_value(""))
                ("add_model", "Store the PMML model FILE in the database --db as the model of the reads --read_type names (pe, se, pb or ont: model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml; pe if not given), replacing the one there. Several models, comma-separated, with as many read types (--add_model pe.xml,se.xml --read_type pe,se), are stored at once. Each model is checked first. database.protal is rewritten once, with its other parts copied as they are, not recompressed.", cxxopts::value<std::string>()->default_value(""))
                ("add_tables", "Store tables in the database --db by their file names, replacing those there: foreign_rates.tsv (scripts/foreign_rates.py: per gene copy, the reads of a tiled scan of the full reference's marker genes and how many came from other species), congener_gaps.tsv or strain_alleles.tsv (what --build writes) or species_priors.tsv (gtdb_to_protal_db.py --priors_only: the priors with the species' genome sizes, for a database converted before 2026-10-08). Several, comma-separated. Each table is read and checked first. database.protal is rewritten once, with its other parts copied as they are; a database of separate files gets the files copied beside them.", cxxopts::value<std::string>()->default_value(""))
                ("write_species_neighbours", "With --db FOLDER of separate files (a converted release before any build: reference.fna, reference.map, internal_taxonomy.dmp): write every species' nearest congeners by their references' marker genes (species_neighbours.tsv's table, which --build also stores in the database) to this file, without building. build_gtdb_database.py chooses the species its training database leaves out by it (species_clouds.tsv): a species complex is left out or kept whole.", cxxopts::value<std::string>()->default_value(""))
                ("full_reference", "With --build: the marker genes of all genomes (not only the representatives'), to check which k-mers are unique and to estimate how fast each gene diverges within species (gene_conservation.tsv, see --gene_conservation)", cxxopts::value<std::string>()->default_value(""))
                ("suspect_copy_distance", "With --build: a species' copy of a gene within this k-mer distance (about the share of bases that differ) of a copy of a species of another genus (or family, order, class, phylum, domain), and 0.02 farther from its nearest congener's copy or without one, is suspect: contamination or a transferred gene. The suspect copies go into the database (suspect_copies.tsv) and a run leaves their records out (see --keep_suspect_copies); every near pair across genera is reported in gene_incongruence.tsv beside the database. 0: no such scan.", cxxopts::value<double>()->default_value("0.02"))
                ("strain_alleles", "With --build: up to this many alleles per species and gene from the other genomes' copies in --full_reference, stored as edits of the representative's copy (strain_alleles.tsv), for the alignment scores (--no_allele_scores) and the profiler's 'alleles' features. Chosen farthest first among up to 16 distinct ones sampled by hash, each nearer the representative than the gene's nearest congener's copy (congener_gaps.tsv) and within 0.1 of it. 0: none.", cxxopts::value<size_t>()->default_value("4"))
                ("allele_genome_share", "With --build: the share of each species' genomes that may give alleles (--strain_alleles), by a hash of the genome's accession (the full reference's second header word, which the converter writes since 2026-10-08). build_gtdb_database.py passes 0.5 and simulates strains only from the other genomes, so that no simulated strain is its species' own allele. Below 1 it needs a full reference that names its genomes.", cxxopts::value<double>()->default_value("1"))
                ("reference", "Set of reference sequences to build the internal alignment database from", cxxopts::value<std::string>()->default_value(""))
                ("build_gene_subset", "With --build: a file of the gene ids (the numbers of reference.map / gene2geneid.tsv, one per line, # lines are comments) to build the database from; the other genes of --reference stay in reference.fna but get no k-mers, no unique k-mer row, no conservation factor and no suspect copies, and their copies in --full_reference are skipped. Refused when the folder has gene_neighbours.tsv (counted over every gene): derive a folder of the subset with scripts/mini_db/gtdb_to_protal_db.py --from_db --genes, which also leaves the other genes' sequences out (docs/databases.md, reduced marker sets).", cxxopts::value<std::string>()->default_value(""))
                ("preload_genomes_off", "Do not preload complete reference library (reference.fna and reference.map in protal index folder) and instead do dynamic loading. This usually decreases performance but saves memory. Needs the database as separate files with an uncompressed reference.fna (not reference.fna.zst or database.protal; see --unpack_db).")

                ("profile_truth", "Truth files, one per sample (comma-separated; a map's PROFILE_TRUTH column does the same). Each line names a species present, by GTDB lineage in any tab-separated field (d__...;s__Genus species, as simulate_metagenomes --protal_metafile writes) or by internal taxid in the first field. Profiles are annotated with TP/FP (<profile>.truth_annotated).", cxxopts::value<std::string>()->default_value(""))
                ("benchmark_alignment", "Benchmark alignment part of protal based on true taxonomic id and gene id supplied in the read header. Header must fulfill the formatting >taxid_geneid... with the regex: >[0-9]+_[0-9]+([^0-9]+.*)*")
                ("benchmark_alignment_output", "Benchmark alignment output. Output is appended to the file.", cxxopts::value<std::string>());


        // General options
        options.add_options("General")
                ("v,version", "Output version information")
                ("h,help", "Print help.")
                ("full_help", "Get help for developer options.")
                ("map_help", "Get help how to format the map file.")
                ("verbose", "Have verbose program output")
                ("t,threads", "Specify number of threads to use: for alignment (which also compresses the SAM), database loading and profiling.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_THREADS)))
                ("force", "Force redo alignment even if sam files exists. With --profile_only: write the profiles and strain outputs even where an earlier run left its own (without it, protal stops).");

        // Not shown: the arguments that follow no option. With --profile_only they are SAM files (an unquoted pattern,
        // which the shell expands into several arguments: the option takes the first), with --build the reference.
        options.add_options("_positional")
                ("positional", "Arguments that follow no option", cxxopts::value<std::vector<std::string>>());

        return options;
    }

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
        size_t strain_alleles = 4;
        double allele_genome_share = 1.0;
        bool no_phasing = false;
        std::string strain_spill;  // --strain_spill: the folder of the strain stage's spill files, or none
        bool fastalign = false;
        bool profile_only = false;
        bool force = false;
        bool verbose = false;
        bool mapq_debug_out = false;
        bool whole_read_alignment = false;
        bool no_alignment_screen = false;
        size_t long_read_budget = 0;
        bool sequential_load = false;
        bool taxon_statistics = false;
        bool full_sam_header = false;
        bool write_unmapped_reads = false;
        bool serial_index_passes = false;
        bool profile_ahead = false;
        bool no_unknown_share = false;
        size_t index_batch_kb = 1024;

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
        size_t m_strain_alleles = 4;  // --strain_alleles (--build)
        double m_allele_genome_share = 1.0;  // --allele_genome_share (--build)
        bool m_no_phasing = false;
        std::string m_strain_spill;  // --strain_spill
        bool m_fastalign = false;
        bool m_profile_only = false;
        bool m_force = false;
        bool m_verbose = false;

        bool m_mapq_debug_out = false;
        bool m_whole_read_alignment = false;
        bool m_no_alignment_screen = false;
        size_t m_long_read_budget = 0;
        bool m_sequential_load = false;
        bool m_taxon_statistics = false;
        bool m_full_sam_header = false;
        bool m_write_unmapped_reads = false;
        bool m_serial_index_passes = false;
        bool m_profile_ahead = false;
        bool m_no_unknown_share = false;
        size_t m_index_batch_kb = 1024;

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
        static inline const std::string PROTAL_SPECIES_NEIGHBOURS_FILE = species_neighbours::kFileName;
        static inline const std::string PROTAL_CONGENER_GAPS_FILE = congener_gaps::kFileName;
        static inline const std::string PROTAL_STRAIN_ALLELES_FILE = strain_alleles::kFileName;
        static inline const std::string PROTAL_FOREIGN_RATES_FILE = foreign_rates::kFileName;
        static inline const std::string PROTAL_GENE_NEIGHBOURS_FILE = gene_neighbours::kFileName;
        static inline const std::string PROTAL_GENE_POSITIONS_FILE = gene_neighbours::kPositionsFileName;
        static inline const std::string PROTAL_TAXONOMY_FILE = "internal_taxonomy.dmp";
        // The presence models of paired-end and of single-end reads.

        static inline const std::string MAP_SAMPLEID = "#SAMPLEID";
        static inline const std::string MAP_FIRST_READ = "FIRST";
        static inline const std::string MAP_SECOND_READ = "SECOND";
        static inline const std::string MAP_NO_SECOND_READ = "-";  // SECOND of a single-end sample
        static inline const std::string MAP_READ_TYPE = "READ_TYPE";
        // A sample's reads that seeded on taxa but aligned nowhere: an unmapped record each (write) or counted per taxon
        // in the SAM header (count); '-': as --write_unmapped_reads says for every sample.
        static inline const std::string MAP_UNMAPPED_READS = "UNMAPPED_READS";
        static inline const std::string MAP_UNMAPPED_WRITE = "write";
        static inline const std::string MAP_UNMAPPED_COUNT = "count";
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
                m_long_read_budget(d.long_read_budget),
                m_sequential_load(d.sequential_load),
                m_taxon_statistics(d.taxon_statistics),
                m_full_sam_header(d.full_sam_header),
                m_write_unmapped_reads(d.write_unmapped_reads),
                m_serial_index_passes(d.serial_index_passes),
                m_profile_ahead(d.profile_ahead),
                m_no_unknown_share(d.no_unknown_share),
                m_index_batch_kb(d.index_batch_kb),
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
                m_model(std::move(d.model)),
                m_model_se(std::move(d.model_se)),
                m_model_pb(std::move(d.model_pb)),
                m_model_ont(std::move(d.model_ont)),
                m_read_type(std::move(d.read_type)),
                m_add_model(std::move(d.add_model)),
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
            m_strain_spill = std::move(d.strain_spill);
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

        std::string ToString() const {

            std::string first_list_str = m_first_list.empty() ? "" : m_first_list.front();
            for (auto i = 1; i < m_first_list.size(); i++) first_list_str += ", " + m_first_list[i];
            std::string second_list_str = m_second_list.empty() ? "" : m_second_list.front();
            for (auto i = 1; i < m_second_list.size(); i++) second_list_str += ", " + m_second_list[i];
            std::string sam_list_str = m_sam_list.empty() ? "" : m_sam_list.front();
            for (auto i = 1; i < m_sam_list.size(); i++) sam_list_str += ", " + m_sam_list[i];
            std::string prefix_list_str = m_prefix_list.empty() ? "" : m_prefix_list.front();
            for (auto i = 1; i < m_prefix_list.size(); i++) prefix_list_str += ", " + m_prefix_list[i];
            std::string profile_list_str = m_profile_list.empty() ? "" : m_profile_list.front();
            for (auto i = 1; i < m_profile_list.size(); i++) profile_list_str += ", " + m_profile_list[i];
            // A value that reads of some kinds have their own of (ReadTypeInfo): the option's, then
            // the own one of each kind in the run that differs, as it is used (ONT: -a 0.85).
            auto per_read_type = [this](double value, auto const& of) {
                std::string own;
                for (auto const& info : kReadTypes) {
                    if (std::find(m_read_types.begin(), m_read_types.end(), info.type) == m_read_types.end()) continue;
                    if (of(info.type) != value) own += (own.empty() ? "" : ", ") + info.name + " reads: " + std::to_string(of(info.type));
                }
                return std::to_string(value) + (own.empty() ? "" : " (" + own + ")");
            };
            bool const long_reads = std::any_of(m_read_types.begin(), m_read_types.end(), IsLongReadType);

            std::ostringstream result_str;
            result_str << "------ General ------" << std::string(30, '-') << '\n';
            result_str << "build:               " << std::to_string(m_build) << '\n';
            if (m_build || m_compress_db) {
                result_str << "compress database:   " << (m_compress || m_compress_db ? "zstd level " + std::to_string(m_compress_level) +
                        (m_compress_window_log ? ", window 2^" + std::to_string(m_compress_window_log) : "") +
                        (m_compress_frame_mb ? ", " + std::to_string(m_compress_frame_mb) + " MB frames" : ", one frame") +
                        (WriteBundle() ? ", single file " + db::kFileName : ", separate files") : "no") << '\n';
            }
            result_str << "no strains:          " << std::to_string(m_no_strains) << '\n';
            result_str << "threads:             " << std::to_string(m_threads) << '\n';
            result_str << "-------- I/O --------" << std::string(30, '-') << '\n';
            result_str << "first:               " << (first_list_str.length() > 50 ? std::to_string(m_first_list.size()) + " files" : first_list_str) << '\n';
            result_str << "second:              " << (second_list_str.length() > 50 ? std::to_string(m_second_list.size()) + " files" : second_list_str) << '\n';
            if (!m_build && !m_read_types.empty()) {
                result_str << "read types:          ";
                for (auto const& info : kReadTypes) {
                    result_str << (info.type == kReadTypes.front().type ? "" : ", ")
                               << std::count(m_read_types.begin(), m_read_types.end(), info.type) << ' ' << info.name;
                }
                result_str << " sample(s)" << '\n';
            }
            result_str << "db path:             " << m_database_path << '\n';
            if (!m_build) {
                result_str << "database:            " << (m_bundle ? "single file " + m_bundle->Path() : "separate files in " + m_location.dir)
                           << (m_location.unused_bundle.empty() ? "" : " (not " + m_location.unused_bundle + ")") << '\n';
            }
            result_str << "sequence file:       " << m_sequence_file << '\n';
            result_str << "sam file:            " << (sam_list_str.length() > 50 ? std::to_string(m_sam_list.size()) + " files" : sam_list_str) << '\n';
            result_str << "profile file:        " << (profile_list_str.length() > 50 ? std::to_string(m_profile_list.size()) + " files" : profile_list_str) << '\n';
            result_str << "output prefix:       " << (prefix_list_str.length() > 50 ? std::to_string(m_prefix_list.size()) + " files" : prefix_list_str) << '\n';
            result_str << "Output dir:          " << m_output_dir << '\n';
            result_str << "preload genomes:     " << (m_preload_genomes ? "yes" : "no") << '\n';
            result_str << "----- Alignment -----" << std::string(30, '-') << '\n';
            result_str << "align top:           " << std::to_string(m_align_top) << '\n';
            result_str << "adaptive candidates: " << std::to_string(m_adaptive_candidates) << '\n';
            result_str << "allele scores:       " << (m_no_allele_scores ? "off (--no_allele_scores)" : "the database's strain alleles, if it has them") << '\n';
            result_str << "max key ubiquity:    " << std::to_string(m_max_key_ubiquity) << '\n';
            result_str << "max seed size:       " << std::to_string(m_max_seed_size) << '\n';
            result_str << "max score ani:       " << per_read_type(m_max_score_ani, [this](ReadType type) { return GetMaxScoreAni(type); }) << '\n';
            result_str << "x-drop:              " << std::to_string(m_x_drop) << (long_reads ? " (short reads; long reads: none)" : "") << '\n';
            result_str << "short reads aligned: " << (m_whole_read_alignment ? "as a whole" : "from their anchors") << '\n';
            result_str << "alignment screen:    " << (m_no_alignment_screen ? "off (--no_alignment_screen)" : "k-mers shared with the window before WFA2") << '\n';
            if (long_reads) result_str << "long read budget:    " << (m_long_read_budget ? std::to_string(m_long_read_budget) + " edits past the best candidate (--long_read_budget)" : "the ANI floor's for every candidate") << '\n';
            result_str << "SAM header lists:    " << (m_full_sam_header ? "every gene" : "the genes aligned to") << '\n';
            result_str << "unaligned reads:     " << (WriteUnmappedReads() ? "an unmapped record each" : "counted per taxon in the SAM header")
                       << (m_unmapped_reads_list.empty() || m_full_sam_header ? "" : " (per sample: the map's UNMAPPED_READS)") << '\n';
            result_str << "fastalign:           " << std::to_string(m_fastalign) << '\n';
            result_str << "max out:             " << std::to_string(m_max_out) << '\n';
            result_str << "------ Strains ------" << std::string(30, '-') << '\n';
            result_str << "snp min cov:         " << std::to_string(m_snp_min_cov) << '\n';
            result_str << "snp min phred sum:   " << std::to_string(m_snp_min_phred_sum) << '\n';
            result_str << "snp min mean qual:   " << std::to_string(m_snp_min_mean_qual) << '\n';
            result_str << "snp min af:          " << per_read_type(m_snp_min_af, [this](ReadType type) { return GetSNPMinAF(type); }) << '\n';
            result_str << "snp require strand:  " << (m_snp_require_strand ? "yes" : "no (--snp_no_strand)") << '\n';
            result_str << "snp max alleles:     " << std::to_string(m_snp_max_alleles) << '\n';
            result_str << "run qcmsa:           " << (m_run_qcmsa ? "yes" : "no") << '\n';
            if (m_run_qcmsa && !m_qcmsa_args.empty())
                result_str << "qcmsa extra args:    " << m_qcmsa_args << '\n';
            result_str << "msa species:         " << Utils::join(m_msa_species, ",") << '\n';
            result_str << "msa min hcov:        " << std::to_string(m_msa_min_hcov) << '\n';
            result_str << "msa min depth:       " << std::to_string(m_msa_min_depth) << '\n';
            result_str << "msa identity margin: " << std::to_string(m_msa_identity_margin) << '\n';
            result_str << "msa knob:            " << std::to_string(GetMSAKnob()) << (m_msa_knob ? "" : " (--knob)") << '\n';
            result_str << "phasing:             " << (m_no_phasing ? "no" : "long-read samples") << '\n';
            result_str << "foreign genes:       " << (m_keep_foreign_genes ? "kept" : "left out of depth and MSAs") << '\n';
            result_str << "suspect gene copies: " << (m_keep_suspect_copies ? "kept" : "left out of the evidence") << '\n';
            result_str << "---- Dev Options ----" << std::string(30, '-') << '\n';
            result_str << "verbose:             " << (m_verbose ? "yes" : "no") << '\n';
            result_str << "benchmark alignment: " << (m_benchmark_alignment ? "yes" : "no") << '\n';
            result_str << "profile truth:       " << (HasProfileTruths() ? "yes" : "no") << '\n';
            result_str << "^ output:            " << m_benchmark_alignment_output << '\n';
            result_str << "---------------------" << std::string(30, '-') << '\n';
            return result_str.str();
        }

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

        // The tables --add_tables may store, by file name.
        static std::vector<std::string> AddableTables() {
            return { PROTAL_CONGENER_GAPS_FILE, PROTAL_FOREIGN_RATES_FILE, PROTAL_SPECIES_PRIORS_FILE, PROTAL_STRAIN_ALLELES_FILE };
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

        // The models --add_model stores, each with the reads it is for: its files, comma-separated, and --read_type's
        // tokens at the same places (a single model: pe if --read_type is not given). Checked in PrepareAndCheckValidity;
        // a read type that is no token is paired-end here.
        std::vector<std::pair<std::string, ReadType>> AddModels() const {
            auto split = [](std::string const& text) {
                std::vector<std::string> items;
                std::string item;
                std::istringstream is(text);
                while (std::getline(is, item, ',')) items.push_back(item);
                return items;
            };
            auto const files = split(m_add_model);
            auto const types = split(m_read_type);
            std::vector<std::pair<std::string, ReadType>> models;
            for (size_t i = 0; i < files.size(); i++) {
                models.emplace_back(files[i], i < types.size() ? ReadTypeFromToken(types[i]).value_or(ReadType::Paired)
                                                               : ReadType::Paired);
            }
            return models;
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

        // Whether a long-read sample's strain MSA row is split into its strains' (Haplotypes.h).
        bool Phasing() const {
            return !m_no_phasing;
        }

        // --strain_spill: the folder of the strain stage's spill files, or empty (in memory).
        std::string const& StrainSpill() const {
            return m_strain_spill;
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

        std::string GetFirstFile(int index) const {
            if (index > m_first_list.size()) {
                std::cerr << "Cannot access index " << index << " of first files (Length: " << m_first_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_first_list[index];
        }

        std::string GetSecondFile(int index) const {
            if (index > m_second_list.size()) {
                std::cerr << "Cannot access index " << index << " of second files (Length: " << m_second_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_second_list[index];
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

        // --long_read_budget: edits past a segment's best hit its other candidates may cost (0: the ANI floor's budget).
        size_t GetLongReadBudget() const {
            return m_long_read_budget;
        }

        // --sequential_load: the database's parts loaded one after another (RunProtal::Run).
        bool SequentialLoad() const {
            return m_sequential_load;
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

        void PrintHelp(bool show_dev=false) {
            // print groups in desired order, skip groups that start with '_' (hidden)
            auto opt = CxxOptions();
            std::vector<std::string> groups = { "I/O", "Profiling", "Strains", "qcmsa", "General" };

            if (show_dev) {
                groups.push_back("Alignment");
                groups.push_back("DevOptions");
            }

            // cxxopts::Options::help can accept a vector of groups to print in that order
            std::cout << opt.help(groups) << std::endl;
        }

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
SAM and PROFILE are optional and default to <PREFIX>.sam.zst (--sam_format) and <PREFIX>.profile,
in OUTPUT_DIR. Every sample needs its own SAM and PROFILE file; protal stops if two samples share
one. A sample whose SAM exists is not aligned again (unless --force), and a SAM given as an
absolute path may lie anywhere, e.g. in another run's folder: a map of existing SAMs (protal_map_utils
merge writes one from the runs' maps) builds strain MSAs over the samples of several runs.)" << std::endl;
        }


        static bool CreateDir(std::string dir_path) {
            // trim whitespace (including CR) from both ends to avoid false negatives on existing dirs
            auto trim_ws = [](std::string &s) {
                auto not_space = [](unsigned char c) { return !std::isspace(c); };
                s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
                s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
            };
            trim_ws(dir_path);
            std::filesystem::path p(dir_path);

            if (!std::filesystem::exists(p)) {
                std::error_code ec;
                if (!std::filesystem::create_directories(p, ec)) {
                    std::cerr << "Failed to create SAM output directory: " << p << "\nError: " << ec.message() << std::endl;
                    return false;
                }
            } else if (!std::filesystem::is_directory(p)) {
                std::cerr << "SAM output path exists but is not a directory: " << p << std::endl;
                return false;
            }
            return true;
        }

        // The ending of the SAM files protal names, by --sam_format (zst, gz, sam); empty if unknown.
        static std::string SamEnding(std::string const& format) {
            if (format == "zst") return ".sam.zst";
            if (format == "gz") return ".sam.gz";
            if (format == "sam") return ".sam";
            return "";
        }
        static inline const std::string kDefaultSamEnding = ".sam.zst";

        // sam_ending: of the SAM names given to samples without a SAM column (<prefix><sam_ending>).
        static bool LoadFromMap(std::string map_path, std::string& output_dir, std::string& strain_output_dir, std::string& misc_output_dir,
                                std::vector<std::string>& prefix_list,
                                std::vector<std::string>& first_list, std::vector<std::string>& second_list,
                                std::vector<std::string>& sam_list, std::vector<std::string>& profile_list,
                                std::vector<std::string>& samplenames_list, std::vector<std::string>& profile_truth_list,
                                std::vector<std::string>& read_type_list, std::vector<std::string>& unmapped_reads_list,
                                std::string const& sam_ending = kDefaultSamEnding) {
            using namespace std::filesystem;

            if (!std::filesystem::exists(map_path)) {
                std::cerr << "Map file " << map_path << " does not exist." << std::endl;
                return false;
            }
            std::ifstream is(map_path, std::ios::in);

            std::string global_output_dir = "";
            std::string map_output_dir = "";
            std::string sam_output_dir = "";
            std::string profile_output_dir = "";

            std::string input_dir = "";

            LineSplitter splitter;

            int sample_id_column = 0;
            int first_column = -1;
            int second_column = -1;
            int prefix_column = -1;
            int sam_column = -1;
            int profile_column = -1;
            int profile_truth_column = -1;
            int header_to_species_column = -1;
            int read_type_column = -1;
            int unmapped_reads_column = -1;

            bool header = true;
            size_t line_num = 0;
            std::vector<std::string> bad_sample_ids;
            for (std::string line; std::getline(is, line);) {
                line_num++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                splitter.Split(line);
                auto& tokens = splitter.Tokens();
                
                
                // Header
                if (line.size() >= 1 && line.compare(0, 1, "#") == 0) {
                    if (!header) {
                        std::cerr << "Line " << line_num << ": Did not expect header line but line starts with #" << std::endl;
                        return false;
                    }

                    // Check if variable definition or header
                    if (!tokens.empty() && tokens[0] == MAP_SAMPLEID) {
                        // Column headers
                        splitter.Split(line);

                        for (auto i = 1; i < splitter.Tokens().size(); i++) {
                            auto& token = splitter.Tokens()[i];
                            if (token == MAP_SAMPLEID) continue;
                            if (token == MAP_FIRST_READ) {
                                if (first_column != -1) {
                                    std::cerr << "Column '" << MAP_FIRST_READ << "' is defined twice" << std::endl;
                                    return false;
                                }
                                first_column = i;
                            }
                            if (token == MAP_SECOND_READ) {
                                if (second_column != -1) {
                                    std::cerr << "Column '" << MAP_SECOND_READ << "' is defined twice" << std::endl;
                                    return false;
                                }
                                second_column = i;
                            }
                            if (token == MAP_PREFIX) {
                                if (prefix_column != -1) {
                                    std::cerr << "Column '" << MAP_PREFIX << "' is defined twice" << std::endl;
                                    return false;
                                }
                                prefix_column = i;
                            }
                            if (token == MAP_SAM) {
                                if (sam_column != -1) {
                                    std::cerr << "Column '" << MAP_SAM << "' is defined twice" << std::endl;
                                    return false;
                                }
                                sam_column = i;
                            }
                            if (token == MAP_PROFILE) {
                                if (profile_column != -1) {
                                    std::cerr << "Column '" << MAP_PROFILE << "' is defined twice" << std::endl;
                                    return false;
                                }
                                profile_column = i;
                            }
                            if (token == MAP_PROFILE_TRUTH) {
                                if (profile_truth_column != -1) {
                                    std::cerr << "Column '" << MAP_PROFILE_TRUTH << "' is defined twice" << std::endl;
                                    return false;
                                }
                                profile_truth_column = i;
                            }
                            if (token == MAP_READ_TYPE) {
                                if (read_type_column != -1) {
                                    std::cerr << "Column '" << MAP_READ_TYPE << "' is defined twice" << std::endl;
                                    return false;
                                }
                                read_type_column = i;
                            }
                            if (token == MAP_UNMAPPED_READS) {
                                if (unmapped_reads_column != -1) {
                                    std::cerr << "Column '" << MAP_UNMAPPED_READS << "' is defined twice" << std::endl;
                                    return false;
                                }
                                unmapped_reads_column = i;
                            }
                            if (token == MAP_HEADER_TO_SPECIES) {
                                if (header_to_species_column != -1) {
                                    std::cerr << "Column '" << MAP_HEADER_TO_SPECIES << "' is defined twice" << std::endl;
                                    return false;
                                }
                                header_to_species_column = i;
                            }
                        }
                        header = false;

                        if (!output_dir.empty()) {
                            global_output_dir = output_dir;
                        } else if (!map_output_dir.empty()) {
                            global_output_dir = map_output_dir;
                        } else {
                            std::cerr << "Line " << line_num << ": Output directory not defined. Please define " << MAP_VAR_OUTPUT_DIR << " before the header line or provide an output directory via --outdir." << std::endl;
                            return false;
                        }
                        output_dir = global_output_dir;

                        // A folder the map names (#SAM_OUTPUT_DIR, ...) lies in the output folder unless it is absolute,
                        // one it does not name is <output folder>/<default>; the output folder (-o, else #OUTPUT_DIR) is
                        // relative to the folder protal runs in. Up to 2026-10-08 a relative output folder was put in front
                        // of a default folder twice (-o out wrote out/out/alignments).
                        auto in_output = [&global_output_dir](std::string const& named, std::string const& default_name) {
                            path const dir = named.empty() ? path(default_name) : path(named);
                            return (dir.is_absolute() ? dir : path(global_output_dir) / dir).string();
                        };
                        sam_output_dir = in_output(sam_output_dir, MAP_VAR_DEFAULT_SAM_OUTPUT_DIR);
                        CreateDir(sam_output_dir);

                        strain_output_dir = in_output(strain_output_dir, MAP_VAR_DEFAULT_STRAIN_OUTPUT_DIR);
                        CreateDir(strain_output_dir);

                        profile_output_dir = in_output(profile_output_dir, MAP_VAR_DEFAULT_PROFILE_OUTPUT_DIR);
                        CreateDir(profile_output_dir);

                        misc_output_dir = in_output(misc_output_dir, MAP_VAR_DEFAULT_MISC_OUTPUT_DIR);
                        CreateDir(misc_output_dir);


                        continue;
                    }

                    // Header variable definitions
                    if (!tokens.empty() && tokens[0] == MAP_VAR_INPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_INPUT_DIR << std::endl;
                            return false;
                        }
                        input_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_STRAIN_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_STRAIN_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        strain_output_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_SAM_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_SAM_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        sam_output_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_PROFILE_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_PROFILE_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        profile_output_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        map_output_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_MISC_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_MISC_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        misc_output_dir = tokens[1];
                    }


                } else {

                    // Mapping file content. Without a SECOND column, all samples are single-end.
                    if (prefix_column == -1 || first_column == -1) {
                        std::cerr << "The columns must be specified: " << MAP_PREFIX << ", " << MAP_FIRST_READ
                                  << " (and " << MAP_SECOND_READ << " for paired-end reads)" << std::endl;
                        return false;
                    }

                    if (header) {
                        std::cerr << "Line " << line_num << ": Expected header line starting with #" << std::endl;
                        return false;
                    }

                    // Every column the header declares needs a value in every row: a missing cell would
                    // otherwise shift the per-sample lists and give one sample another's output files.
                    std::pair<int, std::string const*> const columns[] = {
                            { 0, &MAP_SAMPLEID }, { prefix_column, &MAP_PREFIX }, { first_column, &MAP_FIRST_READ }, { second_column, &MAP_SECOND_READ },
                            { sam_column, &MAP_SAM }, { profile_column, &MAP_PROFILE }, { profile_truth_column, &MAP_PROFILE_TRUTH },
                            { read_type_column, &MAP_READ_TYPE }, { unmapped_reads_column, &MAP_UNMAPPED_READS } };
                    for (auto const& [column, name] : columns) {
                        if (column == -1) continue;
                        if (static_cast<size_t>(column) >= tokens.size() || tokens[column].empty()) {
                            std::cerr << "Line " << line_num << ": no value in column " << column + 1 << " (" << *name
                                      << "); every row needs a value for each column of the header. "
                                         "Columns must be tab separated." << std::endl;
                            return false;
                        }
                    }
                    auto sample_id = tokens[0];  // #SAMPLEID
                    // It names files (misc/<sample>_runtime.tsv): every row is checked, then all problems reported.
                    if (auto const problem = SampleIdProblem(sample_id); !problem.empty()) {
                        bad_sample_ids.push_back("Line " + std::to_string(line_num) + ": sample ID '" + sample_id +
                                                 "' cannot name a file safely on Linux and macOS: " + problem);
                    }
                    auto prefix_path = path(global_output_dir).append(tokens[prefix_column]);
                    auto first_path = path(input_dir).append(tokens[first_column]);
                    // No second file ('-', or no SECOND column): single-end reads.
                    bool const single_end = second_column == -1 || tokens[second_column] == MAP_NO_SECOND_READ;
                    auto second_path = single_end ? path() : path(input_dir).append(tokens[second_column]);

                    // mandatory
                    samplenames_list.emplace_back(sample_id);
                    prefix_list.emplace_back(prefix_path);
                    first_list.emplace_back(first_path);
                    second_list.emplace_back(second_path);

                    // optional
                    if (sam_column != -1 && sam_column < tokens.size()) {
                        auto sam_path = path(sam_output_dir).append(tokens[sam_column]);
                        sam_list.emplace_back(sam_path);
                    }
                    if (profile_column != -1 && profile_column < tokens.size()) {
//                        std::cout << "Profile paths" << std::endl;
//                        std::cout << profile_output_dir << " " << tokens[profile_column] << std::endl;
                        auto profile_path = path(profile_output_dir).append(tokens[profile_column]);
                        profile_list.emplace_back(profile_path);
                    }
                    if (profile_truth_column != -1 && profile_truth_column < tokens.size()) {
                        auto profile_truth_path = tokens[profile_truth_column];
                        profile_truth_list.emplace_back(profile_truth_path);
                    }
                    if (read_type_column != -1) {
                        read_type_list.emplace_back(tokens[read_type_column]);
                    }
                    if (unmapped_reads_column != -1) {
                        auto const& value = tokens[unmapped_reads_column];
                        if (value != MAP_UNMAPPED_WRITE && value != MAP_UNMAPPED_COUNT && value != MAP_NO_SECOND_READ) {
                            std::cerr << "Line " << line_num << ": " << MAP_UNMAPPED_READS << " is '" << value << "'; expected "
                                      << MAP_UNMAPPED_WRITE << ", " << MAP_UNMAPPED_COUNT << " or -" << std::endl;
                            return false;
                        }
                        unmapped_reads_list.emplace_back(value);
                    }
                }
            }


            if (!bad_sample_ids.empty()) {
                for (auto const& problem : bad_sample_ids) std::cerr << problem << std::endl;
                std::cerr << "Rename these samples in the #SAMPLEID column: sample IDs name files and the rows of the strain "
                             "MSAs (no '/' or ':', control characters or a leading '-', UTF-8, at most "
                          << kMaxSampleIdBytes << " bytes)." << std::endl;
                return false;
            }

            // Fill in gaps
            if (profile_list.empty()) {
                if (prefix_list.empty()) {
                    std::cout << "column PREFIX must be specified if column PROFILE is not specified." << std::endl;
                }
                for (auto prefix : prefix_list) {
                    profile_list.emplace_back(prefix + ".profile");
                }
            }
            if (sam_list.empty()) {
                if (prefix_list.empty()) {
                    std::cout << "column PREFIX must be specified if column SAM is not specified." << std::endl;
                }
                for (auto prefix : prefix_list) {
                    sam_list.emplace_back(prefix + sam_ending);
                }
            }


//            std::cout << "strain_output_dir " << strain_output_dir << std::endl;
//            std::cout << "sam_output_dir " << sam_output_dir << std::endl;
//            std::cout << "profile_output_dir " << sam_output_dir << std::endl;
//            std::cout << "input_dir " << input_dir << std::endl;
//
//            std::cout << "\nPrefixes" << std::endl;
//            for (auto& prefix : prefix_list) {
//                std::cout << prefix << std::endl;
//            }
//
//            std::cout << "\nFirsts" << std::endl;
//            for (auto& first : first_list) {
//                std::cout << first << std::endl;
//            }
//
//            std::cout << "\nSeconds" << std::endl;
//            for (auto& second : second_list) {
//                std::cout << second << std::endl;
//            }
//
//            std::cout << "\nSams" << std::endl;
//            for (auto& sam : sam_list) {
//                std::cout << sam << std::endl;
//            }
//
//            std::cout << "\nProfiles" << std::endl;
//            for (auto& profile : profile_list) {
//                std::cout << profile << std::endl;
//            }


            is.close();
            return true;
        }

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
        // or qualities), noted in note_log for samples of this run. With --profile_only, from
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
                    if (!given && single && in_run && i < m_first_list.size()) {
                        auto const guess = GuessReadType(SampleReads(m_first_list[i]), MAX_SHORT_READ_LENGTH);
                        if (guess && guess->type != ReadType::Single) {
                            m_read_types[i] = guess->type;
                            if (note_log) {
                                note_log->emplace_back("Sample " + (i < m_sampleid_list.size() ? m_sampleid_list[i] : m_prefix_list[i]) + ": " + ReadTypeName(guess->type) +
                                                       " reads (" + guess->evidence + " in " + m_first_list[i] + "), aligned and profiled as such; "
                                                       "--read_type (or a map's READ_TYPE) sets the kind of reads");
                            }
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

        bool PrepareAndCheckValidity(bool force_read_check=false) {
            std::vector<std::string> error_log;
            std::vector<std::string> warning_log;
            std::vector<std::string> note_log = m_input_notes;
            error_log = m_input_errors;  // from reading the arguments: --profile_only patterns without a SAM

            ResolveDatabase(error_log);
            bool const add_model = !m_add_model.empty();
            bool const add_tables = !m_add_tables.empty();
            bool const write_neighbours = !m_write_species_neighbours.empty();
            bool const db_mode = m_compress_db || m_decompress_db || m_unpack_db || add_model || add_tables || write_neighbours;
            if (!m_build && !db_mode) ResolveReadTypes(warning_log, &note_log);
            if (int(m_compress_db) + int(m_decompress_db) + int(m_unpack_db) + int(add_model) + int(add_tables) + int(write_neighbours) +
                int(m_build) > 1) {
                error_log.emplace_back("--build, --compress_db, --decompress_db, --unpack_db, --add_model, --add_tables and "
                                       "--write_species_neighbours cannot be combined (--build compresses unless --no_compress)");
            }
            if (write_neighbours && m_bundle) {
                error_log.emplace_back("--write_species_neighbours needs a folder of separate files (reference.fna, reference.map, "
                                       "internal_taxonomy.dmp), not a single-file database: --unpack_db it first");
            }
            if (add_tables) {
                auto const allowed = AddableTables();
                std::set<std::string> names;
                for (auto const& file : AddTables()) {
                    std::string const name = std::filesystem::path(file).filename().string();
                    if (std::find(allowed.begin(), allowed.end(), name) == allowed.end()) {
                        std::string list;
                        for (auto const& table : allowed) list += (list.empty() ? "" : ", ") + table;
                        error_log.emplace_back("--add_tables: " + file + " is none of the tables it stores (" + list + ")");
                    } else if (!names.insert(name).second) {
                        error_log.emplace_back("--add_tables names " + name + " twice");
                    }
                    if (!std::filesystem::is_regular_file(file)) error_log.emplace_back("--add_tables: " + file + " is not a file");
                }
            }
            if (add_model) {
                // --add_model FILE[,FILE...] with --read_type TYPE[,TYPE...]: one read type per model, each once.
                size_t const models = static_cast<size_t>(std::count(m_add_model.begin(), m_add_model.end(), ',')) + 1;
                size_t const types = m_read_type.empty() ? 0 : static_cast<size_t>(std::count(m_read_type.begin(), m_read_type.end(), ',')) + 1;
                if (models > 1 && types != models) {
                    error_log.emplace_back("--add_model gives " + std::to_string(models) + " models: --read_type must give as many "
                                           "read types, comma-separated, one per model (e.g. --add_model pe.xml,se.xml "
                                           "--read_type pe,se), not '" + m_read_type + "'");
                } else if (types > models) {
                    error_log.emplace_back("--read_type gives " + std::to_string(types) + " read types for --add_model's " +
                                           std::to_string(models) + " model");
                }
                std::istringstream type_tokens(m_read_type);
                std::string token;
                std::set<std::string> seen;
                while (!m_read_type.empty() && std::getline(type_tokens, token, ',')) {
                    if (!ReadTypeFromToken(token)) {
                        error_log.emplace_back("--read_type must be one of " + ReadTypeTokens() + ", not '" + token + "'");
                    } else if (!seen.insert(token).second) {
                        error_log.emplace_back("--read_type names " + token + " twice for --add_model");
                    }
                }
                for (auto const& [file, _] : AddModels()) {
                    if (!std::filesystem::is_regular_file(file)) error_log.emplace_back("--add_model: " + file + " is not a file");
                }
            }
            if (m_unpack_db && !m_bundle && error_log.empty()) {
                error_log.emplace_back("--unpack_db needs a single-file database, but " + m_database_path + " holds separate files" +
                                       (m_location.unused_bundle.empty() ? "" : "; to unpack " + m_location.unused_bundle + ", give it as --db"));
            }
            if (m_bundle) {
                std::pair<std::string, std::string> const required[] = {
                        {PROTAL_INDEX_FILE, "index"}, {PROTAL_SEQUENCE_FILE, "reference"},
                        {PROTAL_SEQUENCE_MAP_FILE, "sequence map"}, {PROTAL_TAXONOMY_FILE, "taxonomy"}};
                for (auto const& [name, what] : required) {
                    if (!m_bundle->Find(name)) error_log.emplace_back("The " + what + " " + name + " is not in " + m_bundle->Path());
                }
                if (!m_preload_genomes && !db_mode) error_log.emplace_back(PreloadOffNeedsFilesMessage());
            } else if (!m_unpack_db && m_location.missing.empty() && (m_build || m_location.bundle.empty())) {  // a folder: not missing, nor a single file that failed to open
                size_t const before = error_log.size();
                if (!std::filesystem::exists(ResolvedSequenceFile())) {
                    error_log.emplace_back("Sequence file does not exist: " + GetSequenceFile() + " (nor " +
                                           GetSequenceFile() + zstd::kExtension + ")");
                } else if (!m_build && !m_preload_genomes && !db_mode && zstd::IsCompressed(ResolvedSequenceFile())) {
                    error_log.emplace_back("--preload_genomes_off needs an uncompressed reference, but " + ResolvedSequenceFile() +
                                           " is compressed. Drop --preload_genomes_off, or decompress it with: zstd -d " +
                                           ResolvedSequenceFile());
                }
                if (!std::filesystem::exists(GetSequenceMapFile())) {
                    error_log.emplace_back("Sequence map file does not exist: " + GetSequenceMapFile());
                }
                if (!std::filesystem::exists(GetInternalTaxonomyFile())) {
                    error_log.emplace_back("Taxonomy file does not exist: " + GetInternalTaxonomyFile());
                }
                // --write_species_neighbours reads the references and the taxonomy only: a converted release has no index yet.
                if (!m_build && !write_neighbours && !std::filesystem::exists(ResolvedIndexFile())) {
                    error_log.emplace_back("Index file does not exist: " + GetIndexFile() + " (nor " + GetIndexFile() +
                                           zstd::kExtension + ")");
                }
                std::string const packed = !m_location.bundle.empty() ? m_location.bundle : m_location.unused_bundle;
                if (m_build && error_log.size() > before && !packed.empty()) {
                    error_log.emplace_back(packed + " holds a built database; to rebuild it, first unpack its files with: " +
                                           ShellWord(ProgramName()) + " --unpack_db --db " + ShellWord(packed));
                }
            }
            if (((m_build && WriteBundle()) || (m_compress_db && m_bundle_db)) && m_compress_frame_mb == 0) {
                error_log.emplace_back("--compress_frame_mb 0 (one frame) needs --no_bundle: a single-file database is made of frames");
            }
            if ((m_build && m_compress) || m_compress_db) {
                if (m_compress_level < 1 || m_compress_level > ZSTD_maxCLevel()) {
                    error_log.emplace_back("--compress_level must be between 1 and " + std::to_string(ZSTD_maxCLevel()));
                }
                auto const window = ZSTD_cParam_getBounds(ZSTD_c_windowLog);
                if (m_compress_window_log != 0 &&
                    (m_compress_window_log < window.lowerBound || m_compress_window_log > window.upperBound)) {
                    error_log.emplace_back("--compress_window_log must be 0 or between " + std::to_string(window.lowerBound) +
                                           " and " + std::to_string(window.upperBound));
                }
                if (m_compress_frame_mb < 0 || m_compress_frame_mb > 4095) {
                    error_log.emplace_back("--compress_frame_mb must be between 0 and 4095");
                }
            }
            if (db_mode) {
                // Only the database files are needed.
            } else if (m_build) {
                // Either file may be zstd-compressed, or have a .zst sibling instead.
                if (m_sequence_file.empty() || !std::filesystem::exists(zstd::Resolve(m_sequence_file))) {
                    error_log.emplace_back("--reference does not exist: '" + m_sequence_file + "'");
                }
                if (!m_full_sequence_file.empty() && !std::filesystem::exists(zstd::Resolve(m_full_sequence_file))) {
                    error_log.emplace_back("--full_reference does not exist: " + m_full_sequence_file);
                }
                if (m_full_sequence_file.empty()) {
                    // Without it no unique_kmers.tsv would be written, and a query could detect nothing.
                    warning_log.emplace_back("no --full_reference: unique k-mers are checked against --reference only");
                    m_full_sequence_file = m_sequence_file;
                }
            } else if (m_bundle || (m_location.missing.empty() && m_location.bundle.empty())) {  // not missing, nor a single file that failed to open
                // Nothing is profiled with --no_profile or --write_species_neighbours: no unique k-mers or models needed.
                if (!m_no_profile && !write_neighbours && !UniqueKmersFileExists()) {
                    error_log.emplace_back("Unique k-mer file does not exist: " + UniqueKmersDbFile().Name() +
                                           " (without it every taxon fails the model; rebuild the database with --build)");
                }
                // The model of each kind of reads the samples have (of paired-end reads without samples).
                bool const any_sample = std::any_of(kReadTypes.begin(), kReadTypes.end(), [this](ReadTypeInfo const& t) { return AnySample(t.type); });
                for (auto const& info : kReadTypes) {
                    auto const type = info.type;
                    if (m_no_profile || write_neighbours || !(AnySample(type) || (type == ReadType::Paired && !any_sample)) ||
                        ModelDbFile(type).Exists()) continue;
                    if (!ModelNamed(type).empty()) {
                        error_log.emplace_back("Model file does not exist: " + ModelDbFile(type).Name());
                        continue;
                    }
                    std::string const db = m_bundle ? m_bundle->Path() : m_database_path;
                    error_log.emplace_back("The database has no model for --read_type " + info.token + " (" + info.name + " reads): " +
                                           ModelDbFile(type).Name() + " does not exist" + (type == ReadType::Paired ? " (nor model.xml)" : "") +
                                           ". Give one with " + info.model_option + ", or store one in the database with: " +
                                           ShellWord(ProgramName()) + " --add_model MODEL.xml --read_type " + info.token + " --db " + ShellWord(db) +
                                           (type == ReadType::Paired ? "" : " (the model of paired-end reads does not fit " + info.name + " reads)"));
                }
            }
            if (!(m_depth_identity_margin >= 0)) {
                error_log.emplace_back("--depth_identity_margin must be 0 or more");
            }
            if (!(m_msa_identity_margin >= 0)) {
                error_log.emplace_back("--msa_identity_margin must be 0 or more");
            }
            if (!m_gene_conservation.empty() && m_gene_conservation != "none" && m_gene_conservation != "db" &&
                !std::filesystem::is_regular_file(m_gene_conservation)) {
                error_log.emplace_back("--gene_conservation does not exist: " + m_gene_conservation + " (none, db or a file)");
            }
            if (!(m_knob >= 0 && m_knob <= 1)) {
                error_log.emplace_back("--knob must be between 0 and 1 (a probability)");
            }
            if (m_fdr && !(*m_fdr >= 0 && *m_fdr < 1)) {
                error_log.emplace_back("--fdr must be at least 0 and below 1 (an expected share of false calls; 0: off)");
            }
            if (m_fdr && *m_fdr > 0 && m_knob_given) {
                error_log.emplace_back("--fdr and --knob exclude each other (--fdr 0 turns the calibrated calls off)");
            }
            if (m_msa_min_depth == 0) {
                error_log.emplace_back("--msa_min_depth must be at least 1");
            }
            if (m_msa_knob && !(*m_msa_knob >= 0 && *m_msa_knob <= 1)) {
                error_log.emplace_back("--msa_knob must be between 0 and 1 (a probability)");
            }
            if (!m_profile_truth_list.empty()) {
                if (m_profile_truth_list.size() != m_profile_list.size()) {
                    error_log.emplace_back("Truth files (--profile_truth or a PROFILE_TRUTH column) must name one file per sample: " +
                                           std::to_string(m_profile_truth_list.size()) + " given for " +
                                           std::to_string(m_profile_list.size()) + " samples");
                }
                for (auto const& truth : m_profile_truth_list) {
                    if (!std::filesystem::exists(truth)) {
                        error_log.emplace_back("Truth file does not exist: '" + truth + "'");
                    }
                }
            }

            // Length of files
            bool valid_lengths1 =
                    m_first_list.size() == m_second_list.size() &&
                    m_first_list.size() == m_prefix_list.size();

            bool valid_lengths2 =
                    m_prefix_list.size() == m_sam_list.size();


            if (!valid_lengths1 && !m_first_list.empty()) {
                std::string error =
                        "You must provide equal amounts of items in options -1, -2 (unless all reads are single-end) and --prefix. "
                        "Provided are -1 (" +
                        std::to_string(m_first_list.size()) +
                        "), -2 (" +
                        std::to_string(m_second_list.size()) +
                        "), and --prefix (" +
                        std::to_string(m_prefix_list.size()) +
                        ")";
                error_log.emplace_back(error);
            }

            if (!valid_lengths2 && m_profile_only) {
                std::string error =
                        "You must provide equal amounts of items in options --prefix and --profile_only";
                error_log.emplace_back(error);
            }

            // In profile-only mode there are no read files, only sam files.
            if (!m_profile_only &&
                (m_first_list.size() != m_second_list.size() || m_first_list.size() != m_sam_list.size())) {
                std::cerr << "First:  " << m_first_list.size() << std::endl;
                std::cerr << "Second: " << m_second_list.size() << std::endl;
                std::cerr << "Sam:    " << m_sam_list.size() << std::endl;
                std::cerr << "lists different sizes" << std::endl;
                exit(9);
            }

            if (m_profile_only) {
                for (auto& sam : m_sam_list) {
                    if (!std::filesystem::exists(sam)) {
                        error_log.emplace_back("--profile_only sam file does not exist: " + sam);
                    }
                }
            }

            // --strain_spill: a folder protal can write its spill files into (made if missing).
            if (!m_strain_spill.empty() && !m_build && !db_mode) {
                if (m_no_strains) {
                    warning_log.emplace_back("--strain_spill is not used with --no_strains");
                } else {
                    std::error_code ec;
                    std::filesystem::create_directories(m_strain_spill, ec);
                    if (!std::filesystem::is_directory(m_strain_spill, ec) || ::access(m_strain_spill.c_str(), W_OK) != 0) {
                        error_log.emplace_back("--strain_spill " + m_strain_spill + " is not a folder protal can write to");
                    }
                }
            }

            // --profile_only does not overwrite an earlier run's results (its rerun into the same folder also added a
            // second set of rows to misc/<sample>_runtime.tsv): the profiles it would write, and the list of species of
            // the strain output folder. --force writes them again.
            if (m_profile_only && !m_force && !m_no_profile) {
                std::vector<std::string> existing;
                for (auto const& profile : m_profile_list) {
                    if (std::filesystem::exists(profile)) existing.push_back(profile);
                }
                if (!m_no_strains && std::filesystem::exists(GetStrainSpeciesListOutput())) {
                    existing.push_back(GetStrainSpeciesListOutput());
                }
                if (!existing.empty()) {
                    std::string listed;
                    for (size_t i = 0; i < existing.size() && i < 3; i++) listed += (i ? ", " : "") + existing[i];
                    if (existing.size() > 3) listed += ", ... (" + std::to_string(existing.size()) + " files)";
                    error_log.emplace_back("--profile_only would overwrite the results of an earlier run: " + listed +
                                           ". Give another output folder (-o), or --force to write them again");
                }
            }

            // Samples sharing an output file would overwrite, or corrupt, each other's output.
            auto no_shared_files = [&error_log](std::vector<std::string> const& paths, std::string const& what) {
                std::map<std::string, size_t> seen;
                for (size_t i = 0; i < paths.size(); i++) {
                    std::error_code ec;
                    auto canonical = std::filesystem::weakly_canonical(paths[i], ec);
                    auto key = ec ? std::filesystem::path(paths[i]).lexically_normal().string() : canonical.string();
                    auto [it, fresh] = seen.emplace(key, i);
                    if (!fresh) {
                        error_log.emplace_back("samples " + std::to_string(it->second + 1) + " and " + std::to_string(i + 1) +
                                               " would both use the " + what + " file " + paths[i] +
                                               "; give each sample its own prefix (--prefix, or PREFIX in a map)");
                    }
                }
            };
            no_shared_files(m_sam_list, "SAM");
            no_shared_files(m_profile_list, "profile");

            // Sample IDs name files (misc/<sample>_runtime.tsv): a map's are checked as it is read (LoadFromMap), these
            // are --prefix's and --profile_only's.
            if (!m_build && !db_mode) {
                for (size_t i = 0; i < m_sampleid_list.size(); i++) {
                    auto const problem = SampleIdProblem(m_sampleid_list[i]);
                    if (problem.empty()) continue;
                    error_log.emplace_back("sample ID '" + m_sampleid_list[i] + "' (sample " + std::to_string(i + 1) +
                                           ") cannot name a file safely on Linux and macOS: " + problem +
                                           "; rename it (--prefix, or #SAMPLEID in a map)");
                }
            }

            // Sample IDs name the rows of the strain MSAs: qcmsa matches rows to .meta.tsv by name,
            // and IQ-TREE cuts a name at its first space.
            if (!m_build && !m_no_profile && !m_no_strains) {
                std::map<std::string, size_t> seen;
                for (size_t i = 0; i < m_sampleid_list.size(); i++) {
                    auto const& id = m_sampleid_list[i];
                    if (std::any_of(id.begin(), id.end(), [](unsigned char c) { return std::isspace(c); })) {
                        error_log.emplace_back("sample ID '" + id + "' (sample " + std::to_string(i + 1) +
                                               ") contains whitespace, which the strain MSAs cannot hold; rename it "
                                               "(#SAMPLEID in a map) or pass --no_strains");
                    }
                    auto [it, fresh] = seen.emplace(id, i);
                    if (!fresh) {
                        error_log.emplace_back("samples " + std::to_string(it->second + 1) + " and " + std::to_string(i + 1) +
                                               " share the sample ID '" + id + "', so their strain MSA rows would be "
                                               "confused; give each its own (#SAMPLEID in a map) or pass --no_strains");
                    }
                }
            }

            // Check files
            for (auto i = 0; i < m_first_list.size(); i++) {

                auto first = m_first_list[i];
                auto second = m_second_list[i];
                auto sam = m_sam_list[i];

                auto first_exists = std::filesystem::exists(first);
                auto second_exists = second.empty() || std::filesystem::exists(second);  // none: single-end
                auto sam_exists = std::filesystem::exists(sam);

                if (!first_exists || !second_exists) {
                    if (!force_read_check && sam_exists) {
                        if (!first_exists) {
                            std::string warning = "-1 file does not exist: " + first + " ( But .sam file does )";
                            warning_log.emplace_back(warning);
                        }
                        if (!second_exists) {
                            std::string warning = "-2 file does not exist: " + second + " ( But .sam file does )";
                            warning_log.emplace_back(warning);
                        }
                    } else {
                        if (!first_exists) {
                            std::string error = "-1 file does not exist: " + first;
                            error_log.emplace_back(error);
                        }
                        if (!second_exists) {
                            std::string error = "-2 file does not exist: " + second;
                            error_log.emplace_back(error);
                        }
                    }
                }
                // A file that exists but cannot be read (permissions, a directory) would read as empty.
                if (force_read_check || !sam_exists) {
                    auto check = [&](char const* flag, std::string const& path) {
                        if (path.empty() || !std::filesystem::exists(path)) return;
                        std::string const problem = Unreadable(path);
                        if (!problem.empty()) error_log.emplace_back(std::string(flag) + " file cannot be read: " + path + " (" + problem + ")");
                    };
                    check("-1", first);
                    check("-2", second);
                }
            }

            // The kinds of reads given: only paired-end reads come in two files, and long reads given as
            // short single-end ones would be aligned as short reads.
            if (!m_build && !db_mode) {
                for (size_t i = 0; i < m_read_type_list.size(); i++) {
                    auto const& token = m_read_type_list[i];
                    if (token.empty()) continue;  // by the second file
                    std::string const sample = i < m_sampleid_list.size() ? m_sampleid_list[i] : std::to_string(i + 1);
                    bool const two_files = i < m_second_list.size() && !m_second_list[i].empty();
                    auto const type = ReadTypeFromToken(token);
                    if (!type) {
                        error_log.emplace_back("The read type of sample " + sample + " is '" + token + "' (--read_type or READ_TYPE): give one of " +
                                               ReadTypeTokens());
                    } else if (!m_profile_only && (*type == ReadType::Paired) != two_files) {
                        error_log.emplace_back("Sample " + sample + " has " + ReadTypeName(*type) + " reads (" + token + "), which come in " +
                                               (*type == ReadType::Paired ? "two files, but it has no second read file" :
                                                                            "one file, but it has a second read file: " + m_second_list[i]));
                    }
                }
                for (auto i : m_range) {
                    if (GetReadType(i) != ReadType::Single || i >= m_first_list.size()) continue;
                    size_t const longest = LongestRead(m_first_list[i], 100);
                    if (longest > MAX_SHORT_READ_LENGTH) {
                        error_log.emplace_back("Sample " + GetSampleId(i) + " has reads of up to " + std::to_string(longest) +
                                               " bp (first 100 of " + m_first_list[i] + "), too long for short reads: give "
                                               "--read_type pb or ont (or in the map's READ_TYPE column) for PacBio or ONT reads");
                    }
                }
            }

            for (auto& line : note_log) {
                std::cerr << "Note: " << line << std::endl;
            }
            for (auto& line : warning_log) {
                std::cerr << "Warning: " << line << std::endl;
            }

            if (!error_log.empty()) {
                for (auto& line : error_log) {
                    std::cerr << "Error: " << line << std::endl;
                }
                return false;
            }
            return true;
        }

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

        static Options OptionsFromArguments(int argc, char *argv[]) {
            auto cxx_options = CxxOptions();


            if (argc <= 1) {
                return { true, false, false , false };
            }

            cxx_options.parse_positional({ "positional" });

            // cxxopts throws on unknown options and on values it cannot convert.
            // Turn that into a readable message instead of an uncaught exception.
            cxxopts::ParseResult result;
            try {
                result = cxx_options.parse(argc, argv);
            } catch (const cxxopts::OptionException& e) {
                std::cerr << "Error parsing options: " << e.what() << std::endl;
                std::cerr << "Run 'protal --help' for the available options ('protal --full_help' "
                             "also lists the developer options)." << std::endl;
                exit(2);
            }

            bool show_help = result.count("help");
            bool show_help_dev = result.count("full_help");
            bool show_map_help = result.count("map_help");
            bool show_version = result.count("version");

            if (show_help || show_version || show_map_help || show_help_dev) {
                return { show_help, show_map_help, show_version, show_help_dev };
            }

            bool no_strains = result.count("no_strains");
            bool no_mate_guidance = result.count("no_mate_guidance");
            bool no_gene_neighbours = result.count("no_gene_neighbours");
            bool keep_foreign_genes = result.count("keep_foreign_genes");
            bool keep_suspect_copies = result.count("keep_suspect_copies");
            bool no_phasing = result.count("no_phasing");
            bool build = result.count("build");
            bool preload_genomes_off = result.count("preload_genomes_off");
            bool benchmark_alignment = result.count("benchmark_alignment");
            bool fastalign = false;//result.count("fastalign");
            bool no_profile = result.count("no_profile");
            bool verbose = result.count("verbose");

            size_t threads = result["threads"].as<size_t>();
            size_t align_top = result["align_top"].as<size_t>();
            size_t x_drop = result["x_drop"].as<size_t>();
            size_t max_key_ubiquity = result["max_key_ubiquity"].as<size_t>();
            size_t min_successful_lookups = result["min_successful_lookups"].as<size_t>();
            size_t max_seed_size = result["max_seed_size"].as<size_t>();
            double max_score_ani = result["max_score_ani"].as<double>();
            size_t msa_min_hcov = result["msa_min_hcov"].as<size_t>();
            size_t msa_min_depth = result["msa_min_depth"].as<size_t>();
            size_t max_out = result["max_out"].as<size_t>();
            auto msa_species_arg = result["msa_species"].as<std::string>();
            std::vector<std::string> msa_species;
            if (!msa_species_arg.empty()) {
                auto tokens = Utils::split(msa_species_arg, ",");
                for (auto& token : tokens) {
                    size_t start = 0;
                    while (start < token.size() && std::isspace(static_cast<unsigned char>(token[start]))) start++;
                    size_t end = token.size();
                    while (end > start && std::isspace(static_cast<unsigned char>(token[end - 1]))) end--;
                    auto trimmed = token.substr(start, end - start);
                    if (!trimmed.empty()) msa_species.emplace_back(std::move(trimmed));
                }
            }

            size_t snp_min_cov        = result["snp_min_cov"].as<size_t>();
            size_t snp_min_phred_sum  = result["snp_min_phred_sum"].as<size_t>();
            double snp_min_af         = result["snp_min_af"].as<double>();
            size_t snp_min_mean_qual  = result["snp_min_mean_qual"].as<size_t>();
            bool   snp_require_strand = !result.count("snp_no_strand");
            size_t snp_max_alleles    = result["snp_max_alleles"].as<size_t>();
            bool   run_qcmsa          = !result.count("no_qcmsa");  // on by default
            std::string qcmsa_script  = result["qcmsa_script"].as<std::string>();
            std::string qcmsa_args    = result["qcmsa_args"].as<std::string>();


            auto reference = result["reference"].as<std::string>();
            auto full_reference = result["full_reference"].as<std::string>();
            auto build_gene_subset = result["build_gene_subset"].as<std::string>();

            // Arguments that follow no option. With --profile_only (and no map) they are SAM files: an unquoted pattern
            // (--profile_only runs/*/alignments/*.sam.zst) reaches protal as one argument per file, of which the option
            // takes the first, and the others were once ignored without a word. With --build, the reference. Elsewhere
            // nothing takes them, so they stop protal (e.g. an unquoted -1 *_R1.fq would align the first file only).
            std::vector<std::string> positional = result.count("positional") ? result["positional"].as<std::vector<std::string>>()
                                                                              : std::vector<std::string>{};
            std::vector<std::string> extra_sams;
            if (!positional.empty()) {
                if (result.count("profile_only") && !result.count("map")) {
                    extra_sams = positional;
                } else if (build && !result.count("reference") && positional.size() == 1) {
                    reference = positional.front();
                } else {
                    std::string listed;
                    for (auto const& arg : positional) listed += (listed.empty() ? "" : " ") + arg;
                    std::cerr << "Error: unexpected argument(s): " << listed << "\nAn option takes one value: give several files "
                                 "comma-separated (-1 a_1.fq,b_1.fq). Only --profile_only takes several SAM files as arguments, or "
                                 "a pattern such as 'runs/*/alignments/*.sam.zst'." << std::endl;
                    exit(2);
                }
            }

            auto map_file = result.count("map") ? result["map"].as<std::string>() : "";
            auto first = result.count("first") ? result["first"].as<std::string>() : "";
            auto second = result.count("second") ? result["second"].as<std::string>() : "";
            auto sam_in = result.count("profile_only") ? result["profile_only"].as<std::string>() : "";
            auto universal_prefix = result.count("prefix") ? result["prefix"].as<std::string>() : "";
            auto output_dir = result.count("outdir") ? result["outdir"].as<std::string>() : "";
            std::string const sam_format = result["sam_format"].as<std::string>();
            std::string const sam_ending = SamEnding(sam_format);
            if (sam_ending.empty()) {
                std::cerr << "Error: --sam_format must be zst, gz or sam, not '" << sam_format << "'" << std::endl;
                exit(2);
            }

            auto range_arg = result.count("map_range") ? result["map_range"].as<std::string>() : "";
            std::vector<size_t> range;


            std::vector<std::string> first_list;
            std::vector<std::string> second_list;
            std::vector<std::string> prefix_list;
            std::vector<std::string> sam_list;
            std::vector<std::string> profile_list;
            std::vector<std::string> samplenames_list;
            std::vector<std::string> profile_truth_list;
            std::vector<std::string> read_type_list;
            std::vector<std::string> unmapped_reads_list;
            std::vector<std::string> input_errors;  // reported by PrepareAndCheckValidity
            std::vector<std::string> input_notes;

            std::string strain_output_dir = "";
            std::string misc_output_dir = "";

            std::vector<uint8_t> build_gene_mask;
            if (!build_gene_subset.empty()) {
                if (!std::filesystem::exists(build_gene_subset)) {
                    std::cerr << "Gene subset file " << build_gene_subset << " does not exist." << std::endl;
                    exit(9);
                }
                std::ifstream subset_stream(build_gene_subset);
                if (!subset_stream.good()) {
                    std::cerr << "Failed to open gene subset file " << build_gene_subset << std::endl;
                    exit(9);
                }

                build_gene_mask.assign(1, 0);
                auto trim = [](std::string &s) {
                    size_t start = 0;
                    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) start++;
                    size_t end = s.size();
                    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) end--;
                    s = s.substr(start, end - start);
                };
                std::string line;
                while (std::getline(subset_stream, line)) {
                    trim(line);
                    if (line.empty() || line[0] == '#') continue;
                    if (!std::all_of(line.begin(), line.end(), [](unsigned char c){ return std::isdigit(c); })) {
                        std::cerr << "Invalid gene id in subset file: " << line << std::endl;
                        exit(9);
                    }
                    size_t gene_id = std::stoul(line);
                    if (gene_id < 1) {
                        std::cerr << "Gene id out of range (>=1): " << gene_id << std::endl;
                        exit(9);
                    }
                    if (gene_id >= build_gene_mask.size()) {
                        build_gene_mask.resize(gene_id + 1, 0);
                    }
                    build_gene_mask[gene_id] = 1;
                }
            }
            
            if (!map_file.empty()) {
                if (!std::filesystem::exists(map_file)) {
                    std::cerr << "Map file " << map_file << " does not exist." << std::endl;
                    exit(9);
                }
                if (!LoadFromMap(map_file, output_dir, strain_output_dir, misc_output_dir, prefix_list, first_list, second_list, sam_list, profile_list, samplenames_list, profile_truth_list, read_type_list, unmapped_reads_list, sam_ending)) {
                    std::cerr << "Failed to read map file " << map_file << " (see --map_help)." << std::endl;
                    exit(9);
                }

                if (range_arg != "") {
                    range = ProcessRange(range_arg, first_list.size());
                }
            } else {
                LineSplitter::Split(first, ",", first_list);
                LineSplitter::Split(second, ",", second_list);
                // Without -2, the reads are single-end: no sample has a second file.
                if (second.empty()) second_list.assign(first_list.size(), "");
                LineSplitter::Split(universal_prefix, ",", prefix_list);
                LineSplitter::Split(sam_in, ",", sam_list);
                // --profile_only's items and the arguments an unquoted pattern became, with wildcards expanded.
                if (result.count("profile_only")) {
                    sam_list.insert(sam_list.end(), extra_sams.begin(), extra_sams.end());
                    sam_list = ExpandSamFiles(sam_list, input_errors, input_notes);
                }

                if (strain_output_dir.empty()) {
                    strain_output_dir = std::filesystem::path(output_dir) / std::filesystem::path(MAP_VAR_DEFAULT_STRAIN_OUTPUT_DIR);
                }
                if (misc_output_dir.empty()) {
                    misc_output_dir = std::filesystem::path(output_dir) / std::filesystem::path(MAP_VAR_DEFAULT_MISC_OUTPUT_DIR);
                }
            }


            // The kind of reads of each sample: a map's READ_TYPE column, else --read_type (checked in
            // PrepareAndCheckValidity); with --profile_only, of each SAM.
            if (read_type_list.empty()) {
                read_type_list.assign(std::max(first_list.size(), sam_list.size()), result["read_type"].as<std::string>());
            }
            for (auto& type : read_type_list) {
                if (type == MAP_NO_SECOND_READ) type.clear();  // READ_TYPE '-': by the second file
            }
            for (auto& type : read_type_list) {
                std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c) { return std::tolower(c); });
            }

            auto profile_truth = result.count("profile_truth") ? result["profile_truth"].as<std::string>() : "";
            // One truth file per sample, comma-separated; overrides a map's PROFILE_TRUTH column.
            if (!profile_truth.empty()) {
                profile_truth_list.clear();
                LineSplitter::Split(profile_truth, ",", profile_truth_list);
            }
            bool profile_only = result.count("profile_only");
            bool mapq_debug_output = result.count("mapq_debug_output");
            bool whole_read_alignment = result.count("whole_read_alignment");
            bool no_alignment_screen = result.count("no_alignment_screen");
            size_t long_read_budget = result["long_read_budget"].as<size_t>();
            bool sequential_load = result.count("sequential_load");
            bool taxon_statistics = result.count("taxon_statistics");
            bool full_sam_header = result.count("full_sam_header");
            bool write_unmapped_reads = result.count("write_unmapped_reads");
            bool serial_index_passes = result.count("serial_index_passes");
            bool profile_ahead = result.count("profile_ahead");
            bool no_unknown_share = result.count("no_unknown_share");
            size_t index_batch_kb = std::max<size_t>(result["index_batch_kb"].as<size_t>(), 1);
            bool force = result.count("force");


            // Catch mismatched read lists here with a clear message, instead of running into them
            // further down. A sample has one read file (single-end reads; its second file is empty)
            // or two (paired-end reads).
            // These work on the database only.
            bool const compress_db = result.count("compress_db") > 0;
            bool const decompress_db = result.count("decompress_db") > 0;
            bool const unpack_db = result.count("unpack_db") > 0;
            bool const add_model = !result["add_model"].as<std::string>().empty();
            bool const add_tables = !result["add_tables"].as<std::string>().empty();
            bool const write_neighbours = !result["write_species_neighbours"].as<std::string>().empty();
            if (!build && !profile_only && !compress_db && !decompress_db && !unpack_db && !add_model && !add_tables && !write_neighbours) {
                if (first_list.empty()) {
                    std::cerr << "No input reads given. Provide reads via -1/--first (and -2/--second for paired-end "
                                 "reads), or a map file via --map (see --map_help)." << std::endl;
                    exit(31);
                }
                if (first_list.size() != second_list.size()) {
                    std::cerr << "-1/--first and -2/--second must name the same number of files (paired-end reads), "
                                 "or -2/--second none (single-end reads)." << std::endl;
                    std::cerr << "  -1/--first:  " << first_list.size() << " file(s)" << std::endl;
                    std::cerr << "  -2/--second: " << second_list.size() << " file(s)" << std::endl;
                    exit(31);
                }
            }

            if (profile_only) {
                if (!first_list.empty()) {
                    std::cerr << "Warning: profile only selected but first list is not empty. List is cleared" << std::endl;
                }
                if (std::any_of(second_list.begin(), second_list.end(), [](std::string const& s) { return !s.empty(); })) {
                    std::cerr << "Warning: profile only selected but second list is not empty. List is cleared" << std::endl;
                }
                first_list.clear();
                second_list.clear();

                if (!prefix_list.empty() && prefix_list.size() != sam_list.size()) {
                    std::cerr << "Warning: If profile only is selected, the user can either specify output prefixes for all sam files or for non. Sizes are not equal." << std::endl;
                    exit(34);
                }

                if (sam_list.empty() && input_errors.empty() && map_file.empty()) {
                    input_errors.emplace_back("--profile_only names no SAM file");
                }

                if (prefix_list.empty()) {
                    prefix_list.resize(sam_list.size());
                    for (auto i = 0; i < sam_list.size(); i++) {
                        if (!IsSamName(sam_list[i])) {
                            std::cerr << sam_list[i] << " does not end with .sam, .sam.gz or .sam.zst" << std::endl;
                            exit(35);
                        }
                    }
                    // The samples' names: the SAMs' file names, or their folders where file names repeat. A map's
                    // samples keep its #SAMPLEID.
                    auto const names = SamSampleNames(sam_list);
                    std::string by_folder;  // the samples named by their folders, if any
                    for (auto i = 0; i < sam_list.size(); i++) {
                        std::string const uncompressed = UncompressedSamName(sam_list[i]);  // without .gz or .zst
                        std::string const stem = uncompressed.substr(0, uncompressed.size() - 4);
                        // The outputs go to -o if it is given (the name is joined with it below), else next to the SAM.
                        prefix_list[i] = output_dir.empty() ? stem : names[i];
                        if (names[i] != std::filesystem::path(stem).filename().string() && i < 3) {
                            by_folder += (by_folder.empty() ? "" : ", ") + names[i] + " (" + sam_list[i] + ")";
                        }
                    }
                    if (!by_folder.empty()) {
                        input_notes.emplace_back("SAM file names repeat, so the samples are named by the folders in which their "
                                                 "paths differ: " + by_folder + (sam_list.size() > 3 ? ", ..." : "") +
                                                 "; --prefix (one per SAM) names them otherwise");
                    }
                    samplenames_list = names;
                }
            } else {
                if (prefix_list.empty()) {
                    for (int i = 0; i < first_list.size(); i++) {
                        auto first = std::filesystem::path(first_list[i]);
                        if (second_list[i].empty()) {
                            prefix_list.emplace_back(Utils::ReadFileStem(first.filename().string()));
                            continue;
                        }
                        auto second = std::filesystem::path(second_list[i]);

                        if (first.parent_path() != second.parent_path()) {
                            std::cerr << "Cannot infer output prefix from reads files as they are in different folders:" << std::endl;
                            std::cerr << "  " << first << std::endl;
                            std::cerr << "  " << second << std::endl;
                            std::cerr << "Please provide output prefixes via --prefix" << std::endl;
                            exit(32);
                        }
                        prefix_list.emplace_back(Utils::LongestCommonPrefixTrimmed(first.filename(), second.filename()));
                    }
                }

                if (sam_list.empty()) {
                    for (auto i = 0; i < prefix_list.size(); i++) {
                        std::filesystem::path sam { prefix_list[i] + sam_ending };

                        // If absolute, just push it directly
                        if (sam.is_absolute()) {
                            sam_list.emplace_back(sam);
                        } else {
                            // Join output_dir and the relative filename safely
                            sam_list.emplace_back(std::filesystem::path(output_dir) / sam);
                        }
                    }
                }
            }

            // Generate profile names if not provided through prefix or read-infered prefix
            if (!no_profile && profile_list.empty()) {
                for (auto i = 0; i < prefix_list.size(); i++) {
                    std::filesystem::path profile_file { prefix_list[i] + ".profile"};

                    // If absolute, just push it directly
                    if (profile_file.is_absolute()) {
                        profile_list.emplace_back(profile_file);
                    } else {
                        // Join output_dir and the relative filename safely
                        profile_list.emplace_back(std::filesystem::path(output_dir) / profile_file);
                    }
                }
            }

            auto profile_dir_override = result["profile_dir"].as<std::string>();
            if (!profile_dir_override.empty()) {
                for (auto& profile : profile_list) {
                    profile = (std::filesystem::path(profile_dir_override) / std::filesystem::path(profile).filename()).string();
                }
            }

            auto benchmark_alignment_output_file = result.count("benchmark_alignment_output") ? result["benchmark_alignment_output"].as<std::string>() : "";

            std::string db_path;
            if (result.count("db")) {
                db_path = result["db"].as<std::string>();
            } else {
                auto db_path_env = std::getenv(PROTAL_DB_ENV_VARIABLE.c_str());
                std::cout << "Get DB from environment variable $" << PROTAL_DB_ENV_VARIABLE << std::endl;
                if (db_path_env) db_path = db_path_env;  // else PrepareAndCheckValidity reports that none is given
            }

            if (range.empty()) {
                range.resize(prefix_list.size());
                std::iota(range.begin(), range.end(), 0);
            }

            OptionsData d;
            d.build                    = build;
            d.no_profile               = no_profile;
            d.profile_only             = profile_only;
            d.no_strains               = no_strains;
            d.no_mate_guidance         = no_mate_guidance;
            d.no_gene_neighbours       = no_gene_neighbours;
            d.keep_foreign_genes       = keep_foreign_genes;
            d.no_phasing               = no_phasing;
            d.strain_spill             = result["strain_spill"].as<std::string>();
            d.preload_genomes          = !preload_genomes_off;
            d.benchmark_alignment      = benchmark_alignment;
            d.benchmark_alignment_output = benchmark_alignment_output_file;
            d.show_help                = show_help;
            d.show_help_dev            = show_help_dev;
            d.show_map_help            = show_map_help;
            d.show_version             = show_version;
            d.mapq_debug_out           = mapq_debug_output;
            d.whole_read_alignment     = whole_read_alignment;
            d.no_alignment_screen      = no_alignment_screen;
            d.long_read_budget         = long_read_budget;
            d.sequential_load          = sequential_load;
            d.taxon_statistics         = taxon_statistics;
            d.full_sam_header          = full_sam_header;
            d.write_unmapped_reads     = write_unmapped_reads;
            d.serial_index_passes      = serial_index_passes;
            d.profile_ahead  = profile_ahead;
            d.no_unknown_share         = no_unknown_share;
            d.index_batch_kb           = index_batch_kb;
            d.first_list               = std::move(first_list);
            d.second_list              = std::move(second_list);
            d.samplename_list          = std::move(samplenames_list);
            d.database_path            = db_path;
            d.prefix_list              = std::move(prefix_list);
            d.sequence_file            = reference;
            d.full_sequence_file       = full_reference;
            d.map_file                 = map_file;
            d.strain_output_dir        = strain_output_dir;
            d.misc_output_dir          = misc_output_dir;
            d.output_dir               = output_dir;
            d.threads                  = threads;
            d.align_top                = align_top;
            d.adaptive_candidates      = result["adaptive_candidates"].as<size_t>();
            d.no_allele_scores         = result["no_allele_scores"].as<bool>();
            d.max_out                  = max_out;
            d.max_score_ani            = max_score_ani;
            d.max_score_ani_given      = result.count("max_score_ani") > 0;
            d.msa_min_hcov             = msa_min_hcov;
            d.msa_min_depth            = msa_min_depth;
            d.msa_species              = std::move(msa_species);
            d.snp_min_phred_sum        = snp_min_phred_sum;
            d.snp_min_cov              = snp_min_cov;
            d.snp_min_af               = snp_min_af;
            d.snp_min_af_given         = result.count("snp_min_af") > 0;
            d.snp_min_mean_qual        = snp_min_mean_qual;
            d.snp_require_strand       = snp_require_strand;
            d.run_qcmsa                = run_qcmsa;
            d.qcmsa_script             = std::move(qcmsa_script);
            d.qcmsa_args               = std::move(qcmsa_args);
            d.snp_max_alleles          = snp_max_alleles;
            d.x_drop                   = x_drop;
            d.max_key_ubiquity         = max_key_ubiquity;
            d.min_successful_lookups   = min_successful_lookups;
            d.max_seed_size            = max_seed_size;
            d.fastalign                = fastalign;
            d.sam_list                 = std::move(sam_list);
            d.profile_list             = std::move(profile_list);
            d.profile_truth_list       = std::move(profile_truth_list);
            d.read_type_list           = std::move(read_type_list);
            d.unmapped_reads_list      = std::move(unmapped_reads_list);
            d.input_errors             = std::move(input_errors);
            d.input_notes              = std::move(input_notes);
            d.profile_truth            = profile_truth;
            d.force                    = force;
            d.verbose                  = verbose;
            d.range                    = std::move(range);
            d.build_gene_mask          = std::move(build_gene_mask);
            d.compress                 = !result.count("no_compress");
            d.bundle                   = !result.count("no_bundle");
            d.compress_db              = compress_db;
            d.decompress_db            = decompress_db;
            d.unpack_db                = unpack_db;
            d.unpack_dir               = result["unpack_dir"].as<std::string>();
            d.command_line.assign(argv, argv + argc);
            d.compress_level           = result["compress_level"].as<int>();
            d.compress_window_log      = result["compress_window_log"].as<int>();
            d.compress_frame_mb        = result["compress_frame_mb"].as<int>();
            d.knob                     = result["knob"].as<double>();
            d.knob_given               = result.count("knob") > 0;
            if (result.count("fdr")) d.fdr = result["fdr"].as<double>();
            d.singleton_congener       = result["singleton_congener"].as<size_t>();
            if (result.count("msa_knob")) d.msa_knob = result["msa_knob"].as<double>();
            d.depth_identity_margin    = result["depth_identity_margin"].as<double>();
            d.keep_suspect_copies      = keep_suspect_copies;
            d.suspect_copy_distance    = result["suspect_copy_distance"].as<double>();
            d.strain_alleles           = result["strain_alleles"].as<size_t>();
            d.allele_genome_share      = result["allele_genome_share"].as<double>();
            d.msa_identity_margin      = result["msa_identity_margin"].as<double>();
            d.gene_conservation        = result["gene_conservation"].as<std::string>();
            d.model                    = result["model"].as<std::string>();
            d.model_se                 = result["model_se"].as<std::string>();
            d.model_pb                 = result["model_pb"].as<std::string>();
            d.model_ont                = result["model_ont"].as<std::string>();
            d.read_type                = result["read_type"].as<std::string>();
            d.add_model                = result["add_model"].as<std::string>();
            d.add_tables               = result["add_tables"].as<std::string>();
            d.write_species_neighbours = result["write_species_neighbours"].as<std::string>();

            auto options = Options(std::move(d));

            if (!options.PrepareAndCheckValidity()) {
                std::cerr << "Exit Program" << std::endl;
                exit(30);
            }

            return options;
        }
    };


}
