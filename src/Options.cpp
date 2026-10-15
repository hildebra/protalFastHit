// The options' definitions on the command line (CxxOptions), the map file, the parsing and the checks before a run
// (Options.h declares them). Here rather than in Options.h, which nearly every unit of protal includes
// (docs/claude/2026-10-14-compile-time).

#include <cxxopts.hpp>
#include "Options.h"
#include "FeatureVersions.h"

namespace protal {
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
                ("one_sample", "The files of -1/--first (and -2/--second) are one sample's, its lanes or runs, read one after the other as one sample instead of a sample each; -1's and -2's pair up by position, and a lane whose two files hold different numbers of reads fails the sample. --prefix names it, else the longest beginning the files' names share (without a lane number such as _L001). In a map, a FIRST and SECOND cell lists a sample's files in the same way, comma-separated.")
                ("read_type", "The reads of -1/--first: pe (paired-end reads, with -2/--second), se (single-end reads), pb (PacBio long reads, e.g. HiFi), ont (Oxford Nanopore long reads). Without it, pe with -2/--second; with one read file, its first reads decide, and protal says what it found: se for short reads, for long reads pb or ont by their names (PacBio movie/ZMW, the UUIDs of MinKNOW and dorado), else by their median read quality (pb from Q25 on, ont below; pb without qualities). A map gives it per sample in a READ_TYPE column. With --profile_only, the SAM's header names its reads, and --read_type replaces that (e.g. to profile a SAM with another read type's model); with --add_model, it names the model's reads (pe if not given).", cxxopts::value<std::string>()->default_value(""))
                ("prefix", "Comma separated list of output prefixes (optional). If not specified, output file prefixes are generated from the input file names: the longest common prefix of the two files of paired-end reads (which must then be in the same folder), the file name without its FASTQ/FASTA and compression extensions for single-end reads.", cxxopts::value<std::string>()->default_value(""))
                ("o,outdir", "Overwrites #OUTPUT_DIR in map and needs to be defined if #OUTPUT_DIR is not defined in the map. If not otherwise specified in the map file, sam files, profiles, msas, and other miscellaneous files will be stored in the subfolders to this directory 'alignments', 'profiles', 'strains', and 'misc'. With -1/-2: the folder of the SAM files and profiles, with 'strains' and 'misc' in it (default: the current folder).", cxxopts::value<std::string>())
                
                ("map", "For larger datasets you can define parameters -1, -2, --prefix and -o in a tsv-file.", cxxopts::value<std::string>()->default_value(""))
                ("map_range", "If you specified a map file with --map you can also pass a range to protal to run protal only on a subset. The first entry is 1, the end is inclusive. e.g.: 1-10. If the end open or larger than the number of entries in the map file, the last entry in the map file is selected as end.", cxxopts::value<std::string>()->default_value(""))
                ("profile_only", "Comma separated list of existing sam files to profile without re-running the alignment, e.g. to build strain MSAs over the samples of several runs. An item with a wildcard, such as 'runs/*/alignments/*.sam.zst' (quoted: protal expands it; unquoted, the shell does, which works as well), stands for the .sam, .sam.gz and .sam.zst files it matches, sorted. Read files given via -1/-2 are ignored. Output prefixes are either given via --prefix (one per sam file) or derived from the sam file names (where two files share a name, as aln.sam.zst in a folder per sample, from the folders in which their paths differ); the outputs then go to -o if it is given, else the profiles next to each sam file and 'strains' and 'misc' into the current folder.", cxxopts::value<std::string>()->default_value(""))
                ("sam_format", "Format of the SAM files that protal names (with -1/-2, or a map without a SAM column): zst (zstd-compressed, <prefix>.sam.zst), gz (gzip, .sam.gz) or sam (plain). A map's SAM names keep their own ending.", cxxopts::value<std::string>()->default_value("zst"))
                ("scratch_dir", "A folder for protal's intermediate files, best on a disk of the machine itself (e.g. $TMPDIR on a cluster node) when the outputs go to a network file system. protal makes a folder of its own in it (protal_<host>_<pid>_<time>) and removes it at the end of the run, also when the run fails (not when it is killed). In it: what each sample keeps for the strain MSAs, written once its profile is written and read back one species at a time, so that the memory of many samples (e.g. --profile_only over several runs) does not grow with them (a sample's stays in memory if its file cannot be written; replaces --strain_spill). With --no_strains, the SAM files that protal names (with -1/-2, or a map without a SAM column), removed once their sample is profiled: nothing after the run reads them; without --no_strains they stay in the output folder, for strain MSAs over several runs, and a SAM an earlier run left in the output folder is profiled as before. And the SAM records of --sam_format gz and sam before their header is written, and qcMSA's working files.", cxxopts::value<std::string>()->default_value(""));

        // Alignment / algorithm options
        options.add_options("Alignment")
                ("c,align_top", "After seeding, anchor are sorted by quality passed to alignment. <take_top> specifies how many anchors should be aligned starting with the most promising anchor. auto (default): the database's (build_settings.tsv), else " + std::to_string(DEFAULT_ALIGN_TOP) + "; a run told otherwise than its database is warned, as for every setting --build records (docs/running.md).", cxxopts::value<std::string>()->default_value("auto"))
                ("adaptive_candidates", "A short read whose best alignment is divergent (identity below 0.99) and whose seeds fit more taxa equally well than --align_top took (ZN) tries up to this many more of those anchors, of taxa of the best alignment's genus: a strain's own species that its seeds ranked below its congeners'. 0: never. Needs the taxonomy (it is loaded for profiling). auto (default): the database's, else " + std::to_string(DEFAULT_ADAPTIVE_CANDIDATES) + ".", cxxopts::value<std::string>()->default_value("auto"))
                ("allele_scores", "Short reads: score candidates with the database's strain alleles (the default; see --no_allele_scores), also on a database built with --no_allele_scores (build_settings.tsv), whose runs otherwise follow it.")
                ("weighted_site_shift", "Short reads, with column weights in the database (column_weights.tsv): the site shift of a read protal is unsure about (--no_allele_scores) also takes half a difference off for each of the candidate's substitutions at a column the family's genera change freely (below 3 nats of conservation among the genus references; StrainAlleles.h ShiftOf), so that the columns that reliably differ settle the read. Off by default again since 2026-10-12: at r226 (v24b) it lowered the short-read species' depth by 1-2% with no F1 gain measurable, and on a small real world, with the same samples, it lowered 19-22% of the called species' depths (docs/claude/2026-10-12-r226-v24b). On by default 2026-10-11 to 2026-10-12: a database built then (r226 v24b, build_metadata.tsv: weighted_site_shift on) trained its models with it, so profile it with this option for the calls of its training. ~1% of the alignment's instructions. Neither this nor --no_weighted_site_shift given (auto): as the database's build_settings.tsv records (since 2026-10-15), else off.")
                ("no_weighted_site_shift", "Short reads: the site shift without the column weights' discount, the default (again since 2026-10-12; kept for the command lines of 2026-10-11). Wins over --weighted_site_shift.")
                ("no_site_weights", "Leave the database's column weights (column_weights.tsv, --build --column_weights) unread: the 'weights' features are unknown (-1), the ancestry sites of a species with fewer than three congeners are the nearest congener's differences, as before 2026-10-10, and --weighted_site_shift has nothing to weigh by. For measuring what the weights are worth, not for a production run: a database's models are trained with them.")
                ("no_allele_scores", "Short reads: score each candidate alignment against the reference alone. By default, with a database that has strain alleles (strain_alleles.tsv, --build --strain_alleles), a candidate's score counts the read's differences from the best of its species' known alleles where those explain them (StrainAlleles.h), so a read of a known strain scores on its species as on the strain's own gene (and with --weighted_site_shift discounts its differences at the columns the family's genera change freely); the records keep their alignment to the reference. A database's models are trained with its alleles' scores: use this to measure them, not on a production run. Neither this nor --allele_scores given (auto): as the database's build_settings.tsv records, else on.")
                ("m,max_out", "Maximum alignments that should be outputted", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_OUT)))
                ("no_mate_guidance", "Paired-end reads: do not let a mate that is sure of its alignment guide the other one when they did not align together (to the guiding mate's taxon from the other's own anchor, or on its gene where the fragment can reach, partly if it runs past the gene's end). Without it (auto): as the database's build_settings.tsv records, else on.")
                ("no_gene_neighbours", "Do not use the database's gene neighbours (gene_neighbours.tsv: how often each marker gene end faces which other in a clade's genomes): no looking for a mate past the end of its guiding mate's gene, no pairs of mates on two neighbouring genes, no looking on a long read for the genes next to its genes, adjacent_expected_share and adjacent_unlikely_share 0 and adjacent_support 0.5.")
                ("keep_suspect_copies", "Count reads on the database's suspect gene copies (suspect_copies.tsv: a species' copy of a marker gene that is near-identical to the copy of a species of another genus, family, order, class, phylum or domain while its own congeners' copies are farther, which --build finds; a contaminating contig or a transferred gene) as evidence of the species. By default their records are left out as if the reads had not aligned: every present organism with such a gene puts a perfect read on the copy, a fifth of the false species calls at GTDB r226 (docs/claude/2026-10-03-false-positive-anatomy). Without it (auto): as the database's build_settings.tsv records, else left out; --keep_suspect_copies=off leaves them out on any database.", cxxopts::value<std::string>()->default_value("auto")->implicit_value("on"))
                ("keep_foreign_genes", "Keep foreign genes in their taxon's depth and strain MSAs. A gene is foreign when the genes next to it on its reads (a pair's mates on two genes, a long read's consecutive genes) are mostly unlikely neighbours in its taxon's clade by the database's gene neighbours (4 or more such links judged, more than half of them unlikely): its reads come from another genome. By default (auto) foreign genes are left out of both, unless the database's build_settings.tsv records them kept; --keep_foreign_genes=off leaves them out on any database. .profile.genes.log lists them either way.", cxxopts::value<std::string>()->default_value("auto")->implicit_value("on"))
                ("u,max_key_ubiquity", "Max key ubiquity. Best matching Flexkey count for seed must be lower or equal", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_KEY_UBIQUITY)))
                ("s,max_seed_size", "Max seed size after which seeding is stopped.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_SEED_SIZE)))
                ("w,min_successful_lookups", "If the number of seeds is >=max_seed_size and the number of successful core-mer lookups is >= min_successful_lookups, stop looking for further seeds.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MIN_SUCCESSFUL_LOOKUPS)))
                ("a,max_score_ani", "A max score makes an alignment stop if the alignment diverges too much. This parameter estimates the score for a given ani and is a tradeoff between speed/accuracy. Given, it applies to all read types; auto (default): each read type's floor the database records (build_settings.tsv), else 0.9 and for ONT reads 0.85 (their indels count twice).", cxxopts::value<std::string>()->default_value("auto"))
                ("x,x_drop", "X-drop of the alignment of short reads (WFA2), on top of its adaptive pruning: alignment branches whose score falls this far behind are cut. 0 turns it off. The default changes no 150 bp alignment in tests; small values (50) lose alignments and change MAPQs. Long reads (pb, ont) are aligned without it: over their gene-long windows even 1000 lost alignments.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_X_DROP)))
                ("wfa_traceback", "How WFA2 recovers each alignment's operations: classic (WFA2's high-memory mode, which keeps every wavefront), singletrack (keeps only the match wavefronts and recomputes a gap's steps while tracing back: less memory to touch, and faster), or auto (default): the database's. Both give the same scores, but among equally good alignments they may choose differently, so SAM records and the model's features can differ. --build records its choice (auto: singletrack; classic before protal 0.8.2) in the database's build_settings.tsv and aligns its strain alleles, congener gaps and column weights with it; a run that aligns otherwise than its database was built is warned. A database built before 2026-10-11 was built with classic.", cxxopts::value<std::string>()->default_value("auto"));

        // Profiling options
        options.add_options("Profiling")
                ("no_profile", "Do NOT perform taxonomic profiling, only output alignments.")
                ("knob", "Prediction threshold, 0 to 1: taxa whose model probability is at least this are reported. Lower finds more of the taxa present, higher reports fewer absent ones. How much a change matters depends on the model and the samples, so choose it on data like yours. Default 0.5, or, for a model with knobs by depth (machine_learning_cmdline.py --depth-knobs; build_gtdb_database.py gives them to its PacBio and ONT models), the model's knob for each sample's depth; a --knob given applies to every sample. A model's calibrated calls (machine_learning_cmdline.py --fdr-calls) are used only with --fdr.", cxxopts::value<double>()->default_value("0.5"))
                ("fdr", "Expected share of false calls, 0 to 1, for a model with calibrated calls (machine_learning_cmdline.py --fdr-calls; build_gtdb_database.py --call-mode fdr): each sample reports its highest-scoring taxa as long as the mean of their 1 - probability stays at or below this, the model's scores made probabilities by its calibration and adjusted to the share of present taxa among the sample's candidates (a deep sample, with many more absent candidates, needs higher scores). Off by default: without --fdr the model's knob curve or --knob applies even to a model with calibrated calls (at GTDB r226 they called 0.001-0.007 F1 below the curve for every read type). --fdr F uses them at F; 0: off. Not with --knob.", cxxopts::value<double>())
                ("fragment_count", "How the model's `fragments` feature counts a taxon's fragments: mapq, those the MAPQ and length filters keep (as before 2026-10-12); soft, those plus each filtered one (a read whose best record is on the taxon) by the share the read EM leaves the taxon of the filtered records (em_own_share of the records the filters drop), so that a divergent strain whose reads tie with a congener's reference (MAPQ 0) is not counted as nothing; or auto (default): the database's. --build records its choice (auto: mapq) in build_settings.tsv, and its models must be trained on runs that count so (build_gtdb_database.py --fragment-count); a run that counts otherwise than its database is warned. Only the feature changes: the profile's fragments, depth, the sample's depth and the singleton rule keep the filters' count.", cxxopts::value<std::string>()->default_value("auto"))
                ("singleton_congener", "A taxon of a single fragment beside a species of its genus with at least this many fragments in the sample is not reported, whatever its score, if its read looks like that species': the abundance-weighted assignment of the reads' alternatives leaves it less than half of it, or its identity is below 0.95. 0 (default): no such rule; on GTDB r226 training data the model called none of the taxa the rule would have vetoed, so it is off (docs/claude/2026-10-03-false-positive-anatomy). The training dump's prediction column and the trainer (machine_learning_cmdline.py --singleton-congener) apply the same rule.", cxxopts::value<size_t>()->default_value("0"))
                ("depth_identity_margin", "Reads count towards a species' abundance when their identity is at most this far below that of its best-matching reads (98th percentile). Reads below that, e.g. of a relative the database lacks, still count for detection. The margin is the same on every gene unless --gene_conservation scales it. The default, 0.08, counts the reads of strains up to about 5% from the reference, of which 0.04 dropped up to 40%. 1 lets every read count. auto (default): the database's (build_settings.tsv), else 0.08.", cxxopts::value<std::string>()->default_value("auto"))
                ("gene_conservation", "Scale --depth_identity_margin per gene by how fast each gene diverges within species: db for the database's factors (gene_conservation.tsv, which --build estimates from --full_reference and stores in the database), or a file of them (geneid, factor, species; 1 for a gene of typical conservation). A gene's margin is then 0.03 for read errors plus the rest times its factor (0.08: 0.05 at factor 0.4, 0.10 at 1.4). none (default): the same margin on every gene, which did best summed over three simulated worlds, among them one with many congeners missing from the database. The factors (the file's, else the database's) give the model's conservation features (conserved_fast_depth_ratio, conserved_hit_share) whatever this says.", cxxopts::value<std::string>()->default_value("none"))
                ("model", "PMML model file: an existing path is used as is, otherwise <name> in the database (<name>.xml without an extension). Default: the database's model of each sample's read type: model_pe.xml (or, in older databases, model.xml) for paired-end, model_se.xml for single-end, model_PB.xml for PacBio and model_ONT.xml for ONT samples, trained at the database's sparse-first gate; a sample gated otherwise (--sparse_first) takes the database's model trained at the gate nearest its own where it has models of other gates (model_pe_c1.xml, model_PB_c32.xml, ...). --model replaces all of them unless --model_se, --model_pb or --model_ont is given.", cxxopts::value<std::string>()->default_value(""))
                ("model_se", "PMML model file for single-end samples, given as --model. Default: --model if given, else the database's model_se.xml.", cxxopts::value<std::string>()->default_value(""))
                ("model_pb", "PMML model file for PacBio samples, given as --model. Default: --model if given, else the database's model_PB.xml.", cxxopts::value<std::string>()->default_value(""))
                ("model_ont", "PMML model file for ONT samples, given as --model. Default: --model if given, else the database's model_ONT.xml.", cxxopts::value<std::string>()->default_value(""))
                ("model_kind", "Which kind of the database's presence models profiles the samples: gbm (gradient-boosted trees) or forest (a random forest), one for every read type, or read types and theirs (pe:forest,pb:gbm; a read type not listed: auto), or auto (default): the database's default models (model_pe.xml, ...). A database built by build_gtdb_database.py since 2026-10-13 holds both kinds of every read type at every gate it has models of: its default kind (--model, gbm by default) as model_pe.xml, model_pe_c1.xml, ..., the other as model_pe_forest.xml, model_pe_c1_forest.xml, ... (--add_model --model_kind); a sample takes the model of that kind trained at the gate nearest its own (--sparse_first), the default files where they are of that kind. A run whose read type has no model of the kind stops (exit 2). Not used for the reads whose model --model, --model_se, --model_pb or --model_ont names. With --add_model: the kind of each model, comma-separated as its files ('-': the read type's default model file, as without it), stored under the kind's name (--add_model pe.xml,pe_rf.xml --read_type pe,pe --model_kind -,forest stores model_pe.xml and model_pe_forest.xml); each model is checked to be of its kind.", cxxopts::value<std::string>()->default_value("auto"))
                ("profile_dir", "Override profile output directory. Takes precedence over the directory specified in the map file.", cxxopts::value<std::string>()->default_value(""))
                ("no_unknown_share", "Report the called species' abundances as shares of their summed depth (adding up to 1), as protal did before 2026-10-08. By default, where the database gives genome sizes (species_priors.tsv) and the SAM says how many reads were scanned, the profile ends with a line '?' of the genomes the called species do not explain (the reads they leave over, at their average genome size), and the abundances are shares of all genomes sequenced. <profile>.composition gives the numbers either way.")
                ("ncbi_taxids", "Add two columns to each line of the profile: the species' NCBI taxid and NCBI's name of it, from the database's ncbi_taxids.tsv (per GTDB species the NCBI species taxid most of its genomes have in GTDB's metadata). For most GTDB species NCBI names no species: their taxid is one of NCBI's placeholders ('uncultured Lachnospiraceae bacterium', 'Streptomyces sp. IB2014 016-6'), which NCBI files under the deepest rank it can name, and which several GTDB species may share (the table's named and shared_by columns say which). NA for a species the table lacks; the line '?' gets '?'. A database converted before 2026-10-11 gets the table with gtdb_to_protal_db.py --priors_only and --add_tables.");

        // Strain / SNP options
        options.add_options("Strains")
                ("no_strains", "Stay on species level: do not write strain MSAs or SNP tables. Variants are still called, as the model uses them, so profiles are the same with or without this flag.")
                ("no_phasing", "Write one strain MSA row per long-read sample (PacBio, ONT), the consensus of its reads. By default a sample whose reads show two or more strains of a species (on 3 genes or more) gets a row per strain, <sample>_hap1, _hap2, ..., the most abundant first, each called from its strain's reads: the reads are clustered by their alleles across the genes each one covers, and the blocks of genes that no read links are joined by the strains' shares of the reads where those tell which strain is which. <species>.haplotypes.tsv says how each block was phased.")
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
                ("flex_indels", "Score the flex cells of every seed lookup (the 8 bases either side of the k-mer's 15-base core, compared with the read's) allowing one insertion or deletion in each flank, for reads of these types: read types comma-separated (pe, se, pb, ont), all, or none (default). Without it an indel in a flank shifts the bases beyond it, which then match by chance: one next to the core costs ~6 of the flank's 8 bases, and a relative's cell can take the seed. With it a single indel costs one base, as a mismatch does, and only an exact match still scores 16 (the unique seeds are the same). Meant for ONT reads, whose errors are mostly indels, and for --sparse_first (with --sparse_seeding), where each lookup counts more; --sparse_first_score counts the bases so too. The scan takes ~3-4 times the instructions. It can change seeds, records and the model's features, and the database's models are trained without it. A map's FLEX_INDELS column sets it per sample.", cxxopts::value<std::string>()->default_value("none"))
                ("whole_read_alignment", "Align each short read as a whole into its gene window, as protal did before it aligned from the anchor's exact matches (slower; the results differ in a few alignments). Long reads are always aligned as a whole.")
                ("no_alignment_screen", "Align every candidate with WFA2, without the k-mer screen that refuses a candidate whose read and gene window share too few k-mers for any alignment within the score budget to exist (AlignmentScreen.h). The screen changes no alignment; this is for measuring it.")
                ("no_indel_bound", "ONT reads: screen the candidates by their k-mers, as other reads and as before 2026-10-09, instead of by the bound that refuses a candidate whose read and gene window are too far apart by their indel distance for any alignment within the score budget to exist (AlignmentScreen.h IndelBound; the k-mer screen refuses next to nothing at their identity floor). Neither changes an alignment; this is for measuring the bound.")
                ("no_single_seed_anchors", "With --sparse_first_pick grid: a short read's gene with one seed makes no anchor, as in other runs (two seeds of a gene on one diagonal), instead of anchoring on it (since 2026-10-13: a grid seed fixes the read's place on its gene by itself, and a read seeded from its grid cores alone has about 1/C of its seeds). For comparison.")
                ("direct_keys", "Hold the index's key map as the file has it, 16-bit offsets of all 2^30 cores (3.2 GB whatever the reference), instead of the compact one for the cores with values alone (the default since 2026-10-14): a presence bit per core with ranks (153 MB) and 3 bytes a used core (~0.45 GB at r226, less on a sparse index), built as the index is loaded (with --sparse_seeding, whose statistics read the direct one, converted after them). Every lookup gives the same values either way; a compact lookup reads a second line of the key map, fetched with the read's keys (docs/claude/2026-10-14-sparse-load-key-map).")
                ("compact_keys", "The compact key map, the default since 2026-10-14 (--direct_keys: the direct one); accepted for the runs that ask for it.")
                ("long_read_budget", "Long reads: once a candidate of a read's gene has aligned, the gene's other candidates are aligned with a budget of this many edits (as mismatches) more than the best so far; one that would cost more fails there and counts as a failed candidate (ZF). Saves WFA2 work on far relatives; MAPQ and the listed alternatives of such genes change, so the models should be retrained with it. 0 (default): every candidate gets the ANI floor's budget.", cxxopts::value<size_t>()->default_value("0"))
                ("sequential_load", "Load the database's parts one after another, as protal did before 0.7.6, rather than the index beside the genome preload and the single-threaded tables (taxonomy, models, gene conservation, suspect copies, species priors, gene neighbours) each on a thread of its own. For measuring; the run is the same either way.")
                ("sparse_first", "Look a read (pair; a long read over all its chunks) up by a sample of 1/C of its k-mer cores first (C = 1, 2, 4, 8, 16 or 32; those whose hash falls in 1/C of its range, stage1_probe::Level, or every C-th in read order: --sparse_first_pick) and give it up, before the rest of its seeding and any alignment, unless at least --sparse_first_hits of those lookups have a flex cell sharing --sparse_first_score or more of the 16 flank bases with the read. A read kept is seeded and aligned as without it. At GTDB r226 (score 14) C = 2 kept 9.4% of a gut sample's pairs, with 99.6% of those that have records, and C = 8 18% of a HiFi sample's reads, with 99.2% (docs/claude/2026-10-12-stage1-probe-r226.md). One C for every read type (0: off), read types and their C (pe:1,se:1,pb:8,ont:1; a read type not listed as the database has it), or auto (default): the database's gates. --build records its value as the database's gates (build_settings.tsv; auto: none, off), and its models must be trained on runs gated so (build_gtdb_database.py --sparse-first): reads given up count no failed candidates, which the models' failed_candidate_rate uses. A database may also hold models trained at other gates (model_pe_c1.xml, ...: build_gtdb_database.py --model-gates, protal --add_model --model_gate): a sample is profiled with the model trained at the gate nearest its own (on the scale of log2 C, off next to 1; between two, the lower). A run gated otherwise than any model of its database is warned. A map's SPARSE_FIRST column sets it per sample.", cxxopts::value<std::string>()->default_value("auto"))
                ("sparse_first_score", "With --sparse_first: the flank bases (of 16) a sampled lookup's best flex cell must share with the read to count (default: the database's, else 14; --build records it with the gates).", cxxopts::value<size_t>())
                ("sparse_first_hits", "With --sparse_first: the sampled lookups of a read (both mates of a pair together) that must count for it to be kept (default: the database's, else 1; --build records it with the gates).", cxxopts::value<size_t>())
                ("sparse_first_pick", "With --sparse_first: which 1/C of a read's cores the gate looks up first. hash: those whose hash falls in 1/C of its range (stage1_probe::Level; the only pick before 2026-10-13), the same cores in every read, so gaps fall by chance; stride: every C-th core of each mate (of each chunk of a long read) in read order (stage1_probe::StridePicked), spread evenly over the read, which keeps far more of the reads with records at the same C (locally at C 8, 99.4% against 92.0% of the pairs; docs/claude/2026-10-13-stride-pick) for a few more lookups. At C 1 both look up every core. One pick for every read type, read types and theirs (pe:stride,pb:hash; a read type not listed as the database has it), or auto (default): the database's picks (stride where it records none; hash for a database that records its gates but no pick, built on 2026-10-12/13). --build records its value with the gates (build_settings.tsv; auto: stride, the default since 2026-10-13), and the models must be trained on runs picked so; a run picked otherwise than its database is warned. A map's SPARSE_FIRST_PICK column sets it per sample. --sparse_seeding stands for an index of the hash pick's cores and takes the hash pick only. grid (since 2026-10-13; a database built with grid levels): every core of a read is looked up, and only the database's 1/C of each gene's cores count, those it keeps at C for the most conserved within the species (its strain alleles, else the columns' conservation in the family) and the least shared with other genes, about every C-th along the gene; a read kept is seeded from those alone, as a sparse index of 1/C of the values holds them (--sparse_index, which loads only those), a short read's gene anchored on a single seed too; on the whole index (--sparse_index off, or samples not all grid-picked) the grid only gates, and a read kept is seeded from every core.", cxxopts::value<std::string>()->default_value("auto"))
                ("sparse_seeding", "With --sparse_first: seed the reads kept from the cores of the gate's level alone (1/C of them), as an index holding only those cores would, and log after the index load which share of its values and keys such an index would hold, per C. For measuring what a sparse index (less memory) would cost in records and calls; the whole index is loaded all the same. Short reads at C 1 are seeded as without it.")
                ("sparse_index", "The index of a run whose every sample is gated with the grid pick (--sparse_first C, --sparse_first_pick grid; C 2 or more): auto (default) loads only the entries the smallest of its gates keeps (grid level log2 C, about 1/C of the values, the compact key map their cores alone); on: the same, and an error where it cannot; off: the whole index, on which the grid only gates (as in a run whose samples are not all grid-picked): a read kept is seeded from every core, a gene anchored on two seeds, as in other runs; emulate: the whole index, each sample seeded as the sparse index would seed it (its gate's grid entries alone, a short read's gene anchored on one seed too: the same records as auto, for comparison and for build_gtdb_database.py's --test-grid rows). Needs a database built with grid levels (since 2026-10-13).", cxxopts::value<std::string>()->default_value("auto"))
                ("stage1_probe", "Measure what the first stage of a two-stage protal would keep: every lookup's best flex score and whether its core is in a sparse index of 1/c of the cores (c = 1 to 32), or in a stride pick of every c-th core of a mate in read order, tallied per read (pair) for rules 'at least m lookups scoring T or more' over the reads and those with records, with the share of the lookups a gate would make, and per taxon of the best records; and per mate of a short read the best score of its first lookups. Writes misc/<sample>_stage1_probe.tsv (hash sample), _stage1_probe_stride.tsv (stride pick), _stage1_probe_taxa.tsv and _first_lookups_probe.tsv and a log line; the SAM and the profile are the same as without it (Stage1Probe.h, docs/claude/2026-10-12-stage1-probe).")
                ("taxon_statistics", "Write misc/<taxon>.statistics.tsv for every taxon with reads: its coverage, reads, ANI and MAPQ in each sample, and whether it is reported. Off by default: a sample on a GTDB-sized database has reads on thousands of taxa, and that many small files took 15 s on a network file system. The profile files (<prefix>.profile, .profile.log, .profile.genes.log) hold the same per sample.")
                ("full_sam_header", "List every gene of the database in the SAM header (@SQ), as protal did before; by default only the genes that alignments name are listed.")
                ("write_unmapped_reads", "Write an unmapped record (FLAG 4, no sequence; its ZF tag names the taxa the read seeded on) for every read that seeded on taxa but aligned nowhere, as protal did up to 0.7.6. By default these reads are counted per taxon in the SAM header instead (@CO protal failed candidates of unaligned reads), which is all the profiler takes from them: on a GTDB-sized database they were 95% of a sample's records. --full_sam_header, whose header is written before the reads, writes them too. A map's UNMAPPED_READS column (write or count) chooses per sample.")
                ("serial_index_passes", "With --build: count and place the reference's k-mers and compute the value pointers on one thread, as protal did before these ran in -t threads (slower; the index is the same).")
                ("profile_ahead", "With several samples: profile each sample whose SAM is complete while the next sample's reads are aligned, on a worker with a quarter of the threads beside the alignment's; the profiling stage after the alignment takes the rest on all threads. The profiles are the same either way. Measured on a 6-core laptop it gained nothing (the worker's CPU time came out of the alignment's); it may pay on a node where the profiling stage leaves cores idle.")
                ("index_batch_kb", "With --build: KB of the reference per batch in the parallel index passes (smaller batches are for testing).", cxxopts::value<size_t>()->default_value("1024"))
                ("build", "Build index from reference file with header format ()")
                ("no_compress", "With --build: write the database as separate, uncompressed files (index.prx, reference.fna, ...). By default --build writes the single-file database database.protal (zstd-compressed; see --no_bundle). protal reads every form.")
                ("no_bundle", "With --build or --compress_db: keep the database as separate compressed files (index.prx.zst, reference.fna.zst, reference.map, internal_taxonomy.dmp, build_settings.tsv, unique_kmers.tsv, gene_conservation.tsv, suspect_copies.tsv, species_priors.tsv, ncbi_taxids.tsv, species_neighbours.tsv, congener_gaps.tsv, strain_alleles.tsv, column_weights.tsv, ancestry_sites.bin, foreign_rates.tsv, gene_neighbours.tsv, gene_positions.tsv, the models model_*.xml) instead of packing them into database.protal.")
                ("compress_level", "With --build: zstd compression level (1-22). Higher levels compress more but more slowly (level 19: ~3 MB/s per thread, -t threads are used); decompression speed barely depends on it.", cxxopts::value<int>()->default_value(std::to_string(DEFAULT_COMPRESS_LEVEL)))
                ("compress_window_log", "With --build: zstd long-distance matching window, as log2 bytes (27 = 128 MB, capped at the frame size); finds repeats between distant related sequences. 0 turns it off.", cxxopts::value<int>()->default_value(std::to_string(DEFAULT_COMPRESS_WINDOW_LOG)))
                ("compress_frame_mb", "With --build or --compress_db: size of the independent zstd frames in MB (1-4095). protal loads a database with -t threads, one frame per thread at a time. 0 writes a single frame, which loads with one thread.", cxxopts::value<int>()->default_value(std::to_string(DEFAULT_COMPRESS_FRAME_MB)))
                ("compress_db", "Compress the database folder --db in place, without rebuilding it: index.prx (raw or compressed in an older way) in protal's column format, reference.fna as seekable zstd (see --compress_level, --compress_frame_mb, -t), all packed into database.protal (--no_bundle: kept as index.prx.zst, reference.fna.zst, ...). Everything is read back and compared before the old files are removed. Needs the index in memory. On a single-file database --db: adds the binary gene table (gene_table.bin) a run loads instead of reference.map and unique_kmers.tsv, and the binary forms of its column weights and strain alleles (column_weights.bin, strain_alleles.bin), which a run reads without parsing the text, where it has no current ones.")
                ("decompress_db", "Write the database --db as separate raw files (index.prx, reference.fna, ...) and remove its compressed files (index.prx.zst, reference.fna.zst, or database.protal), e.g. for older protal versions. zstd -d does not give a raw index.prx from protal's column format.")
                ("unpack_db", "Write the files of the single-file database --db into --unpack_dir (default: the folder it is in): index.prx.zst, an uncompressed reference.fna (as --preload_genomes_off needs), reference.map, internal_taxonomy.dmp, build_settings.tsv, unique_kmers.tsv, gene_conservation.tsv, suspect_copies.tsv, species_priors.tsv, ncbi_taxids.tsv, species_neighbours.tsv, congener_gaps.tsv, strain_alleles.tsv, column_weights.tsv, ancestry_sites.bin, foreign_rates.tsv, gene_neighbours.tsv and gene_positions.tsv (if it has them) and the models it has (model_pe.xml or model.xml, model_se.xml, model_PB.xml, model_ONT.xml). database.protal is kept; protal uses the separate files when both are there.")
                ("unpack_dir", "With --unpack_db: the folder to write the files into (default: the one database.protal is in).", cxxopts::value<std::string>()->default_value(""))
                ("add_model", "Store the PMML model FILE in the database --db as the model of the reads --read_type names (pe, se, pb or ont: model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml; pe if not given), replacing the one there. Several models, comma-separated, with as many read types (--add_model pe.xml,se.xml --read_type pe,se), are stored at once. Each model is checked first. database.protal is rewritten once, with its other parts copied as they are, not recompressed.", cxxopts::value<std::string>()->default_value(""))
                ("model_gate", "With --add_model: the sparse-first gate (--sparse_first C: 0 off, 1 to 32) the runs each model was trained on were gated at, comma-separated as its files (--add_model pe.xml,pe_c1.xml --read_type pe,pe --model_gate -,1). A model of another gate than the database's own is stored as model_<type>_c<C>.xml (model_pe_c1.xml, model_PB_c32.xml), which runs gated at C, or nearer C than any other of the database's models, are profiled with; '-' or the database's own gate (build_settings.tsv): the read type's model (model_pe.xml, ...). A read type may have one model per gate and kind (--model_kind). Default: every model the database's gate's.", cxxopts::value<std::string>()->default_value(""))
                ("add_tables", "Store tables in the database --db by their file names, replacing those there: foreign_rates.tsv (scripts/foreign_rates.py: per gene copy, the reads of a tiled scan of the full reference's marker genes and how many came from other species), congener_gaps.tsv, strain_alleles.tsv or column_weights.tsv (what --build writes), species_priors.tsv (gtdb_to_protal_db.py --priors_only: the priors with the species' genome sizes, for a database converted before 2026-10-08) or ncbi_taxids.tsv (gtdb_to_protal_db.py --priors_only too: the species' NCBI taxids for --ncbi_taxids, for a database converted before 2026-10-11); or the files a run does not read but a database carries (since 2026-10-13): build_metadata.tsv (build_gtdb_database.py: what it was built from and with, and its models' scores), gene2geneid.tsv and genome2tiid.tsv (gtdb_to_protal_db.py: its genes' and genomes' names), stored as they are. Several, comma-separated. Each table is read and checked first. database.protal is rewritten once, with its other parts copied as they are (in place where it held the table and only later members follow it, as for build_metadata.tsv, the last); a database of separate files gets the files copied beside them. With --add_model: the tables and the models in one write.", cxxopts::value<std::string>()->default_value(""))
                ("write_species_neighbours", "With --db FOLDER of separate files (a converted release before any build: reference.fna, reference.map, internal_taxonomy.dmp): write every species' nearest congeners by their references' marker genes (species_neighbours.tsv's table, which --build also stores in the database) to this file, without building. build_gtdb_database.py chooses the species its training database leaves out by it (species_clouds.tsv): a species complex is left out or kept whole.", cxxopts::value<std::string>()->default_value(""))
                ("full_reference", "With --build: the marker genes of all genomes (not only the representatives'), to check which k-mers are unique and to estimate how fast each gene diverges within species (gene_conservation.tsv, see --gene_conservation)", cxxopts::value<std::string>()->default_value(""))
                ("suspect_copy_distance", "With --build: a species' copy of a gene within this k-mer distance (about the share of bases that differ) of a copy of a species of another genus (or family, order, class, phylum, domain), and 0.02 farther from its nearest congener's copy or without one, is suspect: contamination or a transferred gene. The suspect copies go into the database (suspect_copies.tsv) and a run leaves their records out (see --keep_suspect_copies); every near pair across genera is reported in gene_incongruence.tsv beside the database. 0: no such scan.", cxxopts::value<double>()->default_value("0.02"))
                ("strain_alleles", "With --build: up to this many alleles per species and gene from the other genomes' copies in --full_reference, stored as edits of the representative's copy (strain_alleles.tsv), for the alignment scores (--no_allele_scores) and the profiler's 'alleles' features. Chosen among up to 16 distinct ones sampled by hash for the sample's coverage (greedily: together they explain as much of each sampled allele as a read of it would be scored with; farthest first until 2026-10-09), each nearer the representative than the gene's nearest congener's copy (congener_gaps.tsv) and within 0.1 of it. 0: none. Default 16 since 2026-10-11 (4 before): more of a species' strains are called, at no index cost (docs/claude/2026-10-11-strain-false-negatives).", cxxopts::value<size_t>()->default_value("16"))
                ("allele_genome_share", "With --build: the share of each species' genomes that may give alleles (--strain_alleles), by a hash of the genome's accession (the full reference's second header word, which the converter writes since 2026-10-08). build_gtdb_database.py passes 0.5 and simulates strains only from the other genomes, so that no simulated strain is its species' own allele. Below 1 it needs a full reference that names its genomes.", cxxopts::value<double>()->default_value("1"))
                ("index_alleles", "With --build: the k-mers of the strain alleles at least this far from the representative's copy (edits per base of the allele's range) go into the index under the species as non-unique entries, so that a read of a deep strain, which seeds on few of the representative's k-mers, finds its species and is aligned and scored against the representative as before. The seeds are left out of unique_kmers.tsv's totals. Off by default (1 or more: none) since the r226 v22 build, whose F1 gains with 0.01 (and the other changes of that build) were too small for the entries it adds (about three per edit, +5% on the test worlds).", cxxopts::value<double>()->default_value("1"))
                ("grid_choice", "With --build (a developer option, since 2026-10-13): how the grid levels (--sparse_first_pick grid) choose the core each window of C cores of a gene keeps. scored (default): the one the species' strain alleles (else the columns' conservation) most often leave a hit, less its key's ubiquity, anywhere in the window; even: every C-th core of the gene whatever the scores (a stride on the gene, its cores C apart); phase: every C-th core too, the phase chosen per gene and gate by the scores of the cores it keeps.", cxxopts::value<std::string>()->default_value("scored"))
                ("column_weights", "With --build: how conserved each column of every marker gene is in the species' family (column_weights.tsv, ColumnWeights.h): up to this many genera per family, chosen by hash, vote on the consensus of their references (the first gives the family's reference copy, whose columns the family's copies map onto), and every genus's species vote on their genus reference, averaged over genera; per column the weights -log(1 - p) of the two (half nats, 4 bits) and of the amino acid per codon, with the family's consensus base. For the 'weights' features (mismatches at conserved columns are errors or reads from afar, a strain's mutations avoid them) and the polarised ancestry sites of a species with fewer than three congeners (--no_site_weights at run time). 0: no table.", cxxopts::value<size_t>()->default_value("10"))
                ("column_weights_alignment", "With --build and --column_weights: where the family columns come from. auto (default): msa if the database folder has GTDB's alignments and tree (column_msa/, species_tree.nwk: gtdb_to_protal_db.py writes them from the release's *_msa_marker_genes_reps archive and tree), else aa. msa (since 2026-10-10): no alignment, each copy's residues placed on its row of GTDB-Tk's alignment of the marker, every genus of a family voting through its ancestral sequence (Fitch parsimony on GTDB's tree), the family's ancestral base polarising small genera. aa: the copies translated in frame and their proteins aligned (WFA2, up to 0.3 of the amino acids different within a genus, 0.5 among genera), each aligned amino-acid pair mapping its codon's bases. nt: the bases aligned (up to 0.2 and 0.4), which fails for most copies of GTDB's wide genera (docs/claude/2026-10-10-real-ancestry).", cxxopts::value<std::string>()->default_value("auto"))
                ("no_ancestry_sites", "With --build: do not store every gene copy's ancestry sites (ancestry_sites.bin, AncestrySites.h: where its copy differs from the consensus of its congeners' copies, the ancestry features' sites), computed by default from the species neighbours, congener gaps and column weights the build writes, as a run would compute them on a copy's first record. Without the table a run computes them there, as before 2026-10-11, at ~0.3 M instructions plus ~0.6 M per congener compared for every copy its records touch (half of the profiling of a sample built at r226 v23, docs/claude/2026-10-11-r226-v23); the features are the same.")
                ("no_unique_cores", "With --build: do not store which of every gene copy's index cores are unique (unique_cores.bin, UniqueCores.h: two bits a core, ~0.7 GB at GTDB r226). With the table a run counts each record's ZS, ZU and ZT from the gene's cores its alignment covers with exact k-mers, the same on the whole index, a sparse index or with --sparse_seeding; without it (or with this option in a run) from the seeds of the record's anchor, as before 2026-10-14, which a sparse index finds a quarter of at C 4. --build records which (build_settings.tsv unique_kmers: cores or anchors); a run with this option on a database built with the table is warned: the models' unique k-mer features (ZS, ZU, ZT) were trained on the cores.")
                ("reference", "Set of reference sequences to build the internal alignment database from", cxxopts::value<std::string>()->default_value(""))
                ("build_gene_subset", "With --build: a file of the gene ids (the numbers of reference.map / gene2geneid.tsv, one per line, # lines are comments) to build the database from; the other genes of --reference stay in reference.fna but get no k-mers, no unique k-mer row, no conservation factor and no suspect copies, and their copies in --full_reference are skipped. Refused when the folder has gene_neighbours.tsv (counted over every gene): derive a folder of the subset with scripts/mini_db/gtdb_to_protal_db.py --from_db --genes, which also leaves the other genes' sequences out (docs/databases.md, reduced marker sets).", cxxopts::value<std::string>()->default_value(""))
                ("preload_genomes_off", "Do not preload complete reference library (reference.fna and reference.map in protal index folder) and instead do dynamic loading. This usually decreases performance but saves memory. Needs the database as separate files with an uncompressed reference.fna (not reference.fna.zst or database.protal; see --unpack_db).")

                ("profile_truth", "Truth files, one per sample (comma-separated; a map's PROFILE_TRUTH column does the same). Each line names a species present, by GTDB lineage in any tab-separated field (d__...;s__Genus species, as simulate_metagenomes --protal_metafile writes) or by internal taxid in the first field. Profiles are annotated with TP/FP (<profile>.truth_annotated).", cxxopts::value<std::string>()->default_value(""))
                ("benchmark_alignment", "Benchmark alignment part of protal based on true taxonomic id and gene id supplied in the read header. Header must fulfill the formatting >taxid_geneid... with the regex: >[0-9]+_[0-9]+([^0-9]+.*)*")
                ("benchmark_alignment_output", "Benchmark alignment output. Output is appended to the file.", cxxopts::value<std::string>());


        // General options
        options.add_options("General")
                ("v,version", "Output version information")
                ("feature_versions", "Print the presence model's features with the version of each one's definition in this protal, the date of its last change and a note (TSV: feature, version, changed, note), and exit. The trainer stores them in the models it writes (protal_feature_versions), and a run warns when a model was trained on other definitions.")
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

        // Not shown: options kept so that older command lines still run.
        options.add_options("_kept")
                ("strain_spill", "Up to 2026-10-12 the folder of the strain stage's spill files; now the same as --scratch_dir, which spills there (--scratch_dir wins).", cxxopts::value<std::string>()->default_value(""));

        return options;
    }

    std::string Options::ToString() const {

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
        if (!m_scratch_dir.empty()) result_str << "scratch dir:         " << m_scratch_dir << '\n';
        result_str << "preload genomes:     " << (m_preload_genomes ? "yes" : "no") << '\n';
        result_str << "----- Alignment -----" << std::string(30, '-') << '\n';
        // A run setting not given (auto) follows the database's build_settings.tsv where it records one (the block
        // ApplyRunSettings prints after these options); here protal's default.
        auto const& g = m_run_settings_given;
        auto run_auto = [this](bool is_given) { return is_given || m_build ? "" : " (auto: the database's, else this)"; };
        result_str << "align top:           " << std::to_string(m_align_top) << run_auto(g.align_top.has_value()) << '\n';
        result_str << "adaptive candidates: " << std::to_string(m_adaptive_candidates) << run_auto(g.adaptive_candidates.has_value()) << '\n';
        result_str << "allele scores:       " << (m_no_allele_scores ? "off (--no_allele_scores)" : "the database's strain alleles, if it has them") << run_auto(g.allele_scores.has_value()) << '\n';
        result_str << "weighted site shift: " << (m_weighted_site_shift ? "with the database's column weights, if it has them (--weighted_site_shift)" : "off (--weighted_site_shift turns it on)") << run_auto(g.weighted_site_shift.has_value()) << '\n';
        result_str << "max key ubiquity:    " << std::to_string(m_max_key_ubiquity) << '\n';
        result_str << "max seed size:       " << std::to_string(m_max_seed_size) << '\n';
        result_str << "max score ani:       " << per_read_type(m_max_score_ani, [this](ReadType type) { return GetMaxScoreAni(type); }) << run_auto(m_max_score_ani_given) << '\n';
        result_str << "x-drop:              " << std::to_string(m_x_drop) << (long_reads ? " (short reads; long reads: none)" : "") << '\n';
        result_str << "WFA2 traceback:      " << (m_wfa_traceback_given ? std::string(WfaTracebackName(*m_wfa_traceback_given)) :
                                                  m_build ? "singletrack (auto)" : "auto (the database's)") << '\n';
        result_str << "fragment count:      " << (m_fragment_count_given ? std::string(FragmentCountName(*m_fragment_count_given)) :
                                                  m_build ? "mapq (auto)" : "auto (the database's)") << '\n';
        result_str << "short reads aligned: " << (m_whole_read_alignment ? "as a whole" : "from their anchors") << '\n';
        result_str << "alignment screen:    " << (m_no_alignment_screen ? "off (--no_alignment_screen)" : "k-mers shared with the window before WFA2") << '\n';
        if (long_reads) {
            result_str << "indel bound:         " << (m_no_alignment_screen ? "off (--no_alignment_screen)" :
                                                      m_no_indel_bound ? "off (--no_indel_bound: ONT reads screened by their k-mers)" :
                                                      "ONT reads screened by their indel distance to the window instead") << '\n';
        }
        if (std::any_of(m_flex_indels.begin(), m_flex_indels.end(), [](bool on) { return on; }) || !m_flex_indels_list.empty()) {
            result_str << "flex indels:         " << FlexIndelsText(m_flex_indels) << " (an indel per flank allowed in the lookups)"
                       << (m_flex_indels_list.empty() ? "" : "; per sample: the map's FLEX_INDELS") << '\n';
        }
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
        result_str << "foreign genes:       " << (m_keep_foreign_genes ? "kept" : "left out of depth and MSAs") << run_auto(g.keep_foreign_genes.has_value()) << '\n';
        result_str << "suspect gene copies: " << (m_keep_suspect_copies ? "kept" : "left out of the evidence") << run_auto(g.keep_suspect_copies.has_value()) << '\n';
        if (m_add_model.empty()) {
            result_str << "model kind:          " << (m_model_kind == "auto" ? "auto (the database's default models)" : m_model_kind) << '\n';
        }
        result_str << "---- Dev Options ----" << std::string(30, '-') << '\n';
        result_str << "verbose:             " << (m_verbose ? "yes" : "no") << '\n';
        result_str << "benchmark alignment: " << (m_benchmark_alignment ? "yes" : "no") << '\n';
        result_str << "profile truth:       " << (HasProfileTruths() ? "yes" : "no") << '\n';
        if (m_stage1_probe) result_str << "stage-1 probe:       yes (misc/<sample>_stage1_probe*.tsv)" << '\n';
        if (m_sparse_seeding) result_str << "sparse seeding:      yes (kept reads seeded from the gate's cores alone)" << '\n';
        bool const picks_given = std::any_of(m_sparse_first_pick_given.begin(), m_sparse_first_pick_given.end(),
                                             [](auto const& p) { return p.has_value(); });
        if (std::any_of(m_sparse_first_given.begin(), m_sparse_first_given.end(), [](auto const& g) { return g.has_value(); }) ||
            m_sparse_first_score_given || m_sparse_first_hits_given || !m_sparse_first_list.empty() || picks_given ||
            !m_sparse_first_pick_list.empty()) {
            result_str << "sparse first:        " << SparseFirstGatesText(m_sparse_first_given)
                       << " (read types not listed: " << (m_build ? "off" : "the database's") << "), score "
                       << (m_sparse_first_score_given ? std::to_string(*m_sparse_first_score_given) : "auto") << ", hits "
                       << (m_sparse_first_hits_given ? std::to_string(*m_sparse_first_hits_given) : "auto")
                       << (picks_given ? ", pick " + SparseFirstPicksText(m_sparse_first_pick_given) : std::string())
                       << (m_sparse_first_list.empty() ? "" : "; per sample: the map's SPARSE_FIRST")
                       << (m_sparse_first_pick_list.empty() ? "" : "; per sample: the map's SPARSE_FIRST_PICK") << '\n';
        }
        result_str << "^ output:            " << m_benchmark_alignment_output << '\n';
        result_str << "---------------------" << std::string(30, '-') << '\n';
        return result_str.str();
    }

    void Options::PrintHelp(bool show_dev) {
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

    bool Options::LoadFromMap(std::string map_path, std::string& output_dir, std::string& strain_output_dir, std::string& misc_output_dir,
                              std::vector<std::string>& prefix_list,
                              std::vector<std::string>& first_list, std::vector<std::string>& second_list,
                              std::vector<std::string>& sam_list, std::vector<std::string>& profile_list,
                              std::vector<std::string>& samplenames_list, std::vector<std::string>& profile_truth_list,
                              std::vector<std::string>& read_type_list, std::vector<std::string>& unmapped_reads_list,
                              std::string const& sam_ending, bool* sams_named,
                              std::vector<std::string>* sparse_first_list,
                              std::vector<std::string>* flex_indels_list,
                              std::vector<std::string>* sparse_first_pick_list) {
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
        int sparse_first_column = -1;
        int flex_indels_column = -1;
        int sparse_first_pick_column = -1;

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
                        if (token == MAP_SPARSE_FIRST) {
                            if (sparse_first_column != -1) {
                                std::cerr << "Column '" << MAP_SPARSE_FIRST << "' is defined twice" << std::endl;
                                return false;
                            }
                            sparse_first_column = i;
                        }
                        if (token == MAP_FLEX_INDELS) {
                            if (flex_indels_column != -1) {
                                std::cerr << "Column '" << MAP_FLEX_INDELS << "' is defined twice" << std::endl;
                                return false;
                            }
                            flex_indels_column = i;
                        }
                        if (token == MAP_SPARSE_FIRST_PICK) {
                            if (sparse_first_pick_column != -1) {
                                std::cerr << "Column '" << MAP_SPARSE_FIRST_PICK << "' is defined twice" << std::endl;
                                return false;
                            }
                            sparse_first_pick_column = i;
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
                    // None of them is made here: a folder is made when something is written into it (a SAM's or
                    // profile's folder by its writer, misc/ and strains/ at the start of a run that writes them),
                    // so a map without SAM or PROFILE columns, whose files go to the output folder itself, leaves
                    // no empty alignments/ or profiles/ (it did up to 2026-10-13).
                    sam_output_dir = in_output(sam_output_dir, MAP_VAR_DEFAULT_SAM_OUTPUT_DIR);
                    strain_output_dir = in_output(strain_output_dir, MAP_VAR_DEFAULT_STRAIN_OUTPUT_DIR);
                    profile_output_dir = in_output(profile_output_dir, MAP_VAR_DEFAULT_PROFILE_OUTPUT_DIR);
                    misc_output_dir = in_output(misc_output_dir, MAP_VAR_DEFAULT_MISC_OUTPUT_DIR);


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
                        { read_type_column, &MAP_READ_TYPE }, { unmapped_reads_column, &MAP_UNMAPPED_READS },
                        { sparse_first_column, &MAP_SPARSE_FIRST }, { flex_indels_column, &MAP_FLEX_INDELS },
                        { sparse_first_pick_column, &MAP_SPARSE_FIRST_PICK } };
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
                // A cell may list several files, comma-separated, read one after the other as one sample (its lanes or
                // runs), each in INPUT_DIR; a paired-end sample's FIRST and SECOND files pair up by position.
                auto read_files = [&](int column, std::string& joined) {
                    auto const files = SplitReadFiles(tokens[column]);
                    std::vector<std::string> paths;
                    for (auto const& file : files) {
                        if (file.empty()) {
                            std::cerr << "Line " << line_num << ": an empty file name in column " << column + 1 << " ("
                                      << (column == first_column ? MAP_FIRST_READ : MAP_SECOND_READ) << ": '" << tokens[column]
                                      << "'); a sample's files are separated by single commas" << std::endl;
                            return false;
                        }
                        paths.push_back(path(input_dir).append(file).string());
                    }
                    joined = JoinReadFiles(paths);
                    return true;
                };
                std::string first_path, second_path;
                if (!read_files(first_column, first_path)) return false;
                // No second file ('-', or no SECOND column): single-end reads.
                bool const single_end = second_column == -1 || tokens[second_column] == MAP_NO_SECOND_READ;
                if (!single_end) {
                    if (!read_files(second_column, second_path)) return false;
                    size_t const firsts = SplitReadFiles(first_path).size(), seconds = SplitReadFiles(second_path).size();
                    if (firsts != seconds) {
                        std::cerr << "Line " << line_num << ": " << MAP_FIRST_READ << " lists " << firsts << " file(s) and "
                                  << MAP_SECOND_READ << " " << seconds << "; a sample's first and second files pair up "
                                  << "by position (lanes: L001_R1,L002_R1 and L001_R2,L002_R2)" << std::endl;
                        return false;
                    }
                }

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
                if (sparse_first_column != -1) {
                    // A sample's gate (C), or '-': its read type's (--sparse_first, else the database's).
                    auto const& value = tokens[sparse_first_column];
                    if (value != MAP_NO_SECOND_READ && !ParseSparseFirstLevel(value)) {
                        std::cerr << "Line " << line_num << ": " << MAP_SPARSE_FIRST << " is '" << value << "'; expected "
                                  << "0 (off), 1, 2, 4, 8, 16, 32 or -" << std::endl;
                        return false;
                    }
                    if (sparse_first_list) sparse_first_list->emplace_back(value);
                }
                if (flex_indels_column != -1) {
                    // A sample's --flex_indels (1 on, 0 off), or '-': its read type's.
                    auto const& value = tokens[flex_indels_column];
                    if (value != "1" && value != "0" && value != MAP_NO_SECOND_READ) {
                        std::cerr << "Line " << line_num << ": " << MAP_FLEX_INDELS << " is '" << value << "'; expected 1, 0 or -" << std::endl;
                        return false;
                    }
                    if (flex_indels_list) flex_indels_list->emplace_back(value);
                }
                if (sparse_first_pick_column != -1) {
                    // A sample's gate's pick (hash or stride), or '-': its read type's.
                    auto const& value = tokens[sparse_first_pick_column];
                    if (value != MAP_NO_SECOND_READ && !ParseSparseFirstPick(value)) {
                        std::cerr << "Line " << line_num << ": " << MAP_SPARSE_FIRST_PICK << " is '" << value
                                  << "'; expected hash, stride, grid or -" << std::endl;
                        return false;
                    }
                    if (sparse_first_pick_list) sparse_first_pick_list->emplace_back(value);
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
        if (sams_named) *sams_named = !sam_list.empty();
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

    bool Options::PrepareAndCheckValidity(bool force_read_check) {
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
        // --add_model and --add_tables store their files in one write.
        if (int(m_compress_db) + int(m_decompress_db) + int(m_unpack_db) + int(add_model || add_tables) + int(write_neighbours) +
            int(m_build) > 1) {
            error_log.emplace_back("--build, --compress_db, --decompress_db, --unpack_db, --add_model or --add_tables, and "
                                   "--write_species_neighbours cannot be combined (--build compresses unless --no_compress)");
        }
        if (!add_model) {
            // A run's kinds of models (with --add_model, each model's kind: checked below).
            std::string error;
            if (auto const kinds = ModelKindsOf(m_model_kind, error)) m_model_kinds = *kinds;
            else error_log.emplace_back(error);
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
            // --add_model FILE[,FILE...] with --read_type TYPE[,TYPE...]: one read type per model, each once per gate
            // (--model_gate: one per model, '-' the database's).
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
            size_t const gates = m_model_gate.empty() ? 0 : static_cast<size_t>(std::count(m_model_gate.begin(), m_model_gate.end(), ',')) + 1;
            if (gates > 0 && gates != models) {
                error_log.emplace_back("--model_gate must give a gate per model of --add_model (" + std::to_string(models) +
                                       "), comma-separated ('-': the database's own), not '" + m_model_gate + "'");
            }
            std::istringstream gate_tokens(m_model_gate);
            std::string token;
            std::vector<std::string> gate_of;
            while (!m_model_gate.empty() && std::getline(gate_tokens, token, ',')) {
                if (token != "-" && !ParseSparseFirstLevel(token)) {
                    error_log.emplace_back("--model_gate: a gate is -, 0, 1, 2, 4, 8, 16 or 32, not '" + token + "'");
                }
                gate_of.push_back(token);
            }
            // --model_kind: one kind per model ('-': the read type's default file), or auto (every model so).
            std::vector<std::string> kind_of;
            if (m_model_kind != "auto") {
                kind_of = SplitCommas(m_model_kind);
                if (kind_of.size() != models) {
                    error_log.emplace_back("--model_kind with --add_model gives a kind per model (" + std::to_string(models) +
                                           "), comma-separated ('-': the read type's default model file), not '" + m_model_kind + "'");
                }
                for (auto const& kind : kind_of) {
                    if (kind != "-" && !IsModelKind(kind)) {
                        error_log.emplace_back("--model_kind: a model's kind is -, " + ModelKindNames() + ", not '" + kind + "'");
                    }
                }
            }
            std::istringstream type_tokens(m_read_type);
            std::set<std::string> seen;
            for (size_t i = 0; !m_read_type.empty() && std::getline(type_tokens, token, ','); i++) {
                std::string const gate = i < gate_of.size() ? gate_of[i] : "-";
                std::string const kind = i < kind_of.size() ? kind_of[i] : "-";
                if (!ReadTypeFromToken(token)) {
                    error_log.emplace_back("--read_type must be one of " + ReadTypeTokens() + ", not '" + token + "'");
                } else if (!seen.insert(token + "@" + gate + "@" + kind).second) {
                    error_log.emplace_back("--read_type names " + token + " twice for --add_model" +
                                           (gate == "-" ? std::string() : " at --model_gate " + gate) +
                                           (kind == "-" ? std::string() : " of --model_kind " + kind));
                }
            }
            for (auto const& model : AddModels()) {
                if (!std::filesystem::is_regular_file(model.file)) error_log.emplace_back("--add_model: " + model.file + " is not a file");
            }
        } else if (!m_model_gate.empty()) {
            error_log.emplace_back("--model_gate goes with --add_model");
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
            // The model of each kind of reads the samples have (of paired-end reads without samples); one of the kind
            // --model_kind asks for is looked for once the database's gates are known (RunProtal.h KindModels).
            bool const any_sample = std::any_of(kReadTypes.begin(), kReadTypes.end(), [this](ReadTypeInfo const& t) { return AnySample(t.type); });
            for (auto const& info : kReadTypes) {
                auto const type = info.type;
                if (m_no_profile || write_neighbours || !(AnySample(type) || (type == ReadType::Paired && !any_sample)) ||
                    ModelDbFile(type).Exists() || (ModelNamed(type).empty() && !ModelKind(type).empty())) continue;
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
        if (m_column_weights_alignment != "auto" && m_column_weights_alignment != "msa" && m_column_weights_alignment != "aa" && m_column_weights_alignment != "nt") {
            error_log.emplace_back("--column_weights_alignment must be auto, msa, aa or nt, not " + m_column_weights_alignment);
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

        // --scratch_dir: a folder protal can make its run's scratch folder in (made if missing; ScratchSpace.h).
        if (!m_scratch_dir.empty()) {
            if (m_build || db_mode) {
                warning_log.emplace_back("--scratch_dir is not used with --build or the database commands");
            } else {
                std::error_code ec;
                std::filesystem::create_directories(m_scratch_dir, ec);
                if (!std::filesystem::is_directory(m_scratch_dir, ec) || ::access(m_scratch_dir.c_str(), W_OK) != 0) {
                    error_log.emplace_back("--scratch_dir " + m_scratch_dir + " is not a folder protal can write to");
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

        // Check files: each of a sample's (several: its lanes, read one after the other).
        for (auto i = 0; i < m_first_list.size(); i++) {

            auto const firsts = SplitReadFiles(m_first_list[i]);
            auto const seconds = SplitReadFiles(m_second_list[i]);  // none: single-end
            auto sam = m_sam_list[i];
            auto sam_exists = std::filesystem::exists(sam);

            if (firsts.empty()) {
                error_log.emplace_back("sample " + std::to_string(i + 1) + " has no read file (-1, or FIRST in a map)");
            } else if (!seconds.empty() && seconds.size() != firsts.size()) {
                error_log.emplace_back("sample " + std::to_string(i + 1) + " has " + std::to_string(firsts.size()) +
                                       " first-in-pair and " + std::to_string(seconds.size()) + " second-in-pair read files; "
                                       "a sample's files pair up by position");
            }
            for (auto const& [flag, files] : { std::pair{ "-1", firsts }, std::pair{ "-2", seconds } }) {
                for (auto const& file : files) {
                    if (!std::filesystem::exists(file)) {
                        if (!force_read_check && sam_exists) {
                            warning_log.emplace_back(std::string(flag) + " file does not exist: " + file + " ( But .sam file does )");
                        } else {
                            error_log.emplace_back(std::string(flag) + " file does not exist: " + file);
                        }
                    } else if (force_read_check || !sam_exists) {
                        // A file that exists but cannot be read (permissions, a directory) would read as empty.
                        std::string const problem = Unreadable(file);
                        if (!problem.empty()) error_log.emplace_back(std::string(flag) + " file cannot be read: " + file + " (" + problem + ")");
                    }
                }
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
                if (GetReadType(i) != ReadType::Single || i >= m_first_list.size() || m_first_list[i].empty()) continue;
                auto const first_file = SplitReadFiles(m_first_list[i]).front();  // of several (lanes), the first
                size_t const longest = LongestRead(first_file, 100);
                if (longest > MAX_SHORT_READ_LENGTH) {
                    error_log.emplace_back("Sample " + GetSampleId(i) + " has reads of up to " + std::to_string(longest) +
                                           " bp (first 100 of " + first_file + "), too long for short reads: give "
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

    Options Options::OptionsFromArguments(int argc, char *argv[]) {
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

        if (result.count("feature_versions")) {  // FeatureVersions.h
            feature_versions::WriteTable(std::cout);
            exit(0);
        }

        bool show_help = result.count("help");
        bool show_help_dev = result.count("full_help");
        bool show_map_help = result.count("map_help");
        bool show_version = result.count("version");

        if (show_help || show_version || show_map_help || show_help_dev) {
            return { show_help, show_map_help, show_version, show_help_dev };
        }

        // The run settings --build records (BuildSettings.h RunSettings): each as given, nothing for auto (the database's;
        // RunProtal ApplyRunSettings). A value that is none exits 2.
        build_settings::RunSettings given;
        {
            auto bad = [](std::string const& option, std::string const& expected, std::string const& value) {
                std::cerr << "Error: --" << option << " must be auto or " << expected << ", not '" << value << "'" << std::endl;
                exit(2);
            };
            auto text = [&result](std::string const& option) { return result[option].as<std::string>(); };
            auto count = [&](std::string const& option, std::optional<size_t>& target) {
                if (text(option) == "auto") return;
                target = build_settings::ParseCount(text(option));
                if (!target) bad(option, "a whole number", text(option));
            };
            auto decimal = [&](std::string const& option, std::optional<double>& target) {
                if (text(option) == "auto") return;
                target = build_settings::ParseDecimal(text(option));
                if (!target || *target < 0) bad(option, "a decimal of 0 or more", text(option));
            };
            auto on_off = [&](std::string const& option, std::optional<bool>& target) {
                if (text(option) == "auto") return;
                target = build_settings::ParseOnOff(text(option));
                if (!target) bad(option, "on or off", text(option));
            };
            if (result.count("no_weighted_site_shift")) given.weighted_site_shift = false;
            else if (result.count("weighted_site_shift")) given.weighted_site_shift = true;
            if (result.count("no_allele_scores")) given.allele_scores = false;
            else if (result.count("allele_scores")) given.allele_scores = true;
            count("align_top", given.align_top);
            count("adaptive_candidates", given.adaptive_candidates);
            if (result.count("no_mate_guidance")) given.mate_guidance = std::string(build_settings::kMateGuidanceOff);
            std::optional<double> identity_floor;
            decimal("max_score_ani", identity_floor);
            given.max_score_ani.fill(identity_floor);
            if (result.count("no_unique_cores")) given.unique_kmers = std::string(build_settings::kUniqueAnchors);
            on_off("keep_foreign_genes", given.keep_foreign_genes);
            on_off("keep_suspect_copies", given.keep_suspect_copies);
            decimal("depth_identity_margin", given.depth_identity_margin);
        }

        bool no_strains = result.count("no_strains");
        bool no_mate_guidance = given.mate_guidance == std::string(build_settings::kMateGuidanceOff);
        bool no_gene_neighbours = result.count("no_gene_neighbours");
        bool keep_foreign_genes = given.keep_foreign_genes.value_or(false);
        bool keep_suspect_copies = given.keep_suspect_copies.value_or(false);
        bool no_phasing = result.count("no_phasing");
        bool build = result.count("build");
        bool preload_genomes_off = result.count("preload_genomes_off");
        bool benchmark_alignment = result.count("benchmark_alignment");
        bool fastalign = false;//result.count("fastalign");
        bool no_profile = result.count("no_profile");
        bool verbose = result.count("verbose");

        size_t threads = result["threads"].as<size_t>();
        size_t align_top = given.align_top.value_or(DEFAULT_ALIGN_TOP);
        size_t x_drop = result["x_drop"].as<size_t>();
        std::optional<WfaTraceback> wfa_traceback;
        if (auto const name = result["wfa_traceback"].as<std::string>(); name != "auto") {
            wfa_traceback = ParseWfaTraceback(name);
            if (!wfa_traceback) {
                std::cerr << "Error: --wfa_traceback must be auto, classic or singletrack, not '" << name << "'" << std::endl;
                exit(2);
            }
        }
        std::optional<FragmentCount> fragment_count;
        if (auto const name = result["fragment_count"].as<std::string>(); name != "auto") {
            fragment_count = ParseFragmentCount(name);
            if (!fragment_count) {
                std::cerr << "Error: --fragment_count must be auto, mapq or soft, not '" << name << "'" << std::endl;
                exit(2);
            }
        }
        size_t max_key_ubiquity = result["max_key_ubiquity"].as<size_t>();
        size_t min_successful_lookups = result["min_successful_lookups"].as<size_t>();
        size_t max_seed_size = result["max_seed_size"].as<size_t>();
        double max_score_ani = given.max_score_ani.front().value_or(DEFAULT_MAX_SCORE_ANI);  // -a: every read type's
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
        bool const one_sample = result.count("one_sample") > 0;
        if (one_sample && !map_file.empty()) {
            std::cerr << "--one_sample is for the read files of -1/-2; in a map, a sample's FIRST and SECOND cells list "
                         "its files, comma-separated" << std::endl;
            exit(2);
        }
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
        std::vector<std::string> sparse_first_list;  // a map's SPARSE_FIRST
        std::vector<std::string> flex_indels_list;   // a map's FLEX_INDELS
        std::vector<std::string> sparse_first_pick_list;  // a map's SPARSE_FIRST_PICK
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

        bool sams_named = false;  // by a map's SAM column
        if (!map_file.empty()) {
            if (!std::filesystem::exists(map_file)) {
                std::cerr << "Map file " << map_file << " does not exist." << std::endl;
                exit(9);
            }
            if (!LoadFromMap(map_file, output_dir, strain_output_dir, misc_output_dir, prefix_list, first_list, second_list, sam_list, profile_list, samplenames_list, profile_truth_list, read_type_list, unmapped_reads_list, sam_ending, &sams_named, &sparse_first_list, &flex_indels_list, &sparse_first_pick_list)) {
                std::cerr << "Failed to read map file " << map_file << " (see --map_help)." << std::endl;
                exit(9);
            }

            if (range_arg != "") {
                range = ProcessRange(range_arg, first_list.size());
            }
        } else {
            LineSplitter::Split(first, ",", first_list);
            LineSplitter::Split(second, ",", second_list);
            // --one_sample: the files are one sample's (its lanes or runs), read one after the other.
            if (one_sample && !first_list.empty()) {
                if (!second.empty() && first_list.size() != second_list.size()) {
                    std::cerr << "With --one_sample, -1/--first and -2/--second must name the same number of files "
                                 "(they pair up by position): -1 names " << first_list.size() << ", -2 "
                              << second_list.size() << std::endl;
                    exit(31);
                }
                first_list = { JoinReadFiles(first_list) };
                second_list = { second.empty() ? std::string() : JoinReadFiles(second_list) };
            }
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
        bool no_indel_bound = result.count("no_indel_bound");
        bool no_single_seed_anchors = result.count("no_single_seed_anchors");
        if (result.count("compact_keys") && result.count("direct_keys")) {
            std::cerr << "Error: --compact_keys and --direct_keys ask for different key maps: give one (the compact one is the "
                      << "default)" << std::endl;
            exit(2);
        }
        bool compact_keys = !result.count("direct_keys");
        size_t long_read_budget = result["long_read_budget"].as<size_t>();
        bool sequential_load = result.count("sequential_load");
        bool stage1_probe = result.count("stage1_probe");
        bool sparse_seeding = result.count("sparse_seeding");
        std::string const sparse_index = result["sparse_index"].as<std::string>();
        if (sparse_index != "auto" && sparse_index != "on" && sparse_index != "off" && sparse_index != "emulate") {
            std::cerr << "Error: --sparse_index is auto, on, off or emulate, not '" << sparse_index << "'" << std::endl;
            exit(2);
        }
        auto const sparse_first = ParseSparseFirstGates(result["sparse_first"].as<std::string>());
        if (!sparse_first) {
            std::cerr << "Error: --sparse_first is auto, one of 0 (off), 1, 2, 4, 8, 16 or 32, or read types and theirs "
                      << "(e.g. pe:1,se:1,pb:8,ont:1; read types " << ReadTypeTokens() << "), not '"
                      << result["sparse_first"].as<std::string>() << "'" << std::endl;
            exit(2);
        }
        std::optional<size_t> const sparse_first_score =
            result.count("sparse_first_score") ? std::optional<size_t>(result["sparse_first_score"].as<size_t>()) : std::nullopt;
        std::optional<size_t> const sparse_first_hits =
            result.count("sparse_first_hits") ? std::optional<size_t>(result["sparse_first_hits"].as<size_t>()) : std::nullopt;
        if ((sparse_first_score && *sparse_first_score > 16) || (sparse_first_hits && *sparse_first_hits == 0)) {
            std::cerr << "Error: --sparse_first_score is 0 to 16 (flank bases), --sparse_first_hits 1 or more" << std::endl;
            exit(2);
        }
        auto const sparse_first_pick = ParseSparseFirstPicks(result["sparse_first_pick"].as<std::string>());
        if (!sparse_first_pick) {
            std::cerr << "Error: --sparse_first_pick is auto, hash, stride, grid, or read types and theirs (e.g. pe:stride,pb:hash; "
                      << "read types " << ReadTypeTokens() << "), not '" << result["sparse_first_pick"].as<std::string>() << "'"
                      << std::endl;
            exit(2);
        }
        // The probe measures what a gate would give up: a gate given with it is an error; the database's (auto) are
        // not used then (ApplyBuildSettings).
        if (stage1_probe && std::any_of(sparse_first->begin(), sparse_first->end(), [](auto const& g) { return g && *g > 0; })) {
            std::cerr << "Error: --stage1_probe measures what --sparse_first would give up; give one of them" << std::endl;
            exit(2);
        }
        bool taxon_statistics = result.count("taxon_statistics");
        bool full_sam_header = result.count("full_sam_header");
        bool write_unmapped_reads = result.count("write_unmapped_reads");
        bool serial_index_passes = result.count("serial_index_passes");
        bool profile_ahead = result.count("profile_ahead");
        bool no_unknown_share = result.count("no_unknown_share");
        bool ncbi_taxids = result.count("ncbi_taxids");
        size_t index_batch_kb = std::max<size_t>(result["index_batch_kb"].as<size_t>(), 1);
        auto const flex_indels = ParseFlexIndels(result["flex_indels"].as<std::string>());
        if (!flex_indels) {
            std::cerr << "Error: --flex_indels is none, all, or read types comma-separated (" << ReadTypeTokens() << "), not '"
                      << result["flex_indels"].as<std::string>() << "'" << std::endl;
            exit(2);
        }
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
                    // A sample of several files (--one_sample): the beginning their names share, in one folder.
                    auto const files = SplitReadFiles(first_list[i]), seconds = SplitReadFiles(second_list[i]);
                    if (files.size() > 1) {
                        std::vector<std::string> names;
                        std::set<std::filesystem::path> folders;
                        for (auto const& list : { files, seconds }) {
                            for (auto const& file : list) {
                                names.push_back(Utils::ReadFileStem(std::filesystem::path(file).filename().string()));
                                folders.insert(std::filesystem::path(file).parent_path());
                            }
                        }
                        auto const prefix = LaneSamplePrefix(names);
                        if (folders.size() > 1 || prefix.empty()) {
                            std::cerr << "Cannot name the sample of " << first_list[i] << " after its read files ("
                                      << (folders.size() > 1 ? "they are in different folders" : "their names share no beginning")
                                      << "); name it with --prefix" << std::endl;
                            exit(32);
                        }
                        prefix_list.push_back(prefix);
                        continue;
                    }
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
        // --strain_spill (before 2026-10-12) is --scratch_dir now, which wins.
        d.scratch_dir              = result["scratch_dir"].as<std::string>();
        if (d.scratch_dir.empty()) d.scratch_dir = result["strain_spill"].as<std::string>();
        d.sams_named               = sams_named;
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
        d.no_indel_bound           = no_indel_bound;
        d.no_single_seed_anchors   = no_single_seed_anchors;
        d.compact_keys             = compact_keys;
        d.long_read_budget         = long_read_budget;
        d.sequential_load          = sequential_load;
        d.stage1_probe             = stage1_probe;
        d.sparse_seeding           = sparse_seeding;
        d.sparse_index             = sparse_index;
        d.sparse_first             = *sparse_first;
        d.sparse_first_score       = sparse_first_score;
        d.sparse_first_hits        = sparse_first_hits;
        d.sparse_first_pick        = *sparse_first_pick;
        d.taxon_statistics         = taxon_statistics;
        d.full_sam_header          = full_sam_header;
        d.write_unmapped_reads     = write_unmapped_reads;
        d.serial_index_passes      = serial_index_passes;
        d.profile_ahead  = profile_ahead;
        d.no_unknown_share         = no_unknown_share;
        d.ncbi_taxids              = ncbi_taxids;
        d.index_batch_kb           = index_batch_kb;
        d.flex_indels              = *flex_indels;
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
        d.adaptive_candidates      = given.adaptive_candidates.value_or(DEFAULT_ADAPTIVE_CANDIDATES);
        d.no_allele_scores         = !given.allele_scores.value_or(true);
        d.max_out                  = max_out;
        d.max_score_ani            = max_score_ani;
        d.max_score_ani_given      = given.max_score_ani.front().has_value();
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
        d.wfa_traceback            = wfa_traceback;
        d.fragment_count           = fragment_count;
        d.max_key_ubiquity         = max_key_ubiquity;
        d.min_successful_lookups   = min_successful_lookups;
        d.max_seed_size            = max_seed_size;
        d.fastalign                = fastalign;
        d.sam_list                 = std::move(sam_list);
        d.profile_list             = std::move(profile_list);
        d.profile_truth_list       = std::move(profile_truth_list);
        d.read_type_list           = std::move(read_type_list);
        d.unmapped_reads_list      = std::move(unmapped_reads_list);
        d.sparse_first_list        = std::move(sparse_first_list);
        d.sparse_first_pick_list   = std::move(sparse_first_pick_list);
        d.flex_indels_list         = std::move(flex_indels_list);
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
        d.depth_identity_margin    = given.depth_identity_margin.value_or(0.08);
        d.keep_suspect_copies      = keep_suspect_copies;
        d.suspect_copy_distance    = result["suspect_copy_distance"].as<double>();
        d.strain_alleles           = result["strain_alleles"].as<size_t>();
        d.allele_genome_share      = result["allele_genome_share"].as<double>();
        d.index_alleles            = result["index_alleles"].as<double>();
        // The column weights' options (declared with f566cd0 but never read back until 2026-10-10: --column_weights,
        // --no_site_weights and --weighted_site_shift kept their defaults whatever the command line said).
        d.column_weights           = result["column_weights"].as<size_t>();
        d.column_weights_alignment = result["column_weights_alignment"].as<std::string>();
        d.grid_choice              = result["grid_choice"].as<std::string>();
        if (d.grid_choice != "scored" && d.grid_choice != "even" && d.grid_choice != "phase") {
            std::cerr << "Error: --grid_choice is scored, even or phase, not '" << d.grid_choice << "'" << std::endl;
            exit(2);
        }
        d.no_ancestry_sites        = result.count("no_ancestry_sites") > 0;
        d.no_unique_cores          = result.count("no_unique_cores") > 0;
        d.no_site_weights          = result.count("no_site_weights") > 0;
        // Off by default (it was on 2026-10-11 to 2026-10-12); --no_weighted_site_shift wins.
        d.weighted_site_shift      = given.weighted_site_shift.value_or(false);
        d.run_settings_given       = given;
        d.msa_identity_margin      = result["msa_identity_margin"].as<double>();
        d.gene_conservation        = result["gene_conservation"].as<std::string>();
        d.model                    = result["model"].as<std::string>();
        d.model_se                 = result["model_se"].as<std::string>();
        d.model_pb                 = result["model_pb"].as<std::string>();
        d.model_ont                = result["model_ont"].as<std::string>();
        d.read_type                = result["read_type"].as<std::string>();
        d.add_model                = result["add_model"].as<std::string>();
        d.model_gate               = result["model_gate"].as<std::string>();
        d.model_kind               = result["model_kind"].as<std::string>();
        d.add_tables               = result["add_tables"].as<std::string>();
        d.write_species_neighbours = result["write_species_neighbours"].as<std::string>();

        auto options = Options(std::move(d));

        if (!options.PrepareAndCheckValidity()) {
            std::cerr << "Exit Program" << std::endl;
            exit(30);
        }

        return options;
    }

    void Options::FollowDatabaseRunSettings(build_settings::RunSettings const& recorded) {
        auto const& g = m_run_settings_given;
        if (!g.weighted_site_shift && recorded.weighted_site_shift) m_weighted_site_shift = *recorded.weighted_site_shift;
        if (!g.allele_scores && recorded.allele_scores) m_no_allele_scores = !*recorded.allele_scores;
        if (!g.align_top && recorded.align_top) m_align_top = *recorded.align_top;
        if (!g.adaptive_candidates && recorded.adaptive_candidates) m_adaptive_candidates = *recorded.adaptive_candidates;
        // A mode this protal does not know (a later protal's) leaves its default (ApplyRunSettings warns).
        if (!g.mate_guidance && recorded.mate_guidance && build_settings::KnownMateGuidance(*recorded.mate_guidance)) {
            m_no_mate_guidance = *recorded.mate_guidance == build_settings::kMateGuidanceOff;
        }
        m_database_max_score_ani = recorded.max_score_ani;  // GetMaxScoreAni: where -a is not given
        if (!g.keep_foreign_genes && recorded.keep_foreign_genes) m_keep_foreign_genes = *recorded.keep_foreign_genes;
        if (!g.keep_suspect_copies && recorded.keep_suspect_copies) m_keep_suspect_copies = *recorded.keep_suspect_copies;
        if (!g.depth_identity_margin && recorded.depth_identity_margin) m_depth_identity_margin = *recorded.depth_identity_margin;
        // unique_kmers: a run reads the database's unique_cores.bin where it has one, unless --no_unique_cores.
    }

    build_settings::RunSettings Options::RunSettingsInUse() const {
        build_settings::RunSettings s;
        s.weighted_site_shift = m_weighted_site_shift;
        s.allele_scores = !m_no_allele_scores;
        s.align_top = m_align_top;
        s.adaptive_candidates = m_adaptive_candidates;
        s.mate_guidance = std::string(m_no_mate_guidance ? build_settings::kMateGuidanceOff : build_settings::kMateGuidanceOn);
        for (auto const& info : kReadTypes) s.max_score_ani[static_cast<size_t>(info.type)] = GetMaxScoreAni(info.type);
        s.unique_kmers = std::string(m_no_unique_cores ? build_settings::kUniqueAnchors : build_settings::kUniqueCores);
        s.keep_foreign_genes = m_keep_foreign_genes;
        s.keep_suspect_copies = m_keep_suspect_copies;
        s.depth_identity_margin = m_depth_identity_margin;
        return s;
    }

    void Options::SetSampleGate(size_t index, size_t gate, std::optional<SparseFirstPick> pick) {
        size_t const samples = m_prefix_list.size();
        if (index >= samples) return;
        // --sparse_first (--sparse_first_pick) of its read type, or the map's SPARSE_FIRST (SPARSE_FIRST_PICK), wins.
        size_t const t = static_cast<size_t>(GetReadType(index));
        if (!m_sparse_first_given[t] && !(index < m_sparse_first_list.size() && m_sparse_first_list[index])) {
            if (m_sparse_first_list.size() < samples) m_sparse_first_list.resize(samples);
            m_sparse_first_list[index] = gate;
        }
        if (pick && !m_sparse_first_pick_given[t] && !(index < m_sparse_first_pick_list.size() && m_sparse_first_pick_list[index])) {
            if (m_sparse_first_pick_list.size() < samples) m_sparse_first_pick_list.resize(samples);
            m_sparse_first_pick_list[index] = pick;
        }
    }
}
