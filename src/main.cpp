// protal 0.8.0. In 0.7 (since 0.6.0a):
// - Reads: single-end (se), PacBio (pb) and ONT (ont) reads besides paired-end, each profiled with
//   its own model (--read_type, or a map's READ_TYPE); long reads are aligned per gene, and reads
//   over 65 kb in chunks. Pipes, gzip, BGZF and FASTA input; unusable inputs fail their sample.
// - Speed: input inflated outside the reader lock (ISA-L since 0.7.8, before libdeflate and zlib-ng), short reads aligned
//   from their anchors' exact matches, branch-free and AVX2 syncmer scans, SAMs compressed in the
//   alignment threads (.sam.zst by default), the database loaded and the index built in parallel.
// - Database: one file, database.protal (seekable zstd, index in columns), checked against its
//   reference; --add_model for a read type's model; index format 2 with correct unique flags.
// - Profiles and strains: abundance from a species' own reads; variants called per fragment with
//   strand tests, deletions and unbiased depth; MSAs chosen by --msa_knob.
// - Training: build_gtdb_database.py builds and trains a database from a GTDB release (held-out
//   clades, other strains, an independent test set, one model per read type) and resumes reruns.
// - Failures are reported and give a non-zero exit; outputs are written crash-safe.
// In 0.7.1: mates paired across neighbouring genes (gene_neighbours.tsv); the reads' other candidates
//   (MAPQ, congener fits) as model features; depth margin 0.08; genes' within-species conservation
//   stored (--gene_conservation db); syncmers from 2-bit codes, genes decoded in windows.
// In 0.7.2: model features of the reads' divergence beyond their base qualities and of the conservation of the
//   genes a taxon's reads hit, before and after the MAPQ filter (gene_conservation.tsv, now read by every query);
//   knobs by sample depth in a model (--depth-knobs, for the long-read models); one sample profiled on all threads,
//   SAM records parsed as views, a .sam.zst decompressed on threads of its own; --build reports how the genes differ
//   between congeners (gene_congeners.tsv); long training reads drawn by the collector.
// In 0.7.3: one binary for every CPU (hot functions also as x86-64-v3, chosen at run time); gene neighbours from
//   every genome, in database.protal, with lines of a species' own gene order; foreign genes left out of depth and
//   MSA (--keep_foreign_genes); strain rows from long reads (--no_phasing); flanks of one mismatch aligned without
//   WFA2; a single read file's read type from its first reads; depth counts the bases where a pair's mates overlap
//   once; knobs as a curve over sample depth for every read type; --add_model with several models; the simulator's
//   samples on threads, compressed as written; PacBio HiFi training reads; the GTDB build's deeper design points and
//   faster conversion; --version gives the commit built from, which build_gtdb_database.py checks at its start.
// In 0.7.4: species calls judged against the sample's congeners (relatives and distance features, distance by
//   default; a read EM over the reads' alternatives; the singleton rule; opt-in calls at a target share of false
//   calls); long reads aligned through every link of their chain; k-mer lookups prefetched; --profile_ahead
//   (opt-in); faster read EM and congener sketches with the same outputs; the GTDB build's simulations in one
//   queue (long reads largest first, deep samples in chunks, beside the paired-end points; a simulated sample's
//   genomes on threads), the training database on --scratch, and both collections profiled in one protal run.
// In 0.7.5: against false positives: the sample's depth (sample_log_fragments) and the reads' divergence beyond their
//   qualities as model features, in the default feature set (the depth replaces the knob curve); gene copies
//   near-identical to another genus's found at build (suspect_copies.tsv) and left out of the evidence; species
//   priors from GTDB (species_priors.tsv); a taxon's evidence before the MAPQ filter, the taxa a read seeded on
//   but did not align to (ZF tag), prevalence across a run's samples; a mixed training design. The index's values
//   held packed in query runs (42-bit entries at r226, 35 -> 27 GB; the file unchanged). Profiling and alignment
//   cheaper by the same outputs (coverage as a difference array, 2-3 mismatch flanks without WFA2).
// In 0.7.6: at GTDB scale (r226, real samples, 32 threads) a paired-end run 131 -> 58 s and a HiFi run 104 -> 25 s:
//   a k-mer screen refuses hopeless candidates before WFA2 (exact; --no_alignment_screen), a long read's window ends
//   at its chain's last link, failed candidates counted in vectors by taxid, the gene tables loaded in parallel and
//   from a binary member (gene_table.bin), the start-up's loads side by side (--sequential_load), the per-taxon
//   statistics files only with --taxon_statistics; alignment counts and every stage timed per run. Outputs no longer
//   depend on the SAM's record order (exact sums, ties broken by taxid). Training: species priors opt-in, reduced
//   marker sets through every build phase (genes ranked per domain), in-silico strains of one-genome species, denser
//   training depths, real gene speeds for simulated worlds.
// In 0.7.7: at r226 a paired-end run 58 -> 42 s and a HiFi run 25 -> 15.6 s, the peak memory 38.1 -> 34.2 GB: the
//   flex scan with AVX2; the reads that aligned nowhere counted per taxon in the SAM header instead of an unmapped
//   record each (--write_unmapped_reads writes them); the profiler's SAM read, the preload and the strain stage
//   (species over the threads) on all threads; seeds sorted as 128-bit keys, and a seed tied with the last one taken
//   kept for its own diagonal's anchor; the read EM capped at 100 sweeps (200 before); the index packed key by key as
//   it loads (no 8-byte chunk copies), in fractional-bit slots, its load's tail throttled to the run's own peak;
//   memory logged per stage. .fq.zst read input; the build's index packed during --build, its simulations in zstd;
//   training scenarios of real studies.
// In 0.7.8: at r226 a paired-end run 42 -> 33 s and a HiFi run 15.6 -> 12.3 s, the same outputs: the seed lookup's
//   best cells taken by AVX2 masks, flat taxid -> genome and -> gene tables with the anchors' genes prefetched, the
//   k-mer screen from the genes' packed bytes and the read's strands packed once, only the seeds FindPairs can pair
//   sorted; gzip read and written with ISA-L (replaces libdeflate and zlib-ng's inflate; paired-end runs no longer
//   wait for their gzip input; .gz written ~7x faster, ~16% larger), `just static` builds ISA-L and nasm itself;
//   seeding counters in the log. --add_model replaces the models in place (seconds instead of a whole rewrite).
//   Training: gradient-boosted models by default (--model gbm, 250 rounds, folds fitted side by side), the
//   reference's k-mer uniqueness as a feature, scenario samples at depths of their own, the default feature set
//   and --evaluation basic as the build's defaults.
// In 0.7.9: against the false positives of complex communities, where most are a species the database lacks
//   landing on its nearest congener at a strain's identity: the reads' consistency, the shape of their genes and the
//   database's neighbourhood as default features; the sample's complexity; and the ancestry sites (AncestrySites.h):
//   where a species' gene copy differs from its nearest congener's, whether each read carries the species' base or
//   the congener's (a strain the former, a congener that branched off below the latter), from the gene store and
//   species_neighbours.tsv at run time. Strains: a cohort's MSA rows packed (159 -> 37 MB per dense sample), strain
//   evidence spilled to disk (--strain_spill), MSAs over several runs (--profile_only patterns, merging runs). Builds:
//   the reads behind every model error kept (--error-reads, --share-logs: their SAMs, the ancestry report and one
//   archive), the in-silico strains' dN/dS calibrated on real strains, scenarios of varied richness and depth with the
//   knob chosen by bootstrap; simulate_metagenomes makes Illumina, Ultima and long reads itself (no ART or pbsim3)
//   and streams large samples into protal; resume, suspect copies and the index straight into database.protal in
//   --build; every test suite in CI.
// In 0.8.0: strains told from congeners: --build keeps the species' other genomes as alleles of each marker gene
//   (strain_alleles.tsv), which score a short read's candidates and settle the reads protal is unsure about by the
//   species' polymorphic sites (--no_allele_scores), and the gaps to the congeners' copies (congener_gaps.tsv); a
//   divergent read tries its genus's untried candidates (--adaptive_candidates, ZC); the ancestry sites are the
//   congeners' consensus, chained across indels, with indel sites; 85 default features (gaps, untried, alleles,
//   polymorphic), untested at r226. The profile ends with the share of the genomes its species do not explain ("?",
//   from the genome sizes in species_priors.tsv; <profile>.composition; --no_unknown_share). Speed: ONT candidates
//   refused by an exact indel-distance bound before WFA2 (-45% of the bench's ONT alignment instructions;
//   --no_indel_bound), the seed lookup and syncmer scan in AVX-512 (PROTAL_SIMD; seeding -10% on Zen 4), BGZF input
//   inflated on 1-4 threads, BGZF output the same bytes on any thread; misc/cpu.tsv per stage. Builds: the foreign
//   scan's leak removed (the scan opt-in, --foreign-rates), species complexes held out whole (species_clouds.tsv),
//   --outdir sorted into protal_db/, model_logs/, logs/, work/ (~118 files, not ~11,900), the build on many cores
//   (setup as a graph, the derive on all cores, --profile-ahead, CPU logs), the composition checked against the
//   truth, in-silico strains at congener sites. Fixed: a map run's relative output folder was used twice.
#include <iostream>
#include "RunProtal.h"

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    return protal::Run(argc, argv);
}
