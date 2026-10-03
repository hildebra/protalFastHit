// protal 0.7.5. In 0.7 (since 0.6.0a):
// - Reads: single-end (se), PacBio (pb) and ONT (ont) reads besides paired-end, each profiled with
//   its own model (--read_type, or a map's READ_TYPE); long reads are aligned per gene, and reads
//   over 65 kb in chunks. Pipes, gzip, BGZF and FASTA input; unusable inputs fail their sample.
// - Speed: input inflated outside the reader lock (libdeflate, vendored zlib-ng), short reads aligned
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
#include <iostream>
#include "RunProtal.h"

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    return protal::Run(argc, argv);
}
