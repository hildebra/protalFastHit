// protal 0.7.1. In 0.7 (since 0.6.0a):
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
#include <iostream>
#include "RunProtal.h"

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    return protal::Run(argc, argv);
}
