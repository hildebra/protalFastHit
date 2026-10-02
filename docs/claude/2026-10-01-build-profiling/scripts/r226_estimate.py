#!/usr/bin/env python3
"""The collection's compute at GTDB r226 with build_gtdb_database.py's defaults, from the unit costs measured
on real-sized genomes (README.md, "Unit costs"): CPU seconds per stage, now and with the changes proposed.

The design: the training data (3 read setups x 5 depths x 12 samples, 20-200 species, strains 0.3,0.1; 5 x 12
PacBio and Nanopore samples) and the test set (3 x 6 x 4, 10-300 species, strains 0.5,0.2; 5 x 4 long-read
samples of each type). Genomes per sample as measured (170 training, 240 test); the genome table and contigs
per genome are parameters (download_gtdb.py's defaults give at most 6,000 x 3 + 2,000 = 20,000 genomes of
8,000 species; its genomes.tsv has every downloaded genome's contig count).

usage: r226_estimate.py [--genomes 20000] [--reps 8000] [--contigs 78] [--gtdb-factor 3]
"""

import argparse

# Unit costs, CPU seconds on the laptop of the measurements (README.md).
SCAN_PER_GENOME = 0.0209          # simulate_metagenomes reads a genome for its length, every run
ART_PER_CALL = 0.105              # art_illumina's fixed cost for a 3.6 Mb genome
COPY_PER_CALL = 0.020             # the simulator's plain copy of the genome for ART
ART_PER_MB = 32 / 300             # 32 s per million 150 bp pairs
PB_PER_CONTIG, PB_PER_MB = 0.008, 0.45
ONT_PER_CONTIG, ONT_PER_MB = 0.045, 0.44
LONG_CONTIGS_PER_GENOME = 0.090   # the collector's Python copy of each genome for pbsim3 (holds the GIL)
TEMPLATES_PER_GENOME = 0.030      # reading a genome for templates, bytes operations (estimated)
TEMPL_PER_MB = 0.78               # pbsim3 --strategy templ, PacBio and Nanopore alike (measured, c2)
LENGTHS_PER_GENOME = 0.027        # build_gtdb_database.genome_length, once per build (measured)
NEIGHBOURS_PER_REP = 1.45         # gene_neighbours.py per representative genome
# protal on the tuning world's training database (4 threads; thread-seconds per unit):
PE_PER_MPAIRS = 4 / 0.15          # 150k pairs/s on 4 threads
SE_PER_MREADS = 4 / 0.38
PB_ALIGN_PER_MB = 4 / 11.0
ONT_ALIGN_PER_MB = 4 / 7.0


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--genomes", type=int, default=20000, help="genomes in the genome table (default 20000)")
    p.add_argument("--reps", type=int, default=8000, help="representatives among them (default 8000)")
    p.add_argument("--contigs", type=float, default=78, help="mean contigs per genome (default 78, as measured)")
    p.add_argument("--gtdb-factor", type=float, default=3,
                   help="protal's alignment cost per read at GTDB size over the tuning world's (default 3; unknown)")
    opts = p.parse_args()

    collections = {
        "training": dict(points=15, samples=12, genomes=170, pairs=[1000, 5000, 20000, 100000, 500000] * 3,
                         lengths=[100] * 5 + [150] * 5 + [250] * 5,
                         long_bases=[3e5, 1.5e6, 6e6, 3e7, 1.5e8]),
        "test": dict(points=18, samples=4, genomes=240, pairs=[500, 2000, 10000, 50000, 200000, 1000000] * 3,
                     lengths=[100] * 6 + [150] * 6 + [250] * 6,
                     long_bases=[1.5e5, 1e6, 5e6, 2.5e7, 2.5e8]),
    }
    rows = []
    for name, c in collections.items():
        pe_samples = c["points"] * c["samples"]
        calls = pe_samples * c["genomes"]
        mb = sum(p * 2 * l for p, l in zip(c["pairs"], c["lengths"])) * c["samples"] / 1e6
        long_samples = len(c["long_bases"]) * c["samples"]
        long_calls = long_samples * c["genomes"]
        long_mb = sum(c["long_bases"]) * c["samples"] / 1e6
        pairs = sum(c["pairs"]) * c["samples"] / 1e6
        rows += [
            (name, "simulator: genome lengths (proposed: once per build, in the table)",
             c["points"] * opts.genomes * SCAN_PER_GENOME, opts.genomes * LENGTHS_PER_GENOME / 2),
            (name, "ART and its genome copies", calls * (ART_PER_CALL + COPY_PER_CALL) + mb * ART_PER_MB,
             calls * (ART_PER_CALL + COPY_PER_CALL) + mb * ART_PER_MB),
            (name, "pbsim3 PacBio", long_calls * opts.contigs * PB_PER_CONTIG + long_mb * PB_PER_MB, long_mb * TEMPL_PER_MB),
            (name, "pbsim3 Nanopore", long_calls * opts.contigs * ONT_PER_CONTIG + long_mb * ONT_PER_MB, long_mb * TEMPL_PER_MB),
            (name, "collector: genome copies for pbsim3 (GIL)", 2 * long_calls * LONG_CONTIGS_PER_GENOME,
             2 * long_calls * TEMPLATES_PER_GENOME),
            (name, f"protal alignment (x{opts.gtdb_factor:g} at GTDB size)",
             opts.gtdb_factor * (pairs * (PE_PER_MPAIRS + SE_PER_MREADS) + long_mb * (PB_ALIGN_PER_MB + ONT_ALIGN_PER_MB)),
             opts.gtdb_factor * (pairs * (PE_PER_MPAIRS + SE_PER_MREADS) + long_mb * (PB_ALIGN_PER_MB + ONT_ALIGN_PER_MB))),
        ]
    rows.append(("before the builds", "gene_neighbours.py", opts.reps * NEIGHBOURS_PER_REP, opts.reps * NEIGHBOURS_PER_REP / 10))
    print("| stage | part | now, CPU h | proposed, CPU h |")
    print("|---|---|---:|---:|")
    total_now = total_new = 0
    for stage, part, now, new in rows:
        total_now += now
        total_new += new
        print(f"| {stage} | {part} | {now / 3600:.1f} | {new / 3600:.1f} |")
    print(f"| all | | **{total_now / 3600:.1f}** | **{total_new / 3600:.1f}** |")


if __name__ == "__main__":
    main()
