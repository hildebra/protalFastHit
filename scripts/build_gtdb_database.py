#!/usr/bin/env python3
"""Build and train a complete protal database from an extracted GTDB release.

Usage: python3 scripts/build_gtdb_database.py --gtdb GTDB_DIR --outdir OUT_DIR

The GTDB release must include taxonomy, marker gene FASTAs and extracted whole
genome FASTAs. A genome table can be supplied explicitly when the release uses
an unusual layout.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONVERTER = os.path.join(HERE, "mini_db", "gtdb_to_protal_db.py")
TRAINER = os.path.join(HERE, "random_forest_cmdline.py")
COLLECTOR = os.path.join(HERE, "collect_training_data.py")
ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")


def run(command, log):
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "w") as fh:
        rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc:
        sys.exit(f"Command failed ({rc}); see {log}: {' '.join(command)}")


def make_genome_table(gtdb, release, output):
    taxonomy = {}
    for domain in ("bac120", "ar53"):
        for suffix in (".tsv", ".tsv.gz"):
            path = os.path.join(gtdb, f"{domain}_taxonomy_r{release}{suffix}")
            if os.path.isfile(path):
                import gzip
                opener = gzip.open if path.endswith(".gz") else open
                with opener(path, "rt") as fh:
                    for line in fh:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) >= 2:
                            match = ACCESSION.search(fields[0])
                            if match:
                                taxonomy[match.group(1)] = fields[1]
                break
    genome_dirs = [os.path.join(gtdb, "genomic_files_all", f"gtdb_genomes_all_r{release}"),
                   os.path.join(gtdb, "genomic_files_reps", f"gtdb_genomes_reps_r{release}")]
    paths = {}
    for root in genome_dirs:
        if not os.path.isdir(root):
            continue
        for pattern in ("**/*.fna", "**/*.fna.gz", "**/*.fa", "**/*.fa.gz", "**/*.fasta", "**/*.fasta.gz"):
            for path in glob.iglob(os.path.join(root, pattern), recursive=True):
                m = ACCESSION.search(os.path.basename(path))
                if m and m.group(1) in taxonomy:
                    paths.setdefault(m.group(1), path)
    if not paths:
        sys.exit("No extracted whole genome FASTAs found under genomic_files_all/gtdb_genomes_all_r" + release +
                 " or genomic_files_reps/gtdb_genomes_reps_r" + release + "; extract GTDB genome files or pass --genome-table")
    with open(output, "w") as fh:
        for accession in sorted(paths):
            fh.write(f"{accession}\t{taxonomy[accession]}\t{os.path.abspath(paths[accession])}\n")
    return len(paths)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtdb", required=True, help="extracted GTDB release directory")
    p.add_argument("--outdir", required=True, help="output root; database is written to OUTDIR/protal_db")
    p.add_argument("--release", help="GTDB release number; detected from taxonomy filenames by default")
    p.add_argument("--genome-table", help="optional simulator table: accession, taxonomy, whole genome FASTA path")
    p.add_argument("--protal", default="protal", help="protal executable")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes executable")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--samples", type=int, default=4)
    p.add_argument("--read-pairs", default="5000,20000,100000,500000")
    p.add_argument("--read-setups", default="100:HS20:300:40,150:HS25:350:50,250:MSv3:550:50")
    p.add_argument("--archaea", type=int, default=2)
    p.add_argument("--species-per-sample", default="10-40")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--ntree", type=int, default=512)
    p.add_argument("--maxnodes", type=int, default=128)
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    db = os.path.join(args.outdir, "protal_db")
    os.makedirs(db, exist_ok=True)
    if args.release:
        release = str(args.release).removeprefix("r")
    else:
        found = set()
        for name in os.listdir(args.gtdb):
            m = re.match(r"(?:bac120|ar53)_taxonomy_r(\d+)", name)
            if m:
                found.add(m.group(1))
        if len(found) != 1:
            sys.exit(f"Cannot detect one GTDB release in {args.gtdb}; specify --release")
        release = next(iter(found))
    genome_table = args.genome_table or os.path.join(args.outdir, "genomes.tsv")
    if not args.genome_table:
        count = make_genome_table(args.gtdb, release, genome_table)
        print(f"Found {count} genome FASTAs for training", flush=True)

    run([sys.executable, CONVERTER, "--gtdb", args.gtdb, "--outdir", db, "--release", release],
        os.path.join(args.outdir, "convert.log"))
    build_cmd = [args.protal, "--build", "--no_profile", "-t", str(args.threads), "--db", db,
                 "--reference", os.path.join(db, "reference.fna")]
    if os.path.isfile(os.path.join(db, "full_reference.fna")):
        build_cmd += ["--full_reference", os.path.join(db, "full_reference.fna")]
    run(build_cmd,
        os.path.join(args.outdir, "index_and_package.log"))

    training = os.path.join(args.outdir, "training")
    run([sys.executable, COLLECTOR, "--db", db, "--genome_table", genome_table, "-o", training,
         "--protal", args.protal, "--simulator", args.simulator, "--samples", str(args.samples),
         "--read_pairs", args.read_pairs, "--read_setups", args.read_setups,
         "--archaea", str(args.archaea), "--species_per_sample", args.species_per_sample,
         "--seed", str(args.seed), "-t", str(args.threads)], os.path.join(args.outdir, "training_data.log"))
    prefix = os.path.join(args.outdir, "trained_model")
    run([sys.executable, TRAINER, "--truth-file", os.path.join(training, "training_data.tsv"),
         "--output-prefix", prefix, "--features", "normalized", "--ntree", str(args.ntree),
         "--maxnodes", str(args.maxnodes), "--seed", str(args.seed), "--threads", str(args.threads)],
        os.path.join(args.outdir, "classifier_training.log"))
    # --build packs/removes the component files, so unpack before replacing model.xml.
    run([args.protal, "--unpack_db", "--db", db, "-t", str(args.threads)],
        os.path.join(args.outdir, "unpack.log"))
    shutil.copyfile(prefix + ".xml", os.path.join(db, "model.xml"))
    run([args.protal, "--compress_db", "--db", db, "-t", str(args.threads)],
        os.path.join(args.outdir, "final_package.log"))
    with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
        fh.write(f"gtdb_release\tr{release}\nclassifier_features\tnormalized\n")
    print(f"Ready protal database: {db}", flush=True)


if __name__ == "__main__":
    main()
