#!/usr/bin/env python3
"""Build and train a complete protal database from an extracted GTDB release.

Usage: python3 scripts/build_gtdb_database.py --gtdb GTDB_DIR --outdir OUT_DIR

The GTDB release must include taxonomy, marker gene FASTAs and extracted whole
genome FASTAs. A genome table can be supplied explicitly when the release uses
an unusual layout.

The presence model is trained on metagenomes simulated from these genomes. GTDB
distributes whole genomes of the representatives only (genomic_files_reps), and
with only those every simulated species is the database's own reference, closer
to it than real strains are. Add other genomes of the species (from NCBI, by the
accessions in GTDB's metadata) under genomic_files_all/gtdb_genomes_all_r<R> or
in --genome-table, so that species are simulated from other strains, too.
OUT_DIR/genome_table.txt says how often they are.

OUT_DIR/model_logs/ collects what tells whether the model is good: the training
report (how it does on species it was not trained on, against the previous model
and training procedure), its numbers as JSON, the per-taxon predictions, the
threshold table, the parity check with protal and the genome table summary.
"""
import argparse
import collections
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
PARITY = os.path.join(HERE, "check_model_parity.py")
ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")
sys.path.insert(0, os.path.join(HERE, "mini_db"))
from gtdb_to_protal_db import normalize_accession, read_representatives  # noqa: E402


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


def summarize_genome_table(path, reps):
    """What the simulations can draw: genomes and species by domain, and how often a simulated species is
    not the database's representative genome (the simulator picks a species, then one of its genomes)."""
    genomes = collections.Counter()
    domain = {}
    rep_genomes = 0
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lineage = next((f for f in fields if f.startswith("d__") and ";s__" in f), None)
            if lineage is None:
                continue
            species = lineage.split(";")[-1]
            genomes[species] += 1
            domain[species] = lineage.split(";")[0][3:]
            rep_genomes += reps is not None and normalize_accession(fields[0]) in reps
    lines = [f"genome table {path}: {sum(genomes.values())} genomes of {len(genomes)} species"]
    for d in sorted(set(domain.values())):
        sp = [s for s in genomes if domain[s] == d]
        lines.append(f"  {d}: {len(sp)} species, {sum(genomes[s] for s in sp)} genomes")
    several = sum(1 for n in genomes.values() if n > 1)
    other_strain = sum((n - 1) / n for n in genomes.values()) / max(1, len(genomes))
    lines.append(f"  {several} species have more than one genome; a simulated species is another genome than "
                 f"its representative {100 * other_strain:.1f}% of the time")
    if reps is not None:
        lines.append(f"  {rep_genomes} of the genomes are species representatives (the database's references)")
    if other_strain < 0.2:
        lines.append("  WARNING: nearly all simulated species will be the database's own reference genome, closer to it "
                     "than real strains are. Add non-representative genomes (NCBI accessions from GTDB's metadata) "
                     "to train on real strain divergence.")
    return lines


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
    p.add_argument("--ntree", type=int, default=64)
    p.add_argument("--maxnodes", type=int, default=128)
    p.add_argument("--evaluation", choices=["full", "basic", "none"], default="full",
                   help="how much the trainer evaluates (random_forest_cmdline.py --evaluation)")
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
    logs = os.path.join(args.outdir, "model_logs")
    os.makedirs(logs, exist_ok=True)
    summary = summarize_genome_table(genome_table, read_representatives(args.gtdb, release))
    with open(os.path.join(args.outdir, "genome_table.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print("\n".join(summary), flush=True)

    run([sys.executable, CONVERTER, "--gtdb", args.gtdb, "--outdir", db, "--release", release],
        os.path.join(args.outdir, "convert.log"))
    # --build packs the taxonomy into database.protal; the trainer reads the taxa's domains from it.
    taxonomy = os.path.join(args.outdir, "internal_taxonomy.dmp")
    shutil.copyfile(os.path.join(db, "internal_taxonomy.dmp"), taxonomy)
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
         "--maxnodes", str(args.maxnodes), "--seed", str(args.seed), "--threads", str(args.threads),
         "--taxonomy", taxonomy, "--evaluation", args.evaluation],
        os.path.join(args.outdir, "classifier_training.log"))
    # protal must score as the trainer does, and compute the features as it did during collection.
    run([sys.executable, PARITY, "--db", db, "--model", prefix + ".xml", "--training", training,
         "--protal", args.protal, "-t", str(args.threads)], os.path.join(args.outdir, "parity.log"))
    for name in (prefix + ".report.txt", prefix + ".metrics.json", prefix + ".thresholds.tsv", prefix + ".varimp.tsv",
                 prefix + ".predictions.tsv.gz", os.path.join(training, "parity", "parity.txt"),
                 os.path.join(args.outdir, "training_data.log"), os.path.join(args.outdir, "genome_table.txt")):
        if os.path.isfile(name):
            shutil.copy(name, logs)
    # The trained model replaces the shipped one as the paired-end model (model_pe.xml) in
    # database.protal; --add_model checks it and copies the other parts as they are.
    run([args.protal, "--add_model", prefix + ".xml", "--read_type", "pe", "--db", db, "-t", str(args.threads)],
        os.path.join(args.outdir, "final_package.log"))
    with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
        fh.write(f"gtdb_release\tr{release}\nclassifier_features\tnormalized\nclassifier_trees\t{args.ntree}\n"
                 f"classifier_max_leaves\t{args.maxnodes}\n")
    shutil.copy(os.path.join(db, "build_metadata.tsv"), logs)
    print(f"Ready protal database: {db}", flush=True)
    print(f"Model evaluation: {logs} (start with trained_model.report.txt)", flush=True)


if __name__ == "__main__":
    main()
