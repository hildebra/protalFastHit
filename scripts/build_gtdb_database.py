#!/usr/bin/env python3
"""Build and train a complete protal database from an extracted GTDB release.

Usage: python3 scripts/build_gtdb_database.py --inputs INPUTS_DIR --outdir OUT_DIR
       python3 scripts/build_gtdb_database.py --gtdb GTDB_DIR --outdir OUT_DIR

--inputs is a folder of scripts/download_gtdb.py (run on a node with internet): the
release, the genomes to simulate from and the species pool, in one place that
serves every build of that release. It stands for --gtdb INPUTS/release
--extra-genomes INPUTS/genomes --simulate-species INPUTS/simulation_species.txt.
With --gtdb, the release must include taxonomy, marker gene FASTAs and extracted
whole genome FASTAs; a genome table can be supplied explicitly when the release
uses an unusual layout.

The presence model is trained on metagenomes simulated from these genomes, made
harder in two ways than simulating the database's own references:
- Strains. GTDB distributes whole genomes of the representatives only
  (genomic_files_reps), and with only those every simulated species is the
  database's own reference, closer to it than real strains are.
  scripts/download_gtdb.py picks other genomes of the species from GTDB's
  metadata and downloads them from NCBI, so that species are simulated from
  other strains, too. OUT_DIR/genome_table.txt says how often.
- Species the database lacks. The samples are profiled against a training
  database (OUT_DIR/training_db) that leaves whole clades of every rank
  (--holdout-clades) and --holdout of the other species out
  (OUT_DIR/heldout_species.txt): their reads land on relatives, as those of
  organisms GTDB lacks do in real samples, and the model learns to reject those
  relatives. The report gives false positive and false negative rates by rank.
  The finished database has all species and the model trained so. The training
  database costs a second index build (while the first one runs) and its disk
  space.

One model per read type (--read-types, default pe,se,pb,ont), trained in
parallel: paired-end reads (ART), their first reads alone (single-end), and
PacBio and Nanopore reads of the same communities (pbsim3). Besides the training
data, an independent test set of another design (--test-*: other depths,
community sizes, abundances and strain mixes) is profiled and scored by each
model: cross-validation on the training data cannot show what its design lacks.

OUT_DIR/model_logs/ collects what tells whether the models are good: summary.txt
(TP, FP, TN, FN, sensitivity, specificity, precision and F1 of each model), each read
type's training report (how it does on species and clades it was not trained on,
on the independent test set, false positive and false negative rates by rank,
against the previous model and training procedure), its numbers as JSON, the
per-taxon predictions, the threshold table, the parity check with protal, the
genome table summary and build_metadata.tsv (what the database was built from
and with).
"""
import argparse
import collections
import glob
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CONVERTER = os.path.join(HERE, "mini_db", "gtdb_to_protal_db.py")
TRAINER = os.path.join(HERE, "random_forest_cmdline.py")
COLLECTOR = os.path.join(HERE, "collect_training_data.py")
PARITY = os.path.join(HERE, "check_model_parity.py")
ACCESSION = re.compile(r"(?:RS_|GB_)?(GC[AF]_\d{9}\.\d+)")
sys.path.insert(0, os.path.join(HERE, "mini_db"))
sys.path.insert(0, HERE)
import lineages  # noqa: E402
from gtdb_to_protal_db import normalize_accession, read_representatives  # noqa: E402
from model_pmml import MODEL_FILES, write_placeholder  # noqa: E402
from collect_training_data import TABLES  # noqa: E402


def run(command, log):
    os.makedirs(os.path.dirname(log), exist_ok=True)
    started = time.time()
    with open(log, "w") as fh:
        rc = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc:
        stop(f"Command failed ({rc}); see {log}: {' '.join(command)}")
    return time.time() - started


def make_genome_table(gtdb, release, output, extra_dirs=(), species=None):
    """The simulator's genome table: every genome FASTA of a GTDB species found in the release (and in
    extra_dirs), or only those of the species in the set `species`."""
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
                            if match and (species is None or fields[1].split(";")[-1] in species):
                                taxonomy[match.group(1)] = fields[1]
                break
    genome_dirs = [os.path.join(gtdb, "genomic_files_all", f"gtdb_genomes_all_r{release}"),
                   os.path.join(gtdb, "genomic_files_reps", f"gtdb_genomes_reps_r{release}"), *extra_dirs]
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
                     "than real strains are. Add non-representative genomes (scripts/download_gtdb.py, then "
                     "--inputs) to train on real strain divergence.")
    return lines


CLADE_RANKS = ("phylum", "class", "order", "family", "genus")


def parse_clades(text):
    """--holdout-clades "phylum:2,class:4,family:8" -> {"phylum": 2, ...}; "none" or "" -> {}."""
    clades = {}
    for part in (text or "").split(","):
        part = part.strip()
        if not part or part == "none":
            continue
        rank, _, count = part.partition(":")
        if rank not in CLADE_RANKS or not count.isdigit():
            sys.exit(f"--holdout-clades: expected RANK:COUNT with a rank of {', '.join(CLADE_RANKS)}, got {part!r}")
        clades[rank] = int(count)
    return clades


def choose_holdout(genome_table, taxonomy, fraction, clades, max_share, seed):
    """The species a training database leaves out, as {species: (rank, clade)}; rank "species" for single
    species. First whole clades, clades[rank] of each rank from phylum down: drawn among those with at least
    two species the genome table can simulate (so that samples can have them), with at most max_share of the
    database's species, and in no clade drawn before. Then a random `fraction` of the species the genome
    table can simulate that no clade took, the same fraction in each domain."""
    lineage_by_id, _ = lineages.from_taxonomy(taxonomy)
    db = {lin["species"]: lin for lin in lineage_by_id.values() if "species" in lin}
    pool_by_domain = collections.defaultdict(set)
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage and lineage.split(";")[-1] in db:
                pool_by_domain[lineage.split(";")[0]].add(lineage.split(";")[-1])
    pool = set().union(*pool_by_domain.values()) if pool_by_domain else set()
    rng = random.Random(seed)
    chosen = {}
    for rank in CLADE_RANKS:
        if not clades.get(rank):
            continue
        members = collections.defaultdict(set)
        for species, lin in db.items():
            if rank in lin:
                members[lin[rank]].add(species)
        eligible = sorted(c for c, species in members.items() if len(species & pool) >= 2
                          and len(species) <= max_share * len(db) and not any(s in chosen for s in species))
        picked = rng.sample(eligible, min(clades[rank], len(eligible)))
        if len(picked) < clades[rank]:
            print(f"WARNING: only {len(eligible)} {rank} clades can be held out (two or more species to simulate, at "
                  f"most {max_share:.0%} of the species); holding out {len(picked)} instead of {clades[rank]}", flush=True)
        for clade in sorted(picked):
            for species in members[clade]:
                chosen[species] = (rank, clade)
    for domain in sorted(pool_by_domain):
        species = sorted(s for s in pool_by_domain[domain] if s not in chosen)
        for s in rng.sample(species, round(fraction * len(species))):
            chosen[s] = ("species", s)
    return chosen


def read_holdout(path):
    """heldout_species.txt: species, and optionally the rank it was held out at and the clade (a species
    alone: rank "species") -> {species: (rank, clade)}."""
    chosen = {}
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if not fields[0] or fields[0].startswith("#"):
                continue
            species = fields[0] if fields[0].startswith("s__") else "s__" + fields[0]
            rank = fields[1] if len(fields) > 1 and fields[1] else "species"
            chosen[species] = (rank, fields[2] if len(fields) > 2 and fields[2] else species)
    return chosen


def pool_species(genome_table):
    """The species a genome table can simulate."""
    species = set()
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                species.add(lineage.split(";")[-1])
    return species


def describe_holdout(chosen, pool):
    """Lines saying what the training database leaves out."""
    lines = []
    by_rank = collections.defaultdict(lambda: collections.defaultdict(set))
    for species, (rank, clade) in chosen.items():
        by_rank[rank][clade].add(species)
    for rank in (*CLADE_RANKS, "species"):
        if rank not in by_rank:
            continue
        species = set().union(*by_rank[rank].values())
        simulated = len(species & pool)
        if rank == "species":
            lines.append(f"  {len(species)} single species ({simulated} to simulate)")
        else:
            names = ", ".join(f"{c} ({len(s)} species, {len(s & pool)} to simulate)" for c, s in sorted(by_rank[rank].items()))
            lines.append(f"  {len(by_rank[rank])} {rank} clades, {len(species)} species ({simulated} to simulate): {names}")
    return lines


def provenance(args, release, genome_table, heldout, n_heldout, read_types, prefixes):
    """build_metadata.tsv: what the database was built from and with, so that two builds can be compared."""
    def output(command):
        try:
            return subprocess.run(command, capture_output=True, text=True, timeout=60).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    version = output([args.protal, "--version"]).splitlines()
    commit = output(["git", "-C", HERE, "rev-parse", "HEAD"])
    if commit and output(["git", "-C", HERE, "status", "--porcelain", "--", "."]):
        commit += " (scripts changed since)"
    genomes, species = 0, set()
    with open(genome_table) as fh:
        for line in fh:
            lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
            if lineage:
                genomes += 1
                species.add(lineage.split(";")[-1])
    clade_counts = collections.Counter(rank for rank, _ in set(read_holdout(heldout).values())) if n_heldout else {}
    rows = [("gtdb_release", f"r{release}"), ("built", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("protal_version", version[-1] if version else "unknown"), ("protal_binary", args.protal),
            ("scripts_commit", commit or "unknown (not a git checkout)"), ("command", " ".join(sys.argv)),
            ("seed", args.seed), ("genome_table", f"{genomes} genomes of {len(species)} species"),
            ("classifier_features", "normalized"), ("classifier_trees", args.ntree),
            ("classifier_max_leaves", args.maxnodes), ("classifier_training_species_left_out", n_heldout)]
    rows += [(f"classifier_training_{rank}_clades_left_out", clade_counts[rank]) for rank in CLADE_RANKS
             if clade_counts.get(rank)]
    rows += [("classifier_training_samples", f"{args.samples} per design point"),
             ("classifier_training_design", f"read pairs {args.read_pairs}; read setups {args.read_setups}; "
                                            f"species per sample {args.species_per_sample}; strains "
                                            f"{args.strains_per_species or 'one'}; abundance {args.abundance or 'default'}"),
             ("classifier_read_types", ",".join(read_types))]
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            continue
        species_cv = metrics.get("evaluation", {}).get("species", {})
        test = metrics.get("test", {}).get("this one", {})
        rows.append((f"model_{t}", f"species held out F1 {species_cv.get('F1')}, FP per sample "
                                   f"{species_cv.get('FP_per_sample')}; independent test F1 {test.get('F1')}, FP per "
                                   f"sample {test.get('FP_per_sample')}"))
    return rows


def summary_lines(read_types, prefixes, db):
    """model_logs/summary.txt: for each read type's model, TP, FP, TN, FN and the rates, with species held
    out in training (cross-validation) and on the independent test set, from its .metrics.json."""
    header = ("read type", "evaluated on", "taxa", "TP", "FP", "TN", "FN", "sensitivity", "specificity",
              "precision", "F1", "FP/sample")
    rows = []
    for t in read_types:
        try:
            with open(prefixes[t] + ".metrics.json") as fh:
                metrics = json.load(fh)
        except (OSError, ValueError):
            rows.append((t, "no metrics (training failed?)") + ("",) * (len(header) - 2))
            continue
        for label, m in (("species held out", metrics.get("evaluation", {}).get("species")),
                         ("independent test set", metrics.get("test", {}).get("this one"))):
            if not m:
                continue
            tp, fn, fp = m["present"] - m["FN"], m["FN"], m["FP"]
            tn = m["taxa"] - m["present"] - fp
            rate = lambda v: "-" if v is None else f"{v:.4f}"
            rows.append((t, label, str(m["taxa"]), str(tp), str(fp), str(tn), str(fn), rate(m.get("sensitivity")),
                         rate(tn / (tn + fp) if tn + fp else None), rate(m.get("precision")), rate(m.get("F1")),
                         "-" if m.get("FP_per_sample") is None else f"{m['FP_per_sample']:.2f}"))
    widths = [max(len(str(r[i])) for r in [header, *rows]) for i in range(len(header))]
    table = ["  ".join(str(v).ljust(w) for v, w in zip(r, widths)).rstrip() for r in [header, *rows]]
    return [f"Presence models of {db}, taxa scored at knob 0.5: TP present and called, FP absent and called, TN "
            "absent and not called, FN present and not called. species held out: each taxon scored by forests "
            "that did not see its species; independent test set: samples of another design, scored by the "
            "final model. Details: trained_model*.report.txt", ""] + table


def build_command(protal, db, threads, *extra):
    command = [protal, "--build", "--no_profile", "-t", str(threads), "--db", db,
               "--reference", os.path.join(db, "reference.fna"), *extra]
    if os.path.isfile(os.path.join(db, "full_reference.fna")):
        command += ["--full_reference", os.path.join(db, "full_reference.fna")]
    return command


class Background:
    """A command run in the background (its output to log); finish() waits for it and stops the script if it
    failed. Stopped when the script stops before."""
    running = []

    def __init__(self, command, log):
        os.makedirs(os.path.dirname(log), exist_ok=True)
        self.command, self.log, self.started = command, log, time.time()
        self.fh = open(log, "w")
        self.process = subprocess.Popen(command, stdout=self.fh, stderr=subprocess.STDOUT)
        Background.running.append(self)

    def finish(self):
        rc = self.process.wait()
        self.fh.close()
        Background.running.remove(self)
        if rc:
            stop(f"Command failed ({rc}); see {self.log}: {' '.join(self.command)}")
        return time.time() - self.started


def stop(message):
    for job in list(Background.running):
        job.process.terminate()
        job.process.wait()
    sys.exit(message)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", help="folder of scripts/download_gtdb.py: the release, the genomes to simulate from "
                                    "and the species pool (instead of --gtdb, --extra-genomes, --simulate-species)")
    p.add_argument("--gtdb", help="extracted GTDB release directory")
    p.add_argument("--outdir", required=True, help="output root; database is written to OUTDIR/protal_db")
    p.add_argument("--release", help="GTDB release number; detected from taxonomy filenames by default")
    p.add_argument("--genome-table", help="optional simulator table: accession, taxonomy, whole genome FASTA path")
    p.add_argument("--extra-genomes", action="append", default=[],
                   help="folder of more whole genomes of GTDB species, found by the accession in their file names "
                        "(e.g. the NCBI genomes of download_gtdb.py); repeatable")
    p.add_argument("--holdout", type=float, default=0.2,
                   help="fraction of the species (of those the genome table can simulate and no clade of "
                        "--holdout-clades took) left out of a separate training database (default 0.2): their "
                        "reads land on relatives, as those of species GTDB lacks do in real samples, and they are "
                        "the false positives the model must learn to reject. 0 with --holdout-clades none trains "
                        "on the database itself")
    p.add_argument("--holdout-clades", default="phylum:2,class:4,order:6,family:8,genus:12",
                   help="whole clades left out of the training database too, RANK:COUNT for ranks phylum, class, "
                        "order, family, genus (default phylum:2,class:4,order:6,family:8,genus:12; none for "
                        "species only): reads of organisms whose genus, family, order, class or phylum the "
                        "database lacks land on ever more distant relatives. The report gives false positive and "
                        "false negative rates by rank")
    p.add_argument("--holdout-max-share", type=float, default=0.02,
                   help="a held-out clade has at most this share of the database's species (default 0.02)")
    p.add_argument("--novel-clades-per-sample", type=int, default=1,
                   help="species of held-out clades in each sample, per rank held out (default 1; "
                        "collect_training_data.py --novel_clades)")
    p.add_argument("--holdout-species", help="file of the species to leave out (optionally with their rank and "
                                             "clade, as heldout_species.txt), instead of choosing them")
    p.add_argument("--one-build-at-a-time", action="store_true",
                   help="build the finished database after the model is trained, not while the training data are "
                        "collected (the default, which needs the memory of two builds, or of one build and the "
                        "collection's protal runs, at once)")
    p.add_argument("--training-db-level", type=int, default=3,
                   help="zstd level of the training database (default 3; the finished database uses protal's default)")
    p.add_argument("--simulate-species",
                   help="file of the species to simulate from (e.g. simulation_species.txt of download_gtdb.py: "
                        "species with other strains, and some without): the simulator draws species uniformly, so "
                        "among all ~130,000 GTDB species the few with downloaded strains would hardly be drawn")
    p.add_argument("--protal", default="protal", help="protal executable")
    p.add_argument("--simulator", default="simulate_metagenomes", help="simulate_metagenomes executable")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--samples", type=int, default=12,
                   help="samples per design point (default 12; on a GTDB-like world the model still improved "
                        "from 60 to 120 samples)")
    p.add_argument("--read-pairs", default="1000,5000,20000,100000,500000",
                   help="read pairs per sample, one design point each (default 1000,5000,20000,100000,500000: "
                        "without the shallowest, a model missed 8%% of the present taxa of samples of 1000 read "
                        "pairs)")
    p.add_argument("--read-setups", default="100:HS20:300:40,150:HSXt:350:50,250:MSv3:550:50",
                   help="LENGTH:ART_PROFILE:FRAGMENT_MEAN:FRAGMENT_SD, one design point each (HSXt: HiSeq X, the "
                        "closest of ART's profiles to NovaSeq; file=R1.txt+R2.txt: profiles art_profiler_illumina "
                        "made from real reads)")
    p.add_argument("--archaea", type=int, default=2)
    p.add_argument("--species-per-sample", default="20-200",
                   help="species per sample, drawn per sample (default 20-200: real gut samples hold 100-300 GTDB "
                        "species with a long tail of rare ones)")
    p.add_argument("--strains-per-species", default="0.3,0.1",
                   help="probabilities of a second, third, ... strain of a species in a sample (default 0.3,0.1): "
                        "real samples often mix strains, which changes the allele-frequency features")
    p.add_argument("--abundance", default="",
                   help="abundance model of the training samples: lognormal:SIGMA, powerlaw:ALPHA or negbin:R:P "
                        "(default: the simulator's, Poisson-lognormal with sigma 1.3)")
    p.add_argument("--read-types", default="pe,se,pb,ont",
                   help="the read types to train a model for (default pe,se,pb,ont): se from the paired-end "
                        "samples' first reads, pb and ont from long reads of the same communities (pbsim3). The "
                        "models are trained in parallel; read types left out keep placeholders")
    p.add_argument("--long-read-bases", default="300000,1500000,6000000,30000000,150000000",
                   help="bases per long-read sample, one design point each (collect_training_data.py)")
    p.add_argument("--pb-setup", default="errhmm:ERRHMM-SEQUEL:15000:3000:0.999",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of PacBio reads")
    p.add_argument("--ont-setup", default="qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36",
                   help="pbsim3 METHOD:MODEL:LENGTH_MEAN:LENGTH_SD:ACCURACY_MEAN of Nanopore reads")
    p.add_argument("--pbsim", default="pbsim", help="pbsim3 executable, for pb and ont")
    p.add_argument("--pbsim-models", help="folder of pbsim3's .model files (default: found next to the executable)")
    p.add_argument("--test-samples", type=int, default=4,
                   help="samples per design point of the independent test set (default 4; 0: none). The test set "
                        "has another design than the training data (--test-*), so that the report shows what "
                        "cross-validation on the training data cannot")
    p.add_argument("--test-read-pairs", default="500,2000,10000,50000,200000,1000000")
    p.add_argument("--test-species-per-sample", default="10-300")
    p.add_argument("--test-abundance", default="lognormal:2.0", help="(default lognormal:2.0: more uneven than training)")
    p.add_argument("--test-strains-per-species", default="0.5,0.2")
    p.add_argument("--test-long-read-bases", default="150000,1000000,5000000,25000000,250000000")
    p.add_argument("--congeners", type=int, default=0,
                   help="species of one genus in every sample of a design point (collect_training_data.py "
                        "--congeners): relatives that share a sample, as in real samples")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--ntree", type=int, default=64)
    p.add_argument("--maxnodes", type=int, default=128)
    p.add_argument("--no-placeholder-models", action="store_true",
                   help="leave the se, pb and ont models out of the database (a run with such reads then stops "
                        "with an error) instead of placeholders that report no species until trained ones replace "
                        "them (placeholder_models.py)")
    p.add_argument("--evaluation", choices=["full", "basic", "none"], default="full",
                   help="how much the trainer evaluates (random_forest_cmdline.py --evaluation)")
    args = p.parse_args()
    read_types = [t.strip() for t in args.read_types.split(",") if t.strip()]
    if not read_types or any(t not in TABLES for t in read_types):
        p.error(f"--read-types: a comma-separated list of {', '.join(TABLES)}, got {args.read_types!r}")
    if any(t in ("pb", "ont") for t in read_types) and not (shutil.which(args.pbsim) or os.path.isfile(args.pbsim)):
        p.error(f"pb and ont reads are simulated with pbsim3, and {args.pbsim} is not there: install it (e.g. "
                "micromamba install -c conda-forge -c bioconda pbsim3), pass --pbsim, or leave pb and ont out of "
                "--read-types")
    if args.inputs:
        if args.gtdb:
            p.error("give --inputs or --gtdb, not both")
        state_file = os.path.join(args.inputs, "download.json")
        if not os.path.isfile(state_file):
            p.error(f"{args.inputs} has no download.json: is it a folder of download_gtdb.py?")
        with open(state_file) as fh:
            state = json.load(fh)
        args.gtdb = os.path.join(args.inputs, "release")
        args.release = args.release or state["release"]["number"]
        if os.path.isdir(os.path.join(args.inputs, "genomes")):
            args.extra_genomes.append(os.path.join(args.inputs, "genomes"))
        if not args.simulate_species and os.path.isfile(os.path.join(args.inputs, "simulation_species.txt")):
            args.simulate_species = os.path.join(args.inputs, "simulation_species.txt")
    elif not args.gtdb:
        p.error("give --inputs (a folder of download_gtdb.py) or --gtdb")
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
        pool = None
        if args.simulate_species:
            with open(args.simulate_species) as fh:
                pool = {n if n.startswith("s__") else "s__" + n
                        for n in (line.rstrip("\n").split("\t")[0].strip() for line in fh) if n and not n.startswith("#")}
        count = make_genome_table(args.gtdb, release, genome_table, args.extra_genomes, pool)
        print(f"Found {count} genome FASTAs for training", flush=True)
    elif args.extra_genomes or args.simulate_species:
        sys.exit("--extra-genomes and --simulate-species shape the genome table this script makes; "
                 "apply them to --genome-table instead")
    logs = os.path.join(args.outdir, "model_logs")
    os.makedirs(logs, exist_ok=True)
    summary = summarize_genome_table(genome_table, read_representatives(args.gtdb, release))
    with open(os.path.join(args.outdir, "genome_table.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print("\n".join(summary), flush=True)

    run([sys.executable, CONVERTER, "--gtdb", args.gtdb, "--outdir", db, "--release", release, "-t", str(args.threads)],
        os.path.join(args.outdir, "convert.log"))
    # --build packs the taxonomy into database.protal; the collector and the trainer read it (domains,
    # representative genomes). The training database has the same taxonomy.
    taxonomy = os.path.join(args.outdir, "internal_taxonomy.dmp")
    shutil.copyfile(os.path.join(db, "internal_taxonomy.dmp"), taxonomy)
    # One model per read type; only paired-end reads can be simulated and aligned here, so the others get
    # placeholders (protal warns when it loads one), packed by --build like model_pe.xml.
    if not args.no_placeholder_models:
        for read_type in ("se", "pb", "ont"):
            write_placeholder(os.path.join(db, MODEL_FILES[read_type]), read_type)

    # The training database leaves some species out: the model then sees reads of species the database
    # lacks, which land on relatives, and reads of whole families, classes and phyla it lacks, which land on
    # distant ones. It is made from the converted files before --build packs them. heldout_species.txt: the
    # species, the rank they were held out at and the clade.
    heldout = os.path.join(args.outdir, "heldout_species.txt")
    clades = parse_clades(args.holdout_clades)
    if args.holdout_species:
        shutil.copyfile(args.holdout_species, heldout)
    elif args.holdout > 0 or clades:
        chosen = choose_holdout(genome_table, taxonomy, args.holdout, clades, args.holdout_max_share, args.seed)
        with open(heldout, "w") as fh:
            fh.write("".join(f"{s}\t{rank}\t{clade}\n" for s, (rank, clade) in sorted(chosen.items())))
    elif os.path.exists(heldout):
        os.remove(heldout)
    training_db = db
    n_heldout, holdout_lines = 0, []
    if os.path.exists(heldout):
        chosen = read_holdout(heldout)
        n_heldout = len(chosen)
        holdout_lines = [f"training database: {n_heldout} species left out ({heldout})"] + \
            describe_holdout(chosen, pool_species(genome_table))
        with open(os.path.join(logs, "holdout.txt"), "w") as fh:
            fh.write("\n".join(holdout_lines) + "\n")
        print("\n".join(holdout_lines), flush=True)
        training_db = os.path.join(args.outdir, "training_db")
        run([sys.executable, CONVERTER, "--from_db", db, "--exclude_species", heldout, "--outdir", training_db],
            os.path.join(args.outdir, "training_db.log"))

    # The finished database is needed only for --add_model at the end: it is built in the background from the
    # start, while the training database is built and the training data collected, unless one build at a time.
    final_log = os.path.join(args.outdir, "index_and_package.log")
    final_build = None
    if training_db != db and not args.one_build_at_a_time:
        final_build = Background(build_command(args.protal, db, args.threads), final_log)
        print(f"Building {db} in the background ({final_log})", flush=True)
    elif training_db == db:
        print(f"Built {db} in {run(build_command(args.protal, db, args.threads), final_log):.0f} s", flush=True)
    if training_db != db:
        # Read only for the training samples and the parity check: zstd level 3 packs it in a fraction of the
        # time of level 19 (which half of a build spent on), and loads as fast.
        seconds = run(build_command(args.protal, training_db, args.threads, "--compress_level", str(args.training_db_level)),
                      os.path.join(args.outdir, "training_db_index.log"))
        print(f"Built {training_db} in {seconds:.0f} s", flush=True)

    # Training data of every read type (pe, se from its first reads, pb and ont from long reads of the same
    # communities), then an independent test set of another design, both against the training database.
    def collect_command(out, samples, read_pairs, species, abundance, strains, long_bases, seed):
        command = [sys.executable, COLLECTOR, "--db", training_db, "--genome_table", genome_table, "-o", out,
                   "--protal", args.protal, "--simulator", args.simulator, "--samples", str(samples),
                   "--read_pairs", read_pairs, "--read_setups", args.read_setups, "--archaea", str(args.archaea),
                   "--species_per_sample", species, "--seed", str(seed), "-t", str(args.threads),
                   "--taxonomy", taxonomy, "--congeners", str(args.congeners), "--read_types", ",".join(read_types),
                   "--long_read_bases", long_bases, "--pb_setup", args.pb_setup, "--ont_setup", args.ont_setup,
                   "--pbsim", args.pbsim]
        command += ["--abundance", abundance] if abundance else []
        command += ["--strains_per_species", strains] if strains else []
        command += ["--pbsim_models", args.pbsim_models] if args.pbsim_models else []
        if n_heldout:
            command += ["--novel_species", heldout, "--novel_clades", str(args.novel_clades_per_sample)]
        return command

    training = os.path.join(args.outdir, "training")
    seconds = run(collect_command(training, args.samples, args.read_pairs, args.species_per_sample, args.abundance,
                                  args.strains_per_species, args.long_read_bases, args.seed),
                  os.path.join(args.outdir, "training_data.log"))
    print(f"Collected the training data ({', '.join(read_types)}) in {seconds:.0f} s", flush=True)
    test = os.path.join(args.outdir, "test")
    if args.test_samples > 0:
        seconds = run(collect_command(test, args.test_samples, args.test_read_pairs, args.test_species_per_sample,
                                      args.test_abundance, args.test_strains_per_species, args.test_long_read_bases,
                                      args.seed + 1000), os.path.join(args.outdir, "test_data.log"))
        print(f"Collected the independent test set in {seconds:.0f} s", flush=True)

    # One model per read type, trained in parallel.
    prefixes = {t: os.path.join(args.outdir, "trained_model" + ("" if t == "pe" else "_" + t)) for t in read_types}
    trainer_threads = max(1, args.threads // len(read_types))
    trainers = {}
    for t in read_types:
        command = [sys.executable, TRAINER, "--truth-file", os.path.join(training, TABLES[t]),
                   "--output-prefix", prefixes[t], "--features", "normalized", "--ntree", str(args.ntree),
                   "--maxnodes", str(args.maxnodes), "--seed", str(args.seed), "--threads", str(trainer_threads),
                   "--taxonomy", taxonomy, "--evaluation", args.evaluation]
        if args.test_samples > 0 and os.path.isfile(os.path.join(test, TABLES[t])):
            command += ["--test-file", os.path.join(test, TABLES[t])]
        trainers[t] = Background(command, os.path.join(args.outdir, "classifier_training" + ("" if t == "pe" else "_" + t) + ".log"))
    seconds = max(job.finish() for job in trainers.values())
    print(f"Trained the {', '.join(read_types)} model{'s' if len(read_types) > 1 else ''} in parallel in {seconds:.0f} s",
          flush=True)
    # protal must score as the trainer does, and compute the features as it did during collection.
    for t in read_types:
        run([sys.executable, PARITY, "--db", training_db, "--model", prefixes[t] + ".xml", "--training", training,
             "--read_type", t, "--protal", args.protal, "-t", str(args.threads)],
            os.path.join(args.outdir, "parity" + ("" if t == "pe" else "_" + t) + ".log"))
    for t in read_types:
        prefix = prefixes[t]
        parity = os.path.join(training, "parity" if t == "pe" else "parity_" + t, "parity.txt")
        for name in (prefix + ".report.txt", prefix + ".metrics.json", prefix + ".thresholds.tsv", prefix + ".varimp.tsv",
                     prefix + ".predictions.tsv.gz", prefix + ".test_predictions.tsv.gz"):
            if os.path.isfile(name):
                shutil.copy(name, logs)
        if os.path.isfile(parity):
            shutil.copy(parity, os.path.join(logs, "parity.txt" if t == "pe" else f"parity_{t}.txt"))
    for name in (os.path.join(args.outdir, "training_data.log"), os.path.join(args.outdir, "test_data.log"),
                 os.path.join(args.outdir, "genome_table.txt"), heldout):
        if os.path.isfile(name):
            shutil.copy(name, logs)
    if final_build is not None:
        print(f"Built {db} in {final_build.finish():.0f} s (in the background)", flush=True)
    elif training_db != db:
        print(f"Built {db} in {run(build_command(args.protal, db, args.threads), final_log):.0f} s", flush=True)
    # The trained models replace the shipped one and the placeholders in database.protal; --add_model checks
    # each and copies the other parts as they are.
    for t in read_types:
        run([args.protal, "--add_model", prefixes[t] + ".xml", "--read_type", t, "--db", db, "-t", str(args.threads)],
            os.path.join(args.outdir, "final_package" + ("" if t == "pe" else "_" + t) + ".log"))
    with open(os.path.join(db, "build_metadata.tsv"), "w") as fh:
        fh.write("".join(f"{k}\t{v}\n" for k, v in provenance(args, release, genome_table, heldout, n_heldout,
                                                                  read_types, prefixes)))
    shutil.copy(os.path.join(db, "build_metadata.tsv"), logs)
    summary = summary_lines(read_types, prefixes, db)
    with open(os.path.join(logs, "summary.txt"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print("\n".join(summary), flush=True)
    print(f"Ready protal database: {db}", flush=True)
    print(f"Model evaluation: {logs} (start with trained_model.report.txt, and trained_model_<read type>.report.txt)",
          flush=True)


if __name__ == "__main__":
    main()
